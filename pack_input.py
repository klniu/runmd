#!/usr/bin/python
# Filename: molpack.py
# 用于生成packmol运行所需要的输入文件。程序从ini配置文件中读取数据，并生成inp文件。
import argparse
import configparser
import os.path
import moltoolkit


def echoWarning(string):
    '''Give a colorful display for string as a warning'''
    warning = '\033[93m'
    endc = '\033[0m'
    print(warning + "Warning: " + string + endc)

def echoNote(string):
    '''Give a colorful display for string as a note'''
    note = '\033[93m'
    endc = '\033[0m'
    print(note + "Note: " + string + endc)

def echoError(string):
    '''Give a colorful display for string as a fail'''
    fail = '\033[91m'
    endc = '\033[0m'
    print(fail + "Error: " + string + endc)

class MolPack():
    '''输出一个种类分子的堆砌信息.

    输出一个种类分子的堆砌信息，可以大盒子分为小格子，其内放入一定数量的分子，也可以将分子限定在一定的水平线之上。
    '''
    def __init__(self, pack):
        '''初始化输出文件，各项信息。

        Args:
            f: 输出文件，对象类型为file
            pack: 一个包含需要堆砌的分子信息的字典结构，结构如：
                    {'pack1':
                        {'pdb':'m',
                         'box':(0, 0, 0, 10, 10, 10),
                         'fix': True/False
                         'membrane': True/False,
                         'alignNum': (2,3)
                         'takeout': 3
                         'num': 10,
                         'over_x':{10:[23,45],...},
                         'below_x':{10:[23,45],...},
                         'over_x':{10:[23,45],...},
                         'below_x':{10:[23,45],...},
                         'over_x':{10:[23,45],...},
                         'below_x':{10:[23,45],...}},
                    ...}
                    其中:
                    pdb: 分子pdb文件名.
                    box: 盒子尺寸，此外定义的6个坐标为长方体的对顶坐标。e.g. [0,0,0,10,10,10] 
                    fix: 是否固定这个pdb，一般用于已经堆砌完的结果需要继续插入其他pdb的情况，此参数以这个pdb的因有尺寸来设置其盒子尺寸，因此box和num都将失效。当此参数设置时,只有pdb与fix两个参数生效。
                    membrane: True/False, 表示是否是要堆砌单层膜, 如果是的话，表示要在一个单层上分割排列这些分子，分子之间在Z轴上是平行的。
                    alignNum: e.g. (5,4), 如果membrane为真，那么此alignNum表示在x和y轴方向上分子的个数。如下图：
                      y     .   .   .   .
                      y       .   .   .   .
                      y     .   .   .   .
                      y       .   .   .   .
                      y     .   .   .   .
                            x   x   x   x  
                            这是一个x*y=4*5的斜晶，4表示x轴一行排列4个分子，5表示有5列,共20个分子
                            如果设置了membrane为true，那么此值必须设置，且下面num值无效。
                    takeout: 是否從正斜晶裏去掉幾個原子，這個選項應用于膜堆砌時，例如6*5, 6*4之間的數字不能用正斜晶的方法，但可以用此選項去掉幾個，就可以了。只應用于membranePack，對randomPack和LinePack無效。

                    num: 表示分子的个数, 必须为整数，如果membrane为false，此值必须设置。
                    over_x: 使原子位于垂直于规定的x坐标平面之上，结构为{x坐标:[原子序号, ...], ...}， e.g. 使23,45号原子位于垂直于x=10的平面之上可以写为{10:[23,45]}
                    below_x: 使原子位于垂直于规定的x坐标平面之下,结构为{x坐标:[原子序号, ...], ...}, e.g. 使23,45号原子位于垂直于x=10的平面之下可以写为{10:[23,45]}
                    over_y
                    below_y
                    over_z
                    below_z: 以上四值与over_x、below_x类似.
        '''
        self.pack = pack

        # 检查pack内值是否合理
        if not self.validCheck():
            exit(1)

        # 检查并处理水平限制数据
        self.__analysePanel()


        # 获取pdb分子空间尺寸
        if self.pack['membrane']:
            if os.path.exists(self.pack['pdb']):
                # get the smallest and largest coordinates of molecule
                # moltoolkit可以获取分子的空间体积
                mol = moltoolkit.Mol(self.pack['pdb'])
            else:
                echoError('The file', self.pack['pdb'], 'is not exist in current directory')
                exit(1)
            self.pdbSize = mol.most_coordinates
            self.pdbSize.sort()
            echoNote('分子空间尺寸为{0[0]} {0[1]} {0[2]}, 请检查分配空间是否合理'.format(self.pdbSize))

        # 按行排列分子，此時不限定單個分子空間，而限定一行所有分子的空間，這樣可能有更大的靈活性, 默認為False，通過setLine可以更改此值。
        self.line = False

        # 按區排列分子，此時不限定x,y上單個分子空間，而只限定分子高度，在一個固定區域內堆砌所有分子，這樣可能有更大的靈活性, 默認為False，通過setLoose可以更改此值。
        self.loose = False

    def packMol(self):
        '''
        堆砌分子，并返回数据

        Return:
            self.result: packmol输出文件的字符串
        '''
        self.result = ''
        if self.pack['fix']:
            # 固定这个pdb
            self.FixPack()
        elif self.pack['membrane']:
            # 如果按行排列分子
            if self.line:
                self.membraneXLinePack()
            elif self.loose:
                #如果按區域排列 
                self.membraneRandomPack()
            else:
                # 如果以单层膜pack
                self.membranePack()
        else:
            self.normalPack()
        return self.result

    def validCheck(self):
        '''
        检查pack内各项值是否合理

        Return: 合理返回true, 否则返回false
        '''
        #必须要有的键
        keys_needed = ('pdb', 'box')
        for i in keys_needed:
            if not i in self.pack:
                echoError('%s键不存在于分子pack字典结构中' % i)
                return False

        #分子数量
        if self.pack['membrane']:
            if not isinstance(self.pack['alignNum'][0] * self.pack['alignNum'][1], int):
                echoError('分子数量必须为正整数')
                return False
        elif not isinstance(self.pack['num'],int):
            echoError('分子数量必须为正整数')
            return False

        return True

    def __analysePanel(self):
        '''添加使分子处于某一水平线的数据结构至dict，方便输出时调用.

        添加使分子处于某一水平线的数据结构至dict，方便输出时调用，生成之后的数据结构如：
            {'over_x':{line_number:[atoms_index],
             'below_x':{line_number:[atoms_index],
             'over_y':{line_number:[atoms_index],
             'below_y':{line_number:[atoms_index],
             'over_z':{line_number:[atoms_index],
             'below_z':{line_number:[atoms_index]}

        Arributes Changed:
            self.panels

        '''
        self.panels = {}
        #检验panel_type
        panel_types = ['over_x', 'below_x', 'over_y', 'below_y', 'over_z', 'below_z']
        for i in panel_types:
            if not i in tuple(self.pack.keys()):
                raise ValueError("分子堆砌数据结构中必须有%s键" % i)
            #跳过没有定义的水平线限定值
            if len(self.pack[i]) == 0:
                continue
            panel_dict = self.pack[i]
            #检验坐标和原子序号的合理性
            for j, k in panel_dict.items():
                if not isinstance(j, (int, float)):
                    echoError('水平線坐标的值必須為數字,此处为%s' % j)
                    exit(1)
                if len(k) == 0:
                    echoError('当需要将原子位于水平线上或下时，原子序号不能为空')
                    exit(1)
                for l in k:
                    if (not isinstance(l, int)) or l < 0:
                        echoError('原子序号值必須為正整数')
                        exit(1)
            self.panels[i] = panel_dict

    def FixPack(self):
        '''
        固定这个pdb，用pdb原有的空间尺寸为依据，固定pdb

        @Member Changed:
            self.result: 结果的字符串
        '''
        self.result += 'structure ' + self.pack['pdb'] + '\n\tnumber 1\n\tfixed 0. 0. 0. 0. 0. 0.\nend structure\n\n'

    def membranePack(self):
        '''
        在一层上按照正斜晶堆砌分子，分子之间平行，可以设置分子的位置限定，即在x,y,z轴上分子位于水平线上下的位置
           .       .    
               .       .
           .       .    
               .       .
           .       .    
               .       .
           .       .    
               .       .
        这是一个2*8的排列，共16个分子
        '''
        # 要去掉的原子個數
        takeout = self.pack['takeout']
        echoNote('将在{0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]}的盒子内排列 {1[0]}*{1[1]}-{2} 个分子.'.format(self.pack['box'], self.pack['alignNum'], takeout))

        # box size
        boxSize = [self.pack['box'][i + 3] -self.pack['box'][i] for i in range(3)]

        # 分子的空间尺寸
        pdbSize = self.pdbSize

        #以下是排布的算法
        #   p s p s p s p s2 # p为pdb边长，s为分子之间间隔，s2为尾端或第二行的顶端的空隙.y轴的排布与此相同，n是一行上分子的数量，x是边长，这些变量的关系是：
        # y .   .   .   .    # n*p + (n-1)s + 0.5*(p+s) = x
        # y   .   .   .   .  #
        # y .   .   .   .    #
        # y   .   .   .   .  #
        # y .   .   .   .    #
        #   x  x   x   x  
        #  这是一个x*y=4*5的斜晶，4表示x轴一行排列4个分子，5表示有5列,共20个分子, 建议y轴设为偶数

        yCoord = self.__getYCoordIterator(self.pack['box'][1], boxSize[1], pdbSize[1], self.pack['alignNum'][1])
        for (row, y1, y2) in yCoord:
            #是否是偶数行，行从0开始
            even = True
            if row % 2 != 0:
                even = False
            xCoord = self.__getXCoordIterator(self.pack['box'][0], boxSize[0], pdbSize[0], self.pack['alignNum'][0], even)
            y1 = y1 if y1 - 0.2 < self.pack['box'][1] else y1 - 0.2
            y2 = y2 if y2 + 0.2 > self.pack['box'][4] else y2 + 0.2
            for column, x1, x2 in xCoord:
                if takeout > 0 and column == 0:
                    takeout -= 1
                    continue
                # 稍微讓分子空間寬鬆一點，它才能好好的站整齊
                x1 = x1 if x1 - 0.2 < self.pack['box'][0] else x1 - 0.2
                x2 = x2 if x2 + 0.2 > self.pack['box'][3] else x2 + 0.2
                #输出开头
                self.result += 'structure %s\n' % self.pack['pdb']
                #输出分子数量
                self.result += '\tnumber 1\n'
                #输出盒子尺寸: inside box
                self.result += '\tinside box'
                z1 = self.pack['box'][2]
                z2 = self.pack['box'][5]
                for m in x1, y1, z1, x2, y2, z2:
                    if isinstance(m, int):
                        self.result += ' ' + str(round(m)) + '.'
                    elif isinstance(m, float):
                        self.result += ' ' + str(round(m, 1))
                    else:
                        echoError('Error: 只能用数字来表示坐标或尺寸')
                self.result += '\n'

                # 限制水平位置
                self.__setPanel(x1, y1, z1, x2, y2, z2)

                #输出结尾
                self.result += 'end structure\n\n'

    def __setPanel(self, *pdbVol):
        '''
        输出panel部分

        这里存在一个margin值，表示分子距离顶端的距离，比如分子在x方向上的最大坐标为10，最小坐标为5，1号原子如果要位于x=7的平面之上，那么，margin值就要是3; 如果要位于x=7的平面之下，那么，margin值要为2，，over将以最大值-3，below将以最小值+2
        '''
        #输出panel部分
        constraint = (k for k in ('over_x', 'below_x', 'over_y', 'below_y', 'over_z', 'below_z') if k in self.panels)
        for k in constraint:
            # self.panels[k]仍为一dict
           for margin, indices in self.pack[k].items():
               self.result += '\tatoms'
               for m in indices:
                   self.result += ' ' + str(m)
               #如果k的首字母为o，即表示为over_*，因此输出over plane
               if k[0] == 'o':
                   self.result += '\n\t\tover plane '
                   #如果k的尾字母为x，即表示平面为x，因此将输出1. 0. 0.
                   if k[-1] == 'x':
                       self.result += '1. 0. 0. ' + str(round(pdbVol[3] - margin))
                   elif k[-1] == 'y':
                       self.result += '0. 1. 0. ' + str(round(pdbVol[4]  - margin))
                   else:
                       self.result += '0. 0. 1. ' + str(round(pdbVol[5] - margin))
               else:
                   self.result += '\n\t\tbelow plane '
                   if k[-1] == 'x':
                       self.result += '1. 0. 0. ' + str(round(pdbVol[0] + margin))
                   elif k[-1] == 'y':                              
                       self.result += '0. 1. 0. ' + str(round(pdbVol[1] + margin))
                   else:                                           
                       self.result += '0. 0. 1. ' + str(round(pdbVol[2] + margin))
               self.result += '\n\tend atoms\n'

    def __getXCoordIterator(self, initCoord, boxLen, pdbLen, alignNum, even=True):
        '''
        计算在x方向上排布分子的初始与结束坐标，可以定义分子的空间位置

        Arguments:
            initCoord: 初始坐标，如x方向上
            boxLen: 可以排布分子的总长度
            pdbLen: 分子在一方向上的尺寸，如x方向上
            alignNum： 分子在这一方向上排布分子的个数
            even: 是否是偶数行，偶数行在分子开始排布时没有space, 奇数行有1/2的(pdb + space)在前，以形成斜晶形

        Return:
            coord: 行數，初始坐标，結束坐標
        '''
        # 检验盒子尺寸够不够排放分子
        if boxLen - pdbLen * (alignNum + 0.5) < 0:
            echoWarning('x边长{0}不足以容下{1}个尺寸为{2}的分子'.format(boxLen, alignNum, pdbLen))

        # 分子之间间隔
        space = (boxLen - (alignNum + 0.5) * pdbLen) / (alignNum - 0.5)
        #
        if space < 1.5:
            echoWarning('x轴上两分子之间距离%s小于1.5，堆砌时间会很长，且很有可能不是最优化结果' % space)

        #第一个坐标为盒子开始
        coord = initCoord if even else initCoord + 0.5 * (space + pdbLen)
        index = 0
        yield (index, coord, coord + pdbLen)
        for i in range(1, alignNum):
            coord += pdbLen + space
            index += 1
            yield (index, coord, coord + pdbLen)

    def __getYCoordIterator(self, initCoord, boxLen, pdbLen, alignNum):
        '''
        计算在一个方向上排布分子的初始与结束坐标，可以定义分子的空间位置

        Arguments:
            initCoord: 初始坐标，如x方向上
            boxLen: 可以排布分子的总长度
            pdbLen: 分子在一方向上的尺寸，如x方向上
            alignNum： 分子在这一方向上排布分子的个数

        Return:
            coord: 行數，初始坐标，結束坐標
        '''
        # 检验盒子尺寸够不够排放分子
        if boxLen - pdbLen * alignNum / 2  < 0:
            echoWarning('y边长{0}不足以容下{1}个尺寸为{2}的分子'.format(boxLen, alignNum, pdbLen))

        # 單列分子數量
        singleLineNum = alignNum // 2 if alignNum % 2 ==0 else alignNum // 2 + 1
        # 分子之间间隔
        space = (boxLen - (singleLineNum + 0.5) * pdbLen) / (singleLineNum - 0.5)
        # 如果排列的是奇數，表示偶數列的最後一個分子是沒有的，需要把這部分距離添加到space上
        if alignNum % 2 != 0:
            space += (space + pdbLen) / (2 * ((alignNum + 1) / 2 - 1))
        # 计算x轴上分子间距离，以测定y上两个分子的距离，以便输出提醒结果
        xSpace = (self.pack['box'][3] - self.pack['box'][0] - (self.pack['alignNum'][0] + 0.5) * self.pdbSize[0]) / (self.pack['alignNum'][0] - 0.5)
        xySpace = pow(pow(xSpace / 2, 2) + pow(space / 2, 2), 0.5)
        if xySpace < 2.5 and xSpace < 2.5:
            echoWarning('分子四周间隔距离%s小于2.5，堆砌时间会很长，且很有可能不是最优化结果' % xySpace)
        if space < 1.5:
            echoWarning('y轴上两分子之间距离%s小于1.5，堆砌时间会很长，且很有可能不是最优化结果' % space)

        for i in range(0, alignNum):
            # y轴合并排列，如:
            # 3   .  initCoord + 0.5 * (space + pdbLen) + space + len
            # 2 .    initCoord + space + pdbLen
            # 1   .  initCoord + 0.5 * (space + pdbLen)
            # 0 .    initCoord
            # space为一列上两分子之间的距离
            if i % 2 == 0:
                coord = initCoord + i / 2 * (pdbLen + space)
            else:
                coord = initCoord + 0.5 * (space + pdbLen) + (i - 1) / 2 * (pdbLen + space)
            yield (i, coord, coord + pdbLen)

    def setLine(self, line):
        '''
        是否按行排列分子
        按行排列分子，此時不限定單個分子空間，而限定一行所有分子的空間，這樣可能有更大的靈活性, 默認為False，通過setLine可以更改此值。
        '''
        self.line = line

    def membraneXLinePack(self):
        '''
        在一层上按照正斜晶堆砌分子，分子之间平行，可以设置分子的位置限定，即在x,y,z轴上分子位于水平线上下的位置
           .       .    
               .       .
           .       .    
               .       .
           .       .    
               .       .
           .       .    
               .       .
        这是一个2*8的排列，共16个分子
        按行排列分子，此時不限定單個分子空間，而限定一行所有分子的空間，這樣可能有更大的靈活性, 默認為False，通過setLine可以更改此值。
        此時，在x方向上將不能再使用水平線限定
        '''
        echoNote('将在{0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]}的盒子内排列{1[0]}*{1[1]}个分子.'.format(self.pack['box'], self.pack['alignNum']))

        # box size
        boxSize = [self.pack['box'][i + 3] -self.pack['box'][i] for i in range(3)]

        # 分子的空间尺寸
        pdbSize = self.pdbSize
        

        #以下是排布的算法
        #   p s p s p s p s2 # p为pdb边长，s为分子之间间隔，s2为尾端或第二行的顶端的空隙.y轴的排布与此相同，n是一行上分子的数量，x是边长，这些变量的关系是：
        # y .   .   .   .    # n*p + (n-1)s + 0.5*(p+s) = x
        # y   .   .   .   .  #
        # y .   .   .   .    #
        # y   .   .   .   .  #
        # y .   .   .   .    #
        #   x  x   x   x  
        #  这是一个x*y=4*5的斜晶，4表示x轴一行排列4个分子，5表示有5列,共20个分子, 建议y轴设为偶数

        yCoord = self.__getYCoordIterator(self.pack['box'][1], boxSize[1], pdbSize[1], self.pack['alignNum'][1])
        for (row, y1, y2) in yCoord:
            #是否是偶数行，行从0开始
            even = True
            if row % 2 != 0:
                even = False
            x1, x2 = self.__getXLineCoord(self.pack['box'][0], boxSize[0], pdbSize[0], self.pack['alignNum'][0], even)
            #输出开头
            self.result += 'structure %s\n' % self.pack['pdb']
            #输出分子数量
            self.result += '\tnumber ' + str(self.pack['alignNum'][0]) + '\n'
            #输出盒子尺寸: inside box
            self.result += '\tinside box'
            z1 = self.pack['box'][2]
            z2 = self.pack['box'][5]
            for m in x1, y1, z1, x2, y2, z2:
                if isinstance(m, int):
                    self.result += ' ' + str(round(m)) + '.'
                elif isinstance(m, float):
                    self.result += ' ' + str(round(m, 1))
                else:
                    echoError('Error: 只能用数字来表示坐标或尺寸')
            self.result += '\n'

            # 限制水平位置
            for panel in ['over_x', 'below_x']:
                if panel in self.panels:
                    echoError('使用按行排列膜分子時，不能再在x方向上限定分子位置，程序將退出')
                    exit(1)
            self.__setPanel(x1, y1, z1, x2, y2, z2)

            #输出结尾
            self.result += 'end structure\n\n'

    def __getXLineCoord(self, initCoord, boxLen, pdbLen, alignNum, even=True):
        '''
        计算在x方向上去除邊界空白之後這一行分子排布的初始与结束坐标，可以定义分子的空间位置

        Arguments:
            initCoord: 初始坐标，如x方向上
            boxLen: 可以排布分子的总长度
            pdbLen: 分子在一方向上的尺寸，如x方向上
            alignNum： 分子在这一方向上排布分子的个数
            even: 是否是偶数行，偶数行在分子开始排布时没有space, 奇数行有1/2的(pdb + space)在前，以形成斜晶形

        Return:
            coord: 初始坐标，結束坐標
        '''
        # 检验盒子尺寸够不够排放分子
        if boxLen - pdbLen * (alignNum + 0.5) < 0:
            echoWarning('x边长{0}不足以容下{1}个尺寸为{2}的分子'.format(boxLen, alignNum, pdbLen))

        # 分子之间间隔
        space = (boxLen - (alignNum + 0.5) * pdbLen) / (alignNum - 0.5)
        #
        if space < 1.5:
            echoWarning('x轴上两分子之间距离%s小于1.5，堆砌时间会很长，且很有可能不是最优化结果' % space)

        #第一个坐标为盒子开始
        begin_coord = initCoord if even else initCoord + 0.5 * (space + pdbLen)
        end_coord = begin_coord + alignNum * (pdbLen + space)
        return (begin_coord, end_coord)

    def setLoose(self, loose):
        '''
        是否按區域排列分子
        按區域排列分子，此時不限定單個分子空間，而限定一定x y區域所有分子的空間，這樣可能有更大的靈活性, 默認為False，通過setLoose可以更改此值。
        '''
        self.loose= loose

    def membraneRandomPack(self):
        '''
        在一层上隨機堆砌分子，分子之间平行，可以设置分子在Z方向上位置限定，即在z轴上分子位于水平线上下的位置
        按區域排列分子，此時不限定單個分子空間，而限定一定x y區域所有分子的空間，這樣可能有更大的靈活性, 默認為False，通過setLoose可以更改此值。
        此時，在x, y方向上將不能再使用水平線限定
        '''
        echoNote('将在{0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]}的盒子内排列{1[0]}*{1[1]}个分子.'.format(self.pack['box'], self.pack['alignNum']))

        # box size
        #boxSize = [self.pack['box'][i + 3] -self.pack['box'][i] for i in range(3)]

        # 分子的空间尺寸
        #pdbSize = self.pdbSize
        

        #以下是排布的算法
        #   p s p s p s p s2 # p为pdb边长，s为分子之间间隔，s2为尾端或第二行的顶端的空隙.y轴的排布与此相同，n是一行上分子的数量，x是边长，这些变量的关系是：
        # y .   .   .   .    # n*p + (n-1)s + 0.5*(p+s) = x
        # y   .   .   .   .  #
        # y .   .   .   .    #
        # y   .   .   .   .  #
        # y .   .   .   .    #
        #   x  x   x   x  
        #  这是一个x*y=4*5的斜晶，4表示x轴一行排列4个分子，5表示有5列,共20个分子, 建议y轴设为偶数

        x1, x2 = self.pack['box'][0], self.pack['box'][3]
        y1, y2 = self.pack['box'][1], self.pack['box'][4]
        #输出开头
        self.result += 'structure %s\n' % self.pack['pdb']
        #输出分子数量
        self.result += '\tnumber ' + str(self.pack['alignNum'][0] * self.pack['alignNum'][1]) + '\n'
        #输出盒子尺寸: inside box
        self.result += '\tinside box'
        z1 = self.pack['box'][2]
        z2 = self.pack['box'][5]
        for m in x1, y1, z1, x2, y2, z2:
            if isinstance(m, int):
                self.result += ' ' + str(round(m)) + '.'
            elif isinstance(m, float):
                self.result += ' ' + str(round(m, 1))
            else:
                echoError('Error: 只能用数字来表示坐标或尺寸')
        self.result += '\n'

        # 限制水平位置
        for panel in ['over_x', 'below_x', 'over_y', 'below_y']:
            if panel in self.panels:
                echoError('使用随机排列膜分子時，不能再在x, y方向上限定分子位置，程序將退出')
                exit(1)
        self.__setPanel(x1, y1, z1, x2, y2, z2)

        #输出结尾
        self.result += 'end structure\n\n'

    def normalPack(self):
        '''
        在一个盒子内堆砌分子，不应用分子位置限定
        '''
        # 如果是膜堆砌，跳过
        if self.pack['membrane']:
            return ''

        echoNote('将在 %s 内堆砌 %d 个 %s 分子 ' % (self.pack['box'], self.pack['num'], self.pack['pdb']))
        self.result += 'structure {0}\n\tnumber {1}\n\tinside box {2[0]} {2[1]} {2[2]} {2[3]} {2[4]} {2[5]}\nend structure\n\n'.format(self.pack['pdb'], self.pack['num'], self.pack['box'])


class ReadConfig(configparser.ConfigParser):
    '''用于读取ini配置文件，此ini包含了生成inp所需要的分子、原子及盒子信息，其中有些参数必须定义值。'''
    def __init__(self, ini_file):
        #ini_file是一个file类型的值，不是字符串
        #ini_name = ini_file.name
        #初始化父类并将行内注释分隔符定义为#与;
        super(ReadConfig, self).__init__(inline_comment_prefixes=('#', ';'))
        #读取ini文件
        try:
            self.read_file(ini_file)
        except:
            echoError(ini_file.name + ' IO读写错误')
            exit(1)
        #检查[packmol]下有几种分子需要堆砌，即
        if not self.has_section('packmol'):
            raise configparser.NoSectionError('packmol')
        elif not self.has_option('packmol', 'pack_num'):
            raise configparser.NoOptionError('pack_num', 'packmol')
        elif not self['packmol']['pack_num'].isdigit():
            raise ValueError(
                    'ini配置文件中pack_num值不正确，期望为正整数')
        else:
            #需要堆砌的分子种类数量
            self.__pack_num = self.getint('packmol', 'pack_num')
        if not self.has_option('packmol', 'last_pdb'):
            raise configparser.NoOptionError('last_pdb', 'packmol')
        else:
            #最后生成的pdb的文件名
            self.__pack_last_pdb = self.get('packmol', 'last_pdb')

        '''
        定义一个容纳堆砌所有数值的dict，结构为
        {'pack1':
            {'pdb':'m',
             'box':[0, 0, 0, 10, 10, 10],
             'split':[1, 1],
             'pdb_num': 10,
             'over_x_panel':{10:[23,45],...},
             'below_x_panel':{10:[23,45],...},
             'over_x_panel':{10:[23,45],...},
             'below_x_panel':{10:[23,45],...},
             'over_x_panel':{10:[23,45],...},
             'below_x_panel':{10:[23,45],...}},
        ...}
        '''
        self.packs = {}

        #检查[pack*]，(*为数字)下的各项是否有错误
        for i in range(self.__pack_num):
            section = 'pack' + str(i)
            if not self.has_section(section):
                raise configparser.NoSectionError('pack%s' % str(i))
            for j in ['pdb']:
                if not self.has_option(section, j):
                    raise configparser.NoOptionError(j, 'pack%s' % str(i))
            self.packs[section] = {}
            # pdb文件名
            self.packs[section]['pdb'] = self.get(section, 'pdb')
            # 盒子尺寸
            self.packs[section]['box'] = [float(k) for k in self.get(section, 'box', fallback='0 0 0 0 0 0').split()]
            # 是否fix这个pdb
            self.packs[section]['fix'] = True if self.getboolean(section, 'fix', fallback=0) else False
            # 是否是膜堆砌
            self.packs[section]['membrane'] = True if self.getboolean(section, 'membrane', fallback=0) else False
            # 是否減掉幾個原子
            self.packs[section]['takeout'] = self.getint(section, 'takeout', fallback=0)
            # x, y方向上排列的分子数量
            if self.packs[section]['membrane']:
                self.packs[section]['alignNum'] = [int(k) for k in self[section]['alignNum'].split()]
            # 分子数量
            if self.packs[section]['fix']:
                self.packs[section]['num'] = 1
            else:
                self.packs[section]['num'] = self.getint(section, 'num', fallback=0)

            if len(self.packs[section]['box']) != 6:
                raise ValueError('pack%s中box值不正确' % str(i))
            
            if self.packs[section]['fix'] and self.packs[section]['membrane']:
                echoError('fix和membrane不能同时定义')
                exit(1)

            for l in ['over', 'below']:
                for m in ['x', 'y', 'z']:
                    panel_type = l + '_' + m
                    panel_temp = self.get(section, panel_type, fallback='')
                    if len(panel_temp) > 0:
                        #使字符串数字化
                        self.packs[section][panel_type] = eval(panel_temp)
                    else:
                        #如果over_x_panel等定义为空，则字典赋值为空
                        self.packs[section][panel_type] = {}

    def get_pack_num(self):
        '''返回需要堆砌的分子的种类，对应于下面的[pack*]section, *为0,1...'''
        return self.__pack_num

    def get_pack_last_pdb(self):
        '''返回inp中定义的output语句中的pdb, 即packmol 最后要生成的pdb的名字'''
        return self.__pack_last_pdb

    def get_packs(self):
        '''
        返回所有堆砌分子的键值的字典,结构为
        {'pack1':
            {'pdb':'m',
             'box':[0, 0, 0, 10, 10, 10],
             'fix': True/False
             'membrane': True/False
             'alignNum':[1, 1],
             'num': 10,
             'over_x':{10:[23,45],...},
             'below_x':{10:[23,45],...},
             'over_x':{10:[23,45],...},
             'below_x':{10:[23,45],...},
             'over_x':{10:[23,45],...},
             'below_x':{10:[23,45],...}},
        ...}
        '''
        return self.packs


def main():
    arg_parser = argparse.ArgumentParser(
                    description='通过读取ini配置文件，生成packmol所需的配置文件。')
    arg_parser.add_argument('-i', '--input',
                            type=argparse.FileType('r'), required=True, help='读取ini配置文件')
    arg_parser.add_argument('-o', '--output',
                            type=argparse.FileType('w'), required=True, help='生成inp文件')
    arg_parser.add_argument('-l', '--loop', type=int, action='store', default=50000, help='nloop值，指最大迭代次数')
    arg_parser.add_argument('-t', '--tolerance', action='store', help='tolerance值，公差，默认值为2.0，减少公差可以使分子结合更紧密，但过小容易让分子堆叠, 此值不建议小于1.5')
    arg_parser.add_argument('--line', action='store_true', help='在堆砌膜時，不以單個分子而按行排列分子，這樣可以有更大的靈活，堆砌時長過長時可以選擇此選項')
    arg_parser.add_argument('--loose', action='store_true', help='在堆砌膜時，不以單個分子而在x,y區域內隨機排列，這樣可以有更大的靈活，而且此時只考慮分子的總數量x*y，而不是分別考慮x, y上的分子數量，堆砌時長過長時可以選擇此選項')
    args = arg_parser.parse_args()
    config_reader = ReadConfig(args.input)

    f = args.output
    if args.loop:
        loop = args.loop
    tolerance = '2.0'
    if args.tolerance:
        tolerance = args.tolerance

    packs_dict = config_reader.get_packs()
    f.write('tolerance ' + tolerance + ' \noutput ' + config_reader.get_pack_last_pdb() + '\nfiletype pdb\ndiscale 1.5\nnloop ' + str(loop) + '\n\n')
    for i in sorted(list(packs_dict.keys())):
        mol_pack = MolPack(packs_dict[i])
        if args.line:
            mol_pack.setLine(args.line)
        if args.loose:
            mol_pack.setLoose(args.loose)
        f.write(mol_pack.packMol())

if __name__ == '__main__':
    main()
