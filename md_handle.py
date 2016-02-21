#!/usr/bin/python
# -*- coding: UTF-8 -*-
'''@目的:  生成gromacs模拟所需要的命令文件，top文件(此top文件并非生成单个分子的itp，而是将各个itp组合，并添加力场、位置限制及[ system ]信息。)，及mdp文件。
@用法:  md_handle.py 选项 ini配置文件
@选项:
        -t, --top       生成top文件
        -c, --command   生成命令文件
        -m, --mdp       生成mdp文件'''

import configparser, sys, os, argparse
import copy


class ReadConfig(configparser.ConfigParser):
    '''用于读取ini配置文件，此ini包含了分子模拟所需要的参数，其中有些参数必须定义值。
    配置文件分四部分：
    [general]下是通用的一些信息，一般不需要改动
    [input]下是模拟的输入信息
    [md]是关于模拟过程的一些配置
    余下的有一些[npt][nvt][ext]等，与[md]中md_sequence对应，表示具体一个模拟小步骤。
    [md-default]与[em-default]分别表示md与能量最小化的的mdp文件中的默认值。
    
    @Attribute
    self.md: dict, [md]下各键值
    self.mds: dict, {md步骤:{key:value},...}md步骤一般为emsteep, emcg, npt, nvt, ext, npt1等下面各节....
    self.md_sequence: list, md步骤
    self.em_sequence: list, em步骤
    self.ions: dict, {ions name:ions num}
    self.general: dict, [general]下各键值
    self.input: dict, [input]下各键值
    self.residure: dict, {残基名称:分子数量}
    self.gmx_cmds: dict, {命令名称:命令}, 表示单双精度的gmx命令
    '''
    
    def __init__(self, iniPath):
        '''从配置文件读取各节的键值并进行检查与处理

        @Args:
        iniPath: 配置文件的路径
        '''
        if not os.path.exists(iniPath):
           echoFail('配置文件不存在!')
           exit(1)
        # ini文件默认必须具有的section列表
        ini_default_sections = ['general', 'input', 'md', 'em-default', 'md-default']
        # [general]下必须具有的键
        ini_general_keys = ['precision']
        # [input]下必须具有的键
        ini_input_keys = ['pdb', 'top', 'forcefield', 'box', 'itp_files', 'residure_name', 'residure_num', 'residure_itp', 'ions_name', 'ions_num']
        # [md]下必须具有的键
        ini_md_keys = ['em_sequence', 'md_sequence', 'maxwarn', 'pr_residure']

        #初始化父类并将行内注释分隔符定义为#与;
        super(ReadConfig, self).__init__(inline_comment_prefixes=('#', ';'))

        #读取ini文件
        try:
            self.read(iniPath)
        except configparser.DuplicateSectionError:
            echoFail('Error: 重复的Section, 请检查配置文件')
            exit(1)
        except configparser.DuplicateOptionError:
            echoFail('Error: 重复的Key, 请检查配置文件')
            exit(1)
        except:
            echoFail(iniPath, '配置文件读取错误')
            exit(1)
        #检查ini配置文件是否有上述节与键
        for i in ini_default_sections:
            if not self.has_section(i):
                echoFail(iniPath.name + ' 配置文件有错，无[' + i + ']')
                exit(1)
        for i in ['general', 'input', 'md']:
            ini_keys=locals()['ini_' + i + '_keys']
            for j in ini_keys:
                if not self.has_option(i, j):
                    echoFail(iniPath.name + ' 配置文件有错误，[' + i + ']节下无' + j + '键')
                    exit(1)

        #从ini文件中读取各节
        self.general = {k: v for (k, v) in self['general'].items()} 
        self.input ={k: v for (k, v) in self['input'].items()} 
        self.md = {k: v for (k, v) in self['md'].items()} 

        # 初始化残基
        self.__init_residure()
        # 初始化离子
        self.__init_ions()
        # 初始化能量最小化
        self.__init_em()
        # 初始化md
        self.__init_md()
        # 初始化位置限定
        self.__init_position_restraint()
        for i in self.md_sequence:
            if self.mds[i]['md_type'] != 'ext':
                # 初始化控溫分組
                self.__make_tcoupl_group(i)
                # 对mdp配置进行检查
                self.__check_mdp_value_valid(i)
        # 设置单双精度命令
        self.__init_precise()

    def __init_residure(self):
        '''从配置文件读取残基并检查
        
        @Changed
        self.residure: dict, 保存残基的名称及数量
        '''
        #残基名称及数量及itp文件
        residure_name = self.get('input', 'residure_name').split()
        residure_num = self.get('input', 'residure_num').split()
        residure_itp = self.get('input', 'residure_itp').split()
        if len(residure_name) != len(residure_num) != len(residure_itp):
            echoFail('ini配置文件中残基名称与残基数量或残基对应的itp文件个数不对应')
            exit(1)
        self.residure = []
        for i in range(len(residure_name)):
            if not residure_num[i].isdigit():
                echoFail('ini配置文件中残基数量不正确，期望为正整数')
                exit(1)
            self.residure.append((residure_name[i], int(residure_num[i])))

    def __init_ions(self):
        '''从配置文件中读取离子名称及数量

        @Changed
        self.ions: dict, 保存残基的名称及数量
        '''
        #离子名称及数量
        ions_name = self.get('input', 'ions_name').split()
        ions_num = self.get('input', 'ions_num').split()
        if len(ions_name) != len(ions_num):
            echoFail('ini配置文件中离子类型与离子数量个数不对应')
            exit(1)
        self.ions = {}
        for i in range(len(ions_name)):
            if not ions_num[i].isdigit():
                echoFail('ini配置文件中离子数量不正确，期望为正整数')
                exit(1)
            self.ions[ions_name[i]] = int(ions_num[i])

    def __init_em(self):
        '''从配置文件中获取能量最小化各步骤及参数的配置
        [md]下的em_sequence下储存了能量最小化时的步骤，比如emsteep,emcg，这些步骤分别对应一节，那节下保存了进行这个能量最小化的配置，而此处要从[em-default]中提取一些通用的选项，如有重复，则覆盖默认值。

        @Changed
        self.em_sequence: list, 能量最小化的所有步骤
        self.ems: {em:{key:value}}, dict, 记录了所有能量最小化的配置
        '''
        self.em_sequence = self.get('md', 'em_sequence').split()
        self.ems = {}
        for j in self.em_sequence:
            if not self.has_section(j):
                echoFail('无[' + j + ']section')
                exit(1)
            #处理default键值，如果default不为空，则从default对应的值所表示的section复制所有键值.
            default = self.get(j, 'default')
            if len(default) != 0:
                #检验是否存在default所对应的section
                if default in self.sections():
                    self.__copy_key_value(default, j)
            #处理em-default, 将em中未定义的键值均复制来过
            self.__copy_key_value('em-default', j)
            #记录em_sequence中各个md的键值
            self.ems[j] = {k: v for (k, v) in self[j].items()}
 
        #检查mdp键值的合理性
        for i in self.em_sequence:
            #检查键值的存在
            key_exist = ('md_type', 'nodes', 'mdp', 'analysis_dir',
                    'integrator', 'nsteps', 'dt', 'emtol', 'emstep')
            for j in key_exist:
                if not j in self.ems[i].keys():
                    echoFail('[' + i + ']section下无' + j + '键值')
                    exit(1)

    def __init_md(self):
        '''从配置文件中获取实际模拟各步骤及参数的配置
        [md]下的em_sequence下储存了模拟时的步骤，比如npt,nvt，这些步骤分别对应一节，那节下保存了进行这个模拟的配置，而此处要从[md-default]中提取一些通用的选项，如有重复，则覆盖默认值。

        @Changed
        self.md_sequence: list, 模拟的所有步骤
        self.mds: {md:{key, value}}, dict, 记录了模拟过程的所有配置
        '''
        self.md_sequence = self.get('md', 'md_sequence').split()
        self.mds = {}
        for j in self.md_sequence:
            if not self.has_section(j):
                echoFail('无[' + j + ']section')
                exit(1)
            #处理default键值，如果default不为空，则从default对应的值所表示的section复制所有键值.
            default = self.get(j, 'default')
            if len(default) != 0:
                #检验是否存在default所对应的section
                if default in self.sections():
                    self.__copy_key_value(default, j)
            #处理md-default, 将md中未定义的键值均复制来过, 延长的步骤除外
            if self.get(j, 'md_type') != 'ext':
                self.__copy_key_value('md-default', j)
            #记录md_sequence中各个md的键值
            self.mds[j] = {k: v for (k, v) in self[j].items()}

    def __init_position_restraint(self):
        '''初始化位置限定模拟的配置。
        如果一个包含模拟配置的节中定义有pr_residure，表示此步模拟是位置限定性模拟，因此要在self.mds中修改
        @Changed
        self.mds: {md:{key, value}}, dict, 记录了模拟过程的所有配置
        '''
        for i in self.md_sequence:
            # if this section has position_restraint key, it means this section
            # is to do position restraint, so the mdp key 'define=' should be
            # define according pr_residure
            if ('position_restraint' in list(self.mds[i].keys()) and
                self.mds[i]['position_restraint'] == 'yes'):
                self.mds[i]['define'] = '-DPOSRES_' + '_'.join(self.md['pr_residure'].split())

    def __make_tcoupl_group(self, mdSection):
        '''将各个残基分组控温
        
        将影响mdp中的tc-grps, tau_t, ref_t

        此处将水与离子分为一组，其他的各个残基分别分为一组
        @Args:
        mdSection: self.md_sequence中的元素
        '''
        residures = self.get('input', 'residure_name').split()
        if 'SOL' in residures:
            residures.remove('SOL')
            # If ions are added into system, the group should contain  water and all
            if len(self.ions) > 0:
                # This group name need to be generated by make_ndx programme
                water_ions_group = 'SOL_' + '_'.join(sorted(self.ions.keys()))
            else:
                water_ions_group = 'SOL'
            residures.append(water_ions_group)
        self.mds[mdSection]['tc-grps'] = ' '.join(residures)
        self.mds[mdSection]['tau_t'] = ((self.mds[mdSection]['tau_t'] + ' ') * len(residures)).strip()
        self.mds[mdSection]['ref_t'] = ((self.mds[mdSection]['ref_t'] + ' ') * len(residures)).strip()

    def __check_mdp_value_valid(self, mdSection):
        '''检查mdp键值的合理性
        @Args:
        mdSection: self.md_sequence中的元素
        '''
        # The key 'continuation' should be 'yes' and 'gen_vel' should be 'no' in the first course.
        if self.md_sequence.index(mdSection) == 0:
            if 'continuation' in self.mds[mdSection].keys() and self.mds[mdSection]['continuation'] == 'yes':
                echoFail('['+mdSection+']下continuation应为no') 
                exit(1)
            if 'gen_vel' in self.mds[mdSection].keys() and self.mds[mdSection]['gen_vel'] == 'no':
                echoFail('['+mdSection+']下gen_vel应为yes') 
                exit(1)
        # The key 'continuation' should be 'no' and 'gen_vel' should be 'yes' after the first course.
        else:
            if 'continuation' in self.mds[mdSection].keys() and self.mds[mdSection]['continuation'] != 'yes':
                echoFail('['+mdSection+']下continuation应为yes') 
                exit(1)
            if 'gen_vel' in self.mds[mdSection].keys() and self.mds[mdSection]['gen_vel'] != 'no':
                echoFail('['+mdSection+']下gen_vel应为no') 
                exit(1)

        # If there exists a key called 'position_restraint' whose value is 'yes', the value of the keys 'pr_residure' and 'pr_pdb' must exist.
        if 'position_restraint' in self.mds[mdSection].keys() and self.mds[mdSection]['position_restraint'] == 'yes':
            if len(self.md['pr_residure']) == 0 or len(self.md['pr_pdb']) == 0:
                echoFail('You want to use position restraint, but the values of pr_residure and pr_pdb are empty')
                exit(1)

        #检查键值的存在
        key_exist = ('md_type', 'nodes', 'mdp', 'analysis_dir', 'nsteps', 'rlist', 'rcoulomb', 'rvdw', 'fourierspacing', 'tcoupl', 'tc-grps', 'tau_t', 'ref_t')
        npt_key = ('tau_p', 'pcoupl', 'ref_p')
        for j in key_exist:
            if not j in self.mds[mdSection].keys():
                echoFail('[' + mdSection + ']section下无'+j+'键值')
                exit(1)
            if self.mds[mdSection]['md_type'] == 'npt':
                for k in npt_key:
                    if not k in self.mds[mdSection].keys():
                        echoFail('[' + mdSection + ']section下无'+j+'键值')
                        exit(1)

    def __init_precise(self):
        '''设置在不同精度版本下gromacs的命令
        
        @Changed:
        self.gmx_cmds: {command name: cmd}包含命令名称和命令
        '''
        if self.general['precision'] == 'single':
            self.gmx_cmds = {'mdrun':'mdrun', 'genbox':'genbox', 'grompp':'grompp', 'editconf':'editconf', 'g_energy':'g_energy', 'g_density':'g_density', 'g_tune_pme':'g_tune_pme', 'genion':'genion', 'make_ndx':'make_ndx', 'g_rms':'g_rms', 'g_rdf':'g_rdf', 'genrestr':'genrestr', 'tpbconv':'tpbconv', 'trjconv':'trjconv'} 
        elif self.general['precision'] == 'double':
            self.gmx_cmds = {'mdrun':'mdrun_d', 'genbox':'genbox_d', 'grompp':'grompp_d', 'editconf':'editconf_d', 'g_energy':'g_energy_d', 'g_density':'g_density_d', 'g_tune_pme':'g_tune_pme_d', 'genion':'genion_d', 'make_ndx':'make_ndx_d', 'g_rms':'g_rms_d', 'g_rdf':'g_rdf_d', 'genrestr':'genrestr_d', 'tpbconv':'tpbconv_d', 'trjconv':'trjconv_d'} 
        elif self.general['precision'] == 'mpi_single':
            self.gmx_cmds = {'mdrun':'mdrun_mpi', 'genbox':'genbox_mpi', 'grompp':'grompp_mpi', 'editconf':'editconf_mpi', 'g_energy':'g_energy_mpi', 'g_density':'g_density_mpi', 'g_tune_pme':'g_tune_pme_mpi','genion':'genion_mpi', 'make_ndx':'make_ndx_mpi', 'g_rms':'g_rms_mpi', 'g_rdf':'g_rdf_mpi', 'genrestr':'genrestr_mpi', 'tpbconv':'tpbconv_mpi', 'trjconv':'trjconv_mpi'} 
        elif self.general['precision'] == 'mpi_double':
            self.gmx_cmds = {'mdrun':'mdrun_mpi_d', 'genbox':'genbox_mpi_d', 'grompp':'grompp_mpi_d', 'editconf':'editconf_mpi_d', 'g_energy':'g_energy_mpi_d', 'g_density':'g_density_mpi_d', 'g_tune_pme_mpi_d':'g_tune_pme','genion':'genion_mpi_d', 'make_ndx':'make_ndx_mpi_d', 'g_rms':'g_rms_mpi_d', 'g_rdf':'g_rdf_mpi_d', 'genrestr':'genrestr_mpi_d', 'tpbconv':'tpbconv_mpi_d', 'trjconv':'trjconv_mpi_d'} 
        else:
            echoFail('precision值应为single、double、mpi_single、mpi_double中的一個.')
            exit(1)

    def get_md_sequence_list(self):
        '''获取以[md_type,...]列表形式表示的模拟类型集合，比如['npt', 'nvt']'''
        return self.md_sequence

    def get_mds_dict(self):
        '''获取以{md_type:{key:value},...}字典形式表示的模拟类型对应下的键值，例如{'mdnpt':{title:npt}}，表示[mdnpt] section下的title键的值为npt'''
        return self.mds

    def get_em_sequence_list(self):
        '''获取以[em_type,...]列表形式表示的模拟类型集合，比如['emsteep', 'emcg']'''
        return self.em_sequence

    def get_ems_dict(self):
        '''获取以{em_type:{key:value},...}字典形式表示的模拟类型对应下的键值，例如{'emsteep':{title:emsteep}}，表示[emsteep]section下的title键的值为emsteep'''
        return self.ems
    
    def get_general_dict(self):
        '''获取以{key:value,...}字典形式表示的[general]下的键值'''
        return self.general
    
    def get_input_dict(self):
        '''获取以{key:value,...}字典形式表示的[input]下的键值'''
        return self.input
    
    def get_md_dict(self):
        '''获取以{key:value,...}字典形式表示的[md]下的键值'''
        return self.md
    
    def get_residure_list(self):
        '''获取以[(key:value),...]元组列表的形式表示的残基的集合，例如[('DRG',
        20)]表示DRG残基的个数为20'''
        return self.residure
    
    def get_ions_dict(self):
        '''获取以{离子:数量,...}形式表示的离子的集合，例如{'NA+':20}表示NA+的个数为20'''
        return self.ions
    
    def get_gmx_cmds_dict(self):
        '''获取以{命令名称:命令,...}形式表示的gmx命令的集合，例如{'mdrun':'mdrun_d'}表示mdrun的命令为双精度的mdrun_d'''
        return self.gmx_cmds

    def get_gmx_precision(self):
        '''獲取gromacs的命令精度'''
        return self.general['precision']

    def get_mpi_nodes(self):
        '''獲取使用mpi時的節點'''
        return int(self.general['mpi_nodes'])
    
    def __copy_key_value(self, source, target):
        '''将source中所有的键值拷贝至target里，之后再将target原的的键值还原，也就是说从source中拷贝target中未有之键至target'''
        #备份
        backup = {}
        for i in self[target]:
            backup[i] = self[target][i]
        #赋值
        for j in self[source]:
            self[target][j] = self[source][j]
        #还原
        for k in list(backup.keys()):
            self[target][k] = backup[k]


class CommandOut:
    '''输出模拟所需要的命令,及之后结果分析需要的命令,需要使用ReadConfig初始化以读取变量'''

    def __init__(self, read_configer):
        '''从ReadConfig中读取各种配置'''
        #gmx 命令dict，里面定义了单双精度的命令
        self.gmx_cmds = read_configer.get_gmx_cmds_dict()
        # 獲取配置文件中[general]下的鍵值
        self.general = read_configer.get_general_dict()
        # gromacs的精度，一般為single, double, mpi_double, mpi_single
        self.precision = self.general['precision']
        # mpi的節點
        self.mpi_nodes = int(self.general['mpi_nodes'])
        # 是否優化mpi，這將啟用g_tune_pme來優化pme，但可能在模擬過程中出錯，請謹慎操作
        self.opt_mpi = int(self.general['opt_mpi'])
        #[input]section下的键值dict
        self.input = read_configer.get_input_dict()
        #ions列表
        self.ions = read_configer.get_ions_dict()
        #[md]section下的键值
        self.md = read_configer.get_md_dict()
        #md_sequence,包含了md模拟顺序的list,例如mdnpt,mdnvt
        self.md_sequence = read_configer.get_md_sequence_list()
        #对应md_sequence中模拟类型section下的键值，其结构为，其结构为{mdtype:{key:value,...}, ...}
        self.mds = read_configer.get_mds_dict()
        #em_sequence,包含了em模拟顺序的list,例如emsteep,emcg
        self.em_sequence = read_configer.get_em_sequence_list()
        #对应em_sequence中模拟类型section下的键值，其结构为，其结构为{emtype:{key:value,...}, ...}
        self.ems = read_configer.get_ems_dict()
        #残基名称与数量对应dict
        self.residure = read_configer.get_residure_list()
        # 用来记录与跟踪已经生成的类型的文件名, 方便之后调用
        self._md_files = {'gro':'', 'edr':'', 'tpr':'', 'trr':'','cpt':'', 'ndx':''}
        #合并所有从能量最小化开始至模拟结束的md过程，结构为{mdtype:{key:value,...}, ...}
        self.all_mds = {}
        for i in (self.ems, self.mds):
            for j, k in i.items():
                self.all_mds[j] = k

    def __write_commands(self, string, mode='a'):
        '''输出字符串到commands命令文件
        mode 表示打开文件的模式，和python一致，默认值为a，表示追加'''
        with open('commands', mode) as f:
            f.write(string)

    def __write_analysis(self, ass, string, mode='a'):
        '''输出字符串到analysis分析命令文件

        mode 表示打开文件的模式，和python一致，默认值为a，表示追加
        '''
        with open(ass, mode) as f:
            f.write(string)

    def output(self):
        '''输出命令'''
        #新建文件
        self.__write_commands('', 'w')
        #

        # 如果设置了位置限定的残基，就生成位置限定所需要的索引文件
        if len(self.md['pr_residure']) > 0:
            self.__genrestr(self.md['pr_residure'].split())

        #转换pdb
        self.__pdb_convert()

        # 生成全局ndx文件
        self.__write_commands(r'echo -e "q\n"|' + self.gmx_cmds['make_ndx'] + ' -f ' + self._md_files['gro'] + ' -o origin.ndx\n')
        self._md_files['ndx'] = 'origin.ndx'

        #如果需要的话，添加离子
        self.__genion_out()

        #生成能量最小化命令
        for i in self.em_sequence:
            # grompp
            self.__grompp_out(i)
            # mdrun
            self.__mdrun_output(i)

        # 生成模拟过程的命令
        for i in self.md_sequence:
            md_type = self.all_mds[i]['md_type']
            path = self.all_mds[i]['analysis_dir']
            # 如果这一节是需要延长上一步的模拟
            if self.all_mds[i]['md_type'] == 'ext':
                last_md_section = self.__find_not_ext(i)
                md_type = self.all_mds[last_md_section]['md_type']
                # tpbconv
                self.__tpbconv_out(i, self.all_mds[i]['time'])
                # mdrun
                self.__mdrun_output(i)
                # analysis
                self.__analysis_output(i, md_type, path)
            else:
                # grompp
                self.__grompp_out(i)
                # mdrun
                self.__mdrun_output(i)
                # analysis
                self.__analysis_output(i, md_type, path)

    def __find_not_ext(self, mdSection):
        '''找到self.md_sequence中mdSection上一个类型不为ext的模拟步骤。

        如果mdSection是一个ext类型的模拟步骤，表示它是一个延长前面模拟过程的步骤，此程序将找到上一个模拟类型不是ext的模拟步骤并返回。
        @Args:
        mdSection: self.md_sequence的一个成员，其md_type应为ext

        @Return:
        如果存在则返回的是self.md_sequence中的一个值, 否则返回空字符串
        '''
        if self.all_mds[mdSection]['md_type'] == 'ext':
            index = self.md_sequence.index(mdSection)
            for i in range(index-1, -1, -1):
                if self.mds[self.md_sequence[i]]['md_type'] != 'ext':
                    return self.md_sequence[i]
        return ''

    def __pdb_convert(self):
        '''输出pdb转换命令'''
        self.__write_commands(self.gmx_cmds['editconf'] + ' -f ' + self.input['pdb'] + ' -box ' + self.input['box'] + ' -o pdb.gro\n')
        self._md_files['gro'] = 'pdb.gro'

    def __genion_out(self):
        '''输出添加离子的命令'''
        if len(self.ions) == 0:
            return
        for i in sorted(self.ions.keys()):
            # Use emsteep mdp file
            mdp = self.ems['emsteep']['mdp']
            grompp = self.gmx_cmds['grompp'] + ' -f ' + mdp + ' -c ' + self._md_files['gro'] + ' -p ' + self.input['top'] + ' -o ions.tpr'
            genionCmd = 'echo SOL | ' + self.gmx_cmds['genion'] + ' -s ions.tpr -o start.gro -p ' + self.input['top']
            self._md_files['gro'], self._md_files['tpr'] = 'start.gro', 'ions.tpr'
            positive_ions = ('NA', 'K', 'MG', 'CA', 'ZN')
            negative_ions = ('CL')
            if i in positive_ions:
                genionCmd += ' -np ' + str(self.ions[i]) + ' -pname ' + i
            elif i in negative_ions:
                genionCmd += ' -nn ' + str(self.ions[i]) + ' -nname ' + i
            else:
                echoFail('离子名称' + i + '不正确')
                exit(1)
            self.__write_commands(grompp + '\n' + genionCmd + '\n')
        # Generate index file after adding ions
        ndx_syntax = 'r SOL|r ' + '|r '.join(sorted(self.ions.keys()))
        self.__write_commands(r'echo -e "' + ndx_syntax + r'\nq\n"|' + self.gmx_cmds['make_ndx'] + ' -f ' + self._md_files['gro'] + ' -o start.ndx\n')
        self._md_files['ndx'] = 'start.ndx'

    def __genrestr(self, residures):
        '''If need to do position restraint, there must be a porse.itp
        @Args:
        residure: The residures are to restrainted. It is a list.
        '''
        # 所要限制的残基
        pr_residures = self.md['pr_residure'].split()
        # 残基的pdb文件
        pr_pdbs = self.md['pr_pdb'].split()
        for i in range(len(pr_residures)):
            pdb, residure = pr_pdbs[i], pr_residures[i]
            self.__write_commands(r'echo -e "0\nq"|' + self.gmx_cmds['genrestr'] + ' -f ' + pdb + ' -o posre_' + residure + '.itp -fc 1000 1000 1000\n')

    def __tpbconv_out(self, mdSection, time):
        '''此步输出延长模拟上一步的命令'''
        tpr = self._md_files['tpr']
        self.__write_commands(self.gmx_cmds['tpbconv'] + ' -s ' + tpr + ' -extend ' + time +' -o ' + mdSection + '.tpr\n')
        self._md_files['tpr'] = mdSection + '.tpr'

    def __grompp_out(self, mdSection):
        '''输出grompp程序命令，mdSection表示md进行的类型，比如emsteep,emcg,npt等等，其对应于ini文件中的各section标签.'''
        #不带-maxwarn选项的grompp命令
        gromppCmd = self.gmx_cmds['grompp'] + ' -f ' + self.all_mds[mdSection]['mdp'] + ' -c ' + self._md_files['gro'] + ' -p ' + self.input['top'] + ' -o ' + mdSection + '.tpr -n ' + self._md_files['ndx']

        # If previous md type is npt or nvt, there must be a -t option to get checkpoint from it
        if self._md_files['cpt'][:3] in ('npt', 'nvt'):
            gromppCmd += ' -t ' + self._md_files['cpt']

        #检测nodes即cpu个数，及maxwarn值，如果它们为0，其选项-nt及-maxwarn都将忽略，如果为正整数，将使用此选项
        maxwarn = int(self.md['maxwarn'])
        if  maxwarn < 0:
            echoFail('ini配置文件中maxwarn值错误，期望值为0及正整数')
            exit(1)
        elif maxwarn > 0:
            # If there is maxwarn, add -maxwarn option
            gromppCmd += ' -maxwarn ' + str(maxwarn)
        self.__write_commands(gromppCmd + '\n')
        self._md_files['tpr'] = mdSection + '.tpr'

    def __mdrun_output(self, mdSection):
        '''输出md执行程序，mdSection表示md进行的类型，比如emsteep,emcg,npt, npt1等等，其对应于ini文件中的各section标签.'''
        #是否是用mpi
        is_mpi = True if self.precision[:3] == 'mpi' else False
        if is_mpi:
            self.__mdrun_mpi_output(mdSection)
        else:
            self.__mdrun_single_output(mdSection)
    
    def __mdrun_single_output(self, mdSection):
        '''
        输出md执行程序，mdSection表示md进行的类型，比如emsteep,emcg,npt, npt1等等，其对应于ini文件中的各section标签.
        這是一個單機運行的輸出，並行輸出在__mdrun_mpi_output
        '''
        mdCmd = self.gmx_cmds['mdrun'] + ' -deffnm ' + mdSection
        if self.all_mds[mdSection]['md_type'] == 'ext':
            mdCmd += ' -cpi ' + self._md_files['cpt']
        # If CPU number is set, use -nt option
        nodes = int(self.all_mds[mdSection]['nodes'])
        if nodes > 0:
            mdCmd += ' -nt ' + str(nodes)
        mdCmd += '\n'
        self.__write_commands(mdCmd) 
        # Add new files generated
        for i in ('cpt', 'gro', 'edr', 'trr'):
            self._md_files[i] = mdSection + '.' + i

    def __mdrun_mpi_output(self, mdSection):
        '''输出md执行程序，mdSection表示md进行的类型，比如emsteep,emcg,npt, npt1等等，其对应于ini文件中的各section标签.
        這個使用了g_tune_pme程序來檢測pme的合理值，之後再運行mdrun,以達到最高效率,但跳过能量最小化，和不需要優化的情況
        '''
        if self.mpi_nodes > 0:
            if self.all_mds[mdSection]['md_type'] == 'steep' or self.all_mds[mdSection]['md_type'] == 'cg' or self.opt_mpi == 0:
                mdCmd ='mpiexec -np ' + str(self.mpi_nodes) + ' ' + self.gmx_cmds['mdrun']
            else:
                mdCmd =self.gmx_cmds['g_tune_pme'] + ' -np ' + str(self.mpi_nodes) + ' -s ' + self._md_files['tpr'] + ' -launch'
        else:
            echoFail('使用並行程序必須指定一個大於0的節點數')
        mdCmd += ' -deffnm ' + mdSection
        if self.all_mds[mdSection]['md_type'] == 'ext':
            mdCmd += ' -cpi ' + self._md_files['cpt']
        mdCmd += '\n'
        self.__write_commands(mdCmd) 
        # Add new files generated
        for i in ('cpt', 'gro', 'edr', 'trr'):
            self._md_files[i] = mdSection + '.' + i


    def __analysis_output(self, mdSection, mdType, path):
        '''output the analysis commands to file
        
        @Attribute:
            mdSection: md course such as emsteep, npt, npt1
            mdType: md type such as emsteep, nvt, npt
            path: which path does analysis result output
        '''
        # If md type is not npt, nvt, return
        # 不再分析em
        if not mdType in ('npt', 'nvt'):
            return

        ass = mdSection + '_ass'
        self.__write_analysis(ass, '#!/bin/bash', 'w')


        #日志输出
        log = path + '/analysis.log'

        #对分析程序赋值
        for i in ['make_ndx', 'g_energy', 'g_rms', 'g_density', 'g_rdf', 'trjconv']:
            globals()[i] = self.gmx_cmds[i]
        
        # file name for analysis
        gro, edr, tpr, trr = '../' + self._md_files['gro'], '../' + self._md_files['edr'], '../' + self._md_files['tpr'], '../' + self._md_files['trr']
        
        # 生成目录
        self.__write_analysis(ass, '\n[[ -d "' + path + '" ]] || mkdir -pv ' + path)

        # 输出开始语句
        self.__write_analysis(ass, '\necho --------------------' + ass + '-------------------- > ' + log + '\n')

        # Generate ndx file and group some atom or ions; S_SO3, O_SO3 are user defined atom types. OPLS_288 is R4N+
        self.__write_analysis(ass, r'echo -e "t S_SO3\\nt O_SO3\\nt OPLS_288\\na HW1 | a HW2\\na OW\\nq\\n"|' + make_ndx + ' -f ' + tpr + ' -o ' + path + '/system.ndx\n')

        # This is a python script to get the group nr for a group name, it can be used in such way: python3 getnr.py group
        self.__write_analysis(ass, r'''echo -e "#!/usr/bin/env python\nimport sys\n\nstringList = sys.argv[1].split()\ntry:\n\tfound = stringList.index(sys.argv[2])\nexcept ValueError:\n\ttry:\n\t\tfound = stringList.index(sys.argv[2] + '+')\n\texcept ValueError:\n\t\ttry:\n\t\t\tfound = stringList.index(sys.argv[2] + '-')\n\t\texcept ValueError:\n\t\t\texit(1)\nresult = stringList[found-1]\nif result.isdigit():\n\tprint(result)\nelse:\n\texit(1)" > getnr.py''' + '\n')

        # Groups whose nrs are to be gotten
        groupList = ["S_SO3", "O_SO3", "OPLS_288", "HW1_HW2", "OW", "CA2", "NA", "CL", "MG2", "System"] + list(dict(self.residure).keys())
        genNdx = r'groups=$(echo -e "q\\n"|' + make_ndx + ' -n ' + path + '/system.ndx' + ' -o ' + path + '/groups.ndx|grep'
        getNr = ''
        for i in groupList:
            genNdx += ' -e "' + i + '"'
            getNr += i + '=$(python3 getnr.py "$groups" ' + i + ');'
        self.__write_analysis(ass, genNdx + ')\n')
        self.__write_analysis(ass, getNr + '\n')

        #能量分析，分别对应能量，压力，温度，张力，密度
        for energy in ['potential', 'total-energy', 'pressure', 'temperature', '#Surf*SurfTen', 'density', 'Pres-XX', 'Pres-YY', 'Pres-ZZ']:
            self.__write_analysis(ass, 'echo "' + energy + '"| ' + g_energy + ' -f ' + edr + ' -o ' + path + '/' + energy + '.xvg >> ' + log + '\n')

        # density
        self.__write_analysis(ass, 'echo 0 | ' + g_density + ' -f ' +
                trr + ' -n ' + path + '/system.ndx' + ' -s ' + tpr +
                ' -d z -dens mass -o ' + path + '/density_system.xvg >> ' +
                log + '\n')
        meta = [i[0] for i in self.residure] + sorted(self.ions.keys())
        for density in meta:
            self.__write_analysis(ass, 'echo 输出' + density + '密度:\n')
            cmd = (r'reside=$(echo -e "r ' + density + r'\\\\nq\\\\n"|' +
                    make_ndx + ' -f ' + gro + ' -o get_nr.ndx|grep ' +
                    density + ');' + r"if [[ -n $reside ]]; then gr=$(" +
                    r"echo $reside|awk '{print $1}'); else set $?=1; fi;" +
                    'echo -e "' + r'${gr}\\\\nq\\\\n" |' + g_density +
                    ' -f ' + trr + ' -n ' + path + '/system.ndx' + ' -s ' +
                    tpr + ' -d z -dens mass -o ' + path + '/density_' +
                    density + '.xvg >> ' + log + '\n')
            self.__write_analysis(ass, cmd)

        #rms
        for i in ["System"] + list(dict(self.residure).keys()):
            self.__write_analysis(ass, 'echo $' + i + ' $' + i + ' | ' + g_rms + ' -f ' + trr + ' -s ' + tpr +' -o ' + path + '/rmsd-' + i + '-vs-nptstart.xvg >> ' + log + '\n')

        # rdf
        if self.input['forcefield'] == 'oplsaa':
            # for other forcefields
            #self.__write_analysis(ass, r'echo -e "q\\n"|' + make_ndx + ' -f ' + gro + ' -o ' + path + '/rdf.ndx\n')

            for i in ['S_SO3', 'O_SO3', 'OPLS_288']:
                for j in ["HW1_HW2", "OW", "CA2", "NA", "CL", "MG2"]:
                    self.__write_analysis(ass, 'if [[ -n ${' + i + '} && -n ${' + j + '} ]]\nthen\n\techo -e "${' + i + r'}\\n${' + j + r'}\\n"|' + g_rdf  + ' -f ' + trr + ' -n ' + path + '/system.ndx' + ' -o ' + path + '/' + i + '_' + j + '_rdf.xvg\nfi\n')
            self.__write_analysis(ass, 'rm getnr.py\n')


class TopOut:
    '''生成top文件，需要使用ReadConfig初始化值'''
    def __init__(self, read_configer):
        if not isinstance(read_configer, ReadConfig):
            echoFail('参数传递错误，此处期望read_configer是一个ReadConfig类型的对象.')
            exit(1)
        #gmx 命令dict，里面定义了单双精度的命令
        self.gmx_cmds = read_configer.get_gmx_cmds_dict()
        #[input]section下的键值dict
        self.input = read_configer.get_input_dict()
        #ions列表
        self.ions = read_configer.get_ions_dict()
        #残基名称与数量对应的元组列表
        self.residure = read_configer.get_residure_list()
        # [md] items
        self.md = read_configer.get_md_dict()

    def output(self):
        with open(self.input['top'], 'w') as f:
            f.write('#include "' + self.input['forcefield'] + '.ff/forcefield.itp"\n')
            if self.input['forcefield'] == 'oplsaa':
                f.write('#include "opls_user.itp"\n')
            residure_itps = self.input['residure_itp'].split()
            residure_names = self.input['residure_name'].split()
            pr_residures = self.md['pr_residure'].split()
            pr_define = 'POSRES_' + '-'.join(pr_residures)
            for i in range(len(residure_names)):
                # 包含水模型itp时，要包含力场前缀
                if residure_itps[i] in ('spc', 'spce', 'tip3p', 'tip4p', 'tip5p'):
                    f.write('#include "' + self.input['forcefield'] + '.ff/' + residure_itps[i] + '.itp"\n')
                else:
                    f.write('#include "' + residure_itps[i] + '.itp"\n')
                if residure_names[i] in pr_residures:
                    # Position restraint
                    f.write('\n; Include Position restraint file\n#ifdef ' + pr_define + '\n#include "posre_' + residure_names[i] + '.itp"\n#endif\n\n')
            for i in self.input['itp_files'].split():
                f.write('#include "' + self.input['forcefield'] + '.ff/' + i + '.itp"\n')

            # [ system ]
            f.write('\n[ system ]\n; name\nsurfactants in water\n\n[ molecules ]\n; Compound\t#mols\n')
            for res, resnum in self.residure:
                f.write(res + '\t\t' + str(resnum) + '\n')


class MdpOut:
    '''生成模拟所需的mdp文件，需要使用ReadConfig初始化值'''
    def __init__(self, read_configer):
        if not isinstance(read_configer, ReadConfig):
            echoFail('参数传递错误，此处期望read_configer是一个ReadConfig类型的对象.')
            exit(1)
        #md_sequence,包含了md模拟顺序的list,例如mdnpt,mdnvt
        self.md_sequence = read_configer.get_md_sequence_list()
        #对应md_sequence中模拟类型section下的键值，其结构为，其结构为{mdtype:{key:value,...}, ...}
        self.mds = read_configer.get_mds_dict()
        #em_sequence,包含了em模拟顺序的list,例如emsteep,emcg
        self.em_sequence = read_configer.get_em_sequence_list()
        #对应em_sequence中模拟类型section下的键值，其结构为，其结构为{emtype:{key:value,...}, ...}
        self.ems = read_configer.get_ems_dict()
        #合并所有从能量最小化开始至模拟结束的md过程，结构为{mdtype:{key:value,...}, ...}
        self.all_mds = {}
        for i in (self.ems, self.mds):
            for j, k in i.items():
                self.all_mds[j] = k
        #以下键值不是mdp中对应的键值，需要忽略
        self.ignored_keys = ['default', 'md_type', 'nodes', 'mdp', 'analysis_dir', 'position_restraint']

    def output(self):
        for i in self.em_sequence:
            self.mdp_out(i)
        for i in self.md_sequence:
            if self.mds[i]['md_type'] != 'ext':
                self.mdp_out(i)

    def mdp_out(self, mdSection):
        for i in self.ignored_keys:
            if i in self.all_mds[mdSection]:
                self.all_mds[mdSection].pop(i)
        with open(mdSection + ".mdp", 'w') as f:
            mdp_dict=self.all_mds[mdSection]
            keys_list=list(mdp_dict.keys())
            keys_list.sort()
            for j in keys_list:
                f.write(j + "\t=\t" + mdp_dict[j] + '\n')


def main():
    arg_parser = argparse.ArgumentParser(description='通过读取ini配置文件，\
            生成GROMACS模拟所需的命令文件、top文件及mdp文件。')
    arg_parser.add_argument('-c', '--command', action='store_true',
            help='生成命令文件, 默认文件名：commands')
    arg_parser.add_argument('-t', '--top', action='store_true',
            help='生成top文件, 默认文件名：m0.top')
    arg_parser.add_argument('-m', '--mdp', action='store_true',
            help='生成mdp文件，默认文件名: *.mdp，*以模拟所运行的md步骤来表示，比如emsteep.mdp')
    arg_parser.add_argument('-i', '--input', action='store', required=True,
            help='ini配置文件')
    args = arg_parser.parse_args()
    config_reader = ReadConfig(args.input)
    if args.command:
        command= CommandOut(config_reader)
        command.output()
    if args.top:
        top= TopOut(config_reader)
        top.output()
    if args.mdp:
        mdp= MdpOut(config_reader)
        mdp.output()

def echoWarning(string):
    '''Give a colorful display for string as a warning'''
    warning = '\033[93m'
    endc = '\033[0m'
    print(warning + "Warning: " + string + endc)

def echoFail(string):
    '''Give a colorful display for string as a fail'''
    fail = '\033[91m'
    endc = '\033[0m'
    print(fail + "Error: " + string + endc)

if __name__ == '__main__':
    main()
