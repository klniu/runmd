#!/usr/bin/env python
'''
在制作make_ndx索引时，有时需要将一组相同位置的原子索引出来，比如得到3个分子中nr为1的原子的所有索引，此分子中有72个原子，需要在make_ndx中输入a 1 | a 73 | a 145,此程序的作用即是如此。
用户提供三个值即可，依次为，分子数，分子内原子个数，原子在分子中的索引
'''

import sys

molsNum = int(sys.argv[1])
atomsNum = int(sys.argv[2])
atomIndex = eval(sys.argv[3])
for idx in atomIndex:
    result = ''
    for i in range(molsNum):
        result += '|a ' + str(i * atomsNum + idx)
    print(result[1:])
    print('')
