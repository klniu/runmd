#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''删除所有除了指定的索引组之外的索引组

add_group_to_ndx.py "except_group1 except_group2" ndx_file save_file

'''

import sys
import os
from gromacs.fileformats import ndx

if len(sys.argv) != 4:
    print(__doc__)
    exit(1)

groups, ndx_file, save_file = sys.argv[1:7]
groups = groups.split()
if not os.path.exists(ndx_file):
    print('ndx_file', ndx_file, 'does not exist.')
    print(__doc__)
    exit(1)

if len(groups) == 0:
    print("The length of groups list must bigger than 1.")
    exit(1)

ndx_exa = ndx.NDX(ndx_file)
for i in groups:
    if not i in ndx_exa:
        print('res_name', i, 'is not a group of', ndx_file)
        print(__doc__)
        exit(1)

ndx_exa.del_all_groups_except_given(groups, save_file)
print("Delete successfully.")
