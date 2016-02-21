#!/usr/bin/env python
'''指定残基名称，残基所包含的原子数量，要添加的原子在残基中的索引，然后在指定ndx文件中添加一个新的组，组内包含这些指定的残基内的原子.

add_group_to_ndx.py new_group_name res_name atoms_num_of_res index_of_atoms num_mols ndx_file

残基名称必须是ndx文件的一个group.
原子数量为整数
索引如果包含多个原子请用空格分开并用引号包围
'''

import sys
import os
from gromacs.fileformats import ndx

if len(sys.argv) != 7:
    print(__doc__)
    exit(1)

group, res, num, indices_str, num_mols, ndx_file = sys.argv[1:7]
num, num_mols = int(num), int(num_mols)
if not os.path.exists(ndx_file):
    print('ndx_file', ndx_file, 'does not exist.')
    print(__doc__)
    exit(1)

if len(group) <= 1:
    print("The length of group name must bigger than 1.")
    exit(1)

#print("\033[0;33mNote: The index of atoms should be numbered from 0.\033[0m")
indices = [int(i) for i in indices_str.split()]

ndx_exa = ndx.NDX(ndx_file)
if not res in ndx_exa:
    print('res_name', res, 'is not a group of', ndx_file)
    print(__doc__)
    exit(1)

ndx_exa.add_group_of_atoms_in_res(group, res, num, indices, num_mols)
print('Added', group, 'including', len(ndx_exa[group]) // len(indices), 'molecules or', len(ndx_exa[group]), 'atoms successfully.')
