#!/usr/bin/env python2
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Hugh Gao at http://klniu.com/

# Released under the GNU Public Licence, v2 or any higher version
#

"""@package docstring
获取原子的坐标。

通过给定残基名称，残基内原子数目，原子在残基内的索引(从0开始)，计算原子的坐标。

command resname atoms_num_of_res index topology_file
"""

import argparse
from MDAnalysis import Universe


def main():
    arg_parser = argparse.ArgumentParser(description='通过给定残基名称，残基内原子数目，原子在残基内的索引(从0开始)，计算原子的坐标。')
    arg_parser.add_argument('resname', action='store', help='残基名称')
    arg_parser.add_argument('atoms_num', type=int, action='store', help='残基内原子数目')
    arg_parser.add_argument('index', type=int, action='store', help='原子的索引，索引从0开始')
    arg_parser.add_argument('topology_file', action='store', help='拓扑文件，例如gro, pdb')
    args = arg_parser.parse_args()

    resname, atoms_num, index = args.resname, args.atoms_num, args.index

    universe = Universe(args.topology_file)
    atom_groups = universe.selectAtoms("resname " + resname)
    if len(atom_groups) % atoms_num != 0:
        print("拓扑文件内对应残基原子总数不是所给原子数目的整数倍，请给予正确的原子数目。")
        exit(1)

    positions = []
    for i in range(0, len(atom_groups), atoms_num):
        positions.append(atom_groups[i:i + atoms_num][index].position)

    print("The positions of atoms %s is:" % (index))
    for i in positions:
        print(i)

if __name__ == '__main__':
    main()
