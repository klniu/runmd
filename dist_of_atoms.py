#!/usr/bin/env python2
# -*- coding: UTF-8 -*-

# Copyright (c) 2012 Hugh Gao at http://klniu.com/

# Released under the GNU Public Licence, v2 or any higher version
#

"""@package docstring
计算原子之间的距离。

通过给定残基名称，残基内原子数目，两个原子在残基内的索引(从0开始)，计算所有残基内这两个原子之间的直线距离。

command resname atoms_num_of_res index1 index2 topology_file
"""

import argparse
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.distances import dist
from MDAnalysis.core.AtomGroup import AtomGroup


def main():
    arg_parser = argparse.ArgumentParser(description='通过给定残基名称，残基内原子数目，两个原子在残基内的索引(从0开始)，计算所有残基内这两个原子之间的直线距离。')
    arg_parser.add_argument('resname', action='store', help='残基名称')
    arg_parser.add_argument('atoms_num', type=int, action='store', help='残基内原子数目')
    arg_parser.add_argument('index1', type=int, action='store', help='第一个原子的索引，索引从0开始')
    arg_parser.add_argument('index2', type=int, action='store', help='第二个原子的索引，索引从0开始')
    arg_parser.add_argument('topology_file', action='store', help='拓扑文件，例如gro, pdb')
    args = arg_parser.parse_args()

    resname, atoms_num, index1, index2 = args.resname, args.atoms_num, args.index1, args.index2

    universe = Universe(args.topology_file)
    atom_groups = universe.selectAtoms("resname " + resname)
    if len(atom_groups) % atoms_num != 0:
        print("拓扑文件内对应残基原子总数不是所给原子数目的整数倍，请给予正确的原子数目。")
        exit(1)

    atoms1 = []
    atoms2 = []
    for i in range(0, len(atom_groups), atoms_num):
        atoms1.append(atom_groups[i:i + atoms_num][index1])
        atoms2.append(atom_groups[i:i + atoms_num][index2])

    dists = dist(AtomGroup(atoms1), AtomGroup(atoms2))
    print("The distance between atoms %s and %s is:" % (index1, index2))
    for i in dists[2]:
        print(i)
    print("The average distance between atoms %s and %s is:" % (index1, index2))
    print(np.average(dists[2]))

if __name__ == '__main__':
    main()
