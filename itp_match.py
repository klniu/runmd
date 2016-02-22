#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Copyright (C) Hugh Gao
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
##
# @file fragitpjoin.py
# @brief According the indics match of a molecule and its fragment, convert the itp file of fragment to molecule to use.
# @author Hugh Gao, hugh8505@gmail.com
# @version 0.1alpa
# @date 2012-04-21
import os
import moltoolkit


def main():
    import argparse
    parser = argparse.ArgumentParser(description='According the indics match of a molecule and its fragment, convert the itp file of fragment to molecule to use.')
    parser.add_argument('-m', '--molfile', required=True, help='Molecule file')
    parser.add_argument('-f', '--fragmol', required=True, help='The molecule file containing fragment.')
    parser.add_argument('-t', '--top', required=True, help='The topology file of the fragment')
    parser.add_argument('-p', '--pairs', help='Optional. The match pairs of atoms indices in the whole atoms and its fragment, optional. Please give a string whose format is python list. e.g. [(1, 2), (4, 5)]. If you do not assign, the program will calculate it. ')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    itpname = os.path.basename(args.top)
    mol = moltoolkit.Mol(args.molfile)
    if args.pairs:
        matches = eval(args.pairs)
    else:
        matches = mol.get_matches_of_mols(moltoolkit.Mol(args.fragmol))

    if len(matches) == 0:
        print("Error: The small molecule is not a fragment of the molecule.")
        exit(1)
    matches_dict = {}
    for (i, j) in matches:
        matches_dict[j] = i

    # handle topology file
    with open(args.top) as f:
        content = f.readlines()
    for idx in range(len(content)):
        content[idx] = content[idx].strip()
    if not content.index('[ bonds ]') < content.index('[ pairs ]') < content.index('[ angles ]') < content.index('[ dihedrals ]') < content.index('[ exclusions ]'):
        print("The order of sections in topology file is not correct. It should be bond, pairs, angles, dihedrals, exclusions.")
        exit(1)
    # get bonds
    bonds = content[content.index('[ bonds ]'):content.index('[ pairs ]')]
    # get pairs
    pairs = content[content.index('[ pairs ]'):content.index('[ angles ]')]
    # get angles
    angles = content[content.index('[ angles ]'):content.index('[ dihedrals ]')]
    # get dihedrals
    first_dihe = content.index('[ dihedrals ]')
    #second_dihe = content.index(r'[ dihedrals ]', first_dihe + len(r'[ dihedrals ]'), -1)
    second_dihe = content.index(r'[ dihedrals ]', first_dihe + 1, -1)
    dihedrals = content[first_dihe:second_dihe]
    # get second dihedrals
    dihedrals1 = content[second_dihe:content.index('[ exclusions ]')]
    # get exclusions
    exclusions = content[content.index('[ exclusions ]'):]

    # handle bonds
    bond_result = []
    for idx, bond in enumerate(bonds):
        bond_split = bond.split()
        if len(bond_split) == 0 or not bond_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(bond_split[0]), '%%')
        atom1 = matches_dict.get(int(bond_split[1]), '%%')
        if atom0 == '%%' or atom1 == '%%' or atom0 < atom1:
            part1 = str(atom0).rjust(5)
            part2 = str(atom1).rjust(5)
            type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
            type2 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
        else:
            part1 = str(atom1).rjust(5)
            part2 = str(atom0).rjust(5)
            type1 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
            type2 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        # find the second atom from the third column
        part3 = bond[bond.index(bond_split[1], 3) + len(bond_split[1]):]
        bond_result.append(part1 + part2 + part3 + "\t\t;" + type1.ljust(3) + ' ' + type2.ljust(3) + ' ' + itpname + '  ' + bond)
    bond_result.sort()

    # handle pairs
    pair_result = []
    for idx, pair in enumerate(pairs):
        pair_split = pair.split()
        if len(pair_split) == 0 or not pair_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(pair_split[0]), '%%')
        atom1 = matches_dict.get(int(pair_split[1]), '%%')
        if atom0 == '%%' or atom1 == '%%' or atom0 < atom1:
            part1 = str(atom0).rjust(5)
            part2 = str(atom1).rjust(5)
            type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
            type2 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
        else:
            part1 = str(atom1).rjust(5)
            part2 = str(atom0).rjust(5)
            type1 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
            type2 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        # find the second atom from the third column
        part3 = pair[pair.index(pair_split[1], 3) + len(pair_split[1]):]
        pair_result.append(part1 + part2 + part3 + "\t\t;" + type1.ljust(3) + ' ' + type2.ljust(3) + ' ' + itpname + '  ' + pair)
    pair_result.sort()

    # handle angles
    angle_result = []
    for idx, angle in enumerate(angles):
        angle_split = angle.split()
        if len(angle_split) == 0 or not angle_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(angle_split[0]), '%%')
        atom1 = matches_dict.get(int(angle_split[1]), '%%')
        atom2 = matches_dict.get(int(angle_split[2]), '%%')
        if atom0 == '%%' or atom2 == '%%' or atom0 < atom2:
            part1 = str(atom0).rjust(5)
            part3 = str(atom2).rjust(5)
            type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
            type3 = mol.get_atom(atom2).type if isinstance(atom2, int) else "%%"
        else:
            part1 = str(atom2).rjust(5)
            part3 = str(atom0).rjust(5)
            type1 = mol.get_atom(atom2).type if isinstance(atom2, int) else "%%"
            type3 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        part2 = str(atom1).rjust(5)
        type2 = mol.get_atom(atom1)
        # find the second atom from the eight column
        part4 = angle[angle.index(angle_split[2], 8) + len(angle_split[2]):]
        angle_result.append(part1 + part2 + part3 + part4 + "\t\t;" + type1.ljust(3) + ' ' + type2.ljust(3) + ' ' + type3.ljust(3) + ' ' + itpname + '  ' + angle)
    angle_result.sort()

    # handle dihedrals
    dihedral_result = []
    for idx, dihedral in enumerate(dihedrals):
        dihedral_split = dihedral.split()
        if len(dihedral_split) == 0 or not dihedral_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(dihedral_split[0]), '%%')
        atom1 = matches_dict.get(int(dihedral_split[1]), '%%')
        atom2 = matches_dict.get(int(dihedral_split[2]), '%%')
        atom3 = matches_dict.get(int(dihedral_split[3]), '%%')
        part1 = str(atom0).rjust(5)
        part2 = str(atom1).rjust(5)
        part3 = str(atom2).rjust(5)
        part4 = str(atom3).rjust(5)
        type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        type2 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
        type3 = mol.get_atom(atom2).type if isinstance(atom2, int) else "%%"
        type4 = mol.get_atom(atom3).type if isinstance(atom3, int) else "%%"
        # find the second atom from the thirteen column
        part5 = dihedral[dihedral.index(dihedral_split[3], 13) + len(dihedral_split[3]):]
        dihedral_result.append(part1 + part2 + part3 + part4 + part5 + "\t\t;" + type1.ljust(3) + ' ' + type2.ljust(3) + ' ' + type3.ljust(3) + ' ' + type4.ljust(3) + ' ' + itpname + '  ' + dihedral)
    dihedral_result.sort()

    # handle dihedral1s
    dihedral1_result = []
    for idx, dihedral1 in enumerate(dihedrals1):
        dihedral1_split = dihedral1.split()
        if len(dihedral1_split) == 0 or not dihedral1_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(dihedral1_split[0]), '%%')
        atom1 = matches_dict.get(int(dihedral1_split[1]), '%%')
        atom2 = matches_dict.get(int(dihedral1_split[2]), '%%')
        atom3 = matches_dict.get(int(dihedral1_split[3]), '%%')
        part1 = str(atom0).rjust(5)
        part2 = str(atom1).rjust(5)
        part3 = str(atom2).rjust(5)
        part4 = str(atom3).rjust(5)
        type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        type2 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
        type3 = mol.get_atom(atom2).type if isinstance(atom2, int) else "%%"
        type4 = mol.get_atom(atom3).type if isinstance(atom3, int) else "%%"
        # find the second atom from the thirteen column
        part5 = dihedral1[dihedral1.index(dihedral1_split[3], 13) + len(dihedral1_split[3]):]
        dihedral1_result.append(part1 + part2 + part3 + part4 + part5 + "\t\t;" + type1.ljust(3) + ' ' + type2.ljust(3) + ' ' + type3.ljust(3) + ' ' + type4.ljust(3) + ' ' + itpname + '  ' + dihedral1)
    dihedral1_result.sort()

    # handle exclusions
    exclusion_result = []
    for idx, exclusion in enumerate(exclusions):
        exclusion_split = exclusion.split()
        if len(exclusion_split) == 0 or not exclusion_split[0].isdecimal():
            continue
        atom0 = matches_dict.get(int(exclusion_split[0]), '%%')
        atom1 = matches_dict.get(int(exclusion_split[1]), '%%')
        if atom0 == '%%' or atom1 == '%%' or atom0 < atom1:
            part1 = str(atom0).rjust(5)
            part2 = str(atom1).rjust(5)
            type1 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
            type2 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
        else:
            part1 = str(atom1).rjust(5)
            part2 = str(atom0).rjust(5)
            type1 = mol.get_atom(atom1).type if isinstance(atom1, int) else "%%"
            type2 = mol.get_atom(atom0).type if isinstance(atom0, int) else "%%"
        # find the second atom from the third column
        part3 = exclusion[exclusion.index(exclusion_split[1], 3) + len(exclusion_split[1]):]
        exclusion_result.append(part1 + part2 + part3 + "\t\t;" + type1.ljust(3) + '\t' + type2.ljust(3) + '\t' + itpname + '  ' + exclusion)
    exclusion_result.sort()

    # output
    with open(args.output, 'w') as fw:
        # bonds
        fw.write('[ bonds ]\n;  ai   aj  funct   c0         c1\n')
        for i in bond_result:
            fw.write(i + '\n')
        # pairs
        fw.write('[ pairs ]\n;  ai   aj  funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp\n')
        for i in pair_result:
            fw.write(i + '\n')
        # angles
        fw.write('[ angles ]\n;  ai   aj   ak  funct   angle     fc\n')
        for i in angle_result:
            fw.write(i + '\n')
        # dihedral
        fw.write('[ dihedrals ]\n; GROMOS improper dihedrals\n;  ai   aj   ak   al  funct   angle     fc\n')
        for i in dihedral_result:
            fw.write(i + '\n')
        # dihedral1
        fw.write('[ dihedrals ]\n;  ai   aj   ak   al  funct    ph0      cp     mult\n')
        for i in dihedral1_result:
            fw.write(i + '\n')
        # exclusion
        fw.write('[ exclusions ]\n;  ai   aj  funct  ;  GROMOS 1-4 exclusions\n')
        for i in exclusion_result:
            fw.write(i + '\n')

if __name__ == '__main__':
    main()
