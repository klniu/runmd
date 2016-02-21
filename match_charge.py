#!/usr/bin/env python
# -*- coding: UTF-8 -*-
##############################################################################
# @file matcharge.py
#
# @date 2012-04-20 10:28
# @author Xiang Gao
# @email email@Klniu.com
#
# @Licence GPL v2
#
# @brief
#
# @detail
#
##############################################################################
from decimal import Decimal
import argparse
import moltoolkit


def main():
    parser = argparse.ArgumentParser(description='Match the charges between two identical molecule and output the new charge file.\n\nThe format of charge file must be like:\nindex\tcharge\nindex\tcharge\n...', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-1', '--molfile1', required=True, help='The first molecule file')
    parser.add_argument('-2', '--molfile2', required=True, help='The second molecule file')
    parser.add_argument('-c', '--chargefile', required=True, help='Charges file for the first molecule')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    mol = moltoolkit.Mol(args.molfile1)
    match_idx = mol.getSub(args.molfile2)

    if len(match_idx) == 0:
        print("The two molecules is not identical.")
        exit(1)
    mol1 = {}
    # 原子计数
    index = 0
    with open(args.chargefile) as f:
        for line in f:
            charge = line.strip()
            if charge[0] == '#':
                continue
            index += 1
            try:
                mol1[index] = round(Decimal(charge), 3)
            except(Decimal.InvalidOperation):
                print("Fatal Error: The format of the charge file is incorrect.")
                exit(1)
    mol2 = {}
    for (atom1, atom2) in match_idx:
        mol2[atom2] = mol1[atom1]

    with open(args.output, 'w') as f:
        for i in sorted(mol2.keys()):
            f.write(str(mol2[i]).rjust(6) + '\n')

if __name__ == '__main__':
    main()
