#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Copyright (C) Hugh Gao
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
##
# @file itpparser.py
# @brief Parser topology file of gromacs.
# @author Hugh Gao, hugh8505@gmail.com
# @version 0.1
# @date 2012-04-26
'''
Gromacs topology file parser.

The Gromacs topology file like this:

; multiline comments
; multiline comments
; The line begins with a delimiter ; is a comment line.
[ moleculetype ]
; Name   nrexcl
DRG      3
[ atoms ]
;  nr  type  resnr  resid  atom  cgnr  charge    mass    total_charge
    1   CH3    1     DRG      C    1   -0.136  15.0350      ;  0.000
    2   CH0    1     DRG      C    2    0.349  12.0110
; total charge of the molecule:  -0.000
[ bonds ]
;  ai   aj  funct   c0         c1
    1    2    2   0.1530   7.1500e+06		;C3  C3  1.itp  7    9    2   0.1530   7.1500e+06
    2    3    2   0.1530   7.1500e+06		;C3  Car 1_nona_1o.itp  16   17    2   0.1530   7.1500e+06
[ pairs ]
;  ai   aj  funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp
    1    5    1
    1    9    1
[ angles ]
;  ai   aj   ak  funct   angle     fc
    1    2    3    2    109.50   285.00		;C3  C3  Car 1.itp  8    7    9    2    109.50   285.00
    1    2    4    2    109.50   285.00		;C3  C3  C3  1.itp  1    7    9    2    109.50   285.00
[ dihedrals ]
; GROMOS improper dihedrals
;  ai   aj   ak   al  funct   angle     fc
    3   17   15   13    2      0.00   167.36		;2.itp  4   20   16    9    2      0.00   167.36
   17   15   13   11    2      0.00   167.36		;2.itp  20   16    9   18    2      0.00   167.36
[ dihedrals ]
;  ai   aj   ak   al  funct    ph0      cp     mult
    4    2    3    9    1    180.00     1.00    6		;2.itp  1    3    4    7    1    180.00     1.00    6
    1    2    4    5    1      0.00     3.77    3		;1.itp  3    1    7   10    1      0.00     3.77    3
[ exclusions ]
;  ai   aj  funct  ;  GROMOS 1-4 exclusions
    2   10		;1.itp  7   12
    2   11		;1.itp  7   15
'''
import re
from collections import OrderedDict


class ItpParser:
    '''Parser the topology file of gromacs.'''
    def __init__(self, fp):
        '''Initialize.

        @param itp itp file-like object.
        '''
        self._file = fp
        # The delimiter of the comment from beginning of the line or inline.
        self._delimiter = ';'

        # The lines in the sections, exclude comment lines
        self._sections_line = OrderedDict()
        # The section header comment, e.g.  ;  nr  type  resnr  resid  atom  cgnr  charge    mass    total_charge
        self._sections_header = OrderedDict()
        # The items in every lines such as 2, 3, 2, 3.25e+5
        self._sections_items = OrderedDict()
        # The inline comment in every line in section
        self._sections_comments = OrderedDict()
        # The comments in the header of the file
        self._header_comment = []


        self._read()
        self._splitItems()

    def _read(self):
        '''Parser the topology file.
        '''
        sec_re = r'^\[ (?P<header>[^ ]+) \]'
        sec_re_obj = re.compile(sec_re)
        current_sec = None
        header_comment = False
        for lineno, line in enumerate(self._file, start=1):
            # Handle header comments.
            if (lineno == 1 or header_comment) and line.rstrip().startswith(self._delimiter):
                header_comment = True
                self._header_comment.append(line.rstrip())
                continue
            else:
                header_comment = False

            # Handle [ moleculetype ], [ bond ], [ angles ], [ dihedrals ], [ dihedrals ], [ exclusions ]
            match = sec_re_obj.match(line.rstrip())
            if match != None:
                current_sec = match.group('header')
                self._sections_line[current_sec] = []
                sec_lineno = lineno
                continue
            if current_sec == None:
                print("The first line no comments should be a secion like [ header ].")
                exit(1)
            else:
                if line.rstrip().startswith(self._delimiter):
                    # Record the sections header
                    if lineno == sec_lineno + 1:
                        self._sections_header[current_sec] = line.rstrip()
                else:
                    self._sections_line[current_sec].append(line.rstrip())
        print(self._sections_line)
        print(self._sections_header)

    def _splitItems(self):
        '''Split the lines in every section into items.'''
        for sec in self._sections_line.keys():
            self._sections_comments[sec] = []
            self._sections_items[sec] = []
            for line in self._sections_line[sec]:
                comment_idx = line.find(self._delimiter) 
                if comment_idx != -1:
                    self._sections_comments[sec].append(line[comment_idx:])
                else:
                    self._sections_comments[sec].append('')
                content = line[:comment_idx]
                self._sections_items[sec].append(content.split())
        print(self._sections_items)
        print(self._sections_comments)



def __main__():
    import argparse
    parser = argparse.ArgumentParser(description='According the indics match of a molecule and its fragment, convert the itp file of fragment to molecule to use.')
    parser.add_argument('-t', '--top', required=True, help='The topology file of the fragment')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()
    if args.top:
        f = open(args.top)
        itpParser = ItpParser(f)

if __name__ == '__main__':
    __main__()
