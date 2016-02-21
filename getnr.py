#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
Get the index of residure in Gromacs ndx file.

Usuage:
    getnr.py residure ndx_file
'''
import sys
from gromacs.fileformats import ndx

ndx_file = sys.argv[1]
resi = sys.argv[2:]

index = ndx.NDX()
index.read(ndx_file)
idces = []
for i in resi:
    try:
        idces.append(str(index.groups.index(i)))
    except ValueError:
        idces.append(str(-1))
print(' '.join(idces), end='')
