#!/usr/bin/env python
'''
Calculated the hbond of radials

command hbond.ndx system.ndx resname 'radials'
'''

import sys
from gromacs.fileformats import ndx


def main():
    hbond_ndx = ndx.NDX(sys.argv[1])
    system_ndx = ndx.NDX(sys.argv[2])
    resname = sys.argv[3]
    radials = sys.argv[4].split()

    non_exist_radials = list(set(radials).difference(set(system_ndx.groups)))
    if len(non_exist_radials) != 0:
        print('The following radials are not in the %s groups: %s' % (sys.argv[2], non_exist_radials))
    radials = list(set(radials).intersection(set(system_ndx.groups)))

    donors_str = 'donors_hydrogens_' + resname
    acceptors_str = 'acceptors_' + resname

    if not donors_str in hbond_ndx.groups or not acceptors_str in hbond_ndx.groups:
        print("There is no %s or %s in %s groups." % (donors_str, acceptors_str, sys.argv[1]))
        exit(1)

    donors = {}
    acceptors = {}
    hbond_donors = hbond_ndx.get(donors_str)
    hbond_acceptors = hbond_ndx.get(acceptors_str)
    for radial in radials:
        donors[radial] = []
        acceptors[radial] = []
        for index in system_ndx.get(radial):
            if index in hbond_donors:
                donors[radial].append(index)
            if index in hbond_acceptors:
                acceptors[radial].append(index)

    for radial in radials:
        print('There are %s %s radial in %s. They are:\n%s' % (len(donors[radial]), radial, donors_str, donors[radial]))
        print('There are %s %s radial in %s. They are:\n%s' % (len(acceptors[radial]), radial, acceptors_str, acceptors[radial]))


if __name__ == '__main__':
    main()
