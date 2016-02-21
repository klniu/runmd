#!/usr/bin/env python
class Four2RB:
    '''Convert fourier dihedrals to Ryckaert-Bellemans parameters
    
    cite: gromacs manual 4.5.4 p76
    '''

    def __init__(self, v1, v2, v3, v4):
        self.v1 = float(v1)
        self.v2 = float(v2)
        self.v3 = float(v3)
        self.v4 = float(v4)

    def getC0(self):
        '''C0 = V2 + (V1 + V3) / 2'''
        return self.v2 + (self.v1 + self.v3) / 2

    def getC1(self):
        '''C1 = (-V1 + 3 * V3) / 2'''
        return (3 * self.v3 - self.v1) / 2

    def getC2(self):
        '''C2 = -V2 + 4V4'''
        return 4 * self.v4 - self.v2

    def getC3(self):
        '''C3 = -2 * V3'''
        return -2 * self.v3

    def getC4(self):
        '''C4 = -4 * V4'''
        return -4 * self.v4

    def getC5(self):
        '''C5 = 0'''
        return 0

def main():
    import sys
    four2RB = Four2RB(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print('C0:' + str(four2RB.getC0()))
    print('C1:' + str(four2RB.getC1()))
    print('C2:' + str(four2RB.getC2()))
    print('C3:' + str(four2RB.getC3()))
    print('C4:' + str(four2RB.getC4()))
    print('C5:' + str(four2RB.getC5()))

if __name__ == '__main__':
    main()
