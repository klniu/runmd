#!/usr/bin/env python
# -*- coding: UTF-8 -*-

class Xvg:
    '''Parser the xvg file.'''

    def __init__(self, filename):
        self.filename = filename
        self._read_file()

    def _read_file(self):
        '''Parser the xvg file.'''
        self.x = []
        self.y = []
        with open(self.filename) as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('@') or line == '':
                    continue
                x, y = line.split()[:2]
                self.x.append(float(x))
                self.y.append(float(y))
        if len(self.x) != len(self.y):
            raise Exception("The number of datas at x and y axes is not equal.")

    def alter_x(self, func):
        '''Alter the x values according the function func.

        @para
            func: a function
        '''
        self.x = list(map(func, self.x))
        
    def alter_y(self, func):
        '''Alter the y values according the function func.

        @para
            func: a function
        '''
        self.y = list(map(func, self.y))

    @property
    def average_y(self):
        '''Get the average value of the in y axis.'''
        return sum(self.y) / len(self.y)

    @property
    def data(self):
        '''Get the x and y datas.'''
        return [self.x, self.y]

