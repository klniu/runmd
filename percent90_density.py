#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from scipy.interpolate import interp1d
import numpy as np
from gromacs.fileformats import xvg
import sys
import os.path


def read_file(filename):
    return xvg.Xvg(filename).data


def cal_90_percent_value(y_list):
    max_value = max(y_list)
    max_values = [i for i in y_list if max_value - i < 4]
    return sum(max_values) / len(max_values) * 0.9


def interpolate(x_list, y_list):
    x = np.array(x_list)
    y = np.array(y_list)
    return interp1d(x, y)


def main():
    if len(sys.argv) != 2:
        print("Usage: command data_file")
        exit(1)
    if not os.path.exists(sys.argv[1]):
        print("Error: The data files does not exist.")
        exit(1)

    x_list, y_list = read_file(sys.argv[1])
    f = interpolate(x_list, y_list)
    percent = cal_90_percent_value(y_list)
    print("Note: The 90 percent value is %s" % percent)

    print("The x and y you desired maybe as following:\n")
    x_all = np.arange(x_list[0], x_list[-1], 0.0001)
    y_all = f(x_all)
    zeros = np.zeros(y_all.shape)
    zeros.fill(percent)
    diff_all = np.abs(np.subtract(y_all, zeros))
    proper = lambda y: abs(y[1] - percent) < 0.1
    results = list(filter(proper, np.column_stack((x_all, y_all, diff_all))))
    for i, j, diff in results:
        print(i, j)

    # Calculate the best result
    # split into groups
    groups = []
    if len(results) > 0:
        groups.append([])
        groups[-1].append(results[0])
    for i, (x, y, diff) in enumerate(results[1:], 1):
        last_x, last_y = results[i - 1][:2]
        if abs(x - last_x) < 0.1 and abs(y - last_y) < 0.1:
            groups[-1].append([x, y, diff])
        else:
            groups.append([])
            groups[-1].append([x, y, diff])
    print("\nThe best result is:\n")
    for i in groups:
        best = min(i, key=lambda j: j[2])
        print(round(best[0], 4), best[1])


if __name__ == "__main__":
    main()
