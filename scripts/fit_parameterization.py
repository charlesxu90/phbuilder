#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def loadxvg(fname: str, col: list = [0, 1], dt: int = 1, b: int = 0):
    """Loads an .xvg file into a list of lists.
    May also be used to load float columns from files in general.

    Args:
        fname (str): file name.
        col (list, optional): Columns to load. Defaults to [0, 1].
        dt (int, optional): Step size. Defaults to 1.
        b (int, optional): Starting point. Defaults to 0.

    Returns:
        list of lists : contains the columns that were loaded.
    """

    count = -1
    data = [[] for _ in range(len(col))]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        # THIS IS FOR THE dt PART.
        count += 1
        if count % dt != 0:
            continue

        listLine = stringLine.split()
        # AND THIS IS FOR THE b PART.
        if b != 0 and float(listLine[col[0]]) < b:
            continue

        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data

if __name__ == "__main__":

    order = 5  # Order of th fit (= #ncoeffs - 1).

    dVdlInitList = []  # To hold the lambda-coordinate points.
    dVdlMeanList = []  # To hold the mean dV/dl values.

    # Loop through the directories and load the data:
    for val in np.arange(-0.10, 1.11, 0.05):

        path = f"r_{val:.2f}_{1 - val:.2f}/cphmd-dvdl-1.xvg"
        print(path)  # debug.

        dVdlInitList.append(val)
        dVdlMeanList.append(np.mean(loadxvg(path)[1][1000:]))  # drop first ns.

    # Perform the polynomial fit.
    coeffs = np.polyfit(dVdlInitList, dVdlMeanList, deg=order)

    # Print result.
    format = 'dvdl_1 = '
    for coeff in coeffs:
        format += f"{coeff:.3f} "
    print(format)

    # Compare the data point and our fit in a plot.
    plt.scatter(dVdlInitList, dVdlMeanList, label='Data', color='r')

    fit = []
    for i in dVdlInitList:
        value = 0
        for j in range(0, order + 1):
            value += coeffs[::-1][j] * i**j
        fit.append(value)
    plt.plot(dVdlInitList, fit, label=f"{order}-order fit")

    plt.xlabel(r'$\lambda$-coordinate')
    plt.ylabel(r'$dV/d\lambda$')
    plt.legend()
    plt.savefig('fit.png')
