#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

class GroupType:
    def __init__(self):
        self.name = ''
        self.incl = []
        self.nstates = -1
        self.atoms = []
        self.charges = {}
        self.pka = {}
        self.dvdls = {}

    def setName(self, name):
        self.name = name

    def setIncl(self, incl):
        self.incl = incl

    def setNstates(self, nstates):
        self.nstates = nstates

    def setAtoms(self, atoms):
        self.atoms = atoms

    def addCharges(self, charges, state):
        self.charges[state] = charges

    def addpKa(self, pka, state):
        self.pka[state] = pka

    def addDvdl(self, dvdl, state):
        self.dvdls[state] = dvdl

class LambdaGroup:
    def __init__(self, group_name, mdpfile, lambdafile, output):
        self.mdp = mdpfile
        self.lambdaF = lambdafile
        self.outputF = output
        self.group_name = group_name
        self.group = GroupType()
        self.getGroupData()

    def getGroupData(self):
        self.readmdp()
        self.readlambda()

    def getHeaderEnd(self, line, start):
        lline = line.strip().split("=")[0].strip().lower()
        return lline[len(start)+1:]

    def getValue(self, line, vtype):
        rline = line.strip().split("=")[1].strip()
        if vtype == "str":
            return rline
        if vtype == "int":
            return int(rline)
        if vtype == "float":
            return float(rline)
        if vtype == "flist":
            return [float(x) for x in rline.split()]
        if vtype == "slist":
            return [x for x in rline.split()]

    def readmdp(self):
        with open(self.mdp, "r") as f:
            data = f.readlines()

        for line in data:
            if line.strip().startswith("lambda-dynamics-group-type"):
                if (
                    line.split("=")[-1].strip().lower() ==
                        self.group_name.lower()
                ):
                    group_ndx = line.strip()[
                        len("lambda-dynamics-group-type")
                    ]

        for line in data:
            if line.strip().startswith(
                f"lambda-dynamics-group-type{group_ndx}"
            ):
                headerEnd = self.getHeaderEnd(
                    line,
                    f"lambda-dynamics-group-type{group_ndx}"
                )
                if headerEnd == "name":
                    self.group.setName(
                        self.getValue(line, "str")
                    )
                if headerEnd == "n-states":
                    self.group.setNstates(
                        self.getValue(line, "int")
                    )
                if headerEnd.startswith("state"):
                    stateN = int(headerEnd.split("-")[1])
                    stateHeaderEnd = self.getHeaderEnd(
                        headerEnd,
                        f"state-{stateN}"
                    )
                    if stateHeaderEnd == "charges":
                        self.group.addCharges(
                            self.getValue(line, "flist"),
                            stateN
                        )
                    if stateHeaderEnd == "reference-pka":
                        self.group.addpKa(
                            self.getValue(line, "float"),
                            stateN
                        )
                    if stateHeaderEnd == "dvdl-coefficients":
                        self.group.addDvdl(
                            self.getValue(line, "flist"),
                            stateN
                        )
    
    def readlambda(self):
        with open(self.lambdaF, "r") as f:
            inBlock = False
            for line in f:
                if line.strip().startswith("["):
                    if inBlock:
                        inBlock = False
                    elif line.strip()[1:-1].strip() == self.group_name:
                        inBlock = True
                    else:
                        continue
                elif inBlock:
                    header = line.strip().split("=")[0].strip().lower()
                    value = line.strip().split("=")[0].strip().lower()
                    if header == "incl":
                        self.group.setIncl(
                            self.getValue(line, "slist")
                        )
                    if header == "atoms":
                        self.group.setAtoms(
                            self.getValue(line, "slist")
                        )
    
    def generateMDPline(self, key, value):
        keyStr = f"lambda-dynamics-group-type1-{key}"
        if not isinstance(value, list):
            valStr = value
        elif isinstance(value[0], str):
            valStr = " ".join(str(x) for x in value)
        elif isinstance(value[0], int):
            valStr = " ".join(str(x) for x in value)
        else:
            valStr = np.array2string(
                np.array(value),
                precision = 3,
            )[1:-1]
        return f"{keyStr:<55} = {valStr}\n"
    
    def write_mdp(self, file):
        file.write(self.generateMDPline("name", self.group.name))
        file.write(self.generateMDPline("n-states", self.group.nstates))
        file.write(
            self.generateMDPline("state-0-charges", self.group.charges[0])
        )
        for i in range(self.group.nstates):
            file.write(
                self.generateMDPline(
                    f"state-{i+1}-charges",
                    self.group.charges[i+1]
                )
            )
            file.write(
                self.generateMDPline(
                    f"state-{i+1}-reference-pka",
                    self.group.pka[i+1]
                )
            )
            file.write(
                self.generateMDPline(
                    f"state-{i+1}-dvdl-coefficients",
                    self.group.dvdls[i+1]
                )
            )

    def generateLambdaGroupTypesline(self, key, value):
        keyStr = f"{key}"
        if not isinstance(value, list):
            valStr = value
        elif isinstance(value[0], str):
            valStr = " ".join(str(x) for x in value)
        elif isinstance(value[0], int):
            valStr = " ".join(str(x) for x in value)
        else:
            valStr = np.array2string(
                np.array(value),
                precision = 3,
            )[1:-1]
        return f"{keyStr:<10} = {valStr}\n"

    def write_lambda(self, file):
        file.write(f"[ {self.group.name} ]\n")
        file.write(
            self.generateLambdaGroupTypesline("incl", self.group.incl)
        )
        file.write(
            self.generateLambdaGroupTypesline("atoms", self.group.atoms)
        )
        file.write(
            self.generateLambdaGroupTypesline("qqA", self.group.charges[0])
        )
        for i in range(self.group.nstates):
            file.write(
                self.generateLambdaGroupTypesline(
                    f"pKa_{i+1}",
                    self.group.pka[i+1]
                )
            )
            file.write(
                self.generateLambdaGroupTypesline(
                    f"qqB_{i+1}",
                    self.group.charges[i+1]
                )
            )
            file.write(
                self.generateLambdaGroupTypesline(
                    f"dvdl_{i+1}",
                    self.group.dvdls[i+1]
                )
            )
            
    def write_output(self):
        with open(self.outputF, "w") as f:
            f.write("# Input for .mdp file\n\n")
            self.write_mdp(f)
            f.write("\n# Input for lambdagrouptypes.dat file\n\n")
            self.write_lambda(f)

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

    data = np.loadtxt(
        fname,
        comments = ["@", "#"],
        usecols = col
    )[b::dt]
    return data

def fit_param(order, prefix):
    dVdlInitList = []  # To hold the lambda-coordinate points.
    dVdlMeanList = []  # To hold the mean dV/dl values.

    # Loop through the directories and load the data:
    for val in np.arange(-0.10, 1.11, 0.1):

        path = f"{prefix}_{val:.2f}_{1 - val:.2f}/cphmd-dvdl-1-2.xvg"

        dVdlInitList.append(val)
        dVdlMeanList.append(np.mean(loadxvg(path)[10:,1]))  # drop first 10 frames

    # Perform the polynomial fit.
    coeffs = np.polyfit(dVdlInitList, dVdlMeanList, deg=order)
 
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

    return coeffs

def derivative_coeficients(coeffs):
    n = len(coeffs) - 1
    powers = np.arange(n, -1, -1)
    return np.multiply(coeffs, powers)[:-1]

def run_reweighing(order, prefix, nreplicas):
    nbins = 35
    bin_edges = np.linspace(-0.05, 1.05, nbins + 1)
    centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    centers_physical = centers[(centers>0.05) & (centers <1-0.05)]
    hist_total = np.zeros_like(centers_physical)

    for i in range(nreplicas):
        path = f"{prefix}_{i+1}/cphmd-coord-1-2.xvg"
        hist, bins = np.histogram(
            loadxvg(path)[:,1],
            range = (-0.05, 1.05),
            bins = 35,
            density = True,
        )
        centers = 0.5 * (bins[1:] + bins[:-1])
        hist /= hist.sum()

        hist_physical = hist[(centers>0.05) & (centers <1-0.05)]
        hist_total += hist_physical

    hist_total /= nreplicas
        
    R = 8.314e-3  # Molar gas constant
    T = 300  # Temperature

    u = - R * T * np.log(hist_total)

    coeffs = np.polyfit(centers_physical, u, deg = order + 1)

    return derivative_coeficients(coeffs)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='This script will find optimal parameters'
            + ' for parameterised group. It can work in two'
            + ' modes: (i) parameterisation, and (ii)'
            + ' reweighing. In the parameterisation mode'
            + ' the script will fit optimal polinomial to'
            + 'dvdl computed for a set of fixed lambda-s.'
            + ' In the reweighing mode the script will analyse'
            + ' the distributions of lambda-coordinate and'
            + ' adjust the coefficients in order to get flat'
            + ' distributions. As an output, the script will'
            + ' provide entries for lambdagrouptypes.dat and'
            + ' .mdp files.')

    parser.add_argument('-f', '--mdp', required=True, type=str,
                        help='Input .mdp file')
    parser.add_argument('-i', '--prefix', required=False, type=str,
                        default = None,
                        help='Prefix of input folders. ' 
                        + 'Default r for parameterisation, s for reweighing')
    parser.add_argument('-m', '--mode', required=True, type=str,
                        choices=['p', 's'],
                        help='Input topology file')
    parser.add_argument('-nr', '--nreplicas', required = False, type = int,
                        default = 10,
                        help = 'Number of replicas run for reweighing.'
                        + ' Default is 10')
    parser.add_argument('-l', '--lambdafile', required=False, type=str,
                        default = 'lambdagrouptypes.dat',
                        help = 'The path to lambdagrouptypes.dat file. '
                        + 'By default the script will search a file in the '
                        + 'current directory')
    parser.add_argument('-g', '--group', required=True, type=str,
                        help='Name of the group for parameterisation')
    parser.add_argument('-fo', '--fitorder', required=False, type=int,
                        default = 5, 
                        help = "Fitting order. Default value is 5")
    parser.add_argument('-o', '--out', required=False, type=str,
                        default='out.dat',
                        help='Name of the outputfile')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    # set input filenames
    mdpfile = args.mdp
    lambdafile = args.lambdafile

    # get mode and prefixes of folders
    mode = args.mode
    if not args.prefix and mode == 'p':
        prefix = 'r'
    elif not args.prefix and mode == 's':
        prefix = 's'
    else:
        prefix = args.prefix

    # get group name to analyse
    group_name = args.group

    # get fitting order
    fit_order = args.fitorder

    # get nreplicas
    nreplicas = args.nreplicas

    # get output file name
    output = args.out

    if mode == 'p':
        new_dvdl = fit_param(fit_order, prefix)
    else:
        new_dvdl = run_reweighing(fit_order, prefix, nreplicas)

    lambda_group = LambdaGroup(group_name, mdpfile, lambdafile, output)

    if mode == 'p':
        lambda_group.group.addDvdl(new_dvdl.tolist(), 1)
    else:
        lambda_group.group.addDvdl(
            (
                np.array(lambda_group.group.dvdls[1]) + new_dvdl
            ).tolist(), 1
        )
    lambda_group.write_output()
