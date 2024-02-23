#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np

from sympy.parsing.sympy_parser import parse_expr
from sympy import symbols

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Prepare files for titration simulations. \
        The code will create directories, modify .mdp file, \
        run grompp and if needed mdrun for specified pH values\
        and for the desired number of replicas. Also, if some pKa\
        values have to be set to pH, this can be specified by \
        setting the corresponding \
        [lambda-dynamics-group-type1-state-1-reference-pka] to pH.\
        The code will automtically recognize it and set the \
        proper pKa.')

    parser.add_argument('-f', '--mdp', required=True, type=str,
                        help='Input .mdp file')
    parser.add_argument('-c', '--struct', required=True, type=str,
                        help='Input structure file')
    parser.add_argument('-r', '--ref', required=False, type=str,
                        help='Input reference structure')
    parser.add_argument('-p', '--top', required=True, type=str,
                        help='Input topology file')
    parser.add_argument('-n', '--index', required=True, type=str,
                        help='Input index file')
    parser.add_argument('-pH', '--pHrange', required=True, type=str,
                        help='pH values at which we run simulations. \
                        Those can be provided in the form of range "a:b:d" \
                        or as a list "{a,b,c,..,}"')
    parser.add_argument('-nr', '--nrep', required=True, type=int,
                        help='Number of replicas we run for each pH')
    parser.add_argument('-o', '--out', required=False, type=str,
                        default='pH', help='Prefix name')
    parser.add_argument('-gmx', '--gmxpath', required=False, type=str,
                        default='gmx', help='GMXPATH. Default is gmx')
    parser.add_argument('-rc', '--gmxrun', required=False, type=str,
                        default=None, help='Run gmx. Default is None')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    mdpp = os.path.abspath(args.mdp)
    mdpn = os.path.basename(os.path.normpath(mdpp))
    strp = os.path.abspath(args.struct)
    topp = os.path.abspath(args.top)
    ndxp = os.path.abspath(args.index)
    pref = args.out
    n = args.nrep

    if args.ref:
        refp = os.path.abspath(args.ref)

    if args.pHrange.startswith("{"):
        pH = np.array(list(float(x) for x in args.pHrange.strip()[1:-1].split(",")))
    else:
        pH = np.arange(float(args.pHrange.split(":")[0]), float(args.pHrange.split(":")[1]), float(args.pHrange.split(":")[2]))

    for p in pH:
        if os.path.exists("{}_{:.2f}".format(pref, p)):
            os.system("rm -r {}_{:.2f}".format(pref, p))
        os.mkdir("{}_{:.2f}".format(pref, p))
        for i in range(1, n + 1):
            os.mkdir("{}_{:.2f}/r_{}".format(pref, p, i))
            os.system("cp {} {}_{:.2f}/r_{}/.".format(mdpp, pref, p, i))
            # change mdp file
            content = ""
            with open("{}_{:.2f}/r_{}/{}".format(pref, p, i, mdpn), "r") as f:
                for line in f:
                    # change pH
                    if line.startswith("lambda-dynamics-simulation-ph"):
                        newline = line.replace(line, "lambda-dynamics-simulation-ph        = {:.2f}".format(p))
                        content += newline + "\n"
                    # change pKa to pH
                    elif line.startswith("lambda-dynamics-group-type") and all(x in line for x in ["pka", "ph"]):
                        ph = symbols("ph")
                        _key = line.split("=")[0]
                        _other = line[line.index("=") + 1:]
                        try:
                            _comment = _other[_other.index(";") + 1:]
                            pKa_str = _other.split(";")[0]
                        except ValueError:
                            _comment = None
                            pKa_str = _other
                        pKaNew = float(parse_expr(pKa_str, evaluate=False).subs(ph, p))
                        newline = f"{_key}   = {pKaNew} ;{_comment}" if _comment else f"{_key}   = {pKaNew}\n"
                        content += newline
                    else:
                        content += line
            # save changes
            with open("{}_{:.2f}/r_{}/{}".format(pref, p, i, mdpn), "w") as f:
                f.write(content)
            gmxCommand = "{} grompp -p {} -f {} -c {} -n {} -o {}_{:.2f}/r_{}/run.tpr".format(
                args.gmxpath, topp, "{}_{:.2f}/r_{}/{}".format(pref, p, i, mdpn), strp, ndxp, pref, p, i)
            if args.ref:
                gmxCommand += " -r {}".format(refp)
            os.system(gmxCommand)
            if (args.gmxrun):
                os.chdir("{}_{:.2f}/r_{}".format(pref, p, i))
                os.system(args.gmxrun)
                os.chdir("../../")
