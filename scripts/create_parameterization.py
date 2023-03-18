import argparse
import os
import sys
import re
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Prepare files for calibration of 2-state \
        systems. The code will create directories, modify .mdp file, \
        run grompp and if needed mdrun for fixed lambda values \
        in range from -0.1 to 1.1 with 0.05 step. The basic .mdp \
        file is expected to have only two group: titratable and \
        neutralizing buffer. Also position restraints for the \
        buffer and correct atoms of parameterised group should be \
        provided. We expect the user to select those atoms. The \
        recommended strategy is to fix one atom close to the \
        group of atoms that change charge and leave all other \
        system free. This will allow to sample different \
        orientations and conformations of the parameterised \
        molecule, while preserving the distance to neutralising ion.')

    parser.add_argument('-f', '--mdp', required=True, type=str,
                        help='Input .mdp file')
    parser.add_argument('-c', '--struct', required=True, type=str,
                        help='Input structure file')
    parser.add_argument('-r', '--ref', required=True, type=str,
                        help='Input reference structure for position restraining')
    parser.add_argument('-p', '--top', required=True, type=str,
                        help='Input topology file')
    parser.add_argument('-n', '--index', required=True, type=str,
                        help='Input index file')
    parser.add_argument('-o', '--out', required=False, type=str,
                        default='r',
                        help='Prefix name for output folders. Default is r')
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
    print(mdpn)
    strp = os.path.abspath(args.struct)
    topp = os.path.abspath(args.top)
    ndxp = os.path.abspath(args.index)
    pref = args.out
    refp = os.path.abspath(args.ref)

    for l in np.arange(-0.1, 1.15, 0.05):
        folder = f"{pref}_{l:.2f}_{(1-l):.2f}"
        if os.path.exists(f"{folder}"):
            os.system(f"rm -r {folder}")
        os.mkdir(f"{folder}")

        os.system(f"cp {mdpp} {folder}/.")
        # read mdp file
        content = []
        atomSetNaming = {}
        with open(f"{folder}/{mdpn}", "r") as f:
            for line in f:
                # search for lambda-dynamics-atom-setX-name lines
                if re.search("lambda-dynamics-atom-set.-name", line):
                    content.append(line)
                    # Add entry to dictionary
                    atomSetNaming[line.split()[0].split("-")[-2][-1]] = (
                        line.split()[-1]
                    )
                else:
                    content.append(line)

        # modify mdp
        output = ""
        print(atomSetNaming)
        for line in content:
            if not re.search("lambda-dynamics-atom-set.-initial-lambda", line):
                output += line
                continue
            for idx in atomSetNaming:
                if f"lambda-dynamics-atom-set{idx}-initial-lambda" in line:
                    val = l if atomSetNaming[idx] != "BUF" else (1 - l)
                    output += (
                        f"lambda-dynamics-atom-set{idx}-initial-lambda = {val:.2f}\n"
                    )

        # save changes
        with open(f"{folder}/{mdpn}", "w") as f:
            f.write(output)
        gmxCommand = (
            f"{args.gmxpath} grompp -p {topp} -f {f'{folder}/{mdpn}'} "
            + f"-c {strp} -n {ndxp} -o {folder}/run.tpr -r {refp}"
        )
        os.system(gmxCommand)
        if (args.gmxrun):
            os.chdir(f"{folder}")
            os.system(args.gmxrun)
            os.chdir("../")
