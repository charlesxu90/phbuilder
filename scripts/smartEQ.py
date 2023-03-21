#!/bin/python3

import os
import MDAnalysis

from science.utility import gromacs


def extractAndUpdate(base, new):
    """Extract the last lambda-coordinate from a previous simulation (base) and
    use those as starting lambdas for the new (new) simulation.
    Args:
        base (str): old simulation (e.g. NVT).
        new (str): new simulation (e.g. NPT).
    """

    def lastLambdaVal(fname):
        "Retrieves the last lambda coordinate from a lambda file."

        with open(fname) as file:
            for line in file:
                pass
            return float(line.split()[1])

    # Directory stuff
    if not os.path.exists(base):
        os.mkdir(base)

    # Extract and move the lambda-coordinate files.
    gromacs('cphmd -s {}.tpr -e {}.edr -numplot 1'.format(base, base))
    os.system('mv cphmd-coord-*.xvg {}'.format(base))

    values = []
    for num in range(1, len(os.listdir(base)) + 1):
        val = lastLambdaVal('{}/cphmd-coord-{}.xvg'.format(base, num))
        values.append(val)
        print(val)  # debug

    u = MDAnalysis.Universe('EM.gro')
    acidics = list(u.select_atoms('resname ASPT GLUT HSPT').residues.resnames)
    acidics.append('BUF')
    print(acidics)  # debug

    lidx = 0
    for atomset in range(0, len(acidics)):

        if acidics[atomset] == 'HSPT':
            newCoord = str(values[lidx]) + ' ' + str(values[lidx + 1]) + ' ' + str(values[lidx + 2])
            lidx += 3

        if acidics[atomset] in ['ASPT', 'GLUT', 'BUF']:
            newCoord = str(values[lidx])
            lidx += 1

        print(acidics[atomset], newCoord)  # debug

        match   = 'lambda-dynamics-atom-set{}-initial-lambda'.format(atomset + 1)
        replace = match + '               = {}'.format(newCoord)

        # print(match)    # debug
        # print(replace)  # debug

        # Replace the line that contains A with B:
        os.system(f'sed -i \'s/.*{match}.*/{replace}/\' {new}.mdp')

# EM (calibration = True)
gromacs('grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1')
gromacs('mdrun -v -deffnm EM')

# Remove the 'lambda-dynamics-calibration = yes' line from NVT.mdp, NPT.mdp, MD.mdp if it is present.
for fname in ['NVT.mdp', 'NPT.mdp', 'MD.mdp']:
    print(f"Checking/removing 'lambda-dynamics-calibration = yes' from {fname} (if it is present)")
    os.system(f"sed -i  '/lambda-dynamics-calibration                            = yes/d' {fname}")

# NVT (calibration = False, nstout = 1)
gromacs('grompp -f NVT.mdp -c EM.gro -p topol.top -n index.ndx -o NVT.tpr -maxwarn 1')
gromacs('mdrun -v -deffnm NVT -notunepme')
extractAndUpdate('NVT', 'NPT')

# NPT (calibration False, nstout = 1)
gromacs('grompp -f NPT.mdp -c NVT.gro -p topol.top -n index.ndx -o NPT.tpr -maxwarn 1')
gromacs('mdrun -v -deffnm NPT -notunepme')
extractAndUpdate('NPT', 'MD')
