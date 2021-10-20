#!/bin/bash

# Gromacs version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 2
gmx mdrun -v -deffnm EM -c EM.pdb

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r EM.pdb -maxwarn 1
gmx mdrun -v -deffnm NVT -c NVT.pdb -notunepme

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r NVT.pdb -maxwarn 1
gmx mdrun -v -deffnm NPT -c NPT.pdb -notunepme
