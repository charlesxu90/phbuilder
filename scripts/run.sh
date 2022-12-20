#!/bin/bash

# Gromacs version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f MD.mdp -c NPT.gro -p topol.top -n index.ndx -o MD.tpr -maxwarn 1
gmx mdrun -v -deffnm MD -c MD.pdb -x MD.xtc 
