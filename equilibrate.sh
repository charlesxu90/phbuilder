#!/bin/bash

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -o EM.tpr -maxwarn 1
gmx mdrun -v -deffnm EM -c EM.pdb

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -o NVT.tpr -r EM.pdb -maxwarn 1
gmx mdrun -v -deffnm NVT -c NVT.pdb -notunepme

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -o NPT.tpr -r NVT.pdb -maxwarn 1
gmx mdrun -v -deffnm NPT -c NPT.pdb -notunepme
