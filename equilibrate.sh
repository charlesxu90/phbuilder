#!/bin/bash

gmx grompp -f EM.mdp -c phprocessed.pdb -p topol.top -o EM.tpr -r phprocessed.pdb
gmx mdrun -v -deffnm EM -c EM.pdb

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -o NVT.tpr -r EM.pdb
gmx mdrun -v -deffnm NVT -c NVT.pdb

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -o NPT.tpr -r NVT.pdb
gmx mdrun -v -deffnm NPT -c NPT.pdb
