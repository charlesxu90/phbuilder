#!/bin/bash

echo Running EM...
gmx grompp -f EM.mdp -c phprocessed.pdb -p topol.top -o EM.tpr -r phprocessed.pdb >> builder.log 2>&1
gmx mdrun -v -deffnm EM -c EM.pdb >> builder.log 2>&1

echo Running NVT...
gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -o NVT.tpr -r EM.pdb >> builder.log 2>&1
gmx mdrun -v -deffnm NVT -c NVT.pdb >> builder.log 2>&1

echo Running NPT...
gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -o NPT.tpr -r NVT.pdb >> builder.log 2>&1
gmx mdrun -v -deffnm NPT -c NPT.pdb >> builder.log 2>&1
