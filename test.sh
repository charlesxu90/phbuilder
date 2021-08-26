#!/bin/bash

echo 1 | phbuilder gentopol -f proteins/1cvo.pdb

gmx editconf -f phprocessed.pdb -o phprocessed.pdb -bt triclinic -d 2
gmx solvate -cp phprocessed.pdb -p topol.top -o phprocessed.pdb

touch IONS.mdp
gmx grompp -f IONS.mdp -c phprocessed.pdb -p topol.top -o IONS.tpr
echo SOL | gmx genion -s IONS.tpr -o phprocessed.pdb -p topol.top -neutral

phbuilder addbuffers -f phprocessed.pdb -o phprocessed.pdb
# phbuilder genparams -f phprocessed.pdb -ph 4
