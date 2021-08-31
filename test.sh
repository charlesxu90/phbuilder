#!/bin/bash

echo 1 | phbuilder gentopol -f proteins/1cvo.pdb

echo phbuilder : adding BOX...
gmx editconf -f phprocessed.pdb -o phprocessed.pdb -bt triclinic -d 1.5 >> builder.log 2>&1

echo phbuilder : adding Solvent...
gmx solvate -cp phprocessed.pdb -p topol.top -o phprocessed.pdb >> builder.log 2>&1

echo phbuilder : adding Ions...
touch IONS.mdp >> builder.log 2>&1
gmx grompp -f IONS.mdp -c phprocessed.pdb -p topol.top -o IONS.tpr >> builder.log 2>&1
echo SOL | gmx genion -s IONS.tpr -o phprocessed.pdb -p topol.top -neutral >> builder.log 2>&1

phbuilder addbuffers -f phprocessed.pdb -o phprocessed.pdb -nbufs 10
phbuilder genparams -f phprocessed.pdb -ph 4.5 -dwpE 5.0
