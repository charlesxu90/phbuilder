#!/bin/bash

echo 1 | phbuilder gentopol -f proteins/1cvo.pdb

# echo phbuilder : adding BOX...
gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5

# echo phbuilder : adding Solvent...
gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

# echo phbuilder : adding Ions...
touch IONS.mdp
gmx grompp -f IONS.mdp -c solvated.pdb -p topol.top -o IONS.tpr
echo SOL | gmx genion -s IONS.tpr -o ions.pdb -p topol.top -neutral

phbuilder addbuffers -f ions.pdb -o buffers.pdb -nbufs 10
phbuilder genparams -f buffers.pdb -ph 4.0
