#!/bin/bash

echo 1 | phbuilder gentopol -f proteins/1cvo.pdb

gmx editconf -f phprocessed.pdb -o phprocessed.pdb -bt triclinic -d 1
gmx solvate -cp phprocessed.pdb -p topol.top -o phprocessed.pdb

phbuilder addbuffers -f phprocessed.pdb -o phprocessed.pdb

phbuilder genparams -f phprocessed.pdb -ph 4

# touch ions.mdp
# gmx grompp ...
