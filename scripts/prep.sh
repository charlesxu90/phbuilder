#!/bin/bash

echo 1 | phbuilder gentopol -f ../misc/proteins/1cvo.pdb -auto 4.0

gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5
gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

phbuilder neutralize -f solvated.pdb
phbuilder genparams -f phneutral.pdb -ph 4.0
