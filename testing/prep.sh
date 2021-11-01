#!/bin/bash

echo 1 | phbuilder gentopol -f ../proteins/1cvo.pdb -ph 7.0 -auto

gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5

gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

phbuilder neutralize -f solvated.pdb
phbuilder genparams -f phneutral.pdb -ph 7.0
