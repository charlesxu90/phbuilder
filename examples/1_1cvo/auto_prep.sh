#!/bin/bash

# Use the -ph option with pH = 4.0 to automatically guess the initial lambda
# states. In addition, we echo 1 to prevent prompting for pdb2gmx.
echo 1 | phbuilder gentopol -f 1cvo.pdb -ph 4.0

# non-phbuilder stuff:
gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5
gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

# add appropriate ions and buffers for neutralizing the system.
phbuilder neutralize -f solvated.pdb

# generate the .mdp (files and) cpHMD-parameters
# echo 0 to prevent prompting for barostat parameters.
echo 0 | phbuilder genparams -f phneutral.pdb -ph 4.0
