#!/bin/bash

# Prepare the topology for 1cvo.pdb
phbuilder gentopol -f 1cvo.pdb

# non-phbuilder stuff:
gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5
gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

# add appropriate ions and buffers for neutralizing the system.
phbuilder neutralize -f solvated.pdb

# generate the .mdp (files and) cpHMD-parameters
phbuilder genparams -f phneutral.pdb -ph 4.0 -inter
