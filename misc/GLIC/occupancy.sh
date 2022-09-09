#!/bin/bash

source /usr/local/gromacs_constantph/bin/GMXRC

for sim in 4HFI_4 4HFI_7 6ZGD_4 6ZGD_7
do
    cd $sim
    for rep in 01 02 03 04
    do
        cd $rep

        # Extract the .pdb trajectory from MD_conv.xtc
        gmx cphmd -s MD.tpr -e MD.edr -f MD_conv.xtc -fo MD_occ.pdb

        # Remove solvent to save on size
        sed -i '/SOL/d'  MD_occ.pdb

        # Remove membrane to save on size
        sed -i '/POPC/d' MD_occ.pdb

        # Add back the chain Identifiers
        cp ../../repairchain.py .
        python3 repairchain.py

        # Cleanup
        mv MD_occ_chain.pdb MD_occ.pdb
        rm -f cphmd-coord-*-*.xvg repairchain.py \#*\#

        cd ..
    done
    cd ..
done
