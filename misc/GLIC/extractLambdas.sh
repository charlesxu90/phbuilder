#!/bin/bash

# Load gmx-cphmd version
source /usr/local/gromacs_constantph/bin/GMXRC

for sim in 4HFI_4 4HFI_7 6ZGD_4 6ZGD_7
do
    cd $sim
    {
        cd 01
        addChainIdentifiers CA.pdb CA.pdb
        gmx cphmd -s MD.tpr -e MD.edr -numplot 1
        rm -f \#*\#
        cd ..
    } &
    {
        cd 02
        addChainIdentifiers CA.pdb CA.pdb
        gmx cphmd -s MD.tpr -e MD.edr -numplot 1
        rm -f \#*\#
        cd ..
    } &
    {
        cd 03
        addChainIdentifiers CA.pdb CA.pdb
        gmx cphmd -s MD.tpr -e MD.edr -numplot 1
        rm -f \#*\#
        cd ..
    } &
    {
        cd 04
        addChainIdentifiers CA.pdb CA.pdb
        gmx cphmd -s MD.tpr -e MD.edr -numplot 1
        rm -f \#*\#
        cd ..
    } &
    wait
    cd ..
done
