#!/bin/bash

source /usr/local/gromacs_constantph/bin/GMXRC

for sim in 4HFI_4 4HFI_7 6ZGD_4 6ZGD_7
do
    cd $sim
    for rep in 01 02 03 04
    do
        cd $rep
        
        echo 0   | gmx trjconv -s MD.tpr -n index.ndx -f MD.xtc -o A.xtc -pbc nojump -dt 1000 -e 1000000
        echo 1 0 | gmx trjconv -s MD.tpr -n index.ndx -f A.xtc -o B.xtc -center
        echo 0   | gmx trjconv -s MD.tpr -n index.ndx -f B.xtc -o MD_conv.xtc -pbc mol
        rm A.xtc B.xtc

        rm -f \#*\#
        cd ..
    done
    cd ..
done
