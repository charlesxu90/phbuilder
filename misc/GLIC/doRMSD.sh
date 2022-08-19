#!/bin/bash

for sim in 4HFI_4 4HFI_7 6ZGD_4 6ZGD_7
do
    cd $sim
    for rep in 01 02 03 04
    do
        cd $rep
        cp /home/anton/GIT/rex/misc/GLIC/misc/rmsd.ndx . # residues 18-200
        echo A | gmx rmsdist -f MD.xtc -s MD.tpr -n rmsd.ndx -o A.xvg -dt 1000 -e 1000000 &
        echo B | gmx rmsdist -f MD.xtc -s MD.tpr -n rmsd.ndx -o B.xvg -dt 1000 -e 1000000 &
        echo C | gmx rmsdist -f MD.xtc -s MD.tpr -n rmsd.ndx -o C.xvg -dt 1000 -e 1000000 &
        echo D | gmx rmsdist -f MD.xtc -s MD.tpr -n rmsd.ndx -o D.xvg -dt 1000 -e 1000000 &
        echo E | gmx rmsdist -f MD.xtc -s MD.tpr -n rmsd.ndx -o E.xvg -dt 1000 -e 1000000 &
        wait
        cd ..
    done
    cd ..
done
