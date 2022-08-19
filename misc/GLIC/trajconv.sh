for sim in 4HFI_4 4HFI_7 6ZGD_4 6ZGD_7
do
    cd $sim
    for rep in 01 02 03 04
    do
        cd $rep
        echo 1 0 | gmx trjconv -f MD.xtc -o MD_conv.xtc -s MD.tpr -e 1000000 -dt 1000 -pbc mol -center
        rm -f \#*\#
        cd ..
    done
    cd ..
done
