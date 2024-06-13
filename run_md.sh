echo "starting phbuilder"
echo -e "1" > input.data

GROMACS_PATH=~/enzyme-md/gmx_cph/install/bin
phbuilder gentopol -f original.pdb -ph 4.0 < input.data

$GROMACS_PATH/gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.2 
$GROMACS_PATH/gmx solvate -cp box.pdb -p topol.top -o solvated.pdb

phbuilder neutralize -f solvated.pdb

phbuilder genparams -f phneutral.pdb -ph 4.0 -time 100 -temp 374