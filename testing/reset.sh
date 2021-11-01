#!/bin/bash

rm -rf __py* charmm36*
rm -f phset.pdb phprocessed.pdb phrecord.dat phions.pdb phneutral.pdb
rm -f residuetypes.dat box.pdb solvated.pdb pdb2gmxtemp.pdb

rm -f \#*
rm -f *.top *.gro *.itp
rm -f builder.log EM.pdb NVT.pdb NPT.pdb MD.pdb lambda_*.dat
rm -f *.ndx *.mdp *.tpr *.log *.edr *.out step*.pdb *.cpt *.trr *.xtc
