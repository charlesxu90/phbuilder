#!/bin/python3

import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib.pyplot as plt
import os

for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
    for rep in [1, 2, 3, 4]:

        path1 = '{}/{:02d}/CA.pdb'.format(sim, rep)
        path2 = '{}/{:02d}/MD_whole.xtc'.format(sim, rep)
        u = MDAnalysis.Universe(path1, path2)

        for seg in ['A', 'B', 'C', 'D', 'E']:
            R = MDAnalysis.analysis.rms.RMSD(u, select='segid {} and (name C CA N) and resid 18-200'.format(seg))
            R.run()

            t  = [val / 1000.0 for val in R.rmsd.T[1]]
            x1 = R.rmsd.T[2]
            plt.plot(t, x1, label=seg, linewidth=0.5)

        plt.title('{} replica {}'.format(sim, rep))
        plt.ylabel(r'RMSD ($\AA$)')
        plt.xlabel('Time (ns)')
        plt.legend()
        plt.xlim([0, 1000])
        plt.ylim([0, 5])
        plt.grid()
        plt.tight_layout()
        plt.savefig('{}_{}.png'.format(sim, rep))
        plt.clf(); plt.close()

# Make the final image and cleanup
os.system('convert 4HFI_4_*.png +append A.png')
os.system('convert 4HFI_7_*.png +append B.png')
os.system('convert 6ZGD_4_*.png +append C.png')
os.system('convert 6ZGD_7_*.png +append D.png')
os.system('convert A.png B.png C.png D.png -append ECDRMSD.png')
os.system('rm -f A.png B.png C.png D.png 4HFI_*.png 6ZGD_*.png')
