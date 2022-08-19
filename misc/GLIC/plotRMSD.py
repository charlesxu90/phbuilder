#!/bin/python3

from analysis import loadxvg
import matplotlib.pyplot as plt

for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
    for rep in ['01', '02', '03', '04']:

        for chain in ['A', 'B', 'C', 'D', 'E']:
            data = loadxvg('{}/{}/{}.xvg'.format(sim, rep, chain))
            plt.plot([0.001 * x for x in data[0]], [10 * x for x in data[1]], label=chain, linewidth=0.5)

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
