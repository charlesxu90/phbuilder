#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np

import MDAnalysis
import MDAnalysis.analysis.rms

def loadxvg(fname, col=[0, 1], dt=1, b=0):
    """
    This function loads an .xvg file into a list of lists.
    fname: file name (e.g. 'rmsd.xvg')
    col: which columns to load
    dt: skip every dt steps
    b: start from ... in first column
    """
    count = -1
    data = [ [] for _ in range(len(col)) ]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        # This is for the dt part.
        count += 1
        if count % dt != 0:
            continue
        
        listLine = stringLine.split()
        # And this is for the b part.
        if b != 0 and float(listLine[col[0]]) < b:
            continue

        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data

# Handle user input.
def inputOptionHandler(message, options):
    '''Handles user input'''

    valids = []
    msgstring = "{}:".format(message)

    # Loop through the options list and create string for dislay
    for idx in range(0, len(options)):
        msgstring += "\n{}. {}".format(idx, options[idx])
        valids.append(str(idx))
        
    while True:
        print(msgstring)
        val = input("Type a number: ")

        if val in valids:
            print()
            return int(val)

        print("{} is not a valid option, please try again:\n".format(val))

def getLambdaFileIndices(structure, residue):
    """
    Returns an array containing the lambda-file indices for the specified residue.
    Example: MD.pdb, residue: E35
    """

    u     = MDAnalysis.Universe(structure).segments[0].atoms
    array = []

    for res in u.select_atoms('resname ASPT GLUT HSPT').residues:

        resname = str(res.resname)
        for idx in range(0, 3):
            if resname == ['ASPT','GLUT','HSPT'][idx]:
                resname = ['D','E','H'][idx]
                break

        array.append(resname + str(res.resid))

    count  = 1
    factor = len(u.select_atoms('resname ASPT GLUT').residues) + 3 * len(u.select_atoms('resname HSPT').residues)
    
    for idx in range(0, len(array)):
        if array[idx] == residue:
            return [count, count+factor, count+factor*2, count+factor*3, count+factor*4]
        
        if array[idx][0] in ['D','E']:
            count += 1
        elif array[idx][0] == 'H':
            count += 3

################################################################################

for simulation in ['4HFI_4']:
# for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:

    # GET THE CHARGES ON SPECIFIED RESIDUE(S) AND CHAIN(S) VERSUS TIME #########

    for residue in ['D32', 'E35']:
        chain = ['A', 'B', 'C', 'D', 'E']
        array = getLambdaFileIndices('{}/01/CA.pdb'.format(simulation), residue)

        for idx in range(0, len(array)):
            data = loadxvg('{}/01/cphmd-coord-{}.xvg'.format(simulation, array[idx]), dt=1000, b=0)
            t    = [val / 1000.0 for val in data[0]] # ps -> ns
            x    = [1.0 - val for val in data[1]] # deprotonation -> protonation

            plt.plot(t, x, linewidth=0.5, label=chain[idx])

        plt.ylabel('Protonation')
        plt.xlabel('Time (ns)')
        plt.legend()
        plt.tight_layout()
        plt.savefig('panels/{}_{}.png'.format(simulation, residue))
        plt.clf(); plt.close()

    # GET THE RMSD OF SPECIFIED SECTION(S) VERSUS TIME #########################

    if inputOptionHandler('Run rms analysis for {} again?'.format(simulation), ['no', 'yes']):

        u = MDAnalysis.Universe('{}/01/CA.pdb'.format(simulation), '{}/01/MD.xtc'.format(simulation))

        R = MDAnalysis.analysis.rms.RMSD(u, select='(segid A B C D E) and backbone and resid 32-35')
        R.run(step=100) # speed things up a bit

        t  = [val / 1000.0 for val in R.rmsd.T[1]]
        x1 = R.rmsd.T[2]

        plt.plot(t, x1)
        plt.xlabel("time (ns)")
        plt.ylabel(r"RMSD ($\AA$)")
        plt.legend()
        plt.tight_layout()
        plt.savefig('panels/{}_{}-{}.png'.format(simulation, 32, 35))
        plt.clf(); plt.close()

    # MAKE THE MINDIST PLOT FOR SPECIFIED RESIDUE COMBINATION ##################
    
    u  = MDAnalysis.Universe('{}/01/CA.pdb'.format(simulation), '{}/01/MD.xtc'.format(simulation))

    selA = u.select_atoms('segid A and resid 243')
    selB = u.select_atoms('segid A and resid 248')

    t = []
    d = []
    for ts in u.trajectory[::100]:
        com1 = selA.center_of_mass()
        com2 = selB.center_of_mass()
        
        t.append(u.trajectory.time)
        d.append(np.linalg.norm(com1 - com2))

    t = [val / 1000.0 for val in t] # ps -> ns

    plt.plot(t, d)
    # plt.plot(t, d, legend='243-248 chain A')

    plt.xlabel("time (ns)")
    plt.ylabel(r"Minimum distance ($\AA$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('panels/{}_{}-{}_dist.png'.format(simulation, 243, 248))
    plt.clf(); plt.close()
