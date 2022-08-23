#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis
import MDAnalysis.analysis.rms

# MDANALYSIS DOES NOT KNOW WHAT THE BACKBONE OF ASPT GLUT HSPT ARE!!!!!

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

def inputOptionHandler(message, options):
    """
    Function for handling user input.
    message: the string you would like to prompt the user.
    options: a list of strings containing the options.
    """

    valids = []
    msgstring = "{}:".format(message)

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

def chargePlot(sim, rep, res):
    """
    Make the charge plot in time for residue. Does not work for histidines.
    sim: the simulation, e.g. 4HFI_4.
    rep: the replica, e.g. 1.
    res: the residue, e.g. E35.
    """

    chain = ['A', 'B', 'C', 'D', 'E']
    array = getLambdaFileIndices('{}/{:02d}/CA.pdb'.format(sim, rep), res)

    store = []
    for idx in range(0, len(array)):
        data = loadxvg('{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, array[idx]), dt=5000, b=0)
        t    = [val / 1000.0 for val in data[0]] # ps -> ns
        x    = [1.0 - val for val in data[1]] # deprotonation -> protonation
        store.append(x)

        plt.plot(t, x, linewidth=0.5, label=chain[idx])

    # PLOT MEAN PROTONATION
    average = [0] * len(store[0])
    for idx in range(0, len(store[0])):
        average[idx] = store[0] + store[1] + store[2] + store[3] + store[4]
    plt.plot(t, x, linewidth=1.5, label='mean', color='b', linestyle=':')

    plt.ylabel('Protonation')
    plt.xlabel('Time (ns)')
    plt.axis([0, 1000, -0.1, 1.1])
    plt.title('{} protonation {}'.format(sim, res))
    plt.legend()
    plt.tight_layout()
    plt.savefig('panels/proto_{}_{}.png'.format(sim, res))
    plt.clf(); plt.close()

def RMSDPlot(sim, rep, sel):
    """
    Creates RMSD plots for a selection of residues. Load MD_whole.xtc.
    sim: the simulation, e.g. 4HFI_4.
    rep: the replica, e.g. 1.
    sel: selection, e.g. '32-35'
    """

    path1 = '{}/{:02d}/CA.pdb'.format(sim, rep)
    path2 = '{}/{:02d}/MD_whole.xtc'.format(sim, rep)
    u = MDAnalysis.Universe(path1, path2)

    # INDIVIDUAL CHAINS
    chain = ['A', 'B', 'C', 'D', 'E']
    for idx in range(0, len(chain)):
        R = MDAnalysis.analysis.rms.RMSD(u, select='segid {} and resid {}'.format(chain[idx], sel))
        R.run()
        t  = [val / 1000.0 for val in R.rmsd.T[1]]
        x1 = R.rmsd.T[2]
        plt.plot(t, x1, label=chain[idx])

    # ALL CHAINS
    R = MDAnalysis.analysis.rms.RMSD(u, select='(segid A B C D E) and resid {}'.format(sel))
    R.run()
    t  = [val / 1000.0 for val in R.rmsd.T[1]]
    x1 = R.rmsd.T[2]
    plt.plot(t, x1, label='all')

    plt.xlabel("time (ns)")
    plt.ylabel(r"RMSD ($\AA$)")
    plt.xlim(0, 1000)
    plt.title('{} RMSD segment {}'.format(sim, sel))
    plt.legend()
    plt.tight_layout()
    plt.savefig('panels/rmsd_{}_{}.png'.format(sim, sel))
    plt.clf(); plt.close()

def mindistPlot(sim, rep, sel1, sel2, name):
    """
    Creates mindists plots between sel1 and sel2.
    sim: the simulation, e.g. 4HFI_4.
    rep: the replica, e.g. 1.
    sel1: selection1, e.g. 'resid 243'
    sel2: selection2, e.g. 'resid 248'
    """

    path1 = '{}/{:02d}/CA.pdb'.format(sim, rep)
    path2 = '{}/{:02d}/MD.xtc'.format(sim, rep)
    u  = MDAnalysis.Universe(path1, path2)

    chain = ['A', 'B', 'C', 'D', 'E']
    for idx in range(0, len(chain)):

        selA = u.select_atoms('segid {} and {}'.format(chain[idx], sel1))
        selB = u.select_atoms('segid {} and {}'.format(chain[idx], sel2))

        t = []; d = []
        for ts in u.trajectory[::100]:
            com1 = selA.center_of_mass()
            com2 = selB.center_of_mass()

            t.append(u.trajectory.time)
            d.append(np.linalg.norm(com1 - com2))

        t = [val / 1000.0 for val in t] # ps -> ns

        plt.plot(t, d, label='chain {}'.format(chain[idx]))

    plt.xlabel("time (ns)")
    plt.xlim(0, 1000)
    plt.ylabel(r"Minimum distance ($\AA$)")
    plt.title('{} minidst {}'.format(sim, name))
    plt.legend()
    plt.tight_layout()
    plt.savefig('panels/mindst_{}_{}.png'.format(sim, name))
    plt.clf(); plt.close()

chargePlot('4HFI_4', 1, 'D32')
chargePlot('4HFI_4', 1, 'E35')
RMSDPlot('4HFI_4', 1, '32-35')
mindistPlot('4HFI_4', 1, 'resid 243', 'resid 248', 'E243-K248')
