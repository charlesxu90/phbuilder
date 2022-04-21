#!/bin/python3

from phbuilder.structure import Structure

import matplotlib.pyplot as plt
from matplotlib import cm
import os
import numpy as np

biophys = { # also prevost2012
    'ASPT-13'  : 1,
    'ASPT-31'  : 1,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 1,
    'ASPT-86'  : 0,
    'ASPT-88'  : 0,
    'ASPT-91'  : 1,
    'ASPT-97'  : 1,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 1,
    'ASPT-153' : 1,
    'ASPT-154' : 1,
    'ASPT-161' : 1,
    'ASPT-178' : 1,
    'ASPT-185' : 1,
    'GLUT-14'  : 1,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 0,
    'GLUT-69'  : 1,
    'GLUT-75'  : 0,
    'GLUT-82'  : 0,
    'GLUT-104' : 1,
    'GLUT-147' : 1,
    'GLUT-163' : 1,
    'GLUT-177' : 0,
    'GLUT-181' : 1,
    'GLUT-222' : 1,
    'GLUT-243' : 0,
    'GLUT-272' : 1,
    'GLUT-282' : 1,
    'HSPT-127' : 1,
    'HSPT-235' : 1,
    'HSPT-277' : 0
}

nury2010 = { # this is also cheng2010, calimet2013
    'ASPT-13'  : 1,
    'ASPT-31'  : 1,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 1,
    'ASPT-86'  : 0,
    'ASPT-88'  : 0,
    'ASPT-91'  : 1,
    'ASPT-97'  : 1,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 1,
    'ASPT-153' : 1,
    'ASPT-154' : 1,
    'ASPT-161' : 1,
    'ASPT-178' : 1,
    'ASPT-185' : 1,
    'GLUT-14'  : 1,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 0,
    'GLUT-69'  : 0,
    'GLUT-75'  : 0,
    'GLUT-82'  : 0,
    'GLUT-104' : 1,
    'GLUT-147' : 1,
    'GLUT-163' : 1,
    'GLUT-177' : 0,
    'GLUT-181' : 1,
    'GLUT-222' : 1,
    'GLUT-243' : 0,
    'GLUT-272' : 1,
    'GLUT-282' : 1,
    'HSPT-127' : 1,
    'HSPT-235' : 1,
    'HSPT-277' : 0
}

fritsch2011 = {
    'ASPT-13'  : 0,
    'ASPT-31'  : 0,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 0,
    'ASPT-86'  : 0,
    'ASPT-88'  : 0,
    'ASPT-91'  : 0,
    'ASPT-97'  : 0,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 0,
    'ASPT-153' : 0,
    'ASPT-154' : 0,
    'ASPT-161' : 0,
    'ASPT-178' : 0,
    'ASPT-185' : 0,
    'GLUT-14'  : 0,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 0,
    'GLUT-69'  : 0,
    'GLUT-75'  : 0,
    'GLUT-82'  : 0,
    'GLUT-104' : 1,
    'GLUT-147' : 0,
    'GLUT-163' : 0,
    'GLUT-177' : 0,
    'GLUT-181' : 0,
    'GLUT-222' : 1,
    'GLUT-243' : 0,
    'GLUT-272' : 0,
    'GLUT-282' : 0,
    'HSPT-127' : 0,
    'HSPT-235' : 1,
    'HSPT-277' : 1
}

lev2017 = {
    'ASPT-13'  : 1,
    'ASPT-31'  : 1,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 1,
    'ASPT-86'  : 1,
    'ASPT-88'  : 1,
    'ASPT-91'  : 1,
    'ASPT-97'  : 1,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 1,
    'ASPT-153' : 1,
    'ASPT-154' : 1,
    'ASPT-161' : 1,
    'ASPT-178' : 1,
    'ASPT-185' : 1,
    'GLUT-14'  : 1,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 0,
    'GLUT-69'  : 0,
    'GLUT-75'  : 0,
    'GLUT-82'  : 0,
    'GLUT-104' : 1,
    'GLUT-147' : 1,
    'GLUT-163' : 1,
    'GLUT-177' : 0,
    'GLUT-181' : 1,
    'GLUT-222' : 1,
    'GLUT-243' : 0,
    'GLUT-272' : 1,
    'GLUT-282' : 1,
    'HSPT-127' : 0,
    'HSPT-235' : 1,
    'HSPT-277' : 1
}

nemecz2017 = { # also Hu2018
    'ASPT-13'  : 1,
    'ASPT-31'  : 1,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 1,
    'ASPT-86'  : 0,
    'ASPT-88'  : 0,
    'ASPT-91'  : 1,
    'ASPT-97'  : 1,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 1,
    'ASPT-153' : 1,
    'ASPT-154' : 1,
    'ASPT-161' : 1,
    'ASPT-178' : 1,
    'ASPT-185' : 1,
    'GLUT-14'  : 1,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 1,
    'GLUT-69'  : 1,
    'GLUT-75'  : 1,
    'GLUT-82'  : 1,
    'GLUT-104' : 1,
    'GLUT-147' : 1,
    'GLUT-163' : 1,
    'GLUT-177' : 1,
    'GLUT-181' : 1,
    'GLUT-222' : 0,
    'GLUT-243' : 0,
    'GLUT-272' : 1,
    'GLUT-282' : 1,
    'HSPT-127' : 1,
    'HSPT-235' : 1,
    'HSPT-277' : 0
}

ullman = { # unpublished
    'ASPT-13'  : 1,
    'ASPT-31'  : 1,
    'ASPT-32'  : 1,
    'ASPT-49'  : 1,
    'ASPT-55'  : 1,
    'ASPT-86'  : 1,
    'ASPT-88'  : 1,
    'ASPT-91'  : 1,
    'ASPT-97'  : 1,
    'ASPT-115' : 1,
    'ASPT-122' : 1,
    'ASPT-136' : 1,
    'ASPT-145' : 1,
    'ASPT-153' : 1,
    'ASPT-154' : 1,
    'ASPT-161' : 1,
    'ASPT-178' : 1,
    'ASPT-185' : 1,
    'GLUT-14'  : 1,
    'GLUT-26'  : 0,
    'GLUT-35'  : 0,
    'GLUT-67'  : 0,
    'GLUT-69'  : 0,
    'GLUT-75'  : 0,
    'GLUT-82'  : 1,
    'GLUT-104' : 1,
    'GLUT-147' : 0,
    'GLUT-163' : 0,
    'GLUT-177' : 0,
    'GLUT-181' : 1,
    'GLUT-222' : 1,
    'GLUT-243' : 0,
    'GLUT-272' : 0,
    'GLUT-282' : 1,
    'HSPT-127' : 0,
    'HSPT-235' : 1,
    'HSPT-277' : 1
}

def loadxvg(fname, col=[0, 1]):
    data = [ [] for _ in range(len(col)) ]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        listLine = stringLine.split()
        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data

class TwoState:
    '''Holds data for a two-state (ASP or GLU) titratable group.'''

    def __init__(self, idx, resname, resid, chain, t, x):
        self.d_idx      = idx
        self.d_fname    = 'cphmd-coord-{}.xvg'.format(idx)
        self.d_resname  = resname
        self.d_resid    = resid
        self.d_chain    = chain
        self.d_t        = t
        self.d_x        = x

class MultiState:
    '''Holds data for a multistate-state (HIS) titratable group.'''

    def __init__(self, idx, resname, resid, chain, t1, t2, t3, x1, x2, x3):
        self.d_idx      = [idx, idx +1, idx +2]
        self.d_fname    = ['cphmd-coord-{}.xvg'.format(idx), 'cphmd-coord-{}.xvg'.format(idx+1), 'cphmd-coord-{}.xvg'.format(idx+2)]
        self.d_resname  = resname
        self.d_resid    = resid
        self.d_chain    = chain
        self.d_t        = [t1, t2, t3]
        self.d_x        = [x1, x2, x3]

class Buffer:
    '''Holds data for the buffer group.'''
    def __init__(self, idx, t, x, count):
        self.d_idx     = idx
        self.d_fname   = 'cphmd-coord-{}.xvg'.format(idx)
        self.d_resname = 'BUF'
        self.d_t       = t
        self.d_x       = x
        self.d_count   = count

class Analysis:
    '''Provides functionality for analyzing GLIC simulations.'''

    def __init__(self, file, name, pH, nstOut, dump):
        """
        file: filename, e.g. NPT.pdb or MD.pdb
        name: protein name, e.g. 1cvo.pdb
        pH: Simulation pH
        nstOut: nstout for the lambda_xxx.dat files
        dump: number of ns to burn
        """

        # Loaded
        self.d_file   = file
        self.d_name   = name
        self.d_pH     = pH
        self.d_dump   = int((1000 * dump) / (0.002 * nstOut))

        # Load d_residues
        self.d_residues = Structure(self.d_file, 2).d_residues

        # Parse all the data
        idx = 1
        foundBUF = False
        self.d_twoStateList   = []
        self.d_multiStateList = []

        print("Will skip first {} ns (= {} lines)".format(dump, self.d_dump))

        for residue in self.d_residues:

            if residue.d_resname in ['ASPT', 'GLUT']:
                print('Loading {}-{} in chain {}...'.format(residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                xvgdata = loadxvg('cphmd-coord-{}.xvg'.format(idx))

                self.d_twoStateList.append(TwoState(
                    idx, 
                    residue.d_resname, 
                    residue.d_resid,
                    residue.d_chain,
                    xvgdata[0][self.d_dump:], # time
                    xvgdata[1][self.d_dump:]  # coordinates
                    ))

                idx += 1

            if residue.d_resname == 'HSPT':
                print('Loading {}-{} in chain {}...'.format(residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                xvgdata1 = loadxvg('cphmd-coord-{}.xvg'.format(idx))
                xvgdata2 = loadxvg('cphmd-coord-{}.xvg'.format(idx+1))
                xvgdata3 = loadxvg('cphmd-coord-{}.xvg'.format(idx+2))

                self.d_multiStateList.append(MultiState(
                    idx,
                    residue.d_resname,
                    residue.d_resid,
                    residue.d_chain,
                    xvgdata1[0][self.d_dump:],  # file 1 time
                    xvgdata1[1][self.d_dump:],  # file 1 coordinates
                    xvgdata2[0][self.d_dump:],  # file 2 time
                    xvgdata2[1][self.d_dump:],  # file 2 coordinates
                    xvgdata3[0][self.d_dump:],  # file 3 time
                    xvgdata3[1][self.d_dump:])) # file 3 coordinates

                idx += 3

            if residue.d_resname == 'BUF' and not foundBUF:
                print('Loading buffer data')

                xvgdata = loadxvg('cphmd-coord-{}.xvg'.format(idx))

                self.d_buffer = Buffer(
                    idx,
                    xvgdata[0],    # no data dump
                    xvgdata[1],    # no data dump
                    count=185)

                foundBUF = True

        # Some user update (for debug)
        print('Length of twoStateList   = {}'.format(len(self.d_twoStateList)))
        print('Length of multiStateList = {}'.format(len(self.d_multiStateList)))

        # Directory structure
        self.d_dir = "lambdaplots"
        if not os.path.isdir(self.d_dir):
            os.mkdir(self.d_dir)

    def plotASPGLU_traj(self, checkConvergence=False):
        for group in self.d_twoStateList:
            print('plotting {}...'.format(group.d_fname), end='\r')

            plt.plot(group.d_t, group.d_x, linewidth=0.5)

            if checkConvergence:
                plt.plot(group.d_t, self.__movingDeprotonation(group.d_x), label='Moving deprotonation', color='r')

            plt.title('{}-{} in chain {} in {} ({})\npH={}, deprotonation={:.2f}'.format(
                group.d_resname,
                group.d_resid,
                group.d_chain,
                self.d_name,
                group.d_fname,
                self.d_pH,
                self.deprotonation(group.d_x)
            ))

            # Axes and stuff
            plt.ylim(-0.1, 1.1)
            plt.xlabel("Time (ps)")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.ylabel(r"$\lambda$-coordinate")
            plt.grid()
            plt.legend(loc='upper center')

            # Save as .png
            plt.savefig('{}/traj_{}_{}_{}.png'.format(self.d_dir, group.d_chain, group.d_resid, group.d_resname))
            # plt.savefig('{}/traj_{}_{}_{}.pdf'.format(self.d_dir, group.d_chain, group.d_resid, group.d_resname))

            # Clear
            plt.clf(); plt.close()

    def plotASPGLU_histo(self):
        numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
        resPerChain = int(len(self.d_twoStateList)/numChains)

        for ii in range(0, resPerChain):

            # For convenient access
            group = self.d_twoStateList[ii]

            # Gather data
            x = []
            for jj in range(0, numChains):
                x +=          self.d_twoStateList[ii + resPerChain * jj].d_x
                print(len(x), self.d_twoStateList[ii + resPerChain * jj].d_fname)

            print("Plotting histogram for {}-{}".format(group.d_resname, group.d_resid))

            plt.figure(figsize=(8, 6))
            plt.hist(x, density=True, bins=200)

            plt.text(x = 0, y = 12.2, s='Proto\n$q = 0$', ha='center', fontsize=12)
            plt.text(x = 1, y = 12.1, s='Deproto\n$q = -1$', ha='center', fontsize=12)

            plt.title('{}-{} (all chains) in {}\npH = {}, deprotonation = {:.2f}'.format(
                group.d_resname,
                group.d_resid,
                self.d_name,
                self.d_pH,
                self.deprotonation(x)))

            # Axes and stuff
            plt.axis([-0.1, 1.1, -0.1, 12])
            plt.xlabel(r"$\lambda$-coordinate")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.grid()

            # Add green vertical line indicating experimental value
            if self.d_pH == 4.0:
                plt.vlines(x=biophys["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=10, color='r', linewidth=4.0, label="biophysics.se/Prevost2012")
                plt.vlines(x=nury2010["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=8, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013")
                plt.vlines(x=fritsch2011["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=6, color='b', linewidth=4.0, label="Fritsch2011")
                plt.vlines(x=lev2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=4, color='c', linewidth=4.0, label="Lev2017")
                plt.vlines(x=nemecz2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018")
                # plt.vlines(x=ullman["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished)")
            plt.legend(loc='upper center')

            # Save and clear
            plt.savefig('{}/hist_{:03d}-{}.png'.format(self.d_dir, group.d_resid, group.d_resname))
            # plt.savefig('{}/hist_{:03d}-{}.pdf'.format(self.d_dir, group.d_resid, group.d_resname))
            plt.clf(); plt.close()

    def plotHIS_traj(self):
        for group in self.d_multiStateList:
            print('plotting {}...'.format(group.d_fname[0]), end='\r')

            plt.plot(group.d_t[0], group.d_x[0], linewidth=0.5, label='state 1 (double proto) ({})'.format(group.d_fname[0]))
            plt.plot(group.d_t[1], group.d_x[1], linewidth=0.5, label='state 2 (anti) ({})'.format(group.d_fname[1]))
            plt.plot(group.d_t[2], group.d_x[2], linewidth=0.5, label='state 3 (syn) ({})'.format(group.d_fname[2]))

            plt.title('{}-{} in chain {} in {}\n({}, {}, {}), pH={}'.format(
                group.d_resname,
                group.d_resid,
                group.d_chain,
                self.d_name,
                group.d_fname[0],
                group.d_fname[1],
                group.d_fname[2],
                self.d_pH))

            # Axes and stuff
            plt.ylim(-0.1, 1.1)
            plt.xlabel("Time (ps)")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.ylabel(r"$\lambda$-coordinate")
            plt.grid()
            plt.legend()

            plt.savefig('{}/traj_{}_{}_{}.png'.format(self.d_dir, group.d_chain, group.d_resid, group.d_resname))
            # plt.savefig('{}/traj_{}_{}_{}.pdf'.format(self.d_dir, group.d_chain, group.d_resid, group.d_resname))

            # Clear
            plt.clf(); plt.close()

    def plotHIS_histo(self):
        numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
        resPerChain = int(len(self.d_multiStateList)/numChains)

        for ii in range(0, resPerChain):
            # For convenient access
            group = self.d_multiStateList[ii]
            print("Plotting histogram for {}-{}".format(group.d_resname, group.d_resid))

            # Gather data
            x1 = []; x2 = []; x3 = []
            for jj in range(0, numChains):
                x1 += self.d_multiStateList[ii + resPerChain * jj].d_x[0]
                x2 += self.d_multiStateList[ii + resPerChain * jj].d_x[1]
                x3 += self.d_multiStateList[ii + resPerChain * jj].d_x[2]
                print(self.d_multiStateList[ii + resPerChain * jj].d_fname)

            # I - PLOT OLD HISTOGRAM
            ####################################################################

            plt.figure(figsize=(8, 6))
            plt.hist(x1, density=True, bins=200, label='state 1 (double proto.)', histtype='step')
            plt.hist(x2, density=True, bins=200, label='state 2 (anti)', histtype='step')
            plt.hist(x3, density=True, bins=200, label='state 3 (syn)', histtype='step')
            plt.legend()

            plt.title('{}-{} (all chains) in {}, pH = {}'.format(group.d_resname, group.d_resid, self.d_name, self.d_pH))

            plt.axis([-0.1, 1.1, -0.1, 12])
            plt.xlabel(r"$\lambda$-coordinate")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.grid()

            # Add green vertical line indicating experimental value
            if self.d_pH == 4.0:
                plt.vlines(x=biophys["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=10, color='r', linewidth=4.0, label="biophysics.se/Prevost2012")
                plt.vlines(x=nury2010["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=8, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013")
                plt.vlines(x=fritsch2011["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=6, color='b', linewidth=4.0, label="Fritsch2011")
                plt.vlines(x=lev2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=4, color='c', linewidth=4.0, label="Lev2017")
                plt.vlines(x=nemecz2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018")
                # plt.vlines(x=ullman["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished)")
            plt.legend(loc='upper center')

            # Save and clear
            plt.savefig('{}/hist_{:03d}-{}_old.png'.format(self.d_dir, group.d_resid, group.d_resname))
            # plt.savefig('{}/hist_{:03d}-{}_old.pdf'.format(self.d_dir, group.d_resid, group.d_resname))
            plt.clf(); plt.close()

            # II - PLOT NEW HISTOGRAM
            ####################################################################
            # https://stackoverflow.com/questions/8437788/how-to-correctly-generate-a-3d-histogram-using-numpy-or-matplotlib-built-in-func

            # Length needs to be equal
            if len(x1) > len(x2):
                xAmplitudes = x1[:len(x2)]
                yAmplitudes = x2
            elif len(x1) < len(x2):
                xAmplitudes = x1
                yAmplitudes = x2[:len(x1)]
            else:
                xAmplitudes = x1
                yAmplitudes = x2

            x = np.array(xAmplitudes)   #turn x,y data into numpy arrays
            y = np.array(yAmplitudes)

            fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
            ax = fig.add_subplot(111, projection='3d')

            #make histogram stuff - set bins - I choose 20x20 because I have a lot of data
            hist, xedges, yedges = np.histogram2d(x, y, bins=(100, 100), density=True)
            hist = hist.T # GODVERDOMME NOG AAN TOE, dit heeft me echt uren gekost... waarom fixt men dit niet
            xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])

            xpos = xpos.flatten()/2.
            ypos = ypos.flatten()/2.
            zpos = np.zeros_like (xpos)

            dx = xedges [1] - xedges [0]
            dy = yedges [1] - yedges [0]
            dz = hist.flatten()

            cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
            max_height = np.max(dz)   # get range of colorbars so we can normalize
            min_height = np.min(dz)
            # scale each z to [0,1], and get their rgb values
            rgba = [cmap((k-min_height)/max_height) for k in dz]

            ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
            plt.xlabel("$\lambda_1$ (0 = single, 1 = double proto)")
            plt.ylabel("$\lambda_2$ (0 = syn, 1 = anti)")

            ####################################################################

            plt.title('{}-{} (all chains) in {}, pH = {}'.format(group.d_resname, group.d_resid, self.d_name, self.d_pH))

            # Save and clear
            plt.savefig('{}/hist_{:03d}-{}.png'.format(self.d_dir, group.d_resid, group.d_resname))
            # plt.savefig('{}/hist_{:03d}-{}.pdf'.format(self.d_dir, group.d_resid, group.d_resname))
            plt.clf(); plt.close()

    def plotBuffer(self):
        group = self.d_buffer

        # Plot lambda coordinate
        plt.plot(group.d_t, group.d_x, linewidth=0.5, label='buffer $\lambda$-coordinate')

        # Convert lambda coordinate to charge
        charge = [(val - 0.5) for val in group.d_x]

        # Plot charge
        plt.plot(group.d_t, charge, linewidth=0.5, label='buffer charge')

        # Title
        plt.title('Buffer $\lambda$-coordinate and charge ({} particles)'.format(group.d_count))

        # Axes and stuff
        plt.ylim(-0.6, 1.1)
        plt.xlabel("Time (ps)")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
        plt.grid()
        plt.legend()

        # Save as .png
        plt.savefig('{}/buffer.png'.format(self.d_dir))
        # plt.savefig('{}/buffer.pdf'.format(self.d_dir))

        # Clear
        plt.clf(); plt.close()

    def checkConvergence(self):
        numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
        resPerChain = int(len(self.d_twoStateList)/numChains)
        chains = ['A', 'B', 'C', 'D', 'E']

        for ii in range(0, resPerChain):

            # For convenient access
            group = self.d_twoStateList[ii]

            # Make the plot
            for jj in range(0, numChains):
                obj = self.d_twoStateList[ii + resPerChain * jj]
                plt.plot(obj.d_t, self.__movingDeprotonation(obj.d_x), label='chain {}'.format(chains[jj]))

            plt.title('Moving deprotonation for {}-{}'.format(group.d_resname, group.d_resid))

            # Axes and stuff
            plt.ylim(-0.1, 1.1)
            plt.xlabel("Time (ps)")
            plt.ylabel(r"$\lambda$-coordinate")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.grid()
            plt.legend(loc='upper center')

            plt.savefig('{}/conv_{}_{}.png'.format(self.d_dir, group.d_resid, group.d_resname))
            # plt.savefig('{}/conv_{}_{}.pdf'.format(self.d_dir, group.d_resid, group.d_resname))

            # Clear
            plt.clf(); plt.close()

    def __movingDeprotonation(self, xList, cutoff=0.80):
        Av = len(xList) * [0]
        lambda_proto = 1
        lambda_deproto = 0

        for idx in range(0, len(xList)):
            if xList[idx] > cutoff:
                lambda_deproto += 1
            elif xList[idx] < 1 - cutoff:
                lambda_proto += 1

            Av[idx] = float(lambda_deproto) / (lambda_proto + lambda_deproto)

        return Av

    def deprotonation(self, xList, cutoff=0.80):

        lambda_proto   = 0
        lambda_deproto = 0

        for x in xList:
            if x > cutoff:
                lambda_deproto += 1
            if x < 1 - cutoff:
                lambda_proto   += 1

        if lambda_proto + lambda_deproto == 0:
            fraction = 0
        else:
            fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

        return fraction
