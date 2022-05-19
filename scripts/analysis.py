#!/bin/python3

# IMPORT MODULES
import os
import numpy as np
import matplotlib.pyplot as plt
from phbuilder.structure import Structure
from data import biophys, nury2010, fritsch2011, lev2017, nemecz2017

def loadxvg(fname, col=[0, 1], dt=1, b=0):
    count = -1
    data = [ [] for _ in range(len(col)) ]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue

        count += 1
        if count % dt != 0:
            continue

        listLine = stringLine.split()

        if b != 0 and float(listLine[col[0]]) < b:
            continue

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
        self.d_idx      = [idx, idx+1, idx+2]
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

class Replica:
    '''Holds data for one replica'''
    def __init__(self, name, replicaID, twoStateList, multiStateList, buffer):
        self.d_name           = name
        self.d_replicaID      = replicaID
        self.d_twoStateList   = twoStateList
        self.d_multiStateList = multiStateList
        self.d_buffer         = buffer

class ReplicaSet:
    '''Holds the data for one set of GLIC replicas (one structure)'''
    def __init__(self, name, replicaList):
        self.d_name        = name
        self.d_replicaList = replicaList

class GLICSims:
    '''Holds the data for all four GLIC replica sets'''
    def __init__(self, folders, replicas, dt, b):

        self.d_replicaSetList = []
        for ii in folders:
            print(ii)
            os.chdir(ii)

            replicaList = []
            for jj in replicas:
                print(jj)
                os.chdir(jj)

                idx = 1
                foundBUF = False
                twoStateList   = []
                multiStateList = []
                for residue in Structure('CA.pdb', 2).d_residues:

                    if residue.d_resname in ['ASPT', 'GLUT']:
                        print('Loading {}-{} in chain {}...'.format(residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                        xvgdata = loadxvg('cphmd-coord-{}.xvg'.format(idx), dt=dt, b=b)

                        twoStateList.append(TwoState(
                            idx, 
                            residue.d_resname, 
                            residue.d_resid,
                            residue.d_chain,
                            xvgdata[0], # time
                            xvgdata[1]  # coordinates
                            ))

                        idx += 1

                    elif residue.d_resname == 'HSPT':
                        print('Loading {}-{} in chain {}...'.format(residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                        xvgdata1 = loadxvg('cphmd-coord-{}.xvg'.format(idx  ), dt=dt, b=b)
                        xvgdata2 = loadxvg('cphmd-coord-{}.xvg'.format(idx+1), dt=dt, b=b)
                        xvgdata3 = loadxvg('cphmd-coord-{}.xvg'.format(idx+2), dt=dt, b=b)

                        multiStateList.append(MultiState(
                            idx,
                            residue.d_resname,
                            residue.d_resid,
                            residue.d_chain,
                            xvgdata1[0],  # file 1 time
                            xvgdata2[0],  # file 2 time
                            xvgdata3[0],  # file 3 time
                            xvgdata1[1],  # file 1 coordinates
                            xvgdata2[1],  # file 2 coordinates
                            xvgdata3[1])) # file 3 coordinates

                        idx += 3

                    elif residue.d_resname == 'BUF' and not foundBUF:
                        xvgdata  = loadxvg('cphmd-coord-{}.xvg'.format(idx), dt=dt)
                        buffer   = Buffer(idx, xvgdata[0], xvgdata[1], count=185)
                        foundBUF = True

                replicaList.append(Replica(ii, int(jj), twoStateList, multiStateList, buffer))

                os.chdir('..')
            
            self.d_replicaSetList.append(ReplicaSet(ii, replicaList))

            os.chdir('..')


































# # PARAMETERS
# # ns_to_drop = 300   # ns (Correct voor als ik zo overkopieer)
# ns_to_drop = 0     # ns
# window     = 500   # frames to consider for moving average (= 50 ns)

# # FUNCTION FOR MOVING DEPROTONATION
# def movingDeprotonation(tList, xList):
#     Sum = sum(xList[0:window]) # Initialize

#     t = tList[window:]
#     x = len(range(window, len(xList))) * [0]

#     for ii in range(window, len(xList)):
#         x[ii - window] = Sum / float(window)
#         Sum -= xList[ii - window]
#         Sum += xList[ii]

#     return t, x

# # CREATE LAMBDAPLOTS DIR IN MAIN DIR IF NOT ALREADY EXISTING
# if not os.path.isdir('lambdaplots'):
#     os.mkdir('lambdaplots')

# # MAIN LOOP
# for dir in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:

#     os.chdir(dir)

#     # SOME SETTINGS
#     startingPDB = 'CA.pdb'
#     name        = dir
#     pH          = float(dir[-1])

#     # LOAD REPLICAS
#     replicas = []

#     os.chdir('01')
#     replicas.append(analyze2.Analysis(startingPDB, name, pH, 500, ns_to_drop))

#     os.chdir('../02')
#     replicas.append(analyze2.Analysis(startingPDB, name, pH, 500, ns_to_drop))

#     os.chdir('../03')
#     replicas.append(analyze2.Analysis(startingPDB, name, pH, 500, ns_to_drop))

#     os.chdir('../04')
#     replicas.append(analyze2.Analysis(startingPDB, name, pH, 500, ns_to_drop))

#     os.chdir('../..')


#     # MAKE HISTOGRAMS
#     numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
#     resPerChain = int(len(replicas[0].d_twoStateList)/numChains)

#     for ii in range(0, resPerChain):

#         valuesList = []
#         for idx in range(0, len(replicas)):
#             rep = replicas[idx]

#             # FOR CONVENIENT ACCESS
#             group = rep.d_twoStateList[ii]

#             # GATHER DATA
#             for jj in range(0, numChains):
#                 # Load object for a single chain
#                 obj = rep.d_twoStateList[ii + resPerChain * jj]
#                 x = [1.0 - val for val in obj.d_x] # Mirror in vertical x=0.5 axis

#                 # GET HISTOGRAM VALUES, BINS (BUT DO NOT ACTUALLY PLOT)
#                 values, bins = np.histogram(x, density=True, bins=200, range=(-0.1, 1.1))

#                 valuesList.append(values)

#         # COMPUTE MEAN AND STANDARD ERROR
#         meanList  = len(values) * [0] # 200, to hold mean for each bin
#         errorList = len(values) * [0] # 200, to hold error for each bin

#         for ii in range(0, len(values)): # 200

#             # Create list of 20 values
#             temp = [0] * len(valuesList) # 4*5=20
#             for jj in range(0, len(valuesList)): # 4*5=20
#                 temp[jj] = valuesList[jj][ii]

#             meanList[ii] = np.mean(temp)
#             errorList[ii] = np.std(temp)

#         # PLOT MEAN AND SHADED REGION (ERROR)
#         A = []; B = []
#         for idx in range(0, len(meanList)):
#             A.append(meanList[idx] + errorList[idx])
#             B.append(meanList[idx] - errorList[idx])

#         plt.figure(figsize=(8, 6))
#         plt.plot(bins[1:], meanList)
#         # plt.fill_between(bins[1:], A, B, alpha=0.4, edgecolor='#CC4F1B', facecolor='#FF9848')
#         plt.fill_between(bins[1:], A, B, alpha=0.4, color='#1f77b4')

#         # MAKE PLOT MORE NICE
#         plt.text(x = 1, y = 12.2, s='Proto\n$q = 0$', ha='center', fontsize=12)
#         plt.text(x = 0, y = 12.1, s='Deproto\n$q = -1$', ha='center', fontsize=12)
#         plt.title(rep.d_name, fontsize=18)
#         plt.axis([-0.1, 1.1, -0.1, 12])
#         plt.xlabel(r"$\lambda$-coordinate")
#         plt.xticks(ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0]) # because we mirror in vertical x=0.5 axis
#         plt.grid()

#         # ADD EXPERIMENTAL DATA FOR PH=4 CASE
#         if pH == 4.0:
#             plt.vlines(x=biophys["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=10, color='r', linewidth=4.0, label="biophysics.se/Prevost2012")
#             plt.vlines(x=nury2010["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=8, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013")
#             plt.vlines(x=fritsch2011["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=6, color='b', linewidth=4.0, label="Fritsch2011")
#             plt.vlines(x=lev2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=4, color='c', linewidth=4.0, label="Lev2017")
#             plt.vlines(x=nemecz2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018")
#             plt.legend(loc='upper center')

#         # SAVE AND CLEAR
#         plt.tight_layout()
#         plt.savefig('lambdaplots/{}_{:03d}-{}.png'.format(replicas[0].d_name, group.d_resid, group.d_resname))
#         plt.clf(); plt.close()


#     # MAKE CONVERGENCE PLOTS
#     numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
#     resPerChain = int(len(replicas[0].d_twoStateList)/numChains)
#     chains = ['A', 'B', 'C', 'D', 'E']

#     for ii in range(0, resPerChain):

#         for idx in range(0, len(replicas)):
#             rep = replicas[idx]

#             # FOR CONVENIENT ACCESS
#             group = rep.d_twoStateList[ii]

#             # MAKE THE PLOT
#             for jj in range(0, numChains):
#                 # Load object for a single chain
#                 obj = rep.d_twoStateList[ii + resPerChain * jj]
#                 t   = obj.d_t
#                 x   = [1.0 - val for val in obj.d_x] # Mirror in vertical x=0.5 axis

#                 a, b = movingDeprotonation(t, x)
#                 plt.plot(a, b, label='rep {} chain {}'.format(idx + 1, chains[jj]))

#         # MAKE PLOT MORE NICE
#         plt.title(rep.d_name, fontsize=18)
#         plt.ylim(-0.1, 1.1)
#         plt.xlabel("Time (ps)")
#         plt.ylabel("Protonation running average")
#         plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
#         plt.grid()

#         # SAVE AND CLEAR
#         plt.tight_layout()
#         plt.savefig('lambdaplots/{}_{:03d}-{}_conv.png'.format(replicas[0].d_name, group.d_resid, group.d_resname))
#         plt.clf(); plt.close()


#     # DO HISTIDINE ANALYSIS
#     numChains   = 5 # number of numChains five-fold symmetric protein, hardcoded.
#     resPerChain = int(len(replicas[0].d_multiStateList)/numChains)

#     for ii in range(0, resPerChain): # Loop over residues

#         plt.figure(figsize=(8, 6))

#         for qq in [0, 1, 2]: # Loop over lambda groups
#             valuesList = []

#             for idx in range(0, len(replicas)): # Loop over replicas
#                 rep = replicas[idx]

#                 # FOR CONVENIENT ACCESS
#                 group = rep.d_multiStateList[ii]

#                 if qq == 0: # user update
#                     print("Plotting histogram for {}-{}".format(group.d_resname, group.d_resid))

#                 # GATHER DATA
#                 for jj in range(0, numChains): # Loop over chains
#                     # Load object for a single chain
#                     obj = rep.d_multiStateList[ii + resPerChain * jj]
#                     x = [1.0 - val for val in obj.d_x[qq]] # Mirror in vertical x=0.5 axis

#                     if qq == 0: # user update
#                         print(obj.d_fname) 

#                     # GET HISTOGRAM VALUES, BINS
#                     values, bins = np.histogram(x, density=True, bins=200, range=(-0.1, 1.1))
#                     valuesList.append(values)

#             # COMPUTE MEAN AND STANDARD ERROR
#             meanList  = len(values) * [0] # 200, to hold mean for each bin
#             errorList = len(values) * [0] # 200, to hold error for each bin

#             for kk in range(0, len(values)): # 200

#                 # Create list of 20 values
#                 temp = [0] * len(valuesList) # 4*5=20
#                 for pp in range(0, len(valuesList)): # 4*5=20
#                     temp[pp] = valuesList[pp][kk]

#                 meanList[kk] = np.mean(temp)
#                 errorList[kk] = np.std(temp)

#             # PLOT MEAN AND SHADED REGION (ERROR)
#             A = []; B = []
#             for mm in range(0, len(meanList)):
#                 A.append(meanList[mm] + errorList[mm])
#                 B.append(meanList[mm] - errorList[mm])

#             description = ['state 1 (double proto)', 'state 2 (anti)', 'state 3 (syn)']
#             color = ['#1f77b4', '#ff7f0e', '#2ca02c']

#             plt.plot(bins[1:], meanList, color=color[qq], label=description[qq])
#             plt.fill_between(bins[1:], A, B, alpha=0.4, color=color[qq])

#         # MAKE PLOT MORE NICE
#         plt.title(rep.d_name, fontsize=18)
#         plt.axis([-0.1, 1.1, -0.1, 12])
#         plt.xlabel(r"$\lambda$-coordinate")
#         plt.xticks(ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0]) # because we mirror in vertical x=0.5 axis
#         plt.grid()
#         plt.legend(loc='upper center')

#         # ADD EXPERIMENTAL DATA FOR PH=4 CASE
#         # if pH == 4.0:
#         #     plt.vlines(x=biophys["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=10, color='r', linewidth=4.0)
#         #     plt.vlines(x=nury2010["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=8, color='g', linewidth=4.0)
#         #     plt.vlines(x=fritsch2011["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=6, color='b', linewidth=4.0)
#         #     plt.vlines(x=lev2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=4, color='c', linewidth=4.0)
#         #     plt.vlines(x=nemecz2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color = 'm', linewidth=4.0)

#         # SAVE AND CLEAR
#         plt.tight_layout()
#         plt.savefig('lambdaplots/{}_{:03d}-{}.png'.format(replicas[0].d_name, group.d_resid, group.d_resname))
#         plt.clf(); plt.close()





# os.chdir('lambdaplots')
# # MERGE TO CREATE FINAL HISTOGRAM PLOTS
# for res in ['127-HSPT', '235-HSPT', '277-HSPT']:
#     os.system('convert 6ZGD_7_{}.png 4HFI_7_{}.png +append temp1.png'.format(res, res))
#     os.system('convert 6ZGD_4_{}.png 4HFI_4_{}.png +append temp2.png'.format(res, res))
#     os.system('convert temp1.png temp2.png -append hist_{}.png'.format(res))

# for res in ['013-ASPT', '014-GLUT', '026-GLUT', '031-ASPT', '032-ASPT', '035-GLUT', '049-ASPT', '055-ASPT', '067-GLUT', '069-GLUT', '075-GLUT', '082-GLUT', '086-ASPT', '088-ASPT', '091-ASPT', '097-ASPT', '104-GLUT', '115-ASPT', '122-ASPT', '136-ASPT', '145-ASPT', '147-GLUT', '153-ASPT', '154-ASPT', '161-ASPT', '163-GLUT', '177-GLUT', '178-ASPT', '181-GLUT', '185-ASPT', '222-GLUT', '243-GLUT', '272-GLUT', '282-GLUT']:
#     os.system('convert 6ZGD_7_{}.png 4HFI_7_{}.png +append temp1.png'.format(res, res))
#     os.system('convert 6ZGD_4_{}.png 4HFI_4_{}.png +append temp2.png'.format(res, res))
#     os.system('convert temp1.png temp2.png -append hist_{}.png'.format(res))

# # MERGE TO CREATE FINAL CONVERGENCE PLOTS
# for res in ['013-ASPT', '014-GLUT', '026-GLUT', '031-ASPT', '032-ASPT', '035-GLUT', '049-ASPT', '055-ASPT', '067-GLUT', '069-GLUT', '075-GLUT', '082-GLUT', '086-ASPT', '088-ASPT', '091-ASPT', '097-ASPT', '104-GLUT', '115-ASPT', '122-ASPT', '136-ASPT', '145-ASPT', '147-GLUT', '153-ASPT', '154-ASPT', '161-ASPT', '163-GLUT', '177-GLUT', '178-ASPT', '181-GLUT', '185-ASPT', '222-GLUT', '243-GLUT', '272-GLUT', '282-GLUT']:
#     os.system('convert 6ZGD_7_{}_conv.png 4HFI_7_{}_conv.png +append temp1.png'.format(res, res))
#     os.system('convert 6ZGD_4_{}_conv.png 4HFI_4_{}_conv.png +append temp2.png'.format(res, res))
#     os.system('convert temp1.png temp2.png -append conv_{}.png'.format(res))

# # MERGE TO CREATE FINAL PLOTS
# for res in ['013-ASPT', '014-GLUT', '026-GLUT', '031-ASPT', '032-ASPT', '035-GLUT', '049-ASPT', '055-ASPT', '067-GLUT', '069-GLUT', '075-GLUT', '082-GLUT', '086-ASPT', '088-ASPT', '091-ASPT', '097-ASPT', '104-GLUT', '115-ASPT', '122-ASPT', '136-ASPT', '145-ASPT', '147-GLUT', '153-ASPT', '154-ASPT', '161-ASPT', '163-GLUT', '177-GLUT', '178-ASPT', '181-GLUT', '185-ASPT', '222-GLUT', '243-GLUT', '272-GLUT', '282-GLUT']:
#     os.system('convert hist_{}.png conv_{}.png +append final_{}.png'.format(res, res, res))
