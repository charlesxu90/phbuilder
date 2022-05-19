#!/bin/python3

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from phbuilder.structure import Structure

def loadxvg(fname, col=[0, 1]):
    data = [ [] for _ in range(len(col)) ]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        listLine = stringLine.split()
        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data

def loadCol(fileName, col, start=0, stop=0):
    data = []

    try:
        for x, y in enumerate(open(fileName)):
            if start == 0 and stop == 0 and (y.split()[0][0] not in ['@','#','&']):
                data.append(float(y.split()[col-1]))

            elif (x >= start-1) and (x <= stop-1):
                data.append(float(y.split()[col-1]))

    except IndexError:
        pass

    return data

def titrate(lambdaFileName, cutoff=0.80):
    lambda_proto   = 0
    lambda_deproto = 0

    for x in loadCol(lambdaFileName, 2):
        if (x > cutoff):
            lambda_deproto += 1
        if (x < 1 - cutoff):
            lambda_proto   += 1

    if (lambda_proto + lambda_deproto == 0):
        fraction = 0
    else:
        fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return fraction

def compareLambdaFiles(namelist):
    # If you accidentally put a string instead of a list, fix it.
    if (type(namelist) == type("")):
        namelist = [namelist]

    # Define (size of) main figure.
    fig = plt.figure(figsize=(24, 10))

    # Define sub-plots.
    plt1 = fig.add_subplot(2, 4, 1)
    plt2 = fig.add_subplot(2, 4, 2)
    plt3 = fig.add_subplot(2, 4, 3)
    plt4 = fig.add_subplot(2, 4, 4)
    plt5 = fig.add_subplot(2, 4, 5)
    plt6 = fig.add_subplot(2, 4, 6)
    plt7 = fig.add_subplot(2, 4, 7)
    plt8 = fig.add_subplot(2, 4, 8)

    # Get the data and plot.
    for name in namelist:
        time        = loadCol(name, 1)
        lambda_x    = loadCol(name, 2)
        lambda_dvdl = loadCol(name, 3)
        lambda_temp = loadCol(name, 4)
        lambda_vel  = loadCol(name, 5)
        F_coulomb   = loadCol(name, 6)
        F_corr      = loadCol(name, 7)
        F_bias      = loadCol(name, 8)
        F_ph        = loadCol(name, 9)

        plt1.plot(time, lambda_x, linewidth=0.5, label="deprotonation = {:.2f}".format(titrate(name)))
        plt2.plot(time, lambda_temp, linewidth=0.5, label="mean = {:.1f} (K)".format(sum(lambda_temp)/len(lambda_temp)))
        plt3.hist(lambda_vel, density=True)
        plt4.scatter(lambda_x, lambda_dvdl, s=5)
        plt5.scatter(lambda_x, F_coulomb, s=5)
        plt6.scatter(lambda_x, F_corr, s=5)
        plt7.scatter(lambda_x, F_bias, s=5)
        plt8.scatter(lambda_x, F_ph, s=5)

    plt1.set_title("$\lambda$-coordinate vs time")
    plt1.set_xlabel("Time (ps)")
    plt1.set_ylabel("$\lambda$-coordinate")
    plt1.set_ylim(-0.1, 1.1)
    plt1.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt1.legend()

    plt2.set_title("$\lambda$-temperature vs time")
    plt2.set_xlabel("Time (ps)")
    plt2.set_ylabel("$\lambda$-temperature (K)")
    plt2.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    plt2.legend()

    plt3.set_title("$\lambda$-velocity distribution")
    plt3.set_xlabel("$\lambda$-velocity (km/s)")

    plt4.set_title("Force (dV/dl) on $\lambda$-particle")
    plt4.set_xlabel("$\lambda$-coordinate")
    plt4.set_ylabel("dV/dl")
    plt4.set_xlim(-0.1, 1.1)

    plt5.set_title("Coulomb-force on $\lambda$-particle")
    plt5.set_xlabel("$\lambda$-coordinate")
    plt5.set_ylabel("$F_{Coulomb}$")
    plt5.set_xlim(-0.1, 1.1)

    plt6.set_title("Reference-force on $\lambda$-particle")
    plt6.set_xlabel("$\lambda$-coordinate")
    plt6.set_ylabel("$F_{corr}$")
    plt6.set_xlim(-0.1, 1.1)

    plt7.set_title("Bias-force on $\lambda$-particle")
    plt7.set_xlabel("$\lambda$-coordinate")
    plt7.set_ylabel("$F_{bias}$")
    plt7.axis([-0.1, 1.1, -200, 200])

    plt8.set_title("pH-force on $\lambda$-particle")
    plt8.set_xlabel("$\lambda$-coordinate")
    plt8.set_ylabel("$F_{pH}$")
    plt8.set_xlim(-0.1, 1.1)

    # Stuff we do in all subplots we can do in a loop:
    for plot in [plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8]:
        plot.grid()

    fig.legend(namelist, loc="upper center")
    # plt.tight_layout()
    plt.show()

def plotlambda(fileName, plotHSPT=True, plotBUF=False):
    pdb = Structure(fileName, verbosity=3)

    idx = 1
    for residue in pdb.d_residues:
        if residue.d_resname in ['ASPT', 'GLUT']:
            t = loadCol('lambda_{0}.dat'.format(idx), 1)
            x = loadCol("lambda_{0}.dat".format(idx), 2)

            plt.plot(t, x, label='{}-{}'.format(residue.d_resname, residue.d_resid), linewidth=1)
            idx += 1

        elif residue.d_resname == 'HSPT':
            if plotHSPT:
                for val in range(0, 3):
                    t = loadCol('lambda_{0}.dat'.format(idx + val), 1)
                    x = loadCol("lambda_{0}.dat".format(idx + val), 2)

                    plt.plot(t, x, label='{}-{}-{}'.format(residue.d_resname, residue.d_resid, val+1), linewidth=1)

            idx += 3

    if plotBUF:
        t = loadCol("lambda_{0}.dat".format(idx), 1)
        x = loadCol("lambda_{0}.dat".format(idx), 2)

        plt.plot(t, x, label="Buffer", linewidth=1)

    plt.xlabel("Time (ps)")
    plt.ylabel(r"$\lambda$-coordinate")
    plt.ylim(-0.1, 1.1)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))

    plt.legend()
    plt.grid()
    plt.show()

def glicphstates(fileName, pdbName, pH, nstOut, dump=0):
    # EXPERIMENTAL DATA ON PROTONATION STATES AT VARIOUS PH ####################
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

    # DIRECTORY STRUCTURE

    dirname = "lambdaplots"
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        os.system("rm -f {0}/*.png {0}/*.pdf".format(dirname))

    # CREATE LAMBDA PLOT FOR EVERY INDIVIDUAL PROTONATABLE RESIDUE

    pdb = Structure(fileName, verbosity=3)

    print('Writing individual lambda plots...')

    lines = int((1000 * dump) / (0.002 * nstOut))
    print("Will skip first {} ns (= {} lines)".format(dump, lines))

    idx = 1
    # Loop through the residues
    for residue in pdb.d_residues:
        # If it's an ASPT or GLUT...
        if residue.d_resname in ['ASPT', 'GLUT']:
            # User update
            print('processing lambda_{}.dat...'.format(idx), end='\r')
            
            # Load the relevant columns
            t = loadCol('lambda_{0}.dat'.format(idx), 1)[lines:]
            x = loadCol("lambda_{0}.dat".format(idx), 2)[lines:]

            # Create the actual plot
            plt.plot(t, x, linewidth=0.5)

            # Title
            plt.title("{0}-{1} in chain {2} in {3}\npH={4}, nstlambda={5}, deprotonation={6:.2f}".format(
                residue.d_resname,
                residue.d_resid,
                residue.d_chain,
                pdbName,
                pH,
                nstOut,
                titrate("lambda_{}.dat".format(idx))
            ))

            # Axes and stuff
            plt.ylim(-0.1, 1.1)
            plt.xlabel("Time (ps)")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.ylabel(r"$\lambda$-coordinate")
            plt.grid()

            # Save as .png
            plt.savefig("{}/{}_{}-{:03d}.png".format(dirname, residue.d_chain, residue.d_resname, residue.d_resid))

            # Clear: clf = clear the entire current figure. close = closes a window
            plt.clf(); plt.close()

            # Increment the lambda_xxx.dat number.
            idx += 1

        elif residue.d_resname == 'HSPT':
            # We want 3 different lines in the plot, so loop over 3:
            for val in range(0, 3):
                # User update
                print('processing lambda_{}.dat...'.format(idx + val), end='\r')

                # Load the relevant columns
                t = loadCol('lambda_{0}.dat'.format(idx + val), 1)[lines:]
                x = loadCol("lambda_{0}.dat".format(idx + val), 2)[lines:]

                # Create the actual plot
                plt.plot(t, x, linewidth=0.5, label='state {}'.format(val + 1))

            # Title
            plt.title("{}-{} in chain {} in {}\npH={}, nstlambda={}, deprotonation={:.2f}".format(
                residue.d_resname,
                residue.d_resid,
                residue.d_chain,
                pdbName,
                pH,
                nstOut,
                titrate("lambda_{}.dat".format(idx))
            ))

            # Axes and stuff
            plt.ylim(-0.1, 1.1)
            plt.xlabel("Time (ps)")
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
            plt.ylabel(r"$\lambda$-coordinate")
            plt.grid()
            plt.legend()

            # Save as .png
            plt.savefig("{}/{}_{}-{:03d}.png".format(dirname, residue.d_chain, residue.d_resname, residue.d_resid))

            # Clear: clf = clear the entire current figure. close = closes a window
            plt.clf(); plt.close()

            idx += 3

    print('Finished writing individual plots')

    # CREATE HISTOGRAM PLOTS FOR COMBINED PROTO STATE OF ALL FIVE CHAINS

    # GATHER DATA

    numChains = 5
    resPerChain = 43    # 215 / 5

    print('Loading data for writing the histograms...')

    dataList = []
    for ii in range(1, resPerChain + 1):
        data = []
        for jj in range(0, numChains):
            data += loadCol('lambda_{}.dat'.format(ii + resPerChain * jj), 2)[lines:]
        dataList.append(data)

    # PERFORM HISTOGRAM PLOTTING

    idx = 0
    for residue in pdb.d_residues:
        if residue.d_chain == 'A':

            if residue.d_resname in ['ASPT', 'GLUT']:
                print('Plotting {}/{}'.format(idx + 1, resPerChain), end='\r')

                plt.figure(figsize=(8, 6))
                plt.hist(dataList[idx], density=True, bins=200)

                # Title
                plt.title("{}-{} (all chains) in {}\npH={}, nstlambda={}".format(
                    residue.d_resname,
                    residue.d_resid,
                    pdbName,
                    pH,
                    nstOut,
                    ))

                # Axes and stuff
                plt.axis([-0.1, 1.1, -0.1, 12])
                plt.xlabel(r"$\lambda$-coordinate")
                plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
                plt.grid()

                # Add green vertical line indicating experimental value
                if pH == 4.0:
                    plt.vlines(x=biophys["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=12, color='r', linewidth=4.0, label="biophysics.se/Prevost2012 = {}".format(biophys["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=nury2010["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=10, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013 = {}".format(nury2010["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=fritsch2011["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=8, color='b', linewidth=4.0, label="Fritsch2011 = {}".format(fritsch2011["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=lev2017["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=6, color='c', linewidth=4.0, label="Lev2017 = {}".format(lev2017["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=nemecz2017["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=4, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018 = {}".format(nemecz2017["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=ullman["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished) = {}".format(ullman["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                plt.legend()
                
                # Save and clear
                plt.savefig('{}/hist_{}-{:03d}.png'.format(dirname, residue.d_resname, residue.d_resid))
                plt.clf(); plt.close()

                idx += 1

            elif residue.d_resname == 'HSPT':
                plt.figure(figsize=(8, 6))

                # We want 3 different lines in the plot, so loop over 3:                
                for val in range(0, 3):
                    # User update
                    print('Plotting {}/{}'.format(idx + val + 1, resPerChain), end='\r')

                    # Create the acutal plot                    
                    plt.hist(dataList[idx + val], density=True, bins=200, label='state {}'.format(val + 1), histtype='step')

                # Title
                plt.title("{0}-{1} (all chains) in {2}\npH={3}, nstlambda={4}".format(
                    residue.d_resname,
                    residue.d_resid,
                    pdbName,
                    pH,
                    nstOut,
                    ))

                # Axes and stuff
                plt.axis([-0.1, 1.1, -0.1, 12])
                plt.xlabel(r"$\lambda$-coordinate")
                plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
                plt.grid()

                # Add green vertical line indicating experimental value
                if pH == 4.0:
                    plt.vlines(x=biophys["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=12, color='r', linewidth=4.0, label="biophysics.se/Prevost2012 = {}".format(biophys["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=nury2010["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=10, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013 = {}".format(nury2010["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=fritsch2011["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=8, color='b', linewidth=4.0, label="Fritsch2011 = {}".format(fritsch2011["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=lev2017["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=6, color='c', linewidth=4.0, label="Lev2017 = {}".format(lev2017["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=nemecz2017["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=4, color = 'm', linewidth=4.0, label="Nemecz2017/Hu2018 = {}".format(nemecz2017["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                    plt.vlines(x=ullman["{0}-{1}".format(residue.d_resname, residue.d_resid)], ymin=0, ymax=2, color='y', linewidth=4.0, label="Ullman (unpublished) = {}".format(ullman["{0}-{1}".format(residue.d_resname, residue.d_resid)]))
                plt.legend()
                
                # Save and clear
                plt.savefig('{}/hist_{}-{:03d}.png'.format(dirname, residue.d_resname, residue.d_resid))
                plt.clf(); plt.close()

                idx += 3

    print('Finished writing histograms')

def inverseboltzmann(fileName, resid):

    def barrier(l):     
        k = 4.7431      # Coded this for nothing: this is only at 7.5!
        a = 0.0435
        b = 0.0027
        d = 3.75
        s = 0.30
        w = 1000.0
        r = 13.5
        m = 0.2019
    
        A = np.exp(- (l - 1 - b)**2 / (2 * a**2) )
        B = np.exp(- (l + b)**2     / (2 * a**2) )
        C = np.exp(- (l - 0.5)**2   / (2 * s**2) )
        D = 1 - math.erf(r * (l + m))
        E = 1 + math.erf(r * (l - 1 - m))

        return -k * (A + B) + d * C + 0.5 * w * (D + E)

    numChains = 5
    resPerChain = 43

    # GET THE DATA

    pdb = Structure(fileName, verbosity=3)

    idx = 1
    data = []
    for residue in pdb.d_residues:
        if residue.d_chain == 'A':
            if residue.d_resid == resid:
                for jj in range(0, numChains):
                    print("loading lambda_{}.dat...".format(idx + resPerChain * jj))
                    data += loadCol('lambda_{}.dat'.format(idx + resPerChain * jj), 2)[50000:]
                break

            if residue.d_resname in ['ASPT', 'GLUT']:
                idx += 1

            elif residue.d_resname == 'HSPT':
                idx += 3

    # GET THE HISTOGRAM AND CORRESPONDING LAMBDA LISTS
    
    bins1 = [x / 1000. for x in range(-100, 1100, 3)] # 400 bins
    hist, bins2 = np.histogram(data, bins=bins1, density=True)
    bins2 = bins2[1:] # remove first element to make arrays equal size

    # DO THE ACTUAL BOLTZMANN INVERSION

    R = 8.3145  # gas constant (not kb because we have kJ/mol as a unit)
    T = 300     # K

    energyList = []
    for p in hist:
        energyList.append(R * T * -np.log(p))

    energyList = [E / 1000. for E in energyList] # From J to kJ.

    # PLOT AND EYEBALL

    plt.plot(bins2, energyList)
    plt.grid()
    plt.ylabel("Energy (kJ/mol)")
    plt.xlabel(r"$\lambda$-coordinate")
    plt.show()
