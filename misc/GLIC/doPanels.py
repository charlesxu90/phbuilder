#!/bin/python3
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import MDAnalysis.analysis.rms
import os
import subprocess
import pathos.multiprocessing as mp
import pandas
import numpy as np
import scipy.stats as stats
import copy

# Set global font size for figures
matplotlib.rcParams.update({'font.size': 14})


def gromacs(command, stdin=[], basepath='/usr/local/gromacs_constantph'):
    """Python function for handeling calls to GROMACS.

    Args:
        command (string): GROMACS command, e.g. 'make_ndx -f protein.pdb'
        stdin (list, optional): List of input arguments. Defaults to [].
        basepath (str, optional): Base path for GROMACS version to be used. Defaults to 'usr/local/gromacs_constantph'.
    """

    # If we don't pass envvars to subprocess (which happens by default) this will work.
    path_to_gmx = os.path.normpath(basepath + '/' + 'bin/gmx')
    command = "{} {}".format(path_to_gmx, command)

    if stdin:
        xstr = ' << EOF\n'
        for val in stdin:
            xstr += '{}\n'.format(val)
        command += xstr + 'EOF'

    process = subprocess.run(command, shell=True, env={})

    if process.returncode != 0:
        print("Failed to run \"{}\" (exitcode {}).".format(command, process.returncode))
        return 1

    return 0


def loadxvg(fname, col=[0, 1], dt=1, b=0):
    """Loads an .xvg file into a list of lists.
    May also be used to load float columns from files in general.

    Args:
        fname (string): file name.
        col (list, optional): Columns to load. Defaults to [0, 1].
        dt (int, optional): Step size. Defaults to 1.
        b (int, optional): Starting point. Defaults to 0.

    Returns:
        list of lists : contains the columns that were loaded.
    """

    count = -1
    data = [[] for _ in range(len(col))]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        # THIS IS FOR THE dt PART.
        count += 1
        if count % dt != 0:
            continue

        listLine = stringLine.split()
        # AND THIS IS FOR THE b PART.
        if b != 0 and float(listLine[col[0]]) < b:
            continue

        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data


def getLambdaFileIndices(structure, resid):
    """Returns an array containing the lambda-file indices for the specified resid.

    Args:
        structure (string): pdb file name.
        resid (int): residue id.

    Returns:
       list: List of lambda indices.
    """
    u                  = MDAnalysis.Universe(structure)
    numChains          = len(u.segments) - 1
    segmentAatoms      = u.segments[0].atoms
    titratableAtoms    = segmentAatoms.select_atoms('resname ASPT GLUT HSPT')
    titratableResnames = list(titratableAtoms.residues.resnames)
    titratableResids   = list(titratableAtoms.residues.resids)
    targetidx          = titratableResids.index(resid)

    numASPTGLUT        = len(segmentAatoms.select_atoms('resname ASPT GLUT').residues)
    numHSPT            = len(segmentAatoms.select_atoms('resname HSPT').residues)
    factor             = numASPTGLUT + 3 * numHSPT

    count = 1
    for idx in range(0, len(titratableResnames)):

        if idx == targetidx:
            array = []
            for ii in range(0, numChains):
                array.append(count + ii * factor)
            return array

        if titratableResnames[idx] in ['ASPT', 'GLUT']:
            count += 1

        elif titratableResnames[idx] == 'HSPT':
            count += 3


def notExists(fname):
    """Returns True if the file 'panels/fname' does not yet exist.

    Args:
        fname (string): file name.

    Returns:
        bool: boolean.
    """

    path = "panels/{}".format(fname)

    if os.path.exists(path):
        print('{} already exists, not creating it again...'.format(path))
        return False

    return True


def ttestPass(sample1, sample2, alpha=0.05):
    """Returns True if the means of sample1 and sample2 differ SIGNIFICANTLY.
    That is, with a confidence interval of 1 - alpha %. Uses Welch's t-test.

    Args:
        sample1 (list): Some sample.
        sample2 (list): Sample to compare to.
        alpha (float, optional): Significance of the test. Defaults to 0.05.

    Returns:
        bool: Whether or not the two sample means differ signifcantly.
    """

    pvalue = stats.ttest_ind(sample1, sample2, equal_var=False)[1]

    return bool(pvalue < alpha)


class PanelBuilder:
    """Combines multiple analyses and plotting functions and attempts to create
    an entire panel for a residue at once.
    """

    def __init__(self, target, resids, rmsd='', test=False):
        """Create a PanelBuilder object and run the analyses as part of the construction.

        Args:
            target (int): The target residue of interest. E.g. 35.
            resids (list): List of contacts you'd like to check.
            rmsd (str, optional): Selection of residues to check RMSD of. E.g. 'resid 15 to 22'. Defaults to ''.
            test (bool, optional): Is this a test run or not? (faster). Defaults to False.
        """

        self.target   = target
        self.resids   = resids
        self.rmsd     = rmsd
        self.test     = test

        if self.test:
            self.sims = ['4HFI_4']
            self.reps = [1]
        else:
            self.sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
            self.reps = [1, 2, 3, 4]

        # WRAPPER FOR MULTITHREADING
        def __task(sim, rep):

            # CREATE THE CHARGE PLOTS
            if notExists("proto_{}_{}_{}.png".format(sim, rep, self.target)) and self.target not in [127, 235, 277]:
                self.chargePlot(sim, rep)

            # CREATE THE MINIMUM DISTANCE PLOTS
            for resid in self.resids:
                if notExists('mindist_{}_{}_{}-{}.png'.format(sim, rep, self.target, resid)):
                    self.mindistPlot(sim, rep, resid)
                    # self.occupancyBarPlot(sim, resid)

            # CREATE THE RMSD PLOTS
            if self.rmsd != '':
                if notExists('rmsd_{}_{}_{}.png'.format(sim, rep, self.target)):
                    self.RMSDPlot(sim, rep)

        # PREPARE ITERABLES
        items = []
        for sim in self.sims:
            for rep in self.reps:
                items.append((sim, rep))

        # RUN MULTITHREADED
        pool = mp.Pool(processes=mp.cpu_count())
        pool.starmap(__task, items, chunksize=1)

        # CREATE BAR PLOTS
        # (this just makes the plots, so do this after multithreaded run)
        self.occupancyBarPlot()

        # CREATE PANELS (this should be the very last step)
        self.rowCount = 0
        self.createPanel()

    def chargePlot(self, sim, rep):
        """
        Make the charge plot in time for residue. Does not currently work for histidines.
        sim: the simulation, e.g. '4HFI_4'
        rep: the replica, e.g. 1
        res: the residue, e.g. 35
        """
        print('Creating charge plot')

        chain = ['A', 'B', 'C', 'D', 'E']
        array = getLambdaFileIndices('{}/{:02d}/CA.pdb'.format(sim, rep), self.target)

        store = []
        for idx in range(0, len(array)):
            data = loadxvg('{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, array[idx]), dt=5000, b=0)
            t    = [val / 1000.0 for val in data[0]]  # ps -> ns
            x    = [1.0 - val for val in data[1]]     # deprotonation -> protonation
            store.append(x)

            plt.plot(t, x, linewidth=1, label=chain[idx])

        # PLOT MEAN PROTONATION
        # average = [0] * len(store[0])
        # for idx in range(0, len(store[0])):
        #     average[idx] = store[0] + store[1] + store[2] + store[3] + store[4]
        # plt.plot(t, x, linewidth=1.5, label='mean', color='b', linestyle=':')

        plt.ylabel('Protonation')
        plt.xlabel('Time (ns)')
        plt.axis([0, 1000, -0.1, 1.1])
        plt.title('{} rep {} residue {}'.format(sim, rep, self.target), fontsize=16)
        plt.legend()
        plt.tight_layout()
        plt.savefig('panels/proto_{}_{}_{}.png'.format(sim, rep, self.target))
        plt.clf()

    def mindistPlot(self, sim, rep, resid):
        """
        Wrapper for GROMACS mindist. Makes the minimum distance plots in time. Uses MD_conv.xtc.
        sim: the simulation, e.g. '4HFI_4'
        rep: the replica, e.g. 1
        resid: the residue it makes contacts with, e.g. 158c
        """
        print("Creating mindist plot")

        # cutoff distance (nm) (hardcoded)
        cutoff = 0.35

        # Go to the simulation directory
        os.chdir('{}/{:02d}'.format(sim, rep))

        # Process ions / the principal vs complementary identifier.
        # Note: 158c gives the correct distances using these lists, and 158 is
        # in fact complementary to 35 so these chain2 orders are correct.
        # Also, we need to make an exception for NA, CL, as shown below.
        chain1 = ['A', 'B', 'C', 'D', 'E']
        if resid not in ['NA', 'CL']:
            if resid[-1] == 'c':
                chain2 = ['E', 'A', 'B', 'C', 'D']
                temp = int(resid[:-1])
            elif resid[-1] == 'p':
                chain2 = ['B', 'C', 'D', 'E', 'A']
                temp = int(resid[:-1])
            else:
                chain2 = ['A', 'B', 'C', 'D', 'E']
                temp = int(resid)
        else:
            temp = resid
            chain2 = 5 * [temp]

        # Create the index file required for the analysis
        stdin = ['q']
        for chain in ['A', 'B', 'C', 'D', 'E'][::-1]:
            # The last part makes sure we only consider the oxygens of ASPT GLUT
            # or the nitrogens in the ring of HSPT.
            stdin.insert(0, 'r {} & chain {} & a OE1 OE2 OD1 OD2 NE2 ND1'.format(self.target, chain))
            stdin.insert(0, 'r {} & chain {}'.format(temp, chain))

        gromacs('make_ndx -f CA.pdb -o mindist.ndx', stdin=stdin)

        # Call GROMACS mindist using the mindist.ndx we just created:
        for idx in range(0, len(chain1)):
            sel1 = 'r_{}_&_ch{}_&_OE1_OE2_OD1_OD2_NE2_ND1'.format(self.target, chain1[idx])
            if resid == 'NA':
                sel2 = 18  # this group number corresponds to NA
            elif resid == 'CL':
                sel2 = 19  # this group number corresponds to CL
            else:  # business as usual
                sel2 = 'r_{}_&_ch{}'.format(temp, chain2[idx])

            if self.test:
                # Speed things up if this is just a test run.
                gromacs('mindist -s MD.tpr -f MD_conv.xtc -n mindist.ndx -dt 10000', stdin=[sel1, sel2])
            else:
                # Normal case.
                gromacs('mindist -s MD.tpr -f MD_conv.xtc -n mindist.ndx -dt 10', stdin=[sel1, sel2])

            data = loadxvg('mindist.xvg')
            t = [val / 1000.0 for val in data[0]]
            x = data[1]
            plt.plot(t, x, linewidth=0.5, label='{}-{}'.format(chain1[idx], chain2[idx]))

            # Compute occupancy fraction and print to file
            contactCount = 0
            for distance in x:
                if distance < cutoff:
                    contactCount += 1
            fraction = contactCount / len(x)

            with open('../../panels/occ_{}_{}_{}-{}.txt'.format(sim, rep, self.target, resid), 'a+') as file:
                file.write('{:.4s} {:.4s} {:.3f}\n'.format(chain1[idx], chain2[idx], fraction))

        # Cleanup and go back
        os.system('rm -f mindist.ndx \\#*\\#')
        os.chdir('../..')

        plt.xlabel("time (ns)")
        plt.xlim(0, 1000)
        if temp in ['NA', 'CL']:  # Use ylim = 2 nm for ions and 1 nm for rest.
            plt.ylim(0, 2)
        else:
            plt.ylim(0, 1)
        plt.ylabel("Minimum distance (nm)")
        plt.hlines(y=cutoff, xmin=0, xmax=1000, color='black', linestyle=':')
        plt.title('{} rep {} mindist {}-{}'.format(sim, rep, self.target, resid), fontsize=16)
        plt.legend()
        plt.tight_layout()
        plt.savefig('panels/mindist_{}_{}_{}-{}.png'.format(sim, rep, self.target, resid))
        plt.clf()

    def RMSDPlot(self, sim, rep):
        """
        Creates RMSD plots for a selection of residues. Loads MD_conv.xtc.
        sim: the simulation, e.g. '4HFI_4'
        rep: the replica, e.g. 1
        MDAnalysis style selection, e.g. 'resid 32 to 35'
        """
        print("Creating RMSD plot")

        path1 = '{}/{:02d}/CA.pdb'.format(sim, rep)
        path2 = '{}/{:02d}/MD_conv.xtc'.format(sim, rep)
        u = MDAnalysis.Universe(path1, path2)

        if self.test:
            step = 100
        else:
            step = 2

        # INDIVIDUAL CHAINS
        chain = ['A', 'B', 'C', 'D', 'E']
        for idx in range(0, len(chain)):
            R = MDAnalysis.analysis.rms.RMSD(u, select='segid {} and {}'.format(chain[idx], self.rmsd))
            R.run(step=step)
            t  = [val / 1000.0 for val in R.rmsd.T[1]]
            x1 = R.rmsd.T[2]
            plt.plot(t, x1, linewidth=0.5, label=chain[idx])

        # ALL CHAINS
        R = MDAnalysis.analysis.rms.RMSD(u, select='(segid A B C D E) and {}'.format(self.rmsd))
        R.run(step=step)
        t  = [val / 1000.0 for val in R.rmsd.T[1]]
        x1 = R.rmsd.T[2]
        plt.plot(t, x1, linewidth=0.5, label='all', color='black')

        plt.xlabel("time (ns)")
        plt.ylabel(r"RMSD ($\AA$)")
        plt.xlim(0, 1000)
        plt.ylim(0, 5)
        plt.title('{} rep {} RMSD {}'.format(sim, rep, self.rmsd), fontsize=16)
        plt.legend()
        plt.tight_layout()
        plt.savefig('panels/rmsd_{}_{}_{}.png'.format(sim, rep, self.target))
        plt.clf()

    def addResidueLetters(self, residList, fileName='4HFI_4/01/CA.pdb'):
        """Adds residue name letters to a list of resids for more readable (bar) plot labels.

        Args:
            residList (list): List of resid strings. E.g. ['79p', '23', 'NA'].
            fileName (str, optional): Reference structure. Defaults to '4HFI_4/01/CA.pdb'.

        Returns:
            List: The enhanced list of strings.
        """

        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
             'ASPT': 'D', 'GLUT': 'E', 'HSPT': 'H'}

        u = MDAnalysis.Universe(fileName)

        array = []
        for resid in residList:
            # If resid, which is a string (possibly '79p'), contains 'p' or 'c'
            # at the end, trim it off. This is to prevent confusing MDAnalysis.
            if resid[-1] in ['c', 'p']:
                temp = resid[:-1]
            else:
                temp = resid

            # NA is just NA, so we can continue.
            if resid == 'NA':
                array.append('Na+')
                continue

            fullName = u.select_atoms('resid {}'.format(temp)).residues[0].resname
            array.append(d[fullName] + str(resid))

        return array

    def occupancyBarPlot(self, width=0.2):
        print("Making occupancy bar plots")
        # For each target residue (e.g. E35) this function creates four plots:
        # one for 4HFI_4, 4HFI_7, etc. Each of these plots contains multiple bars,
        # reflecting the occupancy of the contact of E35 with the various self.resids.
        # The occupancy value is the mean over the replicas and chains (one file).

        superMeanList = []   # Holds four lists, each corresponding to one sim.
        superSerrList = []   # Holds four lists, each corresponding to one sim.

        #! Construct the superDataList.
        #! This is a dictionary of dictionaries of lists.
        #! Holds all actual data sets, not just the means and stderrors.

        A = {}
        for resid in self.resids:
            A[resid] = []

        superDataList = {}
        for sim in self.sims:
            superDataList[sim] = copy.deepcopy(A)

        # print(superDataList)  # debug

        for sim in self.sims:

            # MAKE BARPLOT FOR A SPECIFIC SIMULATION

            meanList = []
            serrList = []

            for resid in self.resids:

                # GATHER MEANs AND SDEVs

                valueList = []
                for rep in self.reps:
                    fname = 'panels/occ_{}_{}_{}-{}.txt'.format(sim, rep, self.target, resid)
                    df    = pandas.read_csv(fname, header=None, delim_whitespace=True, na_filter=False)
                    occ   = list(df.iloc[:, 2])
                    valueList += occ

                # print(resid, len(valueList)) # debug, should be 20 per resid.

                mean  = np.mean(valueList)
                stder = np.std(valueList) / np.sqrt(len(valueList))

                meanList.append(mean)
                serrList.append(stder)
                superDataList[sim][resid] = valueList

            superMeanList.append(meanList)
            superSerrList.append(serrList)

            # MAKE BARPLOT (occ_sim.png)

            plt.bar(self.resids, meanList)
            plt.errorbar(self.resids, meanList, yerr=serrList, fmt='none', capsize=6, linewidth=2)
            plt.ylim(0, 1.1)
            plt.ylabel('Protonation, Contact occupancy')
            plt.title('{} residue {} contact occupancy'.format(sim, self.target))
            plt.tight_layout()
            plt.savefig('panels/occ_{}.png'.format(sim))
            plt.clf()

        # MAKE SUPER BARPLOT I (occ_xx_full.png)
        labels = self.addResidueLetters(self.resids)
        targetWithLetter = self.addResidueLetters([str(self.target)])[0]

        x = np.arange(len(labels))
        fig, ax = plt.subplots()

        # 6ZGD_7
        mean4 = superMeanList[3]
        serr4 = superSerrList[3]
        ax.bar(     x - width * 1.5, mean4, width, color='c', label='closed, pH 7')
        ax.errorbar(x - width * 1.5, mean4, serr4, color='c', fmt='none', capsize=6, linewidth=2)

        # 6ZGD_4
        mean3 = superMeanList[2]
        serr3 = superSerrList[2]
        ax.bar(     x - width / 2.0, mean3, width, color='r', label='closed, pH 4')
        ax.errorbar(x - width / 2.0, mean3, serr3, color='r', fmt='none', capsize=6, linewidth=2)

        # 4HFI_7
        mean2 = superMeanList[1]
        serr2 = superSerrList[1]
        ax.bar(     x + width / 2.0, mean2, width, color='g', label='open, pH 7')
        ax.errorbar(x + width / 2.0, mean2, serr2, color='g', fmt='none', capsize=6, linewidth=2)

        # 4HFI_4
        mean1 = superMeanList[0]
        serr1 = superSerrList[0]
        ax.bar(     x + width * 1.5, mean1, width, color='b', label='open, pH 4')
        ax.errorbar(x + width * 1.5, mean1, serr1, color='b', fmt='none', capsize=6, linewidth=2)

        ax.set_xticks(x, labels)
        ax.legend()

        plt.ylim(0, 1.1)
        plt.ylabel('Contact occupancy')
        plt.title('Residue {}'.format(targetWithLetter))
        plt.tight_layout()
        plt.savefig('panels/occ_{}_full.png'.format(self.target))
        plt.clf()

        # MAKE SUPER BARPLOT II (occ_xx_half.png)

        x = np.arange(len(labels))
        width *= 2
        fig, ax = plt.subplots()

        # 6ZGD_7
        mean4 = superMeanList[3]
        serr4 = superSerrList[3]
        ax.bar(     x - width / 2, mean4, width, color='c', label='closed, pH 7')
        ax.errorbar(x - width / 2, mean4, serr4, color='c', fmt='none', capsize=6, linewidth=2)

        # 4HFI_4
        mean1 = superMeanList[0]
        serr1 = superSerrList[0]
        ax.bar(     x + width / 2, mean1, width, color='b', label='open, pH 4')
        ax.errorbar(x + width / 2, mean1, serr1, color='b', fmt='none', capsize=6, linewidth=2)

        #! Obtains a list 'passes' containing booleans for which of the self.resids are significant.
        passes = []
        for resid in self.resids:
            sample1 = superDataList['4HFI_4'][resid]
            sample2 = superDataList['6ZGD_7'][resid]
            passes.append(ttestPass(sample1, sample2))
        # print(passes)  # debug

        #! The actual plotting part
        X = []
        Y = []
        for idx in range(0, len(passes)):
            if passes[idx]:
                X.append(x[idx])
                Y.append(max(mean1[idx], mean2[idx]) + 0.1)
        plt.scatter(X, Y, marker="*", color='black', linewidth=2)

        ax.set_xticks(x, labels)
        ax.legend()

        plt.ylim(0, 1.1)
        plt.ylabel('Contact occupancy')
        plt.title('Residue {}'.format(targetWithLetter))
        plt.tight_layout()
        plt.savefig('panels/occ_{}_half.png'.format(self.target))
        plt.clf()

    def createPanel(self):
        """Creates the (temporary) panels."""
        print("Creating panels for {}".format(self.target))

        for rep in self.reps:

            # THE FIRST ROW IS ALWAYS A ROW OF CHARGE PLOTS.
            A = 'panels/proto_6ZGD_7_{}_{}.png'.format(rep, self.target)
            B = 'panels/proto_6ZGD_4_{}_{}.png'.format(rep, self.target)
            C = 'panels/proto_4HFI_7_{}_{}.png'.format(rep, self.target)
            D = 'panels/proto_4HFI_4_{}_{}.png'.format(rep, self.target)
            self.rowCount += 1
            os.system('convert {} {} {} {} +append panels/row_{}.png'.format(A, B, C, D, self.rowCount))

            # THE ROWS BELOW COME FROM MINDIST BETWEEN RESIDS
            for resid in self.resids:
                A = 'panels/mindist_6ZGD_7_{}_{}-{}.png'.format(rep, self.target, resid)
                B = 'panels/mindist_6ZGD_4_{}_{}-{}.png'.format(rep, self.target, resid)
                C = 'panels/mindist_4HFI_7_{}_{}-{}.png'.format(rep, self.target, resid)
                D = 'panels/mindist_4HFI_4_{}_{}-{}.png'.format(rep, self.target, resid)
                self.rowCount += 1
                os.system('convert {} {} {} {} +append panels/row_{}.png'.format(A, B, C, D, self.rowCount))

            # THE FINAL ROW IS OPTIONAL AND COMES FROM RMSD
            if self.rmsd != '':
                A = 'panels/rmsd_6ZGD_7_{}_{}.png'.format(rep, self.target)
                B = 'panels/rmsd_6ZGD_4_{}_{}.png'.format(rep, self.target)
                C = 'panels/rmsd_4HFI_7_{}_{}.png'.format(rep, self.target)
                D = 'panels/rmsd_4HFI_4_{}_{}.png'.format(rep, self.target)
                self.rowCount += 1
                os.system('convert {} {} {} {} +append panels/row_{}.png'.format(A, B, C, D, self.rowCount))

            str = ""
            for num in range(1, self.rowCount + 1):
                str += 'panels/row_{}.png '.format(num)

            os.system('convert {} -append panels/panel_{}_{}.png'.format(str, self.target, rep))
            self.rowCount = 0


if __name__ == "__main__":

    loopF = 'resid 152 to 159'

    # THE NEW ANALYSIS (focus on a subgroup of key-residues)

    # 1 (good)
    PanelBuilder(26, ['79p', '80p', '81p', '155', '156', 'NA'])

    # 2 (good)
    PanelBuilder(35, ['114', '116', '156c', '158c', 'NA'], rmsd=loopF)

    # 3 (good)
    PanelBuilder(32, ['119c', '192', 'NA'])
    PanelBuilder(122, ['116', '119', '192', 'NA'])

    # 4 (good)
    PanelBuilder(243, ['200c', '245c', '248', 'NA'])

    # 5 (good)
    PanelBuilder(222, ['277', 'NA'])
    PanelBuilder(277, ['221', '222', 'NA'])

    # THE OLD ANALYSIS (for making the panels)

    # PanelBuilder(26, ['79p', '105', '155', 'NA'])
    # PanelBuilder(32, ['122', '192', 'NA'])
    # PanelBuilder(35, ['29c', '114', '158c', 'NA'], rmsd=loopF)
    # PanelBuilder(67, ['58', '58p', '62', '64', 'NA'])
    # PanelBuilder(69, ['62p', 'NA'])

    # PanelBuilder(97,  ['48', '50', '99', '95', 'NA'])
    # PanelBuilder(104, ['85', '102', 'NA'])
    # PanelBuilder(122, ['32', '116', '119', '192', 'NA'])
    # PanelBuilder(136, ['62c', '138', '179', 'NA'])

    # PanelBuilder(177, ['44c', '148c', '179', 'NA'])
    # PanelBuilder(178, ['148c', 'NA'])
    # PanelBuilder(181, ['179', 'NA'])
    # PanelBuilder(185, ['127', '183', '187', 'NA'])

    # PanelBuilder(222, ['277', 'NA'])
    # PanelBuilder(235, ['260', '263', 'NA'])
    # PanelBuilder(243, ['200c', '245c', '248', 'NA'])
    # PanelBuilder(277, ['221', '222', 'NA'])
