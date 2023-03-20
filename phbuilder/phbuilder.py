#!/usr/bin/python3

# Import the rest of the modules.
import configparser, os, subprocess

from .user import User
from .structure import Structure
from .mdp import gen_mdp

# Stores the information for a lambda group type.
class LambdaType:
    def __init__(self, groupname, incl, pKa, atoms, qqA, qqB, dvdl):
        self.d_groupname = groupname  # str
        self.d_incl      = incl       # list
        self.d_pKa       = pKa        # list  (previously str)
        self.d_atoms     = atoms      # list
        self.d_qqA       = qqA        # list
        self.d_qqB       = qqB        # list of lists (previously list)
        self.d_dvdl      = dvdl       # list of lists (previously list)

# Main phbuilder object.
class phbuilder(User):
    # Construct phbuilder object (handles input parsed from cmdline).
    def __init__(self, CLI):
        # Handle verbosity
        if CLI.verbosity != None: # If -v flag is set...
            self.d_verbosity = 3
        else:
            self.d_verbosity = 2

        super().__init__(self.d_verbosity)

        # Add universal parameters to the universe (used by all three targets).
        self.d_target = CLI.target
        self.d_file   = CLI.file

        # If we run gentopol...
        if (CLI.target == 'gentopol'):
            self.d_output = CLI.output
            self.ph_auto  = CLI.ph # None if not set, else float

            # Process whether the -list flag was or wasn't set.
            if (CLI.list != None):
                resid = []
                for line in open(CLI.list).readlines():
                    resid.append(line.split()[0])

                self.ph_list_resid = [int(i) for i in resid]

        # If we run neutralize...
        elif (CLI.target == 'neutralize'):
            # Either required or has a default value
            self.d_output  = CLI.output
            self.d_topol   = CLI.topol
            self.d_solname = CLI.solname
            self.d_pname   = CLI.pname
            self.d_nname   = CLI.nname
            self.d_conc    = CLI.conc
            self.d_rmin    = CLI.rmin

            # Optional
            if (CLI.nbufs != None):
                self.ph_nbufs = CLI.nbufs
            
            # this needs to be defined because neutralize calls writeLambda_mdp
            self.ph_cal = False

        # If we run genparams...
        elif (CLI.target == 'genparams'):
            self.d_mdp     = CLI.mdp
            self.d_ndx     = CLI.ndx
            self.ph_ph     = CLI.ph
            self.ph_nstout = CLI.nstout
            self.ph_dwpE   = CLI.dwpE
            self.ph_lmass  = CLI.lmass
            self.ph_ltau   = CLI.ltau

            # Process whether the -inter flag was or wasn't set
            if (CLI.inter != None):
                self.ph_inter = True
            else:
                self.ph_inter = False

            # Process whether the -cal flag was or wasn't set
            if (CLI.cal != None):
                self.ph_cal = True
            else:
                self.ph_cal = False

        # User information.
        self.verbose("Parsed the following input from the command line:")
        self.verbose(vars(CLI))

        self.parseLambdaGroupTypesFile()

    # Parse lambdagrouptypes.dat
    def parseLambdaGroupTypesFile(self):
        # Initialize some entries
        self.ph_lambdaTypes = []
        self.ph_BUF_dvdl    = None

        # Add a lambda residue-type to universe.
        def defineLambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl):

            # Create a temporary LambdaType object.
            NewLambdaType = LambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl)

            # Only add NewLambdaType to ph_lambdaTypes if it does not exist yet.
            alreadyPresent = False
            for entry in self.ph_lambdaTypes:
                if entry.d_groupname == NewLambdaType.d_groupname:
                    self.warning("LambdaType with groupname {} is already defined. Skipping...")
                    alreadyPresent = True
                    break

            if not alreadyPresent:
                self.ph_lambdaTypes.append(NewLambdaType)

        # Internal function to convert string to list of floats.
        def str2floatList(string):
            return [float(val) for val in string.split(' ')]

        # Internal function to convert string to list of strings.
        def str2strList(string):
            return string.split(' ')

        # Get path to the packaging_data/ffield directory.
        # packaging_data is in the same dir as phbuilder.py, so we can do this:
        tail, _ = os.path.split(__file__)
        self.p_ffield = os.path.normpath(tail + '/' + 'ffield')

        # If there is a 'lambdagrouptypes.dat' in the working directory, this
        # overrides the default one (located in the packaging_data/ffield).
        if os.path.exists('lambdagrouptypes.dat'):
            self.warning('Found \'lambdagrouptypes.dat\' in working dir, this overrides the default')
            p_lambdagrouptypes = os.path.normpath(os.getcwd() + '/' + 'lambdagrouptypes.dat')
        else:
            p_lambdagrouptypes = os.path.normpath(self.p_ffield + '/' + 'lambdagrouptypes.dat')
            # Error if lambdagrouptypes.dat is not present in directory.
            if not os.path.isfile(p_lambdagrouptypes):
                self.error('Did not find lambdagrouptypes.dat in {}. Please check your directory and/or update your PHFFIELD environment variable in your ~/.bashrc and reload terminal(s).'.format(self.p_ffield))

        # User update
        self.verbose('path to ffield dir = {}'.format(self.p_ffield))
        self.verbose('path to lambdagrouptypes.dat = {}\n'.format(p_lambdagrouptypes))

        # Do the actual parsing
        parser = configparser.ConfigParser()
        parser.read(p_lambdagrouptypes)

        # Loop through the sections.
        for sect in parser.sections():
            # Parse GROMACS parameters
            if (sect.strip() == "GROMACS"):
                self.d_gmxbasepath = parser.get(sect, 'path')

                # Check if there actually is a GROMACS installation in this path.
                if not os.path.isdir(self.d_gmxbasepath):
                    self.verbose("Default path GROMACS base path {} does not seem to exist, will instead use GMXPH_BASEPATH".format(self.d_gmxbasepath))
                    # Get the GROMACS base path from the environment variable.
                    fromEnvVar = os.getenv('GMXPH_BASEPATH')

                    if fromEnvVar == None: # If empty...
                        self.error("Default GROMACS base path {} does not seem to exist, and GMXPH_BASEPATH is not set. Please update your GMXPH_BASEPATH environment variable.".format(self.d_gmxbasepath))
                    if not os.path.isdir(fromEnvVar): # If not empty but not valid...
                        self.error("GMXPH_BASEPATH was found but the specified path {} does not seem to exist. Please update your GMXPH_BASEPATH environment variable.".format(fromEnvVar))

                    self.d_gmxbasepath = fromEnvVar

                continue

            # Parse force field parameters
            if (sect.strip() == "FORCEFIELD"):
                self.d_modelFF    = parser.get(sect, 'name')
                self.d_modelwater = parser.get(sect, 'water')
                continue

            # Parse buffer parameters
            if (sect.strip() == "BUF"):
                self.ph_BUF_dvdl  = str2floatList(parser.get(sect, 'dvdl'))
                self.ph_BUF_range = str2floatList(parser.get(sect, 'range'))
                continue

            # Parse groupname
            groupname = sect.strip()[0:4]

            # Parse incl
            incl = str2strList(parser.get(sect, 'incl'))

            # Parse atoms
            atoms = str2strList(parser.get(sect, 'atoms'))

            # Parse qqA
            qqA = str2floatList(parser.get(sect, 'qqA'))

            pKa  = []
            qqB  = []
            dvdl = []
            for idx in range(1, 11): # Max 10 multisites
                try:
                    # Parse pKa(s)
                    pKa.append(float(parser.get(sect, 'pKA_{}'.format(idx))))

                    # Parse qqB(s)
                    qqB.append(str2floatList(parser.get(sect, 'qqB_{}'.format(idx))))

                    # Parse dvdl(s)
                    dvdl.append(str2floatList(parser.get(sect, 'dvdl_{}'.format(idx))))
                except:
                    break

            # SANITIZE INPUT

            if (len(groupname) < 2 or len(groupname) > 4):
                self.error("groupname of LambdaType needs to contain between 2 and 4 characters.")

            # Call function that constructs the LambdaType object.
            defineLambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl)

        # USER UPDATE
        self.verbose("gmxpath   = {}".format(self.d_gmxbasepath))
        self.verbose("ffname    = {}".format(self.d_modelFF))
        self.verbose("water     = {}\n".format(self.d_modelwater))

        for obj in self.ph_lambdaTypes:
            self.verbose("groupname = {}".format(obj.d_groupname))
            self.verbose("incl      = {}".format(obj.d_incl))
            self.verbose("pKa       = {}".format(obj.d_pKa))
            self.verbose("atoms     = {}".format(obj.d_atoms))
            self.verbose("qqA       = {}".format(obj.d_qqA))
            self.verbose("qqB       = {}".format(obj.d_qqB))
            self.verbose("dvdl      = {}\n".format(obj.d_dvdl))

        self.verbose("BUF_range = {}".format(self.ph_BUF_range))
        self.verbose("BUF_dvdl  = {}".format(self.ph_BUF_dvdl))

    # Function to encapsulate GROMACS calls
    def gromacs(self, command, stdin=[], terminal=False, logFile='builder.log'):
        # If we do not pass any envvars to subprocess (which happens by default) this will work.
        path_to_gmx = os.path.normpath(self.d_gmxbasepath + '/' + 'bin/gmx')
        command = "{} {}".format(path_to_gmx, command)

        if stdin:
            xstr = ' << EOF\n'
            for val in stdin:
                xstr += '{}\n'.format(val)
            command += xstr + 'EOF'

        self.verbose("Running {} ...".format(command))

        if terminal:
            process = subprocess.run(command, shell=True, env={})
        else:
            with open(logFile, 'a+') as file:
                process = subprocess.run(command, shell=True, stdout=file, stderr=file, env={})
        
        if process.returncode != 0:
            if terminal:
                self.error("Failed to run \"{}\" (exitcode {}).".format(command, process.returncode))

            self.error("Failed to run \"{}\" (exitcode {}). Check your logfile ({})".format(command, process.returncode, logFile))

    # Call correct sub function depending on specified target on cmdline.
    def runner(self):
        if self.d_target == 'gentopol':
            self.update('Running gentopol...')
            self.gentopol()

        elif self.d_target == 'neutralize':
            self.update('Running neutralize...')
            self.neutralize()

        elif self.d_target == 'genparams':
            self.update('Running genparams...')
            self.genparams()

        self.pleaseCite()

    # Prepare topology.
    def gentopol(self):
        # Part I - COPY FORCE FIELD AND RESIDUETYPES.DAT TO WORKING DIR

        # /path/to/ffield/charmm36-mar2019-m6.ff
        p_modelFF = os.path.normpath(self.p_ffield + '/' + self.d_modelFF)

        # /path/to/ffield/residuetypes.dat
        p_residuetypes = os.path.normpath(self.p_ffield + '/' + 'residuetypes.dat')

        # User update
        self.verbose('Force field path stuff:')
        self.verbose('path to ffield dir   = {}'.format(self.p_ffield))
        self.verbose('default ffield       = {}'.format(p_modelFF))
        self.verbose('default residuetypes = {}'.format(p_residuetypes))

        # Copy force field and residuetypes.dat from packaging_data/ffield to
        # working dir IF they are not present. Else, use those in working dir.
        if os.path.exists(self.d_modelFF):
            self.update(f'Found \'{self.d_modelFF}\' in working dir, will not (again) copy from default')
        elif os.path.exists(p_modelFF):
            os.system("cp -r {} .".format(p_modelFF))
        else:
            self.error(f"Did not find {self.d_modelFF} in {self.p_ffield}. Please check your directory and/or update your PHFFIELD environment variable.")

        if os.path.exists('residuetypes.dat'):
            self.update('Found \'residuetypes.dat\' in working dir, will not (again) copy from default')
        elif os.path.exists(p_residuetypes):
            os.system("cp {} .".format(p_residuetypes))
        else:
            self.error(f"Did not find residuetypes.dat in {self.p_ffield}. Please check your directory and/or update your PHFFIELD environment variable.")

        # Remove .ff extension from force field (we need d_modelFF for pdb2gmx).
        tail, head     = os.path.split(p_modelFF)
        self.d_modelFF = os.path.splitext(head)[0]

        # PART II - LOAD DATA

        # Load the .pdb/.gro structure (and implicitly phrecord.dat, if it exists).
        pdb = Structure(self.d_file, self.d_verbosity)

        # Load user-specified list of residues (if any).
        if hasattr(self, 'ph_list_resid'):
            list_resid = self.ph_list_resid

            self.update('Found user-defined list containing residues to consider...')
            self.verbose(list_resid)

        # PART III - LOOP THROUGH THE RESIDUES AND MODIFY

        self.update("Modifying structure file...")

        # Loop through all the residue objects.
        for residue in pdb.d_residues:

            # Get the lambdaType object for which d_incl contains residue.d_resname (if any).
            associatedLambdaType = [lambdaType for lambdaType in self.ph_lambdaTypes if residue.d_resname in lambdaType.d_incl]

            # If this list is not empty, i.e. an associated LambdaType was found...
            if associatedLambdaType:

                # Turn this list, containing one object, into the object.
                associatedLambdaType = associatedLambdaType[0]

                # If -list was set AND the current residue is not in the list, continue to the next residue.
                if hasattr(self, 'ph_list_resid') and (residue.d_resid not in list_resid):
                    continue

                # Store original name here as we'll need it later for a user update.
                origName = residue.d_resname

                # If -ph was set (contains a float), do not ask for every residue individually:
                if not type(self.ph_auto) == float:
                    # List to hold the various options for the user.
                    options = ["Keep current (static) protonation state"]

                    # Is this multistate yes/no? Multistate and 2state need to be treated differently.
                    multistate = len(associatedLambdaType.d_pKa) > 1

                    if multistate:
                        chargeLists = associatedLambdaType.d_qqB
                        for idx in range(0, len(chargeLists)):
                            options.append("Make titratable (rename to {}), start in state {} (q = {:+.2f})".format(associatedLambdaType.d_groupname, idx + 1, sum(chargeLists[idx])))

                        val = self.inputOptionHandler("Choose what to do with residue {}-{} in chain {}".format(residue.d_resname, residue.d_resid, residue.d_chain), options)

                        if val != 0:
                            residue.d_resname = associatedLambdaType.d_groupname
                            residue.d_init = str(val)

                    if not multistate:
                        chargeLists = [associatedLambdaType.d_qqA, associatedLambdaType.d_qqB[0]]
                        for idx in range(0, len(chargeLists)):
                            options.append("Make titratable (rename to {}), start in state {} (q = {:+.2f})".format(associatedLambdaType.d_groupname, idx, sum(chargeLists[idx])))

                        val = self.inputOptionHandler("Choose what to do with residue {}-{} in chain {}".format(residue.d_resname, residue.d_resid, residue.d_chain), options)

                        if val != 0:
                            residue.d_resname = associatedLambdaType.d_groupname
                            residue.d_init = str(val - 1)

                else: # If -ph was set, just protonate everything automatically.
                    residue.d_resname = associatedLambdaType.d_groupname

                    # Additionally, use automaticLambdaInits to set the initial lambda values.
                    init = self.automaticLambdaInits(associatedLambdaType.d_pKa, self.ph_auto)
                    residue.d_init = str(init)

                    self.update("Made residue {:>4s}-{:<4d} in chain {:1s} titratable (changed name to {}, initlambda = {})".format(origName, residue.d_resid, residue.d_chain, residue.d_resname, init))

                continue # Otherwise code below is executed for the same residue we just changed.

            # Get the lambdaType object for which d_groupname = residue.d_resname (if any).
            associatedLambdaType = [lambdaType for lambdaType in self.ph_lambdaTypes if residue.d_resname == lambdaType.d_groupname]

            # If this list is not empty, i.e. an associated LambdaType was found,
            if associatedLambdaType:

                # Turn this list, containing one object, into the object.
                associatedLambdaType = associatedLambdaType[0]

                # If -list was set AND the current residue is not in the list, continue to the next residue.
                if hasattr(self, 'ph_list_resid') and (residue.d_resid not in list_resid):
                    continue

                # If -ph was set (contains a float), do not ask for every residue individually:
                if not type(self.ph_auto) == float:

                    # List to hold the various options for the user.
                    options = []

                    # Is this multistate yes/no? Multistate and 2state need to be treated differently.
                    multistate = len(associatedLambdaType.d_pKa) > 1

                    if multistate:
                        chargeLists = associatedLambdaType.d_qqB
                        for idx in range(0, len(chargeLists)):
                            options.append("Keep titratable, start in state {} (q = {:+.2f})".format(idx + 1, sum(chargeLists[idx])))

                    if not multistate:
                        chargeLists = [associatedLambdaType.d_qqA, associatedLambdaType.d_qqB[0]]
                        for idx in range(0, len(chargeLists)):
                            options.append("Keep titratable, start in state {} (q = {:+.2f})".format(idx, sum(chargeLists[idx])))  

                    for name in associatedLambdaType.d_incl:
                        options.append("Change name to {}".format(name))

                    val = self.inputOptionHandler("Choose what to do with residue {}-{} in chain {} in initial state {}".format(residue.d_resname, residue.d_resid, residue.d_chain, residue.d_init), options)

                    # If user picks one of the "still titratable" options:
                    if val in range(0, len(chargeLists)):
                        residue.d_init = str(val)

                    else: # Else the user picked one of the non-titratable options
                        residue.d_resname = associatedLambdaType.d_incl[val - len(chargeLists)]
                        residue.d_init = ' '

                # Normally, we would not have to do anything. However, if we are
                else: # missing initial lambda values in the record, we need to add them.
                    init = self.automaticLambdaInits(associatedLambdaType.d_pKa, self.ph_auto)
                    residue.d_init = str(init)

                    self.update("Residue {}-{} in chain {} is already titratable (initlambda = {})".format(residue.d_resname, residue.d_resid, residue.d_chain, init))

            # If the residue in question is neither an ASP nor an ASPT, as an extra
            else: # measure we make sure that nothing is present in the init field.
                residue.d_init = ' '

        # Write the modified structure file (input for pdb2gmx).
        pdb.write('phset.pdb')
        self.d_file = 'phset.pdb'

        # PART 3.5 - PRINT A REFERENCE LIST OF LAMBDA COORDINATES

        self.update('writing lambda coordinates reference list (lambdareference.dat)...')
        LambdaTypeDict = {}  # Make dict of names and associated number of coordinate files.

        for obj in self.ph_lambdaTypes:
            if len(obj.d_pKa) > 1:  # if this is a multistate...
                # We have #coordinate (files) = number of qqB lists.
                LambdaTypeDict[obj.d_groupname] = len(obj.d_qqB)
            else:
                # If this is not multisite we only have one coordinate (file).
                LambdaTypeDict[obj.d_groupname] = 1

        with open('lambdareference.dat', 'w+') as file:
            file.write('resname coord# resid chain coordinateFile\n')

            coordCount = 1
            for residue in pdb.d_residues:
                if residue.d_resname in LambdaTypeDict.keys():
                    multiCount = 1
                    for _ in range(0, LambdaTypeDict[residue.d_resname]):
                        file.write(f"{residue.d_resname} {multiCount} {residue.d_resid:4d} {residue.d_chain:1s} cphmd-coord-{coordCount}.xvg\n")
                        coordCount += 1
                        multiCount += 1

        # PART IV - HANDLE UNKNOWN RESIDUE TYPES

        self.update("Checking if every residue type is present in residuetypes.dat...")

        # Load residuetypes.dat into a dictionary
        residueTypes = {}
        for val in open('residuetypes.dat').readlines():
            residueTypes[val.split()[0]] = val.split()[1]

        # Compile lists of the unknown residue(s)(types)
        unknownResidues = []
        unknownResTypeNames = []
        for residue in pdb.d_residues:
            # If it's either not in residuetypes at all, or it is but it is an
            # we add it to the list of unknown residues.
            if residue.d_resname not in residueTypes or residueTypes[residue.d_resname] in ['Ion']:
                unknownResidues.append(residue)

                if residue.d_resname not in unknownResTypeNames:
                    unknownResTypeNames.append(residue.d_resname)

        pathList = []   # Loop through the unknown residue types and add them to 
        skipList = []   # skipList and (the manually specified path to) pathList.
        good = True
        for val in unknownResTypeNames:
            self.warning(f"residue type '{val}' in '{self.d_file}' was not found in residuetypes.dat associated with '{self.d_modelFF}'. phbuilder cannot provide topology information for non-standard residue types. Therefore, we expect the user to provide a separate topology (.itp) file. Additionally, be aware that pdb2gmx cannot add hydrogens to non-standard residues, so you either have to add an entry to the force field .hdb file, or make sure that hydrogens are already present in the structure for this residue type. See also https://manual.gromacs.org/current/how-to/topology.html.")
            path = input("phbuilder : Specify /path/to/some.itp (or enter to ignore) (ions can be ignored): ")

            skipList.append(val)
            pathList.append(path)
            good = False

            self.update("Set path for residue type {} to {}...".format(val, path))

        if good:
            self.update("everything seems OK.")

        else:
            # Create a list knownResidues containing only the known residues
            knownResidues = []
            for residue in pdb.d_residues:
                if residue.d_resname not in skipList:
                    knownResidues.append(residue)

            # Update d_residues in universe with knownResidues
            pdb.d_residues = knownResidues

            # Write temporary .pdb containing only the known residues
            someTempName = 'pdb2gmxtemp.pdb'
            pdb.write(someTempName)
            self.d_file = someTempName

            # Update value of d_output (so that pdb2gmx is called on the temporary 
            # structure), and backup the final output name as specified by user.
            self.d_output_orig = self.d_output
            self.d_output      = someTempName

        # PART IV - RUN PDB2GMX AND ASK USER FOR INPUT ABOUT IT

        # FIX: If 'pdb2gmxtemp.pdb' is empty (e.g. because ALL residuetypes are 
        # unrecognized) pdb2gmx will crash. We therefore skip this part entirely.
        if len(pdb.d_residues) == 0:

            self.update('skipping pdb2gmx step entirely as there are no recognized residues.')
            with open('topol.top', 'w+') as file:
                file.write('; Include forcefield parameters\n')
                file.write('#include \"{}.ff/forcefield.itp\"\n\n'.format(self.d_modelFF))
                file.write('; Include water topology\n')
                file.write('#include \"{}.ff/{}.itp\"\n\n'.format(self.d_modelFF, self.d_modelwater))
                file.write('; Include topology for ions\n')
                file.write('#include \"{}.ff/ions.itp\"\n\n'.format(self.d_modelFF))
                file.write('[ system ]\n; Name\n{}\n\n'.format(pdb.d_title))
                file.write('[ molecules ]\n; Compounds  #mols\n')

        else:
            self.update("\nRecommended pdb2gmx command:")
            self.update("gmx pdb2gmx -f {} -o {} -ff {} -water {} -ignh".format(self.d_file, self.d_output, self.d_modelFF, self.d_modelwater))

            # Ask for input for what to do regarding pdb2gmx
            val = self.inputOptionHandler(
                "Choose whether to", 
                ["Do nothing", "Run", "Add additional flags (https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html)"])

            flags = ""
            if val == 2:
                flags += input("phbuilder : Enter flags: ")

            # Run pdb2gmx:
            if val in [1, 2]:
                self.gromacs("pdb2gmx -f {} -o {} -ff {} -water {} -ignh {}".format(self.d_file, self.d_output, self.d_modelFF, self.d_modelwater, flags), terminal=True)

            # If we do nothing, then return because we do not want to do PART V?
            if val == 0:
                return

        # PART V - MERGE THE .PDB FILES

        # If pathList is not empty, i.e. if we had at least one unknown residue:
        if pathList:
            # Load the structure output of pdb2gmx (update d_residues in universe)
            pdb.read(someTempName)

            # Merge the processed structure of good residues with unknown residues
            mergedResidues = pdb.d_residues + unknownResidues

            # Update d_residues in universe
            pdb.d_residues = mergedResidues

            # Write the final structure
            pdb.write(self.d_output_orig)

        # PART VI - MERGE THE TOPOLOGIES

        # Write manually specified files to topol.top
        def add_mol(itpfname, comment, molname=None, molcount=None):
            # Get the contents of current topol.top.
            topList = []
            with open("topol.top") as file:
                for line in file.readlines():
                    topList.append(line)

            # Add the .itp file (line saying: #include "blabla.itp")
            with open("topol.top", 'w') as file:
                try:
                    for line in range(0, len(topList)):
                        file.write(topList[line])

                        if "[ system ]\n" == topList[line + 1] and itpfname != '':
                            file.write("; {0}\n".format(comment))
                            file.write("#include \"{0}\"\n\n".format(itpfname))

                except IndexError:
                    pass

            # if molcount not present, add it, otherwise do nothing.
                if molname != None and molcount != None and molname not in topList[-1]:
                    file.write("{0}\t\t\t{1}\n".format(molname, molcount))

        # If pathList is not empty, i.e. if we had at least one unknown residue:
        if pathList:
            # Remove temporary .pdb file
            os.remove(someTempName)

            # Loop through the unknown residue types
            for idx in range(0, len(skipList)):

                # For each one, count how many there are
                count = 0
                for residue in pdb.d_residues:
                    if residue.d_resname == skipList[idx]:
                        count += 1

                # Add manually to topol.top
                add_mol(pathList[idx], "Include topology for {}".format(skipList[idx]), skipList[idx], count)

        self.update("Finished generating topology for constant-pH.")

    def countRes(self, Structure, resname):
        count = 0
        for residue in Structure.d_residues:
            if residue.d_resname == resname:
                count += 1
        return count

    # Add appropriate number ions and buffers to make the system neutral.
    def neutralize(self):
        # I - LOAD DATA, PERFORM SOME CHECKS

        # Load the input structure into d_residues.
        pdb = Structure(self.d_file, self.d_verbosity)

        # Perform some basic checks for PBC box, solvent, and titratable groups.
        foundPBC         = hasattr(pdb, 'd_box')
        foundSolvent     = False
        foundTitratables = False

        LambdaTypeNames = []
        for LambdaType in self.ph_lambdaTypes:
            LambdaTypeNames.append(LambdaType.d_groupname)

        for residue in pdb.d_residues:

            if residue.d_resname == self.d_solname:
                foundSolvent = True

            elif residue.d_resname in LambdaTypeNames:
                foundTitratables = True

            if foundSolvent and foundTitratables:
                break

        # Error if no periodic box was found
        if not foundPBC:
            self.error("{} does not have a periodic box! Did you forget to add one?".format(self.d_file))

        # Error if no solvent was found
        if not foundSolvent:
            self.error("{} does not seem to have any solvent with molname {}! Did you forget to add solvent? Alternatively you can change the solvent name by setting -solname.".format(self.d_file, self.d_solname))

        # Error if no titratable residues were found
        # We need to throw an error here and not a warning because grompp in getcpHMDcharge will not compile without at least one lambda group being present.
        if not foundTitratables:
            self.error("{} does not seem to have any titratable groups! Did you forget to run gentopol?".format(self.d_file))

        # II - DO THE IONS PART

        def getcpHMDcharge(name, file, Structure, LambdaTypeNames, constrainCharge):
            # Create empty name.mdp to append the lambda dynamics parameters to.
            open('{}.mdp'.format(name), 'w').close()

            # Set dummy parameters
            self.ph_ph    = 7.0; self.ph_nstout = 500
            self.ph_lmass = 5.0; self.ph_inter  = False
            self.ph_ltau  = 2.0; self.ph_dwpE   = 7.5

            # Write the lambda dynamics parameters to name.mdp
            self.writeLambda_mdp(name, Structure, LambdaTypeNames, constrainCharge)

            # Write index groups to name.ndx
            self.writeLambda_ndx('{}.ndx'.format(name), Structure, LambdaTypeNames, constrainCharge)

            # Execute grompp on this using cpHMD GROMACS version and direct output to charge.log.
            self.update("Getting the net-charge of the system (running grompp, this can take some time for larger systems)...")
            self.gromacs("grompp -f {0}.mdp -c {1} -p {2} -n {0}.ndx -o {0}.tpr".format(name, file, self.d_topol), logFile='{}.log'.format(name))

            # Grep total charge from grompp output message
            QQtotalcpHMD = 0
            for line in open('{}.log'.format(name)).read().splitlines():
                if "non-zero total charge" in line:
                    QQtotalcpHMD = float(line.split()[-1])
                    break

            # Cleanup
            os.remove('{}.mdp'.format(name))
            os.remove('{}.ndx'.format(name))
            os.remove('{}.log'.format(name))

            return QQtotalcpHMD

        # boolean to prevent execution for when the system is already neutral.
        isNeutral = False

        # If we already have ions in the input file, give the user an extra update.
        def getIonConcentration(Structure, Nwater):
            totalIons = self.countRes(Structure, self.d_pname) + self.countRes(Structure, self.d_nname)
            return totalIons / (0.01808 * Nwater)

        Nwater         = self.countRes(pdb, self.d_solname)
        ionConcInInput = getIonConcentration(pdb, Nwater)

        if ionConcInInput:
            self.update("Detected {} {} and {} {} in {}, corresponding to a concentration of {:.3f} mol/L for given solvent volume.".format(
                self.countRes(pdb, self.d_pname),
                self.d_pname,
                self.countRes(pdb, self.d_nname),
                self.d_nname,
                self.d_file,
                ionConcInInput))

        # If we already have ions in the input file, d_conc can never be smaller than what we ALREADY HAVE.
        if self.d_conc and ionConcInInput > self.d_conc:
            self.error('Target ion concentration of {} mol/L is smaller than the current ion concentration in {}. Either remove -conc or increase it.'.format(self.d_conc, self.d_file))

        # Get the cpHMD charge (we set constrainCharge to False because we have no buffers at this point).
        QQtotalcpHMD = getcpHMDcharge('ions', self.d_file, pdb, LambdaTypeNames, constrainCharge=False)

        self.update('Net-charge of system     = {:+.2f}'.format(QQtotalcpHMD))

        np = 0; nn = 0
        if QQtotalcpHMD > 0:
            nn = round(QQtotalcpHMD)
        elif QQtotalcpHMD < 0:
            np = round(-QQtotalcpHMD)
        else:
            self.update('System seems to be neutral already.')
            isNeutral = True

        def fromBoxVol():
            Avo    = 6.02214076 * 10 ** 23
            boxVol = pdb.d_box.d_a * pdb.d_box.d_b * pdb.d_box.d_c
            return int(round(self.d_conc * Avo * boxVol * 10 ** -27))

        # Return the number of ions required to achieve a certain ion concentration.
        def fromSolVol(Structure, concentration):
            Nwater = self.countRes(Structure, self.d_solname)
            return int(round(0.01808 * Nwater * concentration))

        if not (isNeutral and self.d_conc == 0):

            if self.d_conc != 0:
                Nions = fromSolVol(pdb, self.d_conc - ionConcInInput)
                ionsMinRequired = np + nn               

                if Nions < ionsMinRequired:
                    self.error('Target ion concentration of {} mol/L is smaller than the minimum (additional) concentration required to neutralize the system. Either remove -conc or increase it.'.format(self.d_conc))

                # We only want integers so we want to divide an even number.
                if (Nions - np - nn) % 2 != 0:
                    Nions += 1

                factor = int((Nions - np - nn)/2.0)

                np += factor
                nn += factor

                self.update('Specified ion concentration of {} mol/L corresponds to {} (additional) ions for given solvent volume'.format(self.d_conc, Nions))

            self.update('Will add {} positive ({}) and {} negative ({}) ions...'.format(np, self.d_pname, nn, self.d_nname))
            self.update('Total charge to be added = {:+.2f}'.format(np - nn))

            # Run genion to add the appropriate number of ions.
            self.gromacs("genion -s ions.tpr -o phions.pdb -p {} -pname {} -nname {} -np {} -nn {} -rmin {}".format(
                self.d_topol, self.d_pname, self.d_nname, np, nn, self.d_rmin), stdin=['SOL']) # this is always SOL, even if the molname is e.g. HOH...

            self.update('Finished adding ions')

        else:
            # Get the number of ions (that are already making the system neutral from input file).
            # We will need these later for the checking part.
            np = self.countRes(pdb, self.d_pname)
            nn = self.countRes(pdb, self.d_nname)

            # User update and write internal structure to phions.pdb.
            self.update("System is already neutral and -conc was not set. Will not add any ions...")
            pdb.write('phions.pdb')

        # III - DO THE BUFFER PART

        # Load the structure that now contains the ions.
        pdb.read('phions.pdb')

        # Boolean to prevent execution for when the system already has enough buffers.
        hasEnoughBUFs = False

        # Obtain the number of buffers that should be added.
        if not hasattr(self, 'ph_nbufs'):
            # Count number of titratable residues
            titratables = 0
            for residue in pdb.d_residues:
                if residue.d_resname in LambdaTypeNames:
                    titratables += 1

            # This is the equation we use to estimate the required number of buffers (hardcoded).
            q_max = abs(self.ph_BUF_range[0])
            nbufs = int(titratables / (2 * q_max))

            self.update("Will add {} buffer(s) based on guess ( = Tsites / 2q_max )...".format(nbufs))

        else:
            nbufs = self.ph_nbufs
            self.update("Will add {} buffer(s) (user specified)...".format(nbufs))

        nbufspresent = self.countRes(pdb, 'BUF')

        if nbufspresent == nbufs:
            self.warning("The number of buffers in {} is already equal to the request number of buffers. Will not add any buffers...".format(self.d_file))
            hasEnoughBUFs = True
        
        elif nbufspresent > nbufs:
            self.warning("The number of buffers in {} exceeds the requested number of buffers. Will not add any buffers...".format(self.d_file))
            hasEnoughBUFs = True

        if not hasEnoughBUFs:

            # We don't need the cpHMD parameters for adding the buffers, therefore
            # we create a dummy buffers.mdp file (containing nothing).
            open('buffers.mdp', 'w+').close()
            self.gromacs("grompp -f buffers.mdp -c phions.pdb -p {} -o buffers.tpr".format(self.d_topol))

            # Run genion to add the appropriate number of buffers.
            self.gromacs("genion -s buffers.tpr -o {} -p {} -pname BUF -np {} -rmin {}".format(self.d_output, self.d_topol, nbufs, self.d_rmin), stdin=['SOL']) # this is always SOL, even if the molname is e.g. HOH...

            self.update('Finished adding buffers')

            os.remove('buffers.mdp')
            os.remove('buffers.tpr')

        # IV - WRAPUP : CHECK IF EVERYTHING IS NOW CORRECT

        os.remove('ions.tpr')

        # Check if the correct number of ions and buffers are present in the output file.

        self.update('\nChecking whether everything was succesful:\n')

        # Update internal pdb record to phneutral.pdb
        pdb.read(self.d_output)

        countPions = self.countRes(pdb, self.d_pname)
        if countPions != np:
            self.warning('Detected {}/{} required positive ions (ignore this if your input already contained ions).'.format(countPions, np))
        else:
            self.update('Detected correct number of positive ions ({}) - check'.format(np))

        countNions = self.countRes(pdb, self.d_nname)
        if countNions != nn:
            self.warning('Detected {}/{} required negative ions (ignore this if your input already contained ions).'.format(countNions, nn))
        else:
            self.update('Detected correct number of negative ions ({}) - check'.format(nn))

        countBUFs  = self.countRes(pdb, 'BUF')
        if countBUFs != nbufs:
            self.warning('Detected {}/{} required buffers'.format(countBUFs, nbufs))
        else:
            self.update('Detected correct number of buffers ({}) - check'.format(nbufs))

        # Give a user update about the ion concentration in the output file.
        self.update("{} {} and {} {} in {} corresponds to a concentration of {:.3f} mol/L for given solvent volume...".format(
            self.countRes(pdb, self.d_pname),
            self.d_pname,
            self.countRes(pdb, self.d_nname),
            self.d_nname,
            self.d_output,
            getIonConcentration(pdb, self.countRes(pdb, 'SOL')))) # For some reason genion renames whatever solname we had in the input file to 'SOL'...

        # Check if the charge is now neutral.

        QQFinal = getcpHMDcharge('check', self.d_output, pdb, LambdaTypeNames, constrainCharge=True)
        if QQFinal != 0:
            self.warning('Output system ({}) is not neutral at cpHMD (net-charge = {:+.2f})! Something went wrong! Check your log files...'.format(self.d_output, QQFinal))
        else:
            self.update('Output system ({}) is neutral at cpHMD (net-charge = {:+.2f}) - check'.format(self.d_output, QQFinal))

        os.remove('check.tpr')
        self.update("Finished phbuilder neutralize.")

    # Generate the actual lambda dynamics parameters.
    def writeLambda_mdp(self, Type, Structure, LambdaTypeNames, constrainCharge):
        file = open("{}.mdp".format(Type), 'a')

        # Formatting function for adding parameters.
        def addParam(name, value):
            file.write("{:54s} = {:13s}\n".format(name, str(value)))

        # PART 1 - WRITE GENERAL PARAMETERS

        file.write("\n; CONSTANT PH\n")

        addParam('lambda-dynamics', 'yes')
        addParam('lambda-dynamics-simulation-ph', "{:.1f}".format(self.ph_ph))
        addParam('lambda-dynamics-lambda-particle-mass', "{:.1f}".format(self.ph_lmass))
        addParam('lambda-dynamics-tau', "{:.1f}".format(self.ph_ltau))
        addParam('lambda-dynamics-update-nst', self.ph_nstout)

        # If we use charge constraining...
        if constrainCharge:
            addParam('lambda-dynamics-charge-constraints', 'yes')

        # If we are doing EM/NVT/NPT, we do not want the lambdas to move:
        if Type in ['EM', 'NVT', 'NPT'] or self.ph_cal:
            addParam('lambda-dynamics-calibration', 'yes')

        # We need to count how many titratable residues in total we have in the
        # protein. For this we compile a list LambdasFoundinProtein.
        LambdasFoundinProtein = [] # (e.g ASPT ASPT GLUT ASPT GLUT ASPT...)

        # Stores the number of buffer atoms/ions.
        buffersFoundinProtein = 0

        for residue in Structure.d_residues:
            if residue.d_resname in LambdaTypeNames:
                LambdasFoundinProtein.append(residue.d_resname)

            # Also in this loop we count how many buffer ions we have.
            elif constrainCharge and residue.d_resname == 'BUF':
                buffersFoundinProtein += 1

        # The fact that we have a LambdaType in lambdagrouptypes.dat does not mean
        # one of those is also present in the protein. In that case, we want to 
        # prevent counting it, so we compile a subgroup of LambdaType groupnames
        # that are not only in lambdagrouptypes.dat, but ALSO found at least once 
        # in the actual protein.
        LambdaTypeNamesFoundinProtein = list(set(LambdasFoundinProtein)) # (e.g ASPT GLUT)

        # If we not only have multistate LambdaResidueTypes defined in the .dat file,
        # but we also have detected such a LambdaResidueType in the actual protein,
        # we'll need to turn on multistate:
        for LambdaType in self.ph_lambdaTypes:
            if len(LambdaType.d_pKa) > 1 and LambdaType.d_groupname in LambdaTypeNamesFoundinProtein:
                addParam('lambda-dynamics-multistate-constraints', 'yes')

        # If we use charge constraining we also have he BUF residue-type, as well as 
        # one extra lambda group containing all the BUFs.
        if constrainCharge:
            addParam('lambda-dynamics-number-lambda-group-types', len(LambdaTypeNamesFoundinProtein) + 1)
            addParam('lambda-dynamics-number-atom-collections', len(LambdasFoundinProtein) + 1)
        else:
            addParam('lambda-dynamics-number-lambda-group-types', len(LambdaTypeNamesFoundinProtein))
            addParam('lambda-dynamics-number-atom-collections', len(LambdasFoundinProtein))

        file.write('\n')

        # PART 2 - WRITE LAMBDA GROUP TYPES

        # Convert a list to a string
        def to_string(Input, round):
            string = ""
            for element in Input:
                string += "{0:.{arg}f} ".format(element, arg=round)
            return string

        # Writes the lambda group type block
        def writeLambdaGroupTypeBlock(number, name, multistates, qqA, qqB, pKa, dvdl):
            addParam('lambda-dynamics-group-type{}-name'.format(number), name)
            addParam('lambda-dynamics-group-type{}-n-states'.format(number), multistates)
            addParam('lambda-dynamics-group-type{}-state-0-charges'.format(number), to_string(qqA, 3))

            for idx in range(1, multistates + 1):
                # When we have a multistate lambdagrouptype, one of the pKas should 
                # be equal to the simulation-pH. This is done by setting this pKa 
                # to zero in the lambdagrouptypes.dat file.
                pKaNew = pKa[idx-1]
                if multistates > 1 and float(pKaNew) == 0.0:
                    pKaNew = self.ph_ph

                addParam('lambda-dynamics-group-type{}-state-{}-charges'.format(number, idx), to_string(qqB[idx-1], 3))
                addParam('lambda-dynamics-group-type{}-state-{}-reference-pka'.format(number, idx), pKaNew)
                addParam('lambda-dynamics-group-type{}-state-{}-dvdl-coefficients'.format(number, idx), to_string(dvdl[idx-1], 3))

            file.write('\n')

        number = 1
        # We loop over the object itself instead of the d_groupname as we need all
        # the information in the object.
        for LambdaType in self.ph_lambdaTypes:
            # This if-statement prevents writing a block when there are no residues 
            # of this type in the protein.
            if (LambdaType.d_groupname in LambdaTypeNamesFoundinProtein):

                writeLambdaGroupTypeBlock(
                    number,
                    LambdaType.d_groupname,
                    len(LambdaType.d_pKa),
                    LambdaType.d_qqA,
                    LambdaType.d_qqB,
                    LambdaType.d_pKa,
                    LambdaType.d_dvdl
                )

                number += 1

        # If we do charge constraining, we additionally need the block for the buffer.
        if constrainCharge:
            qqA = self.ph_BUF_range[0]
            qqB = self.ph_BUF_range[1]

            writeLambdaGroupTypeBlock(
                number,
                'BUF',
                1,
                [qqA],
                [[qqB]],
                [0],
                [self.ph_BUF_dvdl]
            )

        # PART 3 - WRITE LAMBDA GROUPS

        def writeResBlock(number, name, initList, Edwp):
            addParam('lambda-dynamics-atom-set{}-name'.format(number), name)
            addParam('lambda-dynamics-atom-set{}-index-group-name'.format(number), 'LAMBDA{}'.format(number))
            addParam('lambda-dynamics-atom-set{}-initial-lambda'.format(number), to_string(initList, 1))
            addParam('lambda-dynamics-atom-set{}-barrier'.format(number), Edwp)

            if constrainCharge:
                addParam('lambda-dynamics-atom-set{}-charge-restraint-group-index'.format(number), 1)

            if (name == 'BUF'):
                addParam('lambda-dynamics-atom-set{}-buffer-residue'.format(number), 'yes')
                addParam('lambda-dynamics-atom-set{}-buffer-residue-multiplier'.format(number), buffersFoundinProtein)

            file.write('\n')

        number = 1

        for residue in Structure.d_residues:

            if residue.d_resname in LambdaTypeNamesFoundinProtein:

                LambdaType = [obj for obj in self.ph_lambdaTypes if obj.d_groupname == residue.d_resname][0]
                multistate = len(LambdaType.d_pKa) > 1

                # HANDLE INTERACTIVE EDWP

                # Only do this if -inter was set AND we're writing params for production run (no point for EM/NVT/NPT).
                if self.ph_inter and Type == 'MD':
                    Edwp = input("phbuilder : set bias barrier (kJ/mol) for {}-{} in chain {} (or hit enter for default = {}): ".format(residue.d_resname, residue.d_resid, residue.d_chain, self.ph_dwpE))
                    if Edwp == '':
                        Edwp = self.ph_dwpE
                    else:
                        Edwp = float(Edwp)
                        self.verbose("Setting custom Edwp ({}) for residue {}-{} in chain {} (LAMBDA{})".format(Edwp, residue.d_resname, residue.d_resid, residue.d_chain, number))
                else:
                    Edwp = self.ph_dwpE

                # HANDLE INITIAL LAMBDA FROM PHRECORD.DAT

                init = residue.d_init
                if init == '':
                    init = self.automaticLambdaInits(LambdaType.d_pKa, self.ph_ph)
                    self.warning('Initial lambda value for residue {}-{} in chain {} not found! Check your phrecord.dat file. Will revert to automatic (state {}), but your system might no-longer be neutral.'.format(residue.d_resname, residue.d_resid, residue.d_chain, init))

                if multistate:
                    initList = []
                    for idx in range(0, len(LambdaType.d_pKa)):
                        initList.append(idx == int(init) - 1)

                if not multistate:
                    initList = [int(init)]

                # DO THE WRITING

                writeResBlock(number, residue.d_resname, initList, Edwp)

                number += 1

        if constrainCharge:
            # Write block for the buffer.
            if self.ph_cal:
                writeResBlock(number, 'BUF', initList=[1], Edwp=0.0)
            else:
                writeResBlock(number, 'BUF', initList=[0.5], Edwp=0.0)

        file.close() # MD.mdp

    def writeLambda_ndx(self, fileName, Structure, LambdaTypeNames, constrainCharge):
        # Create list of Lambda group types found in protein.
        LambdasFoundinProtein = []
        for residue in Structure.d_residues:
            if residue.d_resname in LambdaTypeNames:
                LambdasFoundinProtein.append(residue.d_resname)

        file = open(fileName, 'a')

        # Formatting function for writing the index block.
        def writeTheGroup(number, atomIndexList):
            file.write('\n[ LAMBDA{} ]\n'.format(number))
            for index in atomIndexList:
                file.write('{} '.format(index))
            file.write('\n')

        atomCount   = 1
        groupNumber = 1
        bufferIndexList = []

        # Write the lambda index groups for the titratable residues
        for residue in Structure.d_residues:
            # If the residue is titratable
            if residue.d_resname in list(set(LambdasFoundinProtein)):
                # To hold the atom indices corresponding to the titratable atoms            
                atomIndexList = []
                # Corresponding LambdaType object
                LambdaType = [obj for obj in self.ph_lambdaTypes if obj.d_groupname == residue.d_resname][0]

                # Loop through atoms - note that the atoms need to be descending order
                for atom in residue.d_atoms:
                    if atom in LambdaType.d_atoms:
                        atomIndexList.append(atomCount)

                    atomCount += 1

                # Write the lambda index group and increment groupnumber
                writeTheGroup(groupNumber, atomIndexList)
                groupNumber += 1

            elif constrainCharge and residue.d_resname == 'BUF':
                bufferIndexList.append(atomCount)
                atomCount += 1

            else: # Increment atomCount
                for atom in residue.d_atoms:
                    atomCount += 1

        # If we do charge constraining write the lambda index group for the buffer(s)
        if constrainCharge:
            writeTheGroup(groupNumber, bufferIndexList)

        file.close()

    # Generate parameters for MD.mdp.
    def genparams(self):
        # PART I - PREP

        # Load the input structure into d_residues.
        pdb = Structure(self.d_file, self.d_verbosity)

        # List of groupnames of the LambdaTypes specified in lambdagrouptypes.dat.
        LambdaTypeNames = []
        for LambdaType in self.ph_lambdaTypes:
            LambdaTypeNames.append(LambdaType.d_groupname)

        # Check whether we have any titratable residues in the structure, 
        # and also check whether we have any buffers.
        anyTitratables  = False
        constrainCharge = False

        for residue in pdb.d_residues:
            if residue.d_resname in LambdaTypeNames:
                anyTitratables  = True

            elif residue.d_resname == 'BUF':
                constrainCharge = True

            if anyTitratables and constrainCharge:
                break

        if not constrainCharge:
            self.warning("No BUF(s) found. Will not use charge constraining...")

        # If we are doing a calibration, set the BUF charges from 0 to 1 as this is
        # more straightforward for keeping net neutral charge.
        if self.ph_cal:
            self.ph_BUF_range = [1, 0]

        gen_mdp(Type='EM', nsteps=5000, nstxout=0, posRes=self.ph_cal)
        self.update('Wrote a generic EM.mdp file (for energy minimization)...')

        gen_mdp(Type='NVT', nsteps=5000, nstxout=0, posRes=self.ph_cal)
        self.update('Wrote a generic NPT.mdp file (for temperature coupling)...')

        val = self.inputOptionHandler("Simulating a membrane protein? (this will modify some barostat parameters)", ['No', 'Yes'])

        gen_mdp(Type='NPT', nsteps=5000, nstxout=0, membrane=val, posRes=self.ph_cal)
        self.update('Wrote a generic NVT.mdp file (for pressure coupling)...')

        # If no existing .mdp file is specified using the -mdp flag, write a generic one.
        if self.d_mdp == None:
            gen_mdp('MD', 50000, 5000, membrane=val, posRes=self.ph_cal)
            self.update('Wrote a generic MD.mdp file (for production)...')

        # Move this warning here so that we can also build normal MD sims with phbuilder.
        if not anyTitratables:
            self.error("No titratable residues detected! Will not perform any cpHMD related stuff...")

        # PARTS II - WRITE LAMBDA DYNAMICS PARAMETERS IN MDP FILE

        # We do not only MD, but also EM, NVT, and NPT at constant-pH.
        # Therefore loop over these files for PARTS I, II, and III (but not IV):
        for Type in ['EM', 'NVT', 'NPT', 'MD']:
            self.writeLambda_mdp(Type, pdb, LambdaTypeNames, constrainCharge)

        # PART III - WRITE LAMBDA INDEX GROUPS

        # If no .ndx file was specified on the command line, generate our generic one:
        if self.d_ndx == None:
            self.update('No .ndx file was specified. Creating a generic index.ndx file...')
            self.gromacs("make_ndx -f {}".format(self.d_file), stdin=['q'])
            self.d_ndx = 'index.ndx'

        self.writeLambda_ndx(self.d_ndx, pdb, LambdaTypeNames, constrainCharge)

        self.update('Finished phbuilder genparams. Please check the parameters in the generated .mdp files.')

    # Handle user input.
    def inputOptionHandler(self, message, options):

        valids = []
        msgstring = "phbuilder : {}:".format(message)

        # Loop through the options list and create string for display
        for idx in range(0, len(options)):
            msgstring += "\nphbuilder : {}. {}".format(idx, options[idx])
            valids.append(str(idx))

        while True:
            print(msgstring)
            val = input("phbuilder : Type a number: ")

            if val in valids:
                print()
                return int(val)

            print("phbuilder : {} is not a valid option, please try again:\n".format(val))

    # If the user does not specify the -inter flag, all the residues associated
    # with a lambdagrouptype will by default made titratable. The charge states
    # of these residues will be guessed based on an optional pH parameter
    # (default = 7.0) that can be specified for gentopol.
    def automaticLambdaInits(self, pKaList, systempH):
        # multistate case
        if len(pKaList) > 1:

            highest = 0
            for pKa in pKaList:
                if pKa <= systempH and pKa > highest:
                    highest = pKa

            return pKaList.index(highest) + 1

        #2state case
        else:
            if systempH >= pKaList[0]: 
                return 1 # if pH >= pKa, we are in deproto = 1 state.

            return 0 # if pH < pKa, we are in the proto = 0 state.

    def pleaseCite(self):
        self.update('If this software contributed to your research, please cite <doi_phbuilder_paper>.')
