import os, structure, universe, utils

def gentopol():
    # Part I - COPY FORCE FIELD AND RESIDUETYPES.DAT TO WORKING DIR

    d_modelFF    = universe.get('d_modelFF')
    d_modelwater = universe.get('d_modelwater')

    tail, head = os.path.split(d_modelFF)

    os.system("cp -r {} {}/residuetypes.dat .".format(d_modelFF, tail))

    utils.update('Force field path stuff: ')
    utils.update('full-path    = {}'.format(d_modelFF))
    utils.update('tail-path    = {}'.format(tail))
    utils.update('residuetypes = {}/residuetypes.dat'.format(tail))
    utils.update('head-path    = {}'.format(head))

    d_modelFF = head[0:len(head)-3]
    utils.update('ffield name  = {}'.format(d_modelFF))
    
    # PART I - MODIFIY THE STRUCTURE FILE AND WRITE

    utils.update("Modifying structure file...")

    # Load the structure into d_residues.
    structure.load(universe.get('d_file'))

    # List of titratable names, e.g. [ASPT, GLUT].
    lambdaTypeNames = [] 
    
    # List of lists of all the incl names for the LambdaTypes,
    # e.g. [[ASP1, ASPH, ASPP, ASH], [GLU, GLUH, GLUP]].
    residuesToBeConsidered = []

    for lambdaType in universe.get('ph_lambdaTypes'):
        lambdaTypeNames.append(lambdaType.d_groupname)
        residuesToBeConsidered.append(lambdaType.d_incl)

    # Checks whether the resname is in any of the lists of targets.
    # If found, return index of said list. If not found, return None.
    def isTarget(resname):
        for idx in range(0, len(residuesToBeConsidered)):
            if resname in residuesToBeConsidered[idx]:
                return idx
        return None

    # Load residues here as changing by reference does not work for stuff stored in shelve.
    residues = universe.get('d_residues')

    # If a list of specific residues was specified, load it.
    if universe.has('ph_list_resid'):
        list_resid = universe.get('ph_list_resid')
        
        # User udate
        utils.update('Found user-defined list of which residues to consider:')
        utils.update(list_resid)

    # Loop through the residue objects.
    for residue in residues:

        # This will either be None if the residue does not have to be considered 
        # (i.e. it doesn't belong to any incl group), or it will have the index of 
        # the incl grou it belongs to.
        index = isTarget(residue.d_resname)
        
        # If we find a residue that was specified to be titratable...
        if (index != None):
            
            # Original name of the residue before we change it. We store this
            # here because we'll need it later for a user update.
            origName = residue.d_resname

            # If -list was set AND the current resid is not in list, continue to next residue in loop above.
            if universe.has('ph_list_resid') and (residue.d_resid not in list_resid):
                continue

            # If -inter was set, ask for every residue individually:
            if universe.has('ph_inter'):

                # Handle the user input
                takeInput = True

                while (takeInput):

                    val = input("phbuilder : Choose what to do with residue {}-{} in chain {}:\nphbuilder : 0. Keep current (static) protonation state\nphbuilder : 1. Make titratable (change name to {})\nphbuilder : Type a number: ".format(
                        residue.d_resname,
                        residue.d_resid,
                        residue.d_chain,
                        lambdaTypeNames[index]))

                    if (val == '0'):
                        takeInput = False

                    elif (val == '1'):
                        residue.d_resname = lambdaTypeNames[index]
                        takeInput = False

                    else:
                        print("{} is not a valid option, please try again:".format(val))

                    print()
                # End of input handling

            else: # If -inter was not set, just protonate everything automatically.
                residue.d_resname = lambdaTypeNames[index]
                utils.update("Made residue {}-{} in chain {} titratable (changed name to {})".format(origName, residue.d_resid, residue.d_chain, residue.d_resname))

            continue # Otherwise code below is executed for the same residue we just changed.

        # Again loop through the residue objects, but now we look for ones that were already changed.
        if (residue.d_resname in lambdaTypeNames):

            # If -list was set AND the current resid is not in list, continue to next residue in loop.
            if universe.has('ph_list_resid') and (residue.d_resid not in list_resid):
                continue

            # If -inter was set, ask for every residue individually:
            if universe.has('ph_inter'):

                # Handle the user input
                takeInput = True

                while (takeInput):
                    
                    # Some user updates:
                    utils.update("Choose what to do with residue {}-{} in chain {}:".format(residue.d_resname, residue.d_resid, residue.d_chain))
                    utils.update("0. Keep titratable (do not change name)")
                    
                    for idx in range(0, len(residuesToBeConsidered[lambdaTypeNames.index(residue.d_resname)])):
                        utils.update("{}. change name to {} (no longer titratable)".format(idx + 1, residuesToBeConsidered[lambdaTypeNames.index(residue.d_resname)][idx]))

                    val = eval(input("Type a number: "))

                    # Handle incorrect input:
                    if (type(val) != int or val not in range(0, len(residuesToBeConsidered[lambdaTypeNames.index(residue.d_resname)]))):
                        utils.update("{} is not a valid option, please try again:\n".format(val))
                        continue
                    else:
                        takeInput = False

                    # Process the specified option:
                    if (val == 0):
                        continue
                    else:
                        residue.d_resname = residuesToBeConsidered[lambdaTypeNames.index(residue.d_resname)][val-1]

                    print()
                # End of input handling

    # Update d_residues.
    universe.add('d_residues', residues)

    # Write structure.
    structure.write(universe.get('d_output'))

    # PART III - HANDLE UNKNOWN STRUCTURE TYPES

    utils.update("Checking if every residue type is present in residuetypes.dat...")

    # Load residuetypes.dat into list
    residueTypes = []
    for val in open('residuetypes.dat').readlines():
        residueTypes.append(val.split()[0])

    # Compile lists of the unknown residue(s)(types)
    unknownResidues = []
    unknownResTypeNames = []
    for residue in universe.get('d_residues'):
        
        if residue.d_resname not in residueTypes:
            unknownResidues.append(residue)

            if residue.d_resname not in unknownResTypeNames:
                unknownResTypeNames.append(residue.d_resname)

    pathList = []   # Loop through the unknown residue types and add them to 
    skipList = []   # skipList and (the manually specified path to) pathList.
    good = True
    for val in unknownResTypeNames:
        utils.update("WARNING - residue type {} in {} wasn't found in residuetypes.dat associated with {}".format(val, universe.get('d_file'), universe.get('d_modelFF')))
        path = input("phbuilder : Please specify path to .itp file (e.g. /path/to/some.itp): ")

        skipList.append(val)
        pathList.append(path)
        good = False

        utils.update("Set path for residue type {} to {}...".format(val, path))

    if good:
        utils.update("everything seems OK.")

    # If pathList is not empty, i.e. if we had at least one unknown residue
    if pathList:

        # Create a list knownResidues containing only the known residues
        knownResidues = []
        for residue in residues:
            if residue.d_resname not in skipList:
                knownResidues.append(residue)

        # Update d_residues in universe with knownResidues
        universe.add('d_residues', knownResidues)

        # Write temporary .pdb containing only the known residues
        someTempName = 'pdb2gmxtemp.pdb'
        structure.write(someTempName)

        # Update value of d_output (so that pdb2gmx is called on the temporary 
        # structure), and backup the final output name as specified by user.
        universe.add('d_output_orig', universe.get('d_output'))
        universe.add('d_output', someTempName)

    # PART IV - RUN PDB2GMX AND ASK USER FOR INPUT ABOUT IT 

    utils.update("\nRecommended pdb2gmx command:")
    utils.update("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh".format(universe.get('d_output'), d_modelFF, d_modelwater))

    takeInput = True
    while (takeInput):

        # Prompt user for input:
        val = input("phbuilder : Choose whether to:\nphbuilder : 0. Do nothing\nphbuilder : 1. Run\nphbuilder : 2. Add additional flags (https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html)\nphbuilder : Type a number: ")

        # Handle incorrect input:
        if val not in ['0', '1', '2']:
            utils.update("{} is not a valid option, please try again:\n".format(val))
            continue
        else:
            takeInput = False

        # Take additional flags:
        if (val == '2'):
            flags = input("phbuilder : Enter flags: ")
        else:
            flags = ''

        # Run pdb2gmx:
        if (val in ['1', '2']):
            utils.update("Calling GROMACS ({}/gmx)...\n".format(os.environ.get("GMXBIN")))
            os.system("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh {3}".format(universe.get('d_output'), d_modelFF, d_modelwater, flags))

    # PART V - MERGE THE .PDB FILES
    
    # If pathList is not empty, i.e. if we had at least one unknown residue:
    if pathList:
        # Load the structure output of pdb2gmx (update d_residues in universe)
        structure.load(someTempName)

        # Merge the processed structure of good residues with unknown residues
        mergedResidues = universe.get('d_residues') + unknownResidues

        # Update d_residues in universe
        universe.add('d_residues', mergedResidues)

        # Write the final structure
        structure.write(universe.get('d_output_orig'))

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

                    if "[ system ]\n" == topList[line + 1]:
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
            for residue in residues:
                if residue.d_resname == skipList[idx]:
                    count += 1

            # Add manually to topol.top
            add_mol(pathList[idx], "Include topology for {}".format(skipList[idx]), skipList[idx], count)

    utils.update("Finished generating topology for constant-pH.")
