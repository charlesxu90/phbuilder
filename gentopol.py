import os, structure, universe, utils

def gentopol():
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

    # Update d_residues and write output.
    universe.add('d_residues', residues)
    structure.write(universe.get('d_output'))

    # Part II - COPY FORCE FIELD AND RESIDUETYPES.DAT TO WORKING DIR

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

    # Part III - CHECK WHETHER PDB2GMX KNOWS WHAT TO DO WITH THE INPUT

    utils.update("\nChecking if every residue name is present in ./residuetypes.dat...")

    restypes = []
    for val in open('residuetypes.dat').readlines():
        restypes.append(val.split()[0])

    unknownResidues = []
    for residue in universe.get('d_residues'):
        if residue.d_resname not in restypes:
            unknownResidues.append(residue.d_resname)

    good = True
    for val in list(set(unknownResidues)):
        utils.update("WARNING - residue {} in {} was not found in ./residuetypes.dat".format(val, universe.get('d_file')))
        good = False

    if good:
        utils.update("everything seems OK.")
    else:
        utils.update("pdb2gmx will not be able to properly process your structure. You need to either: ")
        utils.update("1. Update the path to your force field in lambdagrouptypes.dat")
        utils.update("2. Modify your force field to accomodate the unknown residue type")
        utils.update("3. Cut out the unknown parts of your structure, run pdb2gmx, and then paste those")
        utils.update("   parts back (whilst manually updating your #includes in topol.top accordingly)")

    utils.update("\nRecommended pdb2gmx command:")
    utils.update("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh".format(universe.get('d_output'), d_modelFF, d_modelwater))

    # PART IV - RUN PDB2GMX AND ASK USER FOR INPUT ABOUT IT 

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
            os.system("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh {3}".format(universe.get('d_output'), d_modelFF, d_modelwater, flags)) # Do this with gmxapi in the future?

            structure.load(universe.get('d_output')) # Update d_residues.
            utils.update("Finished generating topology for constant-pH.")
