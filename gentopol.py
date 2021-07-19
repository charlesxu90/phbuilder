import os, structure, universe, utils

def gentopol():
    # PART I - MODIFIY THE STRUCTURE FILE AND WRITE

    # Load the structure into d_residues.
    structure.load(universe.get('d_file'))

    lambdaTypeNames     = [] # Stores e.g. ASPT, XY
    lambdaTypeBaseNames = [] # Stores e.g. ASP, X
    
    for lambdaType in universe.get('ph_lambdaTypes'):
        lambdaTypeNames.append(lambdaType.d_groupname)
        lambdaTypeBaseNames.append(lambdaType.d_groupname[:-1])

    # Load residues here as references do not work for shelve module.
    residues = universe.get('d_residues')
    
    # If a list of specific residues was specified, load it.
    if universe.has('ph_list_resid'):
        list_resid = universe.get('ph_list_resid')

    for residue in residues:
        
        if (residue.d_resname in lambdaTypeBaseNames):

            # If -list was set AND the current resid is not in list, continue to next residue in loop.
            if universe.has('ph_list_resid') and (residue.d_resid not in list_resid):
                continue

            # If -inter was set, ask for every residue individually:
            if universe.has('ph_inter'):
                takeInput = True
                
                while (takeInput):

                    val = input("phbuilder : Choose what to do with residue {}-{} in chain {}:\nphbuilder : 0. Keep current (static) protonation state\nphbuilder : 1. Make dynamically protonatable (change name to {})\nphbuilder : Type a number: ".format(
                        residue.d_resname,
                        residue.d_resid,
                        residue.d_chain,
                        lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname)]))
                    
                    if (val == '0'):
                        takeInput = False
                    
                    elif (val == '1'):
                        residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname)]
                        takeInput = False
                    
                    else:
                        print("{} is not a valid option, please try again:".format(val))
                    
                    print()

            else: # If -inter was not set, just protonate everything automatically.
                residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname)]

            continue # Otherwise code below is executed for the same residue we just changed.

        if (residue.d_resname in lambdaTypeNames):

            # If -list was set AND the current resid is not in list, continue to next residue in loop.
            if universe.has('ph_list_resid') and (residue.d_resid not in list_resid):
                continue

            # If -inter was set, ask for every residue individually:
            if universe.has('ph_inter'):
                takeInput = True

                while (takeInput):

                    val = input("phbuilder : Choose what to do with residue {}-{} in chain {}:\nphbuilder : 0. Make protonation state static (change name to {})\nphbuilder : 1. Keep dynamically protonatable\nphbuilder : Type a number: ".format(
                        residue.d_resname,
                        residue.d_resid,
                        residue.d_chain,
                        lambdaTypeBaseNames[lambdaTypeNames.index(residue.d_resname)]))

                    if (val == '0'):
                        residue.d_resname = lambdaTypeBaseNames[lambdaTypeNames.index(residue.d_resname)]
                        takeInput = False
                    
                    elif (val == '1'):
                        takeInput = False
                    
                    else:
                        print("{} is not a valid option, please try again:\n".format(val))

                    print()

    # Update d_residues and write output.
    universe.add('d_residues', residues)
    structure.write(universe.get('d_output'))

    # Part II - COPY FORCE FIELD AND RESIDUETYPES.DAT TO WORKING DIR

    d_modelFF    = universe.get('d_modelFF')
    d_modelwater = universe.get('d_modelwater')
    
    tail, head = os.path.split(d_modelFF)

    os.system("cp -r {} {}/residuetypes.dat .".format(d_modelFF, tail))

    utils.update('Force field path stuff: ', 3)
    utils.update('full-path    = {}'.format(d_modelFF), 3)
    utils.update('tail-path    = {}'.format(tail), 3)
    utils.update('residuetypes = {}/residuetypes.dat'.format(tail), 3)
    utils.update('head-path    = {}'.format(head), 3)

    d_modelFF = head[0:len(head)-3]
    utils.update('ffield name  = {}'.format(d_modelFF), 3)

    # Part III - RUN PDB2GMX

    utils.update("Recommended pdb2gmx command:")
    utils.update("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh".format(universe.get('d_output'), d_modelFF, d_modelwater))

    takeInput = True

    while (takeInput):

        val = input("phbuilder : Choose whether to:\nphbuilder : 0. Do nothing\nphbuilder : 1. Run\nphbuilder : 2. Add additional flags (see https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html)\nphbuilder : Type a number: ")

        if val not in ['0', '1', '2']:
            utils.update("{} is not a valid option, please try again:\n".format(val))
            continue
        else:
            takeInput = False

        if (val == '2'):
            flags = input("phbuilder : Enter flags: ")
        else:
            flags = ''

        if (val in ['1', '2']):
            utils.update("Calling GROMACS ({}/gmx)...\n".format(os.environ.get("GMXBIN")))
            os.system("gmx pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh {3}".format(universe.get('d_output'), d_modelFF, d_modelwater, flags))

    # PART IV - WRAPUP

    # Update d_residues.
    structure.load(universe.get('d_output'))
