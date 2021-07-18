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
    
    # If the operationmode was to us a list, load it.
    if (universe.get('ph_mode') == "list"):
        list_resid = universe.get('ph_list_resid')

    for residue in residues:
        
        if (residue.d_resname in lambdaTypeBaseNames):

            if (universe.get('ph_mode') == "all"):
                residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname)]

            if (universe.get('ph_mode') == "list"):
                if (residue.d_resid in list_resid):
                    residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname)]

            if (universe.get('ph_mode') == "interactive"):
                takeInput = True
                
                while (takeInput):

                    val = input("Choose what to do with residue {}-{} in chain {}:\n0. Keep current (static) protonation state\n1. Make dynamically protonatable (change name to {})\n\nType a number: ".format(
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
                        print("{} is not a valid option, please try again:\n".format(val))
            
            continue # Otherwise code below is executed for the same residue we just changed.

        if (residue.d_resname in lambdaTypeNames):

            if (universe.get('ph_mode') == "list"):
                if (residue.d_resid not in list_resid):
                    residue.d_resname = residue.d_resname[:-1]

            if (universe.get('ph_mode') == "interactive"):
                takeInput = True

                while (takeInput):

                    val = input("Choose what to do with residue {}-{} in chain {}:\n0. Make protonation state static (change name to {})\n1. Keep dynamically protonatable\n\nType a number: ".format(
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

    # Update d_residues and write output.
    universe.add('d_residues', residues)
    structure.write(universe.get('d_output'))

    # Part II - COPY FORCE FIELD AND RESIDUETYPES.DAT TO WORKING DIR

    d_modelFF    = universe.get('d_modelFF')
    d_modelwater = universe.get('d_modelwater')
    
    tail, head = os.path.split(d_modelFF)

    os.system("cp -r {} {}/residuetypes.dat .".format(d_modelFF, tail))

    utils.update('generate', 'full-path    = {}'.format(d_modelFF))
    utils.update('generate', 'tail-path    = {}'.format(tail))
    utils.update('generate', 'residuetypes = {}/residuetypes.dat'.format(tail))
    utils.update('generate', 'head-path    = {}'.format(head))

    d_modelFF = head[0:len(head)-3]
    utils.update('generate', 'ffield name  = {}'.format(d_modelFF))

    # Part III - RUN PDB2GMX

    os.system("gmx -nocopyright pdb2gmx -f {0} -o {0} -ff {1} -water {2} -ignh {3}".format(
        universe.get('d_output'), 
        d_modelFF,
        d_modelwater,
        universe.get('d_options')))

    # PART IV - WRAPUP

    # Update d_residues.
    structure.load(universe.get('d_output'))
