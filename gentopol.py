import structure, universe

def gentopol():
    # PART I - MODIFIY THE STRUCTURE FILE AND WRITE ########################################################################################################

    # Load the structure into d_residues.
    structure.load(universe.get('d_file'))

    lambdaTypeNames     = [] # Stores e.g. ASPT
    lambdaTypeBaseNames = [] # Stores e.g. ASP
    
    for lambdaType in universe.get('ph_lambdaTypes'):
        lambdaTypeNames.append(lambdaType.d_groupname)
        lambdaTypeBaseNames.append(lambdaType.d_groupname[0:3])

    # print("lambdaTypeNames", lambdaTypeNames)         # debug
    # print("LambdaTypeBaseNames", lambdaTypeBaseNames) # debug

    residues = universe.get('d_residues')
    
    for residue in residues:
        if (residue.d_resname)[0:3] in lambdaTypeBaseNames:

            if (universe.get('ph_mode') == "all"):
                residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname[0:3])]

            if (universe.get('ph_mode') == "interactive"):
                takeInput = True
                
                while (takeInput):
                
                    val = input("Choose what to do with residue {}-{} in chain {}:\n0. Keep current (static) protonation state\n1. Make dynamically protonatable (resname will be changed to {})\n\nType a number: ".format(
                        residue.d_resname[0:3],
                        residue.d_resid,
                        residue.d_chain,
                        lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname[0:3])]))

                    if (val == '0'):
                        residue.d_resname = residue.d_resname[0:3]
                        takeInput = False
                    
                    elif (val == '1'):
                        residue.d_resname = lambdaTypeNames[lambdaTypeBaseNames.index(residue.d_resname[0:3])]
                        takeInput = False
                    
                    else:
                        print("{} is not a valid option, please try again:\n".format(val))

    # Update d_residues and write output.
    universe.add('d_residues', residues)
    structure.write(universe.get('d_output'))

    # Part II - REGENERATE THE TOPOLOGY ####################################################################################################################

    