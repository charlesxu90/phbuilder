import os, numpy as np, structure, universe, utils

def addbuffers():
    maxChargeOnBuffer = 0.3
    
    # PART I - PREP

    # Load the input structure into d_residues.
    structure.load(universe.get('d_file'))

    # Load d_residues here for speedup as we need it in multiple loops.
    residues = universe.get('d_residues')

    # Perform some basic checks for PBC box, solvent, and titratable groups.
    foundPBC         = universe.has('d_box')
    foundSolvent     = False
    foundTitratables = False

    LambdaTypeNames = []
    for LambdaType in universe.get('ph_lambdaTypes'):
        LambdaTypeNames.append(LambdaType.d_groupname)

    for residue in residues:
        
        if residue.d_resname == universe.get('ph_solname'):
            foundSolvent = True
        
        elif residue.d_resname in LambdaTypeNames:
            foundTitratables = True

        if foundSolvent and foundTitratables:
            break

    # Error if no periodic box was found
    if not foundPBC:
        utils.error("{} doesn't have a periodic box! Did you forget to add one?".format(universe.get('d_file')))

    # Error if no solvent was found
    if not foundSolvent:
        utils.error("{} doesn't seem to have any {} molecules! Did you forget to add solvent?".format(universe.get('d_file'), universe.get('ph_solname')))

    # Error if no titratable residues were found
    if not foundTitratables:
        utils.error("{} doesn't seem to have any titratable groups! Did you forget to run gentopol?".format(universe.get('d_file'), universe.get('ph_solname')))

    # if nbufs wasn't manually set, count the number of titratable residues.
    if not universe.has('ph_nbufs'):
        # Count number of titratable residues
        titratables = 0
        for residue in residues:
            if residue.d_resname in LambdaTypeNames:
                titratables += 1

        nbufs = int(np.ceil(titratables / maxChargeOnBuffer)) # Worst case scenario.
        utils.update("Counted {0} titratable residues, will add ceil({0}/{1}) = {2} buffer(s)...".format(titratables, maxChargeOnBuffer, nbufs))

    else:
        nbufs = universe.get('ph_nbufs')
        utils.update("Will add {} buffer(s) (per user request)...".format(nbufs))

    # PART II - RUN GMX GENION TO REPLACE SOLVENT MOLECULES WITH BUFFER

    # Create dummy mdp file for gmx grompp
    open('dummy.mdp', 'w').close()

    # Run gmx grompp to create a .tpr file for gmx genion
    os.system("gmx grompp -f dummy.mdp -c {} -p {} -o dummy.tpr >> builder.log 2>&1".format(universe.get('d_file'), universe.get('d_topol')))

    # Run gxm genion to replace some solvent molecules with buffers
    os.system("gmx genion -s dummy.tpr -p {} -o {} -pname {} -np {} >> builder.log 2>&1 << EOF\n{}\nEOF".format(
        universe.get('d_topol'),
        universe.get('d_output'),
        'BUF', 
        nbufs,
        universe.get('ph_solname')))

    # PART III - WRAPUP
    
    # Remove dummy files
    os.remove('dummy.tpr'); os.remove('dummy.mdp')
    
    # Update d_residues in universe
    structure.load(universe.get('d_output'))

    # Check whether the required number of buffers was added
    BUFcount = 0
    for residue in universe.get('d_residues'):
        if residue.d_resname == 'BUF':
            BUFcount += 1
    
    if BUFcount == nbufs:
        utils.update("Succesfully added {} buffer molecule(s)...".format(BUFcount))
    else:
        utils.update("Only added {}/{} buffer molecules! Try increasing box size...")
