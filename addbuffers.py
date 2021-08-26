import os, structure, universe, utils

def addbuffers():
    # PART I - PREP

    # Load the input structure into d_residues.
    structure.load(universe.get('d_file'))

    # Load d_residues here for speedup as we need it in multiple loops.
    residues = universe.get('d_residues')

    # Error if no periodic box was found
    if not universe.has('d_box'):
        utils.error("{} doesn't have a periodic box! Did you forget to add one?".format(universe.get('d_file')))

    # Error if no solvent was found
    foundSolvent = False
    for residue in residues:
        if residue.d_resname == universe.get('ph_solname'):
            foundSolvent = True
            break

    if not foundSolvent:
        utils.error("{} doesn't seem to have any {} molecules! Did you forget to add solvent?".format(universe.get('d_file'), universe.get('ph_solname')))

    # if nbufs wasn't manually set, count the number of titratable residues.
    if not universe.has('ph_nbufs'):
        # Compile list of lambdaType groupnames
        lambdaTypeNames = []
        for lambdaType in universe.get('ph_lambdaTypes'):
            lambdaTypeNames.append(lambdaType.d_groupname)

        # Count number of titratable residues
        tits = 0
        for residue in residues:
            if residue.d_resname in lambdaTypeNames:
                tits += 1

        utils.update("Counted {0} titratable residues, will add {0} buffer(s)...".format(tits))

    else:
        tits = universe.get('ph_nbufs')
        utils.update("Will add {} buffer(s)...".format(tits))

    # PART II - RUN GMX GENION TO REPLACE SOLVENT MOLECULES WITH BUFFER

    # Create dummy mdp file for gmx grompp
    open('dummy.mdp', 'w').close()

    # Run gmx grompp to create a .tpr file for gmx genion
    os.system("gmx grompp -f dummy.mdp -c {} -p {} -o dummy.tpr".format(universe.get('d_file'), universe.get('d_topol')))

    # Run gxm genion to replace some solvent molecules with buffers
    os.system("gmx genion -s dummy.tpr -p {} -o {} -pname {} -np {} -rmin {} << EOF\n{}\nEOF".format(
        universe.get('d_topol'),
        universe.get('d_output'),
        'BUF', 
        tits, 
        universe.get('ph_bufmargin'),
        universe.get('ph_solname')))

    # PART III - WRAPUP
    utils.update("Succesfully added {} buffer molecule(s)".format(tits))

    # Remove dummy files
    os.remove('dummy.tpr'); os.remove('dummy.mdp')

    # Update d_residues in universe
    structure.load(universe.get('d_output'))
