import structure, universe

def gentopol():
    # PART I - MODIFIY THE STRUCTURE FILE AND WRITE ########################################################################################################
        
    # Load the structure into d_residues.
    structure.load(universe.get('d_file'))

    # Do stuff to .pdb file here.

    structure.write(universe.get('d_output'))

    # Part II - REGENERATE THE TOPOLOGY ####################################################################################################################

    