import os, structure, universe, utils, mdp
from types import LambdaType

def genparams():
    # PART I - PREP

    # Load the input structure into d_residues.
    structure.load(universe.get('d_file'))

    # Load d_residues here for speedup as we need it in multiple loops.
    residues = universe.get('d_residues')

    # Load ph_lambdaTypes here for speedup as we need it in multiple loops.
    LambdaTypes = universe.get('ph_lambdaTypes')

    # List of groupnames of the LambdaTypes specified in lambdagrouptypes.dat.
    LamdaTypeNames = []
    for LambdaType in LambdaTypes:
        LamdaTypeNames.append(LambdaType.d_groupname)

    # Check whether we have any titratable residues in the structure, and also
    # check whether we have any buffers.
    anyTitratables = False
    restrainCharge = False

    for residue in residues:
        if residue.d_resname in LamdaTypeNames:
            anyTitratables = True

        elif residue.d_resname == 'BUF':
            restrainCharge = True

        if anyTitratables and restrainCharge:
            break

    if not anyTitratables:
        utils.error("No titratable residues detected")

    if not restrainCharge:
        utils.update("No buffer(s) found. Will not use charge restraining...")

    # If no .mdp file was specified on the command line, generate our default one:
    if universe.get('d_mdp') == None:
        utils.update('No .mdp file was specified. Creating a default MD.mdp file...')
        mdp.gen_mdp('MD', 50000, 5000, restrainCharge)
        universe.add('d_mdp', 'MD.mdp')

    # If no .ndx file was specified on the command line, generate our default one:
    if universe.get('d_ndx') == None:
        utils.update('No .ndx file was specified. Creating a default index.ndx file...')
        os.system("gmx make_ndx -f {0} >> builder.log 2>&1 << EOF\nq\nEOF".format(universe.get('d_file')))

    file = open(universe.get('d_mdp'), 'a')

    # Formatting function for adding parameters.
    def addParam(name, value):
            file.write("{:54s} = {:13s}\n".format(name, str(value)))

    # PART 1 - WRITE GENERAL PARAMETERS

    file.write("\n; CONSTANT PH\n")

    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', "{:.1f}".format(universe.get('ph_ph')))
    addParam('lambda-dynamics-lambda-particle-mass', "{:.1f}".format(universe.get('ph_lmass')))
    addParam('lambda-dynamics-tau', "{:.1f}".format(universe.get('ph_ltau')))
    addParam('lambda-dynamics-update-nst', universe.get('ph_nstout'))

    # If we use charge restraining...
    if restrainCharge:
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # We need to count how many titratable residues in total we have in the
    # protein. For this we compile a list LambdasFoundinProtein.
    LambdasFoundinProtein = [] # (e.g ASPT ASPT GLUT ASPT GLUT ASPT...)
    
    # Stores the number of buffer atoms/ions.
    buffersFoundinProtein = 0

    for residue in residues:
        if residue.d_resname in LamdaTypeNames:
            LambdasFoundinProtein.append(residue.d_resname)

        # Also in this loop we count how many buffer ions we have.
        elif restrainCharge and residue.d_resname == 'BUF':
            buffersFoundinProtein += 1

    # The fact that we have a LambdaType in lambdagrouptypes.dat does not mean
    # one of those is also present in the protein. In that case, we want to 
    # prevent counting it, so we compile a subgroup of LambdaType groupnames
    # that are not only in lambdagrouptypes.dat, but ALSO found at least once 
    # in the actual protein.
    LambdaTypeNamesFoundinProtein = list(set(LambdasFoundinProtein)) # (e.g ASPT GLUT)

    # If we use the charge leveling scheme "charge-restraining" (2) we also have
    # the BUF residue-type as well as one extra lambda group containing all the BUFs.
    if restrainCharge:
        addParam('lambda-dynamics-number-lambda-group-types', len(LambdaTypeNamesFoundinProtein) + 1)
        addParam('lambda-dynamics-number-atom-collections', len(LambdasFoundinProtein) + 1)
    else:
        addParam('lambda-dynamics-number-lambda-group-types', len(LambdaTypeNamesFoundinProtein))
        addParam('lambda-dynamics-number-atom-collections', len(LambdasFoundinProtein))

    file.write('\n')

    # PART 2 - WRITE LAMBDA GROUP TYPES
    
    # Convert a list to a string
    def to_string(Input):
        string = ""
        for element in Input:
            string += "{:.3f} ".format(element)
        return string

    # Writes the lambda group type block
    def writeBlock(number, name, dvdl, pKa, barrierE, qqA, qqB):
        addParam('lambda-dynamics-group-type{}-name'.format(number), name)
        addParam('lambda-dynamics-group-type{}-dvdl-coefficients'.format(number), to_string(dvdl))
        addParam('lambda-dynamics-group-type{}-reference-pka'.format(number), pKa)
        addParam('lambda-dynamics-group-type{}-barrier'.format(number), barrierE)
        addParam('lambda-dynamics-group-type{}-charges-state-A'.format(number), to_string(qqA))
        addParam('lambda-dynamics-group-type{}-charges-state-B'.format(number), to_string(qqB))
        file.write('\n')

    number = 1
    # We loop over the object itself instead of the d_groupname as we need all
    # the information in the object.
    for LambdaType in LambdaTypes:
        # This if-statement prevents writing a block when there are no residues 
        # of this type in the protein.
        if (LambdaType.d_groupname in LambdaTypeNamesFoundinProtein):
            writeBlock(number, 
                       LambdaType.d_groupname, 
                       LambdaType.d_dvdl[::-1],
                       LambdaType.d_pKa,
                       universe.get('ph_dwpE'),
                       LambdaType.d_qqA,
                       LambdaType.d_qqB)
            number += 1

    # If we do charge restraining, we additionally need the block for the buffer.
    if (restrainCharge):
        writeBlock(number,
                   'BUF',
                   universe.get('ph_BUF_dvdl')[::-1],
                   0,
                   0,
                   [1.0],
                   [0.0])

    # PART 3 - WRITE LAMBDA GROUPS
   
    def writeResBlock(number, name, indexName):
        addParam('lambda-dynamics-atom-set{}-name'.format(number), name)
        addParam('lambda-dynamics-atom-set{}-index-group-name'.format(number), indexName)
        addParam('lambda-dynamics-atom-set{}-initial-lambda'.format(number), 0.5)

        if restrainCharge:
            addParam('lambda-dynamics-atom-set{}-charge-restraint-group-index'.format(number), 1)

        if (name == 'BUF'):
            addParam('lambda-dynamics-atom-set{}-buffer-residue'.format(number), 'yes')
            addParam('lambda-dynamics-atom-set{}-buffer-residue-multiplier'.format(number), buffersFoundinProtein)

        file.write('\n')

    number = 1
    for groupname in LambdasFoundinProtein:
        writeResBlock(number, groupname, 'LAMBDA{}'.format(number))
        number += 1

    if (restrainCharge):
        writeResBlock(number, 'BUF', 'LAMBDA{}'.format(len(LambdasFoundinProtein) + 1))

    file.close() # MD.mdp

    # PART 4 - WRITE LAMBBDA INDEX GROUPS

    # Append to existing index.ndx
    file = open('index.ndx', 'a')

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
    for residue in residues:
        # If the residue is titratable
        if residue.d_resname in LambdaTypeNamesFoundinProtein:
            # To hold the atom indices corresponding to the titratable atoms            
            atomIndexList = []
            # Corresponding LambdaType object
            LambdaType = [obj for obj in LambdaTypes if obj.d_groupname == residue.d_resname][0]

            # Loop through atoms - note that the atoms need to be descending order
            for atom in residue.d_atoms:
                if atom in LambdaType.d_atoms:
                    atomIndexList.append(atomCount)

                atomCount += 1
            
            # Write the lambda index group and increment groupnumber
            writeTheGroup(groupNumber, atomIndexList)
            groupNumber += 1

        elif restrainCharge and residue.d_resname == 'BUF':
            bufferIndexList.append(atomCount)
            atomCount += 1

        else: # Increment atomCount
            for atom in residue.d_atoms:
                atomCount += 1

    # If we do charge restraining write the lambda index group for the buffer(s)
    if restrainCharge:
        writeTheGroup(groupNumber, bufferIndexList)

    file.close() # index.ndx
