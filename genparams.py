import os, structure, universe, utils, mdp
from types import LambdaType

def genparams():
    # PART I - PREP

    # Load the input structure into d_residues.
    structure.load(universe.get('d_file'))

    # Load d_residues here for speedup as we need it in multiple loops.
    residues = universe.get('d_residues')

    # Check whether we have any titratable residues in the structure, and also
    # check whether we have any buffers.
    anyTitratables = False
    restrainCharge = False
    
    lambdaTypeNamesSpecifed = []
    for LambdaType in universe.get('ph_lambdaTypes'):
        lambdaTypeNamesSpecifed.append(LambdaType.d_groupname)

    for residue in residues:
        if residue.d_resname in lambdaTypeNamesSpecifed:
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

    # PART 1 - WRITE GENERAL PARAMETERS ########################################

    file.write("\n; CONSTANT PH\n")

    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', "{:.1f}".format(universe.get('ph_ph')))
    addParam('lambda-dynamics-lambda-particle-mass', "{:.1f}".format(universe.get('ph_lmass')))
    addParam('lambda-dynamics-tau', "{:.1f}".format(universe.get('ph_ltau')))
    addParam('lambda-dynamics-update-nst', universe.get('ph_nstout'))

    # If we use charge restraining...
    if restrainCharge:
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # Gather a list of the names of the lambda residues in the protein.
    lambdaResidueNameList = []
    for residue in residues:
        if (residue.d_resname in lambdaTypeNamesSpecifed):
            lambdaResidueNameList.append(residue.d_resname)

    # Gather a list of the names of the lambda residue-types IN THE PROTEIN,
    # And make sure the order is the same as in ph_lambdaTypes = LambdaTypeNamesSpecified.
    lambdaResidueTypeList = []
    for obj in universe.get('ph_lambdaTypes'):
        if (obj.d_groupname in set(lambdaResidueNameList)):
            lambdaResidueTypeList.append(obj.d_groupname)

    # If we use the charge leveling scheme "charge-restraining" (2) we also have
    # the BUF residue-type as well as one extra lambda group containing all the BUFs.
    if restrainCharge:
        addParam('lambda-dynamics-number-lambda-group-types', len(lambdaResidueTypeList) + 1)
        addParam('lambda-dynamics-number-atom-collections', len(lambdaResidueNameList) + 1)
    else:
        addParam('lambda-dynamics-number-lambda-group-types', len(lambdaResidueTypeList))
        addParam('lambda-dynamics-number-atom-collections', len(lambdaResidueNameList))

    file.write('\n')

    # # PART 2 - WRITE LAMBDA GROUP TYPES
    
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
    for LambdaType in universe.get('ph_lambdaTypes'):
        # This if-statement prevents writing a block when there are no residues of this type.
        if (LambdaType.d_groupname in lambdaResidueTypeList):
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
                   1.0,
                   0.0)

    # PART 3 - WRITE LAMBDA GROUPS
   
    def writeResBlock(number, name, indexName):
        addParam('lambda-dynamics-atom-set{}-name'.format(number), name)
        addParam('lambda-dynamics-atom-set{}-index-group-name'.format(number), indexName)
        addParam('lambda-dynamics-atom-set{}-initial-lambda'.format(number), 0.5)

        if restrainCharge:
            addParam('lambda-dynamics-atom-set{}-charge-restraint-group-index'.format(number), 1)

        if (name == 'BUF'):
            addParam('lambda-dynamics-atom-set{}-buffer-residue'.format(number), 'yes')
            addParam('lambda-dynamics-atom-set{}-buffer-residue-multiplier'.format(number), universe.get('ph_bufnmol'))

        file.write('\n')

    number = 1
    for name in lambdaResidueNameList:
        writeResBlock(number, name, 'LAMBDA{}'.format(number))
        number += 1

    if (restrainCharge):
        writeResBlock(number, 'BUF', 'LAMBDA{}'.format(len(lambdaResidueNameList) + 1))

    file.close() # MD.mdp

    # PART 4 - WRITE LAMBBDA INDEX GROUPS

    # If we use do charge restraining, we'll need a list of atomIndices of the BUFs:
    if restrainCharge:
        bufferAtomIndexList = []

        atomidx = 1
        for residue in residues:
            for _ in residue.d_atoms:
                if residue.d_resname == 'BUF':
                    bufferAtomIndexList.append(atomidx)
                atomidx += 1

    # Append to existing index.ndx
    file = open('index.ndx', 'a')

    # Formatting file for writing the block.
    def writeTheGroup(number, atomIndexList):
        file.write('\n[ LAMBDA{} ]\n'.format(number))
        for index in atomIndexList:
            file.write('{} '.format(index))
        file.write('\n')

    # ph_lambdaTypes = universe.get('ph_lambdaTypes')

    # atomCount = 1; groupNumber = 1
    # for residue in universe.get('d_residues'):
    #     if residue.d_resname in lambdaResidueTypeList:

    #         atomIndexList = []
    #         obj = ph_lambdaTypes[lambdaTypeNamesSpecifed.index(residue.d_resname)]

    #         for atom in residue.d_atoms:
    #             if atom in obj.d_atoms:
    #                 atomIndexList.append(atomCount)

    #             atomCount += 1

    #         # If we use "charge-coupling" (1), assign the atomIndices of one BUF
    #         # to one protonatable lambda residue (use clever list slicing):            
    #         if (universe.get('ph_QQleveling') == 1):
    #             start = (groupNumber - 1) * len(BUF_qqA)
    #             stop  = start + len(BUF_qqA)
    #             atomIndexList += bufferAtomIndexList[start:stop]                

    #         writeTheGroup(groupNumber, atomIndexList)
    #         groupNumber += 1

    #     else:
    #         for atom in residue.d_atoms:
    #             atomCount += 1

    # If we use "charge-restraining" (2), add everything in bufferAtomIndexList
    # to the last lambda index group:

    # If we do charge restraining, write the buffer index group
    if restrainCharge:
        writeTheGroup(groupNumber, bufferAtomIndexList)

    file.close() # index.ndx
