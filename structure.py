import utility, os
from classes import Residue, Crystal

# Load a .pdb file into d_residues.
def read_pdb(name, d_model=1, d_chain=[]):
    
    utility.add('d_model', d_model)
    utility.add('d_chain', d_chain)

    with open(name) as file:

        correctModel = True
        atomLines    = []

        for line in file.readlines():

            # Only import the specified MODEL number.
            if (line[0:6] == "MODEL "):
                if ("MODEL {:8d}".format(d_model) in line):
                    correctModel = True
                else:
                    correctModel = False

            # Get title.
            elif (line[0:6] == "TITLE "):
                d_title = line[7:80].rstrip(); utility.add('d_title', d_title)

            # Get periodic box information (if any).
            elif (line[0:6] == "CRYST1"):
                # print("'{}'".format(line[ 6:15])); print("'{}'".format(line[15:24]))
                # print("'{}'".format(line[24:33])); print("'{}'".format(line[33:40]))
                # print("'{}'".format(line[40:47])); print("'{}'".format(line[47:54]))
                # print("'{}'".format(line[55:66])); print("'{}'".format(line[66:70]))
                utility.add('d_box', Crystal(float(line[6:15]), float(line[15:24]), float(line[24:33]), float(line[33:40]), float(line[40:47]), float(line[47:54]), line[55:66], int(line[66:70])))

            # If our line is an ATOM,
            elif (line[0:6] == "ATOM  "):
                # and we are currently reading the correct MODEL,
                if (correctModel):
                    # and we want all the chains,
                    if (d_chain == []):
                        # then load the line:
                        atomLines.append(line)
                    # Or, if we want a selection of chains,
                    elif (line[21:22] in d_chain):
                        # load the selection:
                        atomLines.append(line)

    # Loop through the atomLines and create a list of Residue objects.

    d_residues = []
    atoms      = []
    x          = []
    y          = []
    z          = []
    lastLine   = False

    for idx in range(0, len(atomLines)):

        atoms.append(atomLines[idx][12:16])
        x.append(float(atomLines[idx][30:38]))
        y.append(float(atomLines[idx][38:46]))
        z.append(float(atomLines[idx][46:54]))

        try:
            currentResID = int(atomLines[idx][22:26])
            nextResID    = int(atomLines[idx + 1][22:26])
        except IndexError:
            lastLine = True

        if (currentResID != nextResID or lastLine):
            
            currentResName = atomLines[idx][17:21].strip()
            currentChain   = atomLines[idx][21:22]
            
            # Create the Residue object.
            d_residues.append(Residue(atoms, currentResName, currentChain, currentResID, x, y, z))

            # Reset.
            atoms = []
            x     = []
            y     = []
            z     = []

    # Add the list of Residues to utility.
    utility.add('d_residues', d_residues)

def write_pdb(name):
    with open(name, 'w') as file:
        if utility.has('d_title'):
            file.write("TITLE {0}\n".format(utility.get('d_title')))

        if utility.has('d_box'):
            cryst = utility.get('d_box')
            file.write("CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} {:11s}{:>4d}\n".format(cryst.d_a, cryst.d_b, cryst.d_c, cryst.d_alpha, cryst.d_beta, cryst.d_gamma, cryst.d_space, cryst.d_Z))

        file.write("MODEL {:8d}\n".format(utility.get('d_model')))

        atomNumber = 1
        for residue in utility.get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber % 100000, residue.d_atoms[idx], '', residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

def read_gro(name):

    utility.add('d_model', 1)

    atomLines = open(name).read().splitlines()

    # Loop through the atomLines and create a list of Residue objects.

    d_residues = []
    atoms      = []
    x          = []
    y          = []
    z          = []

    for idx in range(0, len(atomLines)):

        # Title.
        if (idx == 0):
            utility.add('d_title', atomLines[idx])
            continue

        # Number of atoms.
        if (idx == 1):
            continue

        # Periodic box information.
        if (idx == len(atomLines) - 1):
            utility.add('d_box', Crystal(10*float(atomLines[idx][0:10]), 10*float(atomLines[idx][10:20]), 10*float(atomLines[idx][20:30]), 90, 90, 90, "P 1", 1))
            continue

        atom = atomLines[idx][11:15].strip()
        if (len(atom) == 3): 
            atom = ' ' + atom
        atoms.append(atom)

        x.append(10 * float(atomLines[idx][20:28]))
        y.append(10 * float(atomLines[idx][28:36]))
        z.append(10 * float(atomLines[idx][36:44]))

        if (idx != len(atomLines) - 2):
            currentResID = int(atomLines[idx][0:5])
            nextResID    = int(atomLines[idx + 1][0:5])

        if (currentResID != nextResID or idx == len(atomLines) - 2):
            currentResName = atomLines[idx][5:10].strip()

            # Create the Residue object.
            d_residues.append(Residue(atoms, currentResName, ' ', currentResID, x, y, z))

            # Reset.
            atoms = []
            x     = []
            y     = []
            z     = []

    # Add the list of Residues to utility.
    utility.add('d_residues', d_residues)

def write_gro(name):
    with open(name, 'w') as file:
        # Title.
        if utility.has('d_title'):
            file.write("{}\n".format(utility.get('d_title').strip()))

        # Total number of atoms.
        total = 0
        for residue in utility.get('d_residues'):
            for _ in residue.d_atoms:
                total += 1
        file.write("{:>5d}\n".format(total))

        # Atoms.
        total = 1
        for residue in utility.get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:>5d}{:5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                    residue.d_resid, residue.d_resname, residue.d_atoms[idx].strip(), 
                    total % 100000, residue.d_x[idx]/10, residue.d_y[idx]/10, residue.d_z[idx]/10))
                total += 1

        # Periodic box.
        if utility.has('d_box'):
            cryst = utility.get('d_box')
            file.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(cryst.d_a/10, cryst.d_b/10, cryst.d_c/10))
        else:
            file.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(0.0, 0.0, 0.0))

def loadstructure(name, model=1, chains=[]):
    extension = os.path.splitext(name)[1]
    
    if (extension == ".pdb"):
        read_pdb(name, d_model=model, d_chain=chains)
    elif (extension == ".gro"):
        read_gro(name)
    else:
        utility.error("loadstructure", "Unknown file format specified. Formats are .pdb or .gro.")

def writestructure(name):
    extension = os.path.splitext(name)[1]

    if (extension == ".pdb"):
        write_pdb(name)
    elif (extension == ".gro"):
        write_gro(name)
    else:
        utility.error("writestructure", "Unknown file format specified. Formats are .pdb or .gro.")
