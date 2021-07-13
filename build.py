#!/bin/python3

import shelve
import os

# UNIVERSE #####################################################################
################################################################################

# Set/update variable to 
def add(varName, value):
    with shelve.open('universe') as shelf:
        shelf[varName] = value

# Check whether universe contains a certain varName
def has(varName):
    with shelve.open('universe') as shelf:
        return varName in shelf

# Retrieve variable from 
def get(varName):
    if has(varName):
        return shelve.open('universe')[varName]

    data = eval(input("couldn't retrieve var \"{0}\" from  Enter manually: ".format(varName)))
    print("add {0} = {1} {2}".format(varName, data, type(data)))
    add(varName, data)
    return data

# Display all variables (name, data, type) stored in the 
def inspect():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in sorted(shelf):
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}] ({3}) {4}".format(item, shelf[item][0], shelf[item][-1], len(shelf[item]), type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1} {2}".format(item, shelf[item], type(shelf[item]), arg=longest))

# UTILS ########################################################################
################################################################################

# Prints an update for the user.
def update(tool, message):
    print("{:18s} : {:s}".format(tool, message))

# Prints a warning for the user.
def warning(tool, message):
    print("{:18s} : WARNING - {:s}".format(tool, message))

def error(tool, message):
    print("{:18s} : ERROR - {:s} quiting...".format(tool, message))
    quit()

# STRUCTURE ####################################################################
################################################################################

class Residue:
    def __init__(self, atoms, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list      holds atom types
        self.d_resname = resname    # string    holds residue name
        self.d_chain   = chain      # string    holds chain name (A, B, etc.)
        self.d_resid   = resid      # int       holds residue number
        self.d_x       = x          # list      holds x-coordinates
        self.d_y       = y          # list      holds y-coordinates
        self.d_z       = z          # list      holds z-coordinates

class Crystal:
    def __init__(self, a, b, c, alpha, beta, gamma, space, Z):
        self.d_a       = a          # Angstroms
        self.d_b       = b          # Angstroms
        self.d_c       = c          # Angstroms
        self.d_alpha   = alpha      # degrees
        self.d_beta    = beta       # degrees
        self.d_gamma   = gamma      # degrees
        self.d_space   = space      # Space group
        self.d_Z       = Z          # Z value

# Load a .pdb file into d_residues.
def read_pdb(name, d_model=1, d_chain=[]):
    
    add('d_model', d_model)
    add('d_chain', d_chain)

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
                d_title = line[7:80].rstrip(); add('d_title', d_title)

            # Get periodic box information (if any).
            elif (line[0:6] == "CRYST1"):
                # print("'{}'".format(line[ 6:15])); print("'{}'".format(line[15:24]))
                # print("'{}'".format(line[24:33])); print("'{}'".format(line[33:40]))
                # print("'{}'".format(line[40:47])); print("'{}'".format(line[47:54]))
                # print("'{}'".format(line[55:66])); print("'{}'".format(line[66:70]))
                add('d_box', Crystal(float(line[6:15]), float(line[15:24]), float(line[24:33]), float(line[33:40]), float(line[40:47]), float(line[47:54]), line[55:66], int(line[66:70])))

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

    # Add the list of Residues to universe.
    add('d_residues', d_residues)

def write_pdb(name):
    with open(name, 'w') as file:
        if has('d_title'):
            file.write("TITLE {0}\n".format(get('d_title')))

        if has('d_box'):
            cryst = get('d_box')
            file.write("CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} {:11s}{:>4d}\n".format(cryst.d_a, cryst.d_b, cryst.d_c, cryst.d_alpha, cryst.d_beta, cryst.d_gamma, cryst.d_space, cryst.d_Z))

        file.write("MODEL {:8d}\n".format(get('d_model')))

        atomNumber = 1
        for residue in get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber % 100000, residue.d_atoms[idx], '', residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

def read_gro(name):

    add('d_model', 1)

    atomLines = open(name).read().splitlines()

    # Loop through the atomLines and create a list of Residue objects.

    d_residues = []
    atoms      = []
    x          = []
    y          = []
    z          = []
    lastLine   = False

    for idx in range(0, len(atomLines)):

        # Title.
        if (idx == 0):
            add('d_title', atomLines[idx])
            continue

        # Number of atoms.
        if (idx == 1):
            continue

        # Periodic box information.
        if (idx == len(atomLines) - 1):
            add('d_box', Crystal(10*float(atomLines[idx][0:10]), 10*float(atomLines[idx][10:20]), 10*float(atomLines[idx][20:30]), 90, 90, 90, "P 1", 1))
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

    # Add the list of Residues to universe.
    add('d_residues', d_residues)

def write_gro(name):
    with open(name, 'w') as file:
        if has('d_title'):
            file.write("{}\n".format(get('d_title').strip()))

        total = 0
        for residue in get('d_residues'):
            for _ in residue.d_atoms:
                total += 1
        file.write("{:>5d}\n".format(total))

        total = 1
        for residue in get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:>5d}{:5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                    residue.d_resid, residue.d_resname, residue.d_atoms[idx].strip(), 
                    total % 100000, residue.d_x[idx]/10, residue.d_y[idx]/10, residue.d_z[idx]/10))
                total += 1

        if has('d_box'):
            cryst = get('d_box')
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
        error("loadstructure", "Unknown file format specified. Formats are .pdb or .gro.")

def writestructure(name):
    extension = os.path.splitext(name)[1]

    if (extension == ".pdb"):
        write_pdb(name)
    elif (extension == ".gro"):
        write_gro(name)
    else:
        error("writestructure", "Unknown file format specified. Formats are .pdb or .gro.")

# MAIN #########################################################################
################################################################################
