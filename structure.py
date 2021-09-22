import os

# Stores the information for a residue.
class Residue:
    def __init__(self, atoms, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms        # list      holds atom types
        self.d_resname = resname      # string    holds residue name
        self.d_chain   = chain        # string    holds chain name (A, B, etc.)
        self.d_resid   = resid        # int       holds residue number
        self.d_x       = x            # list      holds x-coordinates
        self.d_y       = y            # list      holds y-coordinates
        self.d_z       = z            # list      holds z-coordinates

# Stores the information for the periodic box.
class Crystal:
    def __init__(self, a, b, c, alpha, beta, gamma, space, Z):
        self.d_a       = a            # Angstroms
        self.d_b       = b            # Angstroms
        self.d_c       = c            # Angstroms
        self.d_alpha   = alpha        # degrees
        self.d_beta    = beta         # degrees
        self.d_gamma   = gamma        # degrees
        self.d_space   = space        # Space group
        self.d_Z       = Z            # Z va

class Structure:
    def __init__(self, name):
        extension = os.path.splitext(name)[1]
        
        if (extension == ".pdb"):
            self.__read_pdb(name)

        if (extension == ".gro"):
            self.__read_gro(name)

    # Writes d_residues to a structure (pdb/gro) file.
    def write(self, name):
        extension = os.path.splitext(name)[1]

        if (extension == ".pdb"):
            self.__write_pdb(name)

        if (extension == ".gro"):
            self.__write_gro(name)

    # Load a .pdb file into d_residues.
    def __read_pdb(self, name):

        with open(name) as file:

            atomLines = []

            for line in file.readlines():

                # Get title.
                if (line[0:6] == "TITLE "):
                    self.d_title = line[7:80].strip()

                # Get periodic box information (if any).
                elif (line[0:6] == "CRYST1"):
                    self.d_box = Crystal(float(line[6:15]), float(line[15:24]), float(line[24:33]), float(line[33:40]), float(line[40:47]), float(line[47:54]), line[55:66], int(line[66:70]))

                # Get the ATOM lines.
                elif (line[0:6] == "ATOM  "):
                    atomLines.append(line)

                # Only import the first MODEL...
                elif (line[0:6] == "ENDMDL"):
                    break

        # Loop through the atomLines and create a list of Residue objects.

        residues = []
        atoms    = []
        x        = []
        y        = []
        z        = []
        lastLine = False

        for idx in range(0, len(atomLines)):

            atoms.append(atomLines[idx][12:16].strip())
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
                residues.append(Residue(atoms, currentResName, currentChain, currentResID, x, y, z))

                # Reset.
                atoms = []
                x     = []
                y     = []
                z     = []

        # Add the list of Residues to universe.
        self.d_residues = residues

    def __write_pdb(self, name):
        with open(name, 'w') as file:
            if hasattr(self, 'd_title'):
                file.write("TITLE     {0}\n".format(self.d_title))

            if hasattr(self, 'd_box'):
                cryst = self.d_box
                file.write("CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} {:11s}{:>4d}\n".format(cryst.d_a, cryst.d_b, cryst.d_c, cryst.d_alpha, cryst.d_beta, cryst.d_gamma, cryst.d_space, cryst.d_Z))

            file.write("MODEL {:8d}\n".format(1))

            atomNumber = 1
            for residue in self.d_residues:
                for idx in range(0, len(residue.d_atoms)):
                    
                    atom = residue.d_atoms[idx]
                    if len(atom) == 3:
                        atom = ' ' + atom

                    file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber % 100000, atom, '', residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                    atomNumber += 1

            file.write("TER\nENDMDL\n")

    def __read_gro(self, name):

        atomLines = open(name).read().splitlines()

        # Loop through the atomLines and create a list of Residue objects.

        residues = []
        atoms    = []
        x        = []
        y        = []
        z        = []

        for idx in range(0, len(atomLines)):

            # Title.
            if (idx == 0):
                self.d_title = atomLines[idx]
                continue

            # Number of atoms.
            if (idx == 1):
                continue

            # Periodic box information.
            if (idx == len(atomLines) - 1):
                self.d_box = Crystal(10*float(atomLines[idx][0:10]), 10*float(atomLines[idx][10:20]), 10*float(atomLines[idx][20:30]), 90, 90, 90, "P 1", 1)
                continue

            atoms.append(atomLines[idx][11:15].strip())
            x.append(10 * float(atomLines[idx][20:28]))
            y.append(10 * float(atomLines[idx][28:36]))
            z.append(10 * float(atomLines[idx][36:44]))

            if (idx != len(atomLines) - 2):
                currentResID = int(atomLines[idx][0:5])
                nextResID    = int(atomLines[idx + 1][0:5])

            if (currentResID != nextResID or idx == len(atomLines) - 2):
                currentResName = atomLines[idx][5:10].strip()

                # Create the Residue object.
                residues.append(Residue(atoms, currentResName, ' ', currentResID, x, y, z))

                # Reset.
                atoms = []
                x     = []
                y     = []
                z     = []

        # Add the list of Residues to universe.
        self.d_residues = residues

    def __write_gro(self, name):
        with open(name, 'w') as file:
            # Title.
            if hasattr(self, 'd_title'):
                file.write("{}\n".format(self.d_title.strip()))

            # Total number of atoms.
            total = 0
            for residue in self.d_residues:
                for _ in residue.d_atoms:
                    total += 1
            file.write("{:>5d}\n".format(total))

            # Atoms.
            total = 1
            for residue in self.d_residues:
                for idx in range(0, len(residue.d_atoms)):
                    file.write("{:>5d}{:5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(
                        residue.d_resid, residue.d_resname, residue.d_atoms[idx], 
                        total % 100000, residue.d_x[idx]/10, residue.d_y[idx]/10, residue.d_z[idx]/10))
                    total += 1

            # Periodic box.
            if hasattr(self, 'd_box'):
                cryst = self.d_box
                file.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(cryst.d_a/10, cryst.d_b/10, cryst.d_c/10))
            else:
                file.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(0.0, 0.0, 0.0))
