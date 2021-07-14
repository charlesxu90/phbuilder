# Stores the information for a residue.
class Residue:
    def __init__(self, atoms, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list      holds atom types
        self.d_resname = resname    # string    holds residue name
        self.d_chain   = chain      # string    holds chain name (A, B, etc.)
        self.d_resid   = resid      # int       holds residue number
        self.d_x       = x          # list      holds x-coordinates
        self.d_y       = y          # list      holds y-coordinates
        self.d_z       = z          # list      holds z-coordinates


# Stores the information for the periodic box.
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

# Stores the information for a lambda group type.
class LambdaType:
    def __init__(self, groupname, pKa, atoms, qqA, qqB, dvdl):
        self.d_groupname = groupname
        self.d_pKa       = pKa
        self.d_atoms     = atoms
        self.d_qqA       = qqA
        self.d_qqB       = qqB
        self.d_dvdl      = dvdl
