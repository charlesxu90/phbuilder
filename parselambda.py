import utils, universe, configparser

# Stores the information for a lambda group type.
class LambdaType:
    def __init__(self, groupname, incl, pKa, atoms, qqA, qqB, dvdl):
        self.d_groupname = groupname
        self.d_incl      = incl
        self.d_pKa       = pKa
        self.d_atoms     = atoms
        self.d_qqA       = qqA
        self.d_qqB       = qqB
        self.d_dvdl      = dvdl

# Parses lambdagrouptypes.dat.
def parseLambdaGroupTypesFile():

    # Add a lambda residue-type to universe.
    def defineLambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl):

        # Create a temporary LambdaType object.
        NewLambdaType = LambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl)
        
        # If a list of LambdaType objects (ph_lambdaTypes) is already present in the universe,
        # then check whether a LambdaType with the same groupname is already stored.
        if universe.has('ph_lambdaTypes'):
            temp = universe.get('ph_lambdaTypes')

            for entry in temp:
                if entry.d_groupname == NewLambdaType.d_groupname:
                    utils.update("LambdaType with groupname {} is already defined in universe. Skipping...")
                    break
            else:
                # Add th new LambdaType to the list and store it again in universe.
                temp.append(NewLambdaType)
                universe.add('ph_lambdaTypes', temp)
        else: # If no list of LambdaTypes exists, create one.
            universe.add('ph_lambdaTypes', [NewLambdaType])

    # Interal function to convert string to list of floats.
    def str2floatList(string):
        return [float(val) for val in string.split(' ')]

    # Internal function to convert string to list of strings.
    def str2strList(string):
        return string.split(' ')

    parser = configparser.ConfigParser()
    parser.read("lambdagrouptypes.dat") # name is hardcoded.

    # Loop through the sections.
    for sect in parser.sections():

        if (sect.strip() == "FORCEFIELD"):
            universe.add('d_modelFF', parser.get(sect, 'path'))
            universe.add('d_modelwater', parser.get(sect, 'water'))
            continue

        if (sect.strip() == "BUF"):
            universe.add('ph_BUF_dvdl', str2floatList(parser.get(sect, 'dvdl')))
            continue

        groupname = sect.strip()
        pKa       = parser.getfloat(sect, 'pKa')
        incl      = str2strList(parser.get(sect, 'incl'))
        atoms     = str2strList(parser.get(sect, 'atoms'))
        qqA       = str2floatList(parser.get(sect, 'qqA'))
        qqB       = str2floatList(parser.get(sect, 'qqB'))
        dvdl      = str2floatList(parser.get(sect, 'dvdl'))

        # Sanitize input of groupname.
        if (len(groupname) < 2 or len(groupname) > 4):
            utils.error("groupname of LambdaType needs to contain between 2 and 4 characters.")

        # Call function that constructs the LambdaType object and adds it to universe.
        defineLambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl)

    # User update.
    utils.update("ffpath    = {}".format(universe.get('d_modelFF')))
    utils.update("water     = {}".format(universe.get('d_modelwater')))

    for obj in universe.get('ph_lambdaTypes'):
        utils.update("groupname = {}".format(obj.d_groupname))
        utils.update("incl      = {}".format(obj.d_incl))
        utils.update("pKa       = {}".format(obj.d_pKa))
        utils.update("atoms     = {}".format(obj.d_atoms))
        utils.update("qqA       = {}".format(obj.d_qqA))
        utils.update("qqB       = {}".format(obj.d_qqB))
        utils.update("dvdl      = {}\n".format(obj.d_dvdl))

    if (universe.has('ph_BUF_dvdl')):
        utils.update("BUF_dvdl  = {}\n".format(universe.get('ph_BUF_dvdl')))
    else:
        utils.update("dvdl coefficients for buffer were not found in lambdagrouptypes.dat.")
        utils.update("This is fine if you don't plan on using charge restraining/buffers.")
