import utils, universe, configparser

# Stores the information for a lambda group type.
class LambdaType:
    def __init__(self, groupname, incl, pKa, atoms, qqA, qqB, dvdl):
        self.d_groupname = groupname    # str
        self.d_incl      = incl         # list
        self.d_pKa       = pKa          # list  (previously str)
        self.d_atoms     = atoms        # list
        self.d_qqA       = qqA          # list
        self.d_qqB       = qqB          # list of lists (previously list)
        self.d_dvdl      = dvdl         # list of lists (previously list)

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
            universe.add('ph_BUF_range', str2floatList(parser.get(sect, 'range')))
            continue
        
        # Parse groupname
        groupname = sect.strip()[0:4]
        
        # Parse incl
        incl = str2strList(parser.get(sect, 'incl'))

        # Parse atoms
        atoms = str2strList(parser.get(sect, 'atoms'))

        # Parse qqA
        qqA = str2floatList(parser.get(sect, 'qqA'))

        pKa  = []
        qqB  = []
        dvdl = []
        for idx in range(1, 11): # Max 10 multistates
            try:
                # Parse pKa(s)
                pKa.append(parser.get(sect, 'pKA_{}'.format(idx)))

                # Parse qqB(s)
                qqB.append(str2floatList(parser.get(sect, 'qqB_{}'.format(idx))))

                # Parse dvdl(s)
                dvdl.append(str2floatList(parser.get(sect, 'dvdl_{}'.format(idx))))
            except:
                break

        # SANITIZE INPUT

        if (len(groupname) < 2 or len(groupname) > 4):
            utils.error("groupname of LambdaType needs to contain between 2 and 4 characters.")

        # Call function that constructs the LambdaType object and adds it to universe.
        defineLambdaType(groupname, incl, pKa, atoms, qqA, qqB, dvdl)

    # USER UPDATE

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
