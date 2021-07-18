import utils, universe, configparser

# Stores the information for a lambda group type.
class LambdaType:
    def __init__(self, groupname, pKa, atoms, qqA, qqB, dvdl):
        self.d_groupname = groupname
        self.d_pKa       = pKa
        self.d_atoms     = atoms
        self.d_qqA       = qqA
        self.d_qqB       = qqB
        self.d_dvdl      = dvdl

def parseLambdaGroupTypesFile():

    # Add a lambda residue-type to universe.
    def defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl):
        NewLambdaType = LambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)
        if universe.has('ph_lambdaTypes'):
            temp = universe.get('ph_lambdaTypes')
            
            for entry in temp:
                if entry.d_groupname == NewLambdaType.d_groupname:
                    utils.update("LambdaType with groupname {} is already defined in universe. Skipping...")
                    break
            else:
                temp.append(NewLambdaType)
                universe.add('ph_lambdaTypes', temp)
        else:
            universe.add('ph_lambdaTypes', [NewLambdaType])

    # Convert string to list of floats.
    def str2floatList(string):
        return [float(val) for val in string.split(' ')]

    # Convert string to list of strings.
    def str2strList(string):
        return string.split(' ')

    parser = configparser.ConfigParser()
    parser.read("lambdagrouptypes.dat")
    
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
        atoms     = str2strList(parser.get(sect, 'atoms'))
        qqA       = str2floatList(parser.get(sect, 'qqA'))
        qqB       = str2floatList(parser.get(sect, 'qqB'))
        dvdl      = str2floatList(parser.get(sect, 'dvdl'))

        if (len(groupname) < 2 or len(groupname) > 4):
            utils.error("groupname of LambdaType needs to contain between 2 and 4 characters.")

        defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)

    # User update.
    utils.update("ffpath = {}".format(universe.get('d_modelFF')), 3)
    utils.update("water  = {}".format(universe.get('d_modelwater')), 3)

    for obj in universe.get('ph_lambdaTypes'):
        utils.update("groupname = {}".format(obj.d_groupname), 3)
        utils.update("pKa       = {}".format(obj.d_pKa), 3)
        utils.update("atoms     = {}".format(obj.d_atoms), 3)
        utils.update("qqA       = {}".format(obj.d_qqA), 3)
        utils.update("qqB       = {}".format(obj.d_qqB), 3)
        utils.update("dvdl      = {}\n".format(obj.d_dvdl), 3)

    if (universe.has('ph_BUF_dvdl')):
        utils.update("BUF_dvdl  = {}\n".format(universe.get('ph_BUF_dvdl')), 3)
