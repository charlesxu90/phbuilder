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

    utils.pedantic('parseLambdaGroupTypes', 'parsing lambdagrouptypes.dat...\n')

    # Add a lambda residue-type to universe.
    def defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl):
        NewLambdaType = LambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)
        if universe.has('ph_lambdaTypes'):
            temp = universe.get('ph_lambdaTypes')
            
            for entry in temp:
                if entry.d_groupname == NewLambdaType.d_groupname:
                    utils.warning("defineLambdaType", "LambdaType with groupname {} is already defined in ph_lambdaTypes. Skipping...".format(NewLambdaType.d_groupname))
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

        if (len(groupname) != 4):
            utils.error("parseLambdaGroupTypes", "groupname of LambdaType needs to have 4 letters")

        defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)

    if (universe.get('d_verbosity') == 3):
        for obj in universe.get('ph_lambdaTypes'):
            print("groupname = {}".format(obj.d_groupname))
            print("pKa       = {}".format(obj.d_pKa))
            print("atoms     = {}".format(obj.d_atoms))
            print("qqA       = {}".format(obj.d_qqA))
            print("qqB       = {}".format(obj.d_qqB))
            print("dvdl      = {}\n".format(obj.d_dvdl))

        if (universe.has('ph_BUF_dvdl')):
            print("BUF_dvdl  = {}\n".format(universe.get('ph_BUF_dvdl')))
