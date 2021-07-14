import utility
from classes import LambdaType

import configparser

def parseLambdaGroupTypes():

    # Add a lambda residue-type to universe.
    def defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl):
        NewLambdaType = LambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)
        if utility.has('ph_lambdaTypes'):
            temp = utility.get('ph_lambdaTypes')
            
            for entry in temp:
                if entry.d_groupname == NewLambdaType.d_groupname:
                    utility.warning("defineLambdaType", "LambdaType {} is already defined in ph_lambdaTypes. Skipping...".format(NewLambdaType.d_groupname))
                    break
            else:
                temp.append(NewLambdaType)
                utility.add('ph_lambdaTypes', temp)
        else:
            utility.add('ph_lambdaTypes', [NewLambdaType])

    # Convert string to list of floats.
    def str2floatList(string):
        return [float(val) for val in string.split(' ')]

    # Convert string to list of strings.
    def str2strList(string):
        return string.split(' ')

    parser = configparser.ConfigParser()
    parser.read("lambdagrouptypes.dat")
    
    for sect in parser.sections():
        
        if (sect.strip() == "BUF"):
            utility.add('pH_BUF_dvdl', str2floatList(parser.get(sect, 'dvdl')))
            continue

        groupname = sect.strip()
        pKa       = parser.getfloat(sect, 'pKa')
        atoms     = str2strList(parser.get(sect, 'atoms'))
        qqA       = str2floatList(parser.get(sect, 'qqA'))
        qqB       = str2floatList(parser.get(sect, 'qqB'))
        dvdl      = str2floatList(parser.get(sect, 'dvdl'))

        if (len(groupname) != 4):
            utility.error("parseLambdaGroupTypes", "Name of lambdagrouptype should be 4 letters")

        # print(groupname)
        # print(pKa)
        # print(atoms)
        # print(qqA)
        # print(qqB)
        # print(dvdl)

        defineLambdaType(groupname, pKa, atoms, qqA, qqB, dvdl)
