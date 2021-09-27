# PYTHON_ARGCOMPLETE_OK

import argparse, argcomplete

# Parses command line using the python argparse module.
def parsecmdline():
    desc_1 = "phbuilder builds constant-pH simulations for GROMACS."
    desc_2 = "Copyright Anton Jansen 2021."
    desc_3 = "Encapsulates gmx pdb2gmx. Regenerates topology of the specified \
              protein using our modified charmm36 force field. By default, every \
              residue part of the incl groups specified in lambdagrouptyes.dat \
              will be made titratable. This behavior can be modified by specifying \
              the -inter and -list flags."
    desc_4 = "Encapsulated gxm genion. Adds buffer particles."
    desc_5 = "Generates the constant-pH input parameters for the .mdp file."

    parser = argparse.ArgumentParser(prog='phbuilder', description=desc_1, epilog=desc_2)

    subparsers = parser.add_subparsers(required=False)

    parser_1  = subparsers.add_parser('gentopol', help=desc_3)
    required1 = parser_1.add_argument_group("required arguments")

    required1.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='Specify structure file for input (.pdb/.gro).')

    parser_1.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phprocessed.pdb',
                        help='Specify structure file for output (.pdb/.gro).')

    parser_1.add_argument('-inter',
                        required=False,
                        dest='inter',
                        action='store_const',
                        const=1,
                        help='Interactively select which residues to make protonatable.')

    parser_1.add_argument('-list',
                        required=False,
                        dest='list',
                        action='store',
                        help='Specify list of resid(ue)s to be considered.')

    parser_1.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='Set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose.',
                        type=int)

    parser_1.set_defaults(target='gentopol')

    parser_2  = subparsers.add_parser('addbuffers', help=desc_4)
    required2 = parser_2.add_argument_group("required arguments")

    required2.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='Specify structure file for input (.pdb/.gro).')

    required2.add_argument('-p', 
                        required=False,
                        dest='topol',
                        action='store',
                        default='topol.top',
                        help='Specify topology file for input (.top).')

    parser_2.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phbuffers.pdb',
                        help='Specify structure file for output (.pdb/.gro).')

    parser_2.add_argument('-solname',
                        required=False,
                        dest='solname',
                        action='store',
                        default='SOL',
                        help='Specify name of solvent molecule.')

    parser_2.add_argument('-nbufs',
                        required=False,
                        dest='nbufs',
                        action='store',
                        help='Specify number of buffer molecules.',
                        type=int)

    parser_2.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='Set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose.',
                        type=int)

    parser_2.set_defaults(target='addbuffers')

    parser_3  = subparsers.add_parser('genparams', help=desc_5)
    required3 = parser_3.add_argument_group("required arguments")

    required3.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='specify structure file for input (.pdb/.gro).')

    required3.add_argument('-ph', 
                        required=True,
                        dest='ph',
                        action='store',
                        help='specify pH of system',
                        type=float)

    parser_3.add_argument('-mdp',
                        required=False,
                        dest='mdp',
                        action='store',
                        help='Specify .mdp file to write constant-pH parameters to.')

    parser_3.add_argument('-ndx',
                        required=False,
                        dest='ndx',
                        action='store',
                        help='Specify .ndx file to write constant-pH index groups to.')

    parser_3.add_argument('-nstout',
                        required=False,
                        dest='nstout',
                        action='store',
                        default=500,
                        help='Specify output frequency for lambda data.',
                        type=int)

    parser_3.add_argument('-dwpE',
                        required=False,
                        dest='dwpE',
                        action='store',
                        default=7.5,
                        help='Specify bias potential energy.',
                        type=float)

    parser_3.add_argument('-lmass',
                        required=False,
                        dest='lmass',
                        action='store',
                        default=5,
                        help='Specify mass of lambda particle (amu).',
                        type=float)

    parser_3.add_argument('-ltau',
                        required=False,
                        dest='ltau',
                        action='store',
                        default=2,
                        help='Specify tau for lambda dynamics thermostat (ps^-1).',
                        type=float)

    parser_3.add_argument('-inter',
                        required=False,
                        dest='inter',
                        action='store_const',
                        const=1,
                        help='Interactively set barrier energy for different lambda groups.')

    parser_3.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose',
                        type=int)

    parser_3.set_defaults(target='genparams')

    argcomplete.autocomplete(parser)    # Required for autocompleting using argcomplete.
    CLI = parser.parse_args()           # Do the actual parsing.

    return CLI # object containing all the parsed key-value pairs.
