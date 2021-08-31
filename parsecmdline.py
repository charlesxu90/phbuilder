# PYTHON_ARGCOMPLETE_OK

import argparse, argcomplete, universe, utils

def parsecmdline():

    parser = argparse.ArgumentParser(prog='phbuilder', description="pHbuilder builds constant-pH simulations for GROMACS.", epilog="Copyright blabla")

    subparsers = parser.add_subparsers(title='subcommands', description="description of subparsers?", required=False)

    parser_1  = subparsers.add_parser('gentopol', help="blabla about gentopol")
    required1 = parser_1.add_argument_group("required arguments")

    required1.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='specify structure file for input (.pdb/.gro)')

    parser_1.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phprocessed.pdb',
                        help='specify structure file for output (.pdb/.gro)')

    parser_1.add_argument('-inter',
                        required=False,
                        dest='inter',
                        action='store_const',
                        const=1,
                        help='Interactively select which residues to make protonatable')

    parser_1.add_argument('-list',
                        required=False,
                        dest='list',
                        action='store',
                        help='specify list of resid(ue)s to be considered')

    parser_1.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose',
                        type=int)

    parser_1.set_defaults(target='gentopol')

    parser_2  = subparsers.add_parser('addbuffers', help="blabla about addbuffers")
    required2 = parser_2.add_argument_group("required arguments")

    required2.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='specify structure file for input (.pdb/.gro).')

    required2.add_argument('-p', 
                        required=False,
                        dest='topol',
                        action='store',
                        default='topol.top',
                        help='specify topology file for input (.top).')

    parser_2.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phbuffers.pdb',
                        help='specify structure file for output (.pdb/.gro).')

    parser_2.add_argument('-solname',
                        required=False,
                        dest='solname',
                        action='store',
                        default='SOL',
                        help='specify name of solvent molecule.')

    parser_2.add_argument('-nbufs',
                        required=False,
                        dest='nbufs',
                        action='store',
                        help='specify number of buffer molecules.',
                        type=int)

    parser_2.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose.',
                        type=int)

    parser_2.set_defaults(target='addbuffers')

    parser_3  = subparsers.add_parser('genparams', help="blabla about genparams")
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

    # Add universal parameters to the universe (these are used by all three targets).
    universe.add('d_target', CLI.target)
    universe.add('d_verbosity', CLI.verbosity)
    
    # If we run gentopol...
    if (CLI.target == 'gentopol'):
        universe.add('d_file', CLI.file)
        universe.add('d_output', CLI.output)

        # Process whether the -inter flag was or wasn't set.
        if (CLI.inter != None):
            universe.add('ph_inter', True)

        # Process whether the -list flag was or wasn't set.
        if (CLI.list != None):
            resid = []

            for line in open(CLI.list).readlines():
                resid.append(line.split()[0])

            resid = [int(i) for i in resid]

            universe.add('ph_list_resid', resid)

    # If we run addbuffers...
    elif (CLI.target == 'addbuffers'):
        # Either required or has a default value
        universe.add('d_file', CLI.file)
        universe.add('d_topol', CLI.topol)
        universe.add('d_output', CLI.output)
        universe.add('ph_solname', CLI.solname)

        # Optional
        if (CLI.nbufs != None):
            universe.add('ph_nbufs', CLI.nbufs)

    # If we run genparams...
    elif (CLI.target == 'genparams'):
        universe.add('d_file', CLI.file)
        universe.add('d_mdp', CLI.mdp)
        universe.add('d_ndx', CLI.ndx)
        universe.add('ph_ph', CLI.ph)
        universe.add('ph_nstout', CLI.nstout)
        universe.add('ph_dwpE', CLI.dwpE)
        universe.add('ph_lmass', CLI.lmass)
        universe.add('ph_ltau', CLI.ltau)

        # Process whether the -inter flag was or wasn't set
        if (CLI.inter != None):
            universe.add('ph_inter', True)

    # User information.
    utils.update("Parsed the following input from the command line:")
    utils.update(vars(CLI))
