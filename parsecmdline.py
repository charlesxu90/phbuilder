# PYTHON_ARGCOMPLETE_OK

import argparse, argcomplete, universe, utils

def parsecmdline():

    parser = argparse.ArgumentParser(prog='phbuilder', description="pHbuilder builds constant-pH simulations for GROMACS.", epilog="Copyright blabla")

    subparsers = parser.add_subparsers(title='subcommands', description="description of subparsers?", required=False)

    parser_1 = subparsers.add_parser('gentopol', help="blabla about gentopol")
    required = parser_1.add_argument_group("required arguments")

    required.add_argument('-f', 
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

    # parser_2 = subparsers.add_parser('addbuffers', help="blabla about addbuffers")
    # parser_2.set_defaults(target='addbuffers')

    # parser_3 = subparsers.add_parser('genparams', help="blabla about genparams")
    # parser_3.set_defaults(target='genparams')

    argcomplete.autocomplete(parser)    # Required for autocompleting using argcomplete.
    CLI = parser.parse_args()           # Do the actual parsing.

    # Add relevant parameters to the universe.
    universe.add('d_target', CLI.target)
    universe.add('d_file', CLI.file)
    universe.add('d_output', CLI.output)
    universe.add('d_verbosity', CLI.verbosity)

    if (CLI.inter != None):
        universe.add('ph_inter', True)

    if (CLI.list != None):
        resid = []

        for line in open(CLI.list).readlines():
            resid.append(line.split()[0])

        resid = [int(i) for i in resid]

        universe.add('ph_list_resid', resid)

    # User information.
    utils.update("Parsed the following input from the command line:")
    utils.update(vars(CLI))
