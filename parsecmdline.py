# PYTHON_ARGCOMPLETE_OK

import argparse, argcomplete

def parsecmdline():

    parser = argparse.ArgumentParser(prog='phbuilder', description="pHbuilder builds constant-pH simulations for GROMACS.", epilog="Copyright blabla")

    subparsers = parser.add_subparsers(title='subcommands', description="description of subparsers?", required=False)

    parser_1 = subparsers.add_parser('gentopol', help="blabla about gentopol")
    required = parser_1.add_argument_group("required arguments")

    required.add_argument('-f', '--file', 
                        required=True,
                        dest='file',
                        action='store',
                        help='specify structure file for input (.pdb/.gro)')

    required.add_argument('-m', '--mode', 
                        required=True,
                        dest='mode',
                        action='store',
                        choices=['all', 'list', 'interactive', 'none'],
                        help='specify operationmode.')

    parser_1.add_argument('-o', '--output',
                        required=False,
                        dest='output',
                        action='store',
                        default='phprocessed.pdb',
                        help='specify structure file for output (.pdb/.gro)')

    parser_1.add_argument('-r', '--restraincharge', 
                        required=False,
                        dest='restraincharge',
                        action='store',
                        choices=['yes', 'no'],
                        default='yes',
                        help='restrains the charges using buffers')

    parser_1.add_argument('-v', '--verbosity',
                        required=False,
                        dest='verbosity',
                        action='store',
                        default=2,
                        choices=[0, 1, 2, 3],
                        help='set verbosity. 0 : supress all output, 1 only warnings and errors, 2 default, 3 more verbose',
                        type=int)

    parser_1.set_defaults(target='gentopol')

    parser_2 = subparsers.add_parser('addbuffers', help="blabla about addbuffers")
    parser_2.set_defaults(target='addbuffers')

    parser_3 = subparsers.add_parser('genparams', help="blabla about genparams")
    parser_3.set_defaults(target='genparams')

    argcomplete.autocomplete(parser)
    CLI = parser.parse_args()
