# PYTHON_ARGCOMPLETE_OK

import argparse, argcomplete

# Parses command line using the python argparse module.
def parsecmdline():
    desc_1 = "System builder for constant-pH simulations in GROMACS. phbuilder consists of three tools: gentopol, neutralize, and genparams. Each tool performs a specific task for preparing a constant-pH simulation."
    desc_2 = "Please cite [phbuilder_paper]."
    desc_3 = "gentopol encapsulates gmx pdb2gmx and allows you to (re)generate the topology for your system using our modified version of the CHARMM36 force field. This is necessary as some dihedral parameters were modified for titratable residues (ref manuscript 2). gentopol by default allows you to interactively set the initial lambda value (protonation state) for each residue associated with a defined lambdagrouptype. This behavior can be automated by setting the -auto <ph> flag. In this case, every residue associated with a defined lambdagrouptype will automatically be made titratable, and the initial lambda values will be guessed based on the specified ph, together with the pKa defined in the lambdagrouptypes.dat file. Note that you should use the same pH value for genparams."
    desc_4 = "The purpose of this tool is to ensure a charge-neutral system by adding the appropriate number of ions and buffer particles."
    desc_5 = "genparams generates the .mdp files, including all the required constant-pH parameters. genparams requires the existance of a phrecord.dat file for setting the initial lambda values."

    parser = argparse.ArgumentParser(prog='phbuilder', description=desc_1, epilog=desc_2)

    subparsers = parser.add_subparsers(required=False)

    parser_1  = subparsers.add_parser('gentopol', help=desc_3)
    required1 = parser_1.add_argument_group("required arguments")

    required1.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='[<.pdb/.gro>] (required) Specify input structure file.')

    parser_1.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phprocessed.pdb',
                        help='[<.pdb/.gro>] (phprocessed.pdb) Specify output structure file.')

    parser_1.add_argument('-list',
                        required=False,
                        dest='list',
                        action='store',
                        help='[<.txt>] Provide a subset of resid(ue)s to consider. Helpful if you do not want to manually go through many (unimportant) residues.')

    parser_1.add_argument('-auto',
                        required=False,
                        dest='auto',
                        action='store',
                        help='[<real>] Use automatic mode and specify the simulation pH to base guess for initial lambda values on.',
                        type=float)

    parser_1.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store_const',
                        const=1,
                        help='(no) Be more verbose (helpful for debugging).')

    parser_1.set_defaults(target='gentopol')

    parser_2  = subparsers.add_parser('neutralize', help=desc_4)
    required2 = parser_2.add_argument_group("required arguments")

    required2.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='[<.pdb/.gro>] (required) Specify input structure file.')

    required2.add_argument('-p', 
                        required=False,
                        dest='topol',
                        action='store',
                        default='topol.top',
                        help='[<.top>] (topol.top) Specify input topology file.')

    parser_2.add_argument('-o',
                        required=False,
                        dest='output',
                        action='store',
                        default='phneutral.pdb',
                        help='[<.pdb/.gro>] (phneutral.pdb) Specify output structure file.')

    parser_2.add_argument('-solname',
                        required=False,
                        dest='solname',
                        action='store',
                        default='SOL',
                        help=' [<string>] (SOL) Specify solvent name (of which to replace molecules with ions and buffers).')

    parser_2.add_argument('-pname',
                        required=False,
                        dest='pname',
                        action='store',
                        default='NA',
                        help='[<string>] (NA) Specify name of positive ion to use. Analogous to gmx genion.')

    parser_2.add_argument('-nname',
                        required=False,
                        dest='nname',
                        action='store',
                        default='CL',
                        help='[<string>] (CL) Specify name of negative ion to use. Analogous to gmx genion.')

    parser_2.add_argument('-conc',
                        required=False,
                        dest='conc',
                        action='store',
                        default=0.0,
                        help='[<real>] (0.0) Specify ion concentration in mol/L. Analogous to gmx genion but will use the solvent volume for calculating the required number of ions, not the periodic box volume as genion does.',
                        type=float)

    parser_2.add_argument('-nbufs',
                        required=False,
                        dest='nbufs',
                        action='store',
                        help='[<int>] Manually specify the number of buffer particles to add. If this flag is not set, a (more generous than necessarily required) estimate will be made based on the number of titratable sites. Currently N_buf = N_sites / 2q_max with q_max = 0.3.',
                        type=int)

    parser_2.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store_const',
                        const=1,
                        help='(no) Be more verbose (helpful for debugging).')

    parser_2.set_defaults(target='neutralize')

    parser_3  = subparsers.add_parser('genparams', help=desc_5)
    required3 = parser_3.add_argument_group("required arguments")

    required3.add_argument('-f', 
                        required=True,
                        dest='file',
                        action='store',
                        help='[<.pdb/.gro>] (required) Specify input structure file.')

    required3.add_argument('-ph', 
                        required=True,
                        dest='ph',
                        action='store',
                        help='[<real>] (required) Specify simulation pH.',
                        type=float)

    parser_3.add_argument('-mdp',
                        required=False,
                        dest='mdp',
                        action='store',
                        help='[<.mdp>] (MD.mdp) Specify .mdp file for the constant-pH parameters to be appended to. If the specified file does not exist, the .mdp file will be generated from scratch. Note that this only applies to production (MD), for energy minimization (EM) and equilibration (NVT/NPT), the .mdp files will be generated from scratch regardless.')

    parser_3.add_argument('-ndx',
                        required=False,
                        dest='ndx',
                        action='store',
                        help='[<.idx>] (index.ndx) Specify .ndx file for the constant-pH (lambda) groups to be appended to. If the specified file does not exist, the .ndx file will be generated from scratch.')

    parser_3.add_argument('-nstout',
                        required=False,
                        dest='nstout',
                        action='store',
                        default=500,
                        help='[<int>] (500) Specify output frequency for the lambda_xxx.dat files. 500 is large enough for subsequent frames to be uncoupled.',
                        type=int)

    parser_3.add_argument('-dwpE',
                        required=False,
                        dest='dwpE',
                        action='store',
                        default=7.5,
                        help='[<real>] (7.5) Specify default height of bias potential barrier in kJ/mol. 7.5 should be large enough in most cases, but if you observe a lambda coordinate spending a signficant amount of time between physical (i.e. lambda = 0/1) states, you should manually increase (either directly in the .mdp file or by setting the -inter flag).',
                        type=float)

    parser_3.add_argument('-lmass',
                        required=False,
                        dest='lmass',
                        action='store',
                        default=5,
                        help='[<real>] (5.0) Specify mass of the lambda particle(s). The user should probably not touch this.',
                        type=float)

    parser_3.add_argument('-ltau',
                        required=False,
                        dest='ltau',
                        action='store',
                        default=2,
                        help='[<real>] (2.0) Specify thermostat coupling time for the lambda-particles. The user should probably not touch this.',
                        type=float)

    parser_3.add_argument('-inter',
                        required=False,
                        dest='inter',
                        action='store_const',
                        const=1,
                        help='(no) If this flag is set, the user can manually specify the height of the bias potential barrier (in kJ/mol) for every titratable group.')

    parser_3.add_argument('-v',
                        required=False,
                        dest='verbosity',
                        action='store_const',
                        const=1,
                        help='(no) Be more verbose (helpful for debugging).')

    parser_3.set_defaults(target='genparams')

    argcomplete.autocomplete(parser)    # Required for autocompleting using argcomplete.
    CLI = parser.parse_args()           # Do the actual parsing.

    # Prevent Python errors when no subcommand is specified.
    if vars(CLI) == {}:
        parser.print_help()
        quit()

    return CLI # object containing all the parsed key-value pairs.
