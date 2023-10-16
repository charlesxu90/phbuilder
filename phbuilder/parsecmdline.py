# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete
import sys


def parsecmdline():
    """Parses command line using the python argparse module.

    Returns:
        CLI: object containing all the parsed key-value pairs.
    """

    # phbuilder main description.
    desc_1 = "System builder for constant-pH simulations in GROMACS. phbuilder consists of three tools: gentopol, neutralize, and genparams. Each tool performs a specific task for preparing a constant-pH simulation. Functionality for setting up titrations and parameterizations is provided with the help of stand-alone Python scripts, provided on the gitlab. Out of the box, phbuilder comes with the force field and CpHMD topology parameters required for setting up titratable Arg, Lys, Asp, Glu, and His residues in CHARMM36m."

    # Epilogue. Also used by subcommands.
    desc_2 = "phbuilder VERSION 1.2.1. For the user manual visit https://gitlab.com/gromacs-constantph/phbuilder."

    # gentopol main description.
    desc_3 = "Allows you to select which residues to make titratable and which initial lambda (protonation) state they should have. Also (re)generates the system topology using the modified CHARMM36m force field. If you don't want to manually set initial lambda values, you can use the -ph flag to have gentopol automatically choose the appropriate initial lambda values, based on the criterion: pH > pKa means start deprotonated, else start protonated."

    # neutralize main description.
    desc_4 = "Adds the appropriate number of ions to ensure a net-neutral system at t=0, and adds buffer particles in order to maintain a net-neutral system at t>0. The system charge a t=0 depends on the chosen initial lambda (protonation) states. At t>0, protonation states can change dynamically, meaning the resulting charge difference needs to be 'absorbed' by buffer particles. Each buffer particle can absorb up to ±0.5 charge, and charge is distributed evenly across all buffers (-10 system charge and 100 BUF implies every BUF +0.1). The GROMACS CpHMD beta MUST be sourced/loaded for neutralize to work correctly."

    # genparams main description.
    desc_5 = "Generates the CpHMD-specific .mdp and .ndx files. Will write generic EM.mdp EQ.mdp, and MD.mdp files for CHARMM36m and append the CpHMD parameters at the bottom. genparams requires the existence of a phrecord.dat file, which contains the initial lambda values and is created during the gentopol step. Note: if you previously used the auto feature (-ph flag) for gentopol, the pH you specify for genparams should reflect this."

    parser = argparse.ArgumentParser(prog='phbuilder', description=desc_1, epilog=desc_2)

    subparsers = parser.add_subparsers(required=False)

    parser_1 = subparsers.add_parser('gentopol', help=desc_3, epilog=desc_2)
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

    parser_1.add_argument('-ph',
                          required=False,
                          dest='ph',
                          action='store',
                          help='[<float>] Specify intended simulation pH. Will be used together with the macroscopic pKas of the titratable sites to auto set the initial lambdas.',
                          type=float)

    parser_1.add_argument('-v',
                          required=False,
                          dest='verbosity',
                          action='store_const',
                          const=1,
                          help='(no) Be more verbose.')

    parser_1.set_defaults(target='gentopol')

    parser_2 = subparsers.add_parser('neutralize', help=desc_4, epilog=desc_2)
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
                          help=' [<str>] (SOL) Specify solvent name (of which to replace molecules with ions and buffers).')

    parser_2.add_argument('-pname',
                          required=False,
                          dest='pname',
                          action='store',
                          default='NA',
                          help='[<str>] (NA) Specify name of positive ion to use.')

    parser_2.add_argument('-nname',
                          required=False,
                          dest='nname',
                          action='store',
                          default='CL',
                          help='[<str>] (CL) Specify name of negative ion to use.')

    parser_2.add_argument('-conc',
                          required=False,
                          dest='conc',
                          action='store',
                          default=0.0,
                          help='[<float>] (0.0) Specify ion concentration in mol/L. Note: uses solvent volume for calculating the required number of ions, not the periodic box volume as gmx genion does.',
                          type=float)

    parser_2.add_argument('-nbufs',
                          required=False,
                          dest='nbufs',
                          action='store',
                          help="[<int>] Specify number of buffer particles to add. If not set, N_buf = 2N_sites + 1. This ensures enough buffer particles will always be added, although you can likely get away with much less for larger systems.",
                          type=int)

    parser_2.add_argument('-rmin',
                          required=False,
                          dest='rmin',
                          action='store',
                          default=0.60,
                          help='[<float>] (0.6) Set the minimum distance the ions and buffers should be placed from the solute.',
                          type=float)

    parser_2.add_argument('-ignw',
                          required=False,
                          dest='ignw',
                          action='store_const',
                          const=1,
                          help='(no) Ignore all gmx grompp warnings.')

    parser_2.add_argument('-v',
                          required=False,
                          dest='verbosity',
                          action='store_const',
                          const=1,
                          help='(no) Be more verbose.')

    parser_2.set_defaults(target='neutralize')

    parser_3 = subparsers.add_parser('genparams', help=desc_5, epilog=desc_2)
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
                           help='[<float>] (required) Specify simulation pH.',
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
                          help="[<.idx>] (index.ndx) Specify .ndx file to append the CpHMD index groups to. If not set or the specified file does not exist, a generic 'index.ndx' will be created first.")

    parser_3.add_argument('-nstout',
                          required=False,
                          dest='nstout',
                          action='store',
                          default=500,
                          help='[<int>] (500) Specify lambda coordinate output frequency. 500 is large enough for subsequent frames to be uncoupled (with a dt = 0.002).',
                          type=int)

    parser_3.add_argument('-dwpE',
                          required=False,
                          dest='dwpE',
                          action='store',
                          default=7.5,
                          help='[<float>] (7.5) Specify default height of bias potential barrier (kJ/mol).',
                          type=float)

    parser_3.add_argument('-inter',
                          required=False,
                          dest='inter',
                          action='store_const',
                          const=1,
                          help='(no) Interactively set the height of the bias potential barrier (kJ/mol) for every titratable site.')

    parser_3.add_argument('-cal',
                          required=False,
                          dest='cal',
                          action='store_const',
                          const=1,
                          help="(no) Run CpHMD simulation in calibration mode: forces on the lambda coordinates are computed, but their positions won't be updated. This is only used for parameterization purposes.")

    parser_3.add_argument('-v',
                          required=False,
                          dest='verbosity',
                          action='store_const',
                          const=1,
                          help='(no) Be more verbose.')

    parser_3.set_defaults(target='genparams')

    # Required for auto-completing using argcomplete.
    argcomplete.autocomplete(parser)

    # Do the actual parsing.
    CLI = parser.parse_args()

    # Prevent Python errors when no subcommand is specified.
    if vars(CLI) == {}:
        parser.print_help()
        sys.exit(0)

    # Object containing all the parsed key-value pairs.
    return CLI
