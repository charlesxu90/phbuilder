<b>Description</b>
<p>Python-based system builder for constant-pH simulations in GROMACS.</p>

<b>Install instructions</b>
1. If you have a GPU and want to use GPU-acceleration, make sure you first install <a href="https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions">CUDA</a>.
2. Clone <a href="https://bitbucket.org/berkhess/gromacs-constantph/branch/clean-cpHMD-branch">clean-cpHMD-branch</a>.
3. Install using the instructions <a href="https://manual.gromacs.org/documentation/current/install-guide/index.html">here</a>. I personally use the following CMake flags:
`cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_constantph`
4. Clone phbuilder (this) repository.
5. Add phbuilder directory to PYTHONPATH by adding `export PYTHONPATH=$PYTHONPATH:/path/to/phbuilder` to your `~/.bashrc` file (and reload terminal(s)).

<b>phbuilder gentopol</b>

SYNOPSIS

`phbuilder gentopol [-h] -f FILE [-o OUTPUT] [-inter] [-list LIST] [-ph PH] [-v {0,1,2,3}]`

DESCRIPTION

<p>This tool encapsulates pdb2gmx and allows you to (re)generate the topology for your protein using our modified version of the CHARMM36 force field. This is necessary as some dihedral parameters were modified for titratable residues (ref manuscript 2). Gentopol also allows you to interactively set the initial lambda (protonation) state for each residue associated with a defined lambdagrouptype (for this, the `-inter` flag should be set). If the `-inter` flag is not set, the initial lambda values will be guessed based on an optional `-ph` flag that is by default set to `7.0`, together with the pKa defined in the `lambdagrouptypes.dat` file.</p>

OPTIONS

| Flag         | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-o`         | [\<.pdb/.gro>] (phprocessed.pdb) <br /> Specify output structure file. | 
| <nobr>`-inter`</nobr> | (no) <br /> Toggle interactive mode. | 
| `-list`      | [\<.txt>] <br /> Specify subset of residues to treat. | 
| `-ph`        | [\<real>] (7.0) <br /> Simulation pH. If the `-inter` flag is not set, this (together with the pKa specified in `lambdagrouptypes.dat`) will be used to guess the initial lambda state of the titratable residue(s).|
| `-v`         | [\<int>] (2) (phprocessed.pdb) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. | 

---

<b>phbuilder neutralize</b>

SYNOPSIS

`phbuilder neutralize [-h] -f FILE [-p TOPOL] [-o OUTPUT] [-solname SOLNAME] [-pname PNAME] [-nname NNAME] [-conc CONC] [-nbufs NBUFS] [-v {0,1,2,3}]`

DESCRIPTION

<p>The purpose of this tool is to ensure a charge-neutral system by adding the appropriate number of ions and buffer particles.</p>

OPTIONS

| Flag         | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-p`         | [\<.top>] (topol.top) <br /> Specify input topology file. |
| `-o`         | [\<.pdb/.gro>] (phneutral.pdb) <br /> Specify output structure file. |
| `-solname`   | [\<string>] (SOL) <br /> Specify solvent name (of which to replace molecules with ions and buffers). |
| `-pname`     | [\<string>] (NA) <br /> Specify name of positive ion to use. |
| `-nname`     | [\<string>] (CL) <br /> Specify name of negative ion to use. |
| `-conc`      | [\<real>] (0.0) <br /> Specify ion concentration in mol/L (will use solvent volume for calculating the required number of ions, not periodic box volume as `gmx genion` does). |
| `-nbufs`     | [\<int>] <br /> Manually specify the number of buffer particles to add. If this flag is not set, an estimate will be made based on the number of titratable sites. |
| `-v`         | [\<int>] (2) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. |

---

<b>phbuilder genparams</b>

SYNOPSIS

`phbuilder genparams [-h] -f FILE -ph PH [-mdp MDP] [-ndx NDX] [-nstout NSTOUT] [-dwpE DWPE] [-lmass LMASS] [-ltau LTAU] [-inter] [-v {0,1,2,3}]`

DESCRIPTION

`genparams` generates the `.mdp` files, including all the required constant-pH parameters.

OPTIONS

| Flag         | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. |
| `-ph`        | [\<real>] (required) <br /> Specify simulation pH. |
| `-mdp`       | [\<.mdp>] <br /> Specify existing .mdp file for production run to append the constant-pH parameters to, instead of generating a new MD.mdp from scratch. |
| `-ndx`       | [\<.idx>] <br /> Specify existing .ndx file to append the constant-pH groups to, instead of generating a new index.ndx from scratch.
| `-nstout`    | [\<int>] (500) <br /> Specify output frequency (in simulation steps) for the lambda_xxx.dat files. 500 is large enough for subsequent frames to be uncoupled.
| `-dwpE`      | [\<real>] (7.5) <br /> Specify default height of bias potential barrier in kJ/mol. 7.5 should be large enough in most cases, but if you observe a lambda coordinate spending a signficant amount of time between physical (i.e. lambda = 0/1) states, you should manually increase (either directly in the .mdp file or by setting the `-inter` flag).
| `-lmass`    | [\<real>] (5.0) <br /> Specify mass of the lambda particle(s). 5.0 is a good value (Why 5.0? something about autocorrelation times with carboxyl H vibrations?).
| `-ltau`     | [\<real>] (2.0) <br /> Specify thermostat coupling time for the lambda-particles. 2.0 ps^-1 is a good value).
| `-inter`     | (no) <br /> If this flag is set, the user can manually specify the height of the bias potential barrier (in kJ/mol) for every titratable group.
| `-v`         | [\<int>] (2) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. |
