<b>Description</b>
<p>Python-based system builder for constant-pH simulations in GROMACS.</p>

<b>Install instructions</b>
1. If you have a GPU and want to use GPU-acceleration, make sure you first install <a href="https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions">CUDA</a>.
2. Clone <a href="https://bitbucket.org/berkhess/gromacs-constantph/branch/clean-cpHMD-branch">clean-cpHMD-branch</a>.
3. Install using the instructions <a href="https://manual.gromacs.org/documentation/current/install-guide/index.html">here</a>. I personally use the following CMake flags:
`cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_constantph`
4. Clone phbuilder (this) repository.
5. Add phbuilder directory to PYTHONPATH by adding `export PYTHONPATH=$PYTHONPATH:/path/to/phbuilder` to your `~/.bashrc` file (and reload terminal(s)).

<b>Getting started</b>

Limitations: CHARMM36-mar2019 with tip3p water model. We have only parametrized aspartic acid, glutamic acid, histidine, N-terminus, C-terminus.

<b>The `lambdagrouptypes.dat` file</b>


---

<b>phbuilder gentopol</b>

SYNOPSIS

`phbuilder gentopol [-h] -f FILE [-o OUTPUT] [-inter] [-list LIST] [-ph PH] [-v {0,1,2,3}]`

DESCRIPTION

This tool encapsulates pdb2gmx and allows you to (re)generate the topology for your protein using our modified version of the CHARMM36 force field. This is necessary as some dihedral parameters were modified for titratable residues (ref manuscript 2). Gentopol also allows you to interactively set the initial lambda (protonation) state for each residue associated with a defined lambdagrouptype (for this, the `-inter` flag should be set). If the `-inter` flag is not set, the initial lambda values will be guessed based on an optional `-ph` flag that is by default set to `7.0`, together with the pKa defined in the `lambdagrouptypes.dat` file.

OPTIONS

| Flag         | Description    |
|--------------|----------------|
| `-f`         | [<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-o`         | [<.pdb/.gro>] (phprocessed.pdb) <br /> Specify output structure file. | 
| <nobr>`-inter`</nobr> | (no) <br /> Toggle interactive mode. | 
| `-list`      | [<.txt>] <br /> Specify subset of residues to treat. | 
| `-ph`        | [\<real>] (7.0) <br /> Simulation pH. If the `-inter` flag is not set, this (together with the pKa specified in `lambdagrouptypes.dat`) will be used to guess the initial lambda state of the titratable residue(s).|
| `-v`         | [\<int>] (2) (phprocessed.pdb) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. | 

---

<b>phbuilder neutralize</b>

SYNOPSIS

`phbuilder neutralize [-h] -f FILE [-p TOPOL] [-o OUTPUT] [-solname SOLNAME] [-pname PNAME] [-nname NNAME] [-conc CONC] [-nbufs NBUFS] [-v {0,1,2,3}]`

DESCRIPTION

Blabla on this tool.

OPTIONS

| Flag         | Description    |
|--------------|----------------|
| `-f`         | [<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-o`         | [<.pdb/.gro>] (phprocessed.pdb) <br /> Specify output structure file. | 
| <nobr>`-inter`</nobr> | (no) <br /> Toggle interactive mode. | 
| `-list`      | [<.txt>] <br /> Specify subset of residues to treat. | 
| `-ph`        | [\<real>] (7.0) <br /> Simulation pH. If the `-inter` flag is not set, this (together with the pKa specified in `lambdagrouptypes.dat`) will be used to guess the initial lambda state of the titratable residue(s).|
| `-v`         | [\<int>] (2) (phprocessed.pdb) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. | 
