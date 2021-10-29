<b>Description</b>
<p>Python-based system builder for constant-pH simulations in GROMACS.</p>

---

<b>Install instructions</b>
1. If you have a GPU and want to use GPU-acceleration, make sure you first install <a href="https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions">CUDA</a>.
2. Clone <a href="https://bitbucket.org/berkhess/gromacs-constantph/branch/clean-cpHMD-branch">clean-cpHMD-branch</a>.
3. Install using the instructions <a href="https://manual.gromacs.org/documentation/current/install-guide/index.html">here</a>. I personally use the following CMake flags:
`cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_constantph`
4. Clone phbuilder (this) repository.
5. Add phbuilder directory to PYTHONPATH by adding `export PYTHONPATH=$PYTHONPATH:/path/to/phbuilder` to your `~/.bashrc` file (and reload terminal(s)).

---

<b>phbuilder gentopol</b>

SYNOPSIS

`phbuilder gentopol [-h] -f FILE [-o OUTPUT] [-auto] [-list LIST] [-ph PH] [-v {0,1,2,3}]`

DESCRIPTION

This tool encapsulates pdb2gmx and allows you to (re)generate the topology for your system using our modified version of the CHARMM36 force field. This is necessary as some dihedral parameters were modified for titratable residues (ref manuscript 2). gentopol by default allows you to interactively set the initial lambda (protonation) state for each residue associated with a defined lambdagrouptype. This behavior can be automated by setting the `-auto` flag. In this case, every residue associated with a defined lambdagrouptype will automatically be made titratable, and the initial lambda values will be guessed based on an optional `-ph` flag that is by default set to `7.0`, together with the pKa defined in the `lambdagrouptypes.dat` file.

OPTIONS

| Flag__       | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-o`         | [\<.pdb/.gro>] (phprocessed.pdb) <br /> Specify output structure file. | 
| `-auto`      | (no) <br /> Toggle automatic mode. | 
| `-list`      | [\<.txt>] <br /> Specify subset of residues to treat. | 
| `-ph`        | [\<real>] (7.0) <br /> Simulation pH. If the `-auto` flag is set, this (together with the pKa specified in `lambdagrouptypes.dat`) will be used to guess the initial lambda state of the titratable residue(s).|
| `-v`         | [\<int>] (2) (phprocessed.pdb) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. | 

---

<b>phbuilder neutralize</b>

SYNOPSIS

`phbuilder neutralize [-h] -f FILE [-p TOPOL] [-o OUTPUT] [-solname SOLNAME] [-pname PNAME] [-nname NNAME] [-conc CONC] [-nbufs NBUFS] [-v {0,1,2,3}]`

DESCRIPTION

The purpose of this tool is to ensure a charge-neutral system by adding the appropriate number of ions and buffer particles.

OPTIONS

| Flag______   | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-p`         | [\<.top>] (topol.top) <br /> Specify input topology file. |
| `-o`         | [\<.pdb/.gro>] (phneutral.pdb) <br /> Specify output structure file. |
| `-solname`   | [\<string>] (SOL) <br /> Specify solvent name (of which to replace molecules with ions and buffers). |
| `-pname`     | [\<string>] (NA) <br /> Specify name of positive ion to use. Analogous to `gmx genion`.|
| `-nname`     | [\<string>] (CL) <br /> Specify name of negative ion to use. Analogous to `gmx genion`. |
| `-conc`      | [\<real>] (0.0) <br /> Specify ion concentration in mol/L. Analogous to `gmx genion` but will use solvent volume for calculating the required number of ions, not periodic box volume as `gmx genion` does). |
| `-nbufs`     | [\<int>] <br /> Manually specify the number of buffer particles to add. If this flag is not set, a (more generous than necessarily required) estimate will be made based on the number of titratable sites. |
| `-v`         | [\<int>] (2) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. |

---

<b>phbuilder genparams</b>

SYNOPSIS

`phbuilder genparams [-h] -f FILE -ph PH [-mdp MDP] [-ndx NDX] [-nstout NSTOUT] [-dwpE DWPE] [-lmass LMASS] [-ltau LTAU] [-inter] [-v {0,1,2,3}]`

DESCRIPTION

genparams generates the .mdp files, including all the required constant-pH parameters. genparams requires the existance of a record.dat file for setting the initial lambda values.

OPTIONS

| Flag_____    | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. |
| `-ph`        | [\<real>] (required) <br /> Specify simulation pH. |
| `-mdp`       | [\<.mdp>] <br /> Specify existing .mdp file for production run to append the constant-pH parameters to, instead of generating a new MD.mdp from scratch. |
| `-ndx`       | [\<.idx>] <br /> Specify existing .ndx file to append the constant-pH groups to, instead of generating a new index.ndx from scratch.
| `-nstout`    | [\<int>] (500) <br /> Specify output frequency (in simulation steps) for the lambda_xxx.dat files. 500 is large enough for subsequent frames to be uncoupled.
| `-dwpE`      | [\<real>] (7.5) <br /> Specify default height of bias potential barrier in kJ/mol. 7.5 should be large enough in most cases, but if you observe a lambda coordinate spending a signficant amount of time between physical (i.e. lambda = 0/1) states, you should manually increase (either directly in the .mdp file or by setting the `-inter` flag).
| `-lmass`     | [\<real>] (5.0) <br /> Specify mass of the lambda particle(s). 5.0 is a good value.
| `-ltau`      | [\<real>] (2.0) <br /> Specify thermostat coupling time for the lambda-particles. 2.0 ps^-1 is a good value.
| `-inter`     | (no) <br /> If this flag is set, the user can manually specify the height of the bias potential barrier (in kJ/mol) for every titratable group.
| `-v`         | [\<int>] (2) <br /> Verbosity: 0 = no output, 1 = errors and warnings only, 2 = default, 3 = be more verbose. |

---

<b>Basic workflow</b>

1. Prepare your structure file. <br /> This is an important step, and it applies especially to .pdbs that are straight from rcsb.org. Make sure your structure file only contains one MODEL, does not contain alternate location indicators, does not miss certain atoms/residues, etc. etc. It is also strongly recommended that every (non-ion/water) molecule has a chain identifier. If you have only one chain or do not care about this, you can simply set everything to A. One basic check to see if your input file contains mistakes can be to simply run: <br /> `gmx editconf -f input.pdb -o test.pdb` <br /> and see if you get any errors and how test.pdb differs from input.pdb.

2. Check whether your structure file contains any moleculetypes that <i>are</i> part of CHARMM36, but for which `pdb2gmx` <i>cannot</i> generate the topology data. Take as an example the lipid POPC. A CHARMM36 topology for POPC exists in the form of a stand-alone POPC.itp file, but if you supply just POPC.pdb to `gmpdb2gmx`, it won't be able to generate topol.top. If your structure file contains such residues, phbuilder can still be used but you'll be prompted for the path to POPC.itp when gentopol is called.

3. Decide which residues you want to have titratable, and in which protonation state those residues should be at t = 0 (i.e. which initial lambda values they should have). If you do not care about this, you can set the `-auto` flag to have gentopol automatically choose the appropriate initial lambda values based on the system pH and pKa.

3. (Re)generate the topology using: <br /> `phbuilder gentopol -f input.pdb` <br /> Alternatively, you can set the `-auto` flag and run: <br /> `phbuilder gentopol -f input.pdb -auto` <br /> In this case, the initial lambda values will be guessed based on the system ph (optional flag for gentopol by default set to 7.0) together with the pKa specified in lambdagrouptypes.dat.

4. Add a periodix box (if not already present) by e.g.: <br /> `gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5`

5. Add solvent (if not already present) by e.g.: <br /> `gmx solvate -cp box.pdb -p topol.top -o solvated.pdb`

6. Add the appropriate number of positive/negative ions and buffers to ensure a net-neutral system: <br /> `phbuilder neutralize -f solvated.pdb` <br /> Alternatively, if you want a specific ion concentration and/or a specific number of buffer particles , you could do: <br /> `phbuilder neutralize -f solvated.pdb -conc 0.15 -nbufs 20`

7. At this point, if everything went correctly both your structure and topology file(s) should be completed and constitute a net-neutral system when running cpHMD. What is now left is the actual simulation part: energy minimization, equilibration and production using the correct MD parameters.

8. Generate the .mdp files for EM/EQ/MD, including the constant-pH parameters for a specific simulation pH: 
<br /> `phbuilder genparams -f phneutral.pdb -ph 4.0` <br /> This will write the following files: EM.mdp, NVT.mdp, NPT.mdp, MD.mdp, and index.ndx.
9. Check the generated files and modify parameters specific to your system as required. For example if your system explodes upon starting NVT, you might need to adjust the time step for NVT coupling. Also note that by default no position restraints are used for the protein during NVT and NPT coupling.
10. Perform equilibration using a script such as:
```
#!/bin/bash

# GROMACS version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 2
gmx mdrun -v -deffnm EM -c EM.pdb

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r EM.pdb -maxwarn 1
gmx mdrun -v -deffnm NVT -c NVT.pdb -notunepme

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r NVT.pdb -maxwarn 1
gmx mdrun -v -deffnm NPT -c NPT.pdb -notunepme
```

11. Perform production run using a script such as:
```
#!/bin/bash

# GROMACS version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f MD.mdp -c NPT.pdb -p topol.top -n index.ndx -o MD.tpr -maxwarn 1
gmx mdrun -v -deffnm MD -c MD.pdb -x MD.xtc
```

---

To-do

* Improve the currently-written sections.
* Write a section with advanced things you can do.
* Implement the gmx-api for handling GROMACS calls.

---

FAQ

1. No questions here yet.
