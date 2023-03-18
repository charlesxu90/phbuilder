<b>Description</b>

System builder for constant-pH simulations in [GROMACS](https://www.gromacs.org/). phbuilder consists of three tools: gentopol, neutralize, and genparams. Each tool performs a specific task for preparing a constant-pH simulation.

---

<b>Installation instructions</b>

0. This works for Linux and should also work for macOS. If you're on Windows, it is strongly recommended to use [WSL](https://docs.microsoft.com/en-us/windows/wsl/about).

1. If you have a GPU and want to use GPU-acceleration, make sure you first install [CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions).

2. Obtain the [GROMACS constant-pH beta build](https://gitlab.com/gromacs-constantph/).

3. Install using the instructions [here](https://manual.gromacs.org/documentation/current/install-guide/index.html). Suggested CMake command:
    ```
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_constantph
    ```
    The default path phbuilder will look for the GROMACS CpHMD install is `/usr/local/gromacs_constantph`. If GROMACS was installed in a different location, you are required to set the `GMXPH_BASEPATH` environment variable (e.g. in your `~/.bashrc`).
4. Install phbuilder using:
    ```
    pip3 install --index-url https://test.pypi.org/simple/ --no-deps phbuilder
    ```
    This is a testing version of the official PyPI repository (to-be-created later). The installation includes the modified CHARMM36M force field as well as all lambdagrouptypes.dat file, which contains the CpHMD-topology parameters.
5. phbuilder has [argcomplete](https://pypi.org/project/argcomplete/) functionality. To make sure this works, you should run:
    ```
    activate-global-python-argcomplete --user 
    ```
    once (and reload your terminal(s)).

<!-- 5. Clone phbuilder (this) repository.
6. Set the appropriate environment variables in `~/.bashrc` or `~/.zprofile`: <br /> `export PATH=$PATH:/path/to/clone/phbuilder` <br /> `export PHFFIELD=/path/to/clone/ffield` <br /> The ffield dir contains the modified CHARMM36M force field, as well as the lambdagrouptypes.dat file containing CpHMD specific topology data.
7. Make sure that the (base)path to your GROMACS constant-pH build is set correctly in lambdagrouptypes.dat (default = `/usr/local/gromacs_constantph`). -->

<!-- <b>Required Python packages</b>

* argcomplete (requires more than just installing in pip, look into this)
* argparse (should be a simple `pip3 install argparse`)
* configparser (should be a simple `pip3 install configparser`)
* os (should be a simple `pip3 install os`)
* subprocess (should be a simple `pip3 install subprocess`) -->

---

<b>phbuilder gentopol</b>

SYNOPSIS
```
phbuilder gentopol [-h] -f FILE [-o OUTPUT] [-list LIST] [-ph PH] [-v]
```
DESCRIPTION

gentopol encapsulates [gmx pdb2gmx](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html) and allows you to (re)generate the topology for your system using our modified version of the CHARMM36M force field. This is necessary as some dihedral parameters were modified for titratable residues (ref manuscript 2). gentopol by default allows you to interactively set the initial lambda value (protonation state) for each residue associated with a defined lambdagrouptype. This behavior can be automated by setting the `-ph <ph>` flag. In this case, every residue associated with a defined lambdagrouptype will automatically be made titratable, and the initial lambda values will be guessed based on the specified `ph`, together with the pKa defined in the `lambdagrouptypes.dat` file. Note that you should use the same pH value for genparams.

LIMITATIONS

* It is important that your protein(s)/molecule(s) containing the titratable groups is <i>at the top</i> of your .pdb/.gro input file. So first the titratable protein(s), and only then solvent, ions, lipids, etc. (the order of the latter components shouldn't matter).

OPTIONS

| Flag___      | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-o`         | [\<.pdb/.gro>] (phprocessed.pdb) <br /> Specify output structure file. | 
| `-ph`        | [\<real>] <br /> Use automatic mode and specify the simulation pH to base guess for initial lambda values on. |
| `-list`      | [\<.txt>] <br /> Provide a subset of resid(ue)s to consider. Helpful if you do not want to manually go through many (unimportant) residues. |
| `-v`         | (no) <br /> Be more verbose (helpful for debugging). |

---

<b>phbuilder neutralize</b>

SYNOPSIS
```
phbuilder neutralize [-h] -f FILE [-p TOPOL] [-o OUTPUT] [-solname SOLNAME] [-pname PNAME] [-nname NNAME] [-conc CONC] [-nbufs NBUFS] [-v]
```

DESCRIPTION

The purpose of this tool is to ensure a charge-neutral system by adding the appropriate number of ions and buffer particles.

LIMITATIONS

* phbuilder neutralize only keeps track of one type of positive (default NA), and one type of negative (default CL) ion. If you have either no ions or only NA and CL in your input structure, things should work. If you have or want to use a different type, you can use the `-pname` and `-nname` options (see below). If you have or want multiple different types of ions in your system, phbuilder is not guaranteed to work.
* Similar to [gmx genion](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html), phbuilder neutralize neutralizes the system by <i>adding</i> ions to the input structure, not by removing or rebalancing existing ones. This implies the ion concentration in your output files cannot and will not be lower than the ion concentration in your input file.

OPTIONS

| Flag_______  | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. | 
| `-p`         | [\<.top>] (topol.top) <br /> Specify input topology file. |
| `-o`         | [\<.pdb/.gro>] (phneutral.pdb) <br /> Specify output structure file. |
| `-solname`   | [\<string>] (SOL) <br /> Specify solvent name (of which to replace molecules with ions and buffers). |
| `-pname`     | [\<string>] (NA) <br /> Specify name of positive ion to use. Analogous to [gmx genion](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html).|
| `-nname`     | [\<string>] (CL) <br /> Specify name of negative ion to use. Analogous to [gmx genion](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html). |
| `-conc`      | [\<real>] (0.0) <br /> Specify ion concentration in mol/L. Analogous to [gmx genion](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html) but will use the solvent volume for calculating the required number of ions, not the periodic box volume as genion does. |
| `-nbufs`     | [\<int>] <br /> Manually specify the number of buffer particles to add. If this flag is not set, a (more generous than necessarily required) estimate will be made based on the number of titratable sites. Currently N_buf = N_sites / 2q_max with q_max = 0.5.|
| `-v`         | (no) <br /> Be more verbose (helpful for debugging). |

---

<b>phbuilder genparams</b>

SYNOPSIS
```
phbuilder genparams [-h] -f FILE -ph PH [-mdp MDP] [-ndx NDX] [-nstout NSTOUT] [-dwpE DWPE] [-lmass LMASS] [-ltau LTAU] [-inter] [-v]
```
DESCRIPTION

genparams generates the .mdp files, including all the required constant-pH parameters. genparams requires the existence of a phrecord.dat file for setting the initial lambda values.

OPTIONS

| Flag_____    | Description    |
|--------------|----------------|
| `-f`         | [\<.pdb/.gro>] (required) <br /> Specify input structure file. |
| `-ph`        | [\<real>] (required) <br /> Specify simulation pH. |
| `-mdp`       | [\<.mdp>] (MD.mdp) <br /> Specify .mdp file for the constant-pH parameters to be appended to. If the specified file does not exist, the .mdp file will be generated from scratch. Note that this only applies to production (MD), for energy minimization (EM) and equilibration (NVT/NPT), the .mdp files will be generated from scratch regardless. |
| `-ndx`       | [\<.idx>] (index.ndx) <br /> Specify .ndx file for the constant-pH (lambda) groups to be appended to. If the specified file does not exist, the .ndx file will be generated from scratch (using `echo q \| gmx make_ndx -f input.pdb`). |
| `-nstout`    | [\<int>] (500) <br /> Specify output frequency for the lambda_xxx.dat files. 500 is large enough for subsequent frames to be uncoupled.
| `-dwpE`      | [\<real>] (7.5) <br /> Specify default height of bias potential barrier in kJ/mol. 7.5 should be large enough in most cases, but if you observe a lambda coordinate spending a significant amount of time between physical (i.e. lambda = 0/1) states, you should manually increase (either directly in the .mdp file or by setting the `-inter` flag).
| `-lmass`     | [\<real>] (5.0) <br /> Specify mass of the lambda particle(s). The user should probably not touch this.
| `-ltau`      | [\<real>] (2.0) <br /> Specify thermostat coupling time for the lambda-particles. The user should probably not touch this.
| `-inter`     | (no) <br /> If this flag is set, the user can manually specify the height of the bias potential barrier (in kJ/mol) for every titratable group.
| `-cal`       | (no) <br /> If this flag is set, the CpHMD simulation will be run in calibration mode: forces on the lambdas are computed, but they will not be updated. This is used for calibration purposes. |
| `-v`         | (no) <br /> Be more verbose (helpful for debugging). |

---

<b>Basic workflow</b>

1. Prepare your structure file. <br /> This is an important step, and it applies especially to .pdbs that are straight from [rcsb](https://www.rcsb.org/). Make sure your structure file only contains one MODEL, does not contain alternate location indicators, does not miss certain atoms/residues, etc. etc. It is also important that your protein(s)/molecule(s) containing the titratable groups is <i>at the top</i> of your .pdb/.gro input file. So first the titratable protein(s), and only then solvent, ions, lipids, etc. (the order of the latter components shouldn't matter). Furthermore, it is strongly recommended that every (non-ion/water) molecule has a chain identifier. If you have only one chain or do not care about this, you can simply set everything to A. One basic check to see if your input file contains mistakes can be to simply run: 
    ```
    gmx editconf -f input.pdb -o test.pdb
    ```
    and see if you get any errors and how test.pdb differs from input.pdb.

2. Check whether your structure file contains any moleculetypes that <i>are</i> part of CHARMM36M, but for which `gmx pdb2gmx` <i>cannot</i> generate the topology data. Take as an example the lipid POPC. A CHARMM36M topology for POPC exists in the form of a stand-alone POPC.itp file, but if you supply just POPC.pdb to `gmx pdb2gmx`, it won't be able to generate topol.top. If your structure file contains such residues, phbuilder can still be used but you'll be prompted for the path to POPC.itp when gentopol is called.

3. Decide which residues you want to have titratable, and in which protonation state those residues should be at t = 0 (i.e. which initial lambda values they should have). If you do not care about this, you can use the `-ph <ph>` flag to have gentopol automatically choose the appropriate initial lambda values based on the system pH and pKas of the lambdagrouptypes.

4. (Re)generate the topology using:
    ```
    phbuilder gentopol -f input.pdb
    ``` 
    Alternatively, you can set the -ph flag and run:
    ```
    phbuilder gentopol -f input.pdb -ph <ph>
    ```
    In the latter case, the initial lambda values will be guessed based on the system ph together with the pKa specified in lambdagrouptypes.dat.

5. Add a periodic box (if not already present) by e.g. (see [gmx editconf](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html)):
    ```
    gmx editconf -f phprocessed.pdb -o box.pdb -bt cubic -d 1.5
    ```

6. Add solvent (if not already present) by e.g. (see [gmx solvate](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html)):
    ```
    gmx solvate -cp box.pdb -p topol.top -o solvated.pdb
    ```

7. (Optional) [gmx solvate](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html) is a relatively basic program and will by default add solvent simply based on the van der Waals radii of the solute. This can lead to water molecules being generated in locations where there should not be any (e.g. hydrophobic pockets). It is good practice to check this, and if this occurs in your system, we recommend you to utilize the included [clean_after_solvate.py](scripts/clean_after_solvate.py) script.

8. Add the appropriate number of positive/negative ions and buffers to ensure a net-neutral system:
    ```
    phbuilder neutralize -f solvated.pdb
    ```
    Alternatively, if you want a specific ion concentration (in mol/L) and/or a specific number of buffer particles , you could do:
    ```
    phbuilder neutralize -f solvated.pdb -conc 0.15 -nbufs 20
    ```
    Note that phbuilder neutralize neutralizes the system by <i>adding</i> ions to the input structure, not by removing or rebalancing existing ones. This implies the ion concentration in your output files cannot and will not be lower than the ion concentration in your input file.

9. At this point, if everything went correctly both your structure and topology file(s) should be completed and constitute a net-neutral system when running CpHMD. What is now left is the actual simulation part: energy minimization, equilibration and production using the correct MD parameters.

10. Generate the .mdp files for EM/EQ/MD, including the constant-pH parameters for a specific simulation pH: 
    ```
    phbuilder genparams -f phneutral.pdb -ph 4.0 
    ```
    By default, the following files will be written:
    * EM.mdp
    * NVT.mdp
    * NPT.mdp
    * MD.mdp
    * index.ndx

11. Check the generated files and modify parameters specific to your system as required. For example if your system explodes upon starting NVT, you might need to adjust the time step for NVT coupling. Also note that by default no position restraints are used for the protein during NVT and NPT coupling.

12. Perform equilibration using a script such as:
```
#!/bin/bash

# GROMACS version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1
gmx mdrun -v -deffnm EM -c EM.pdb

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr
gmx mdrun -v -deffnm NVT -c NVT.pdb -notunepme

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -n index.ndx -o NPT.tpr
gmx mdrun -v -deffnm NPT -c NPT.pdb -notunepme
```

13. Perform production run using a script such as:
```
#!/bin/bash

# GROMACS version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f MD.mdp -c NPT.pdb -p topol.top -n index.ndx -o MD.tpr
gmx mdrun -v -deffnm MD -c MD.pdb -x MD.xtc
```

---

<b>Performing titrations</b>

Performing a computational titration is helpful for determining the microscopic pKas of titratable sites. After steps 1 to 12 of the basic workflow have been completed, one can use the included [create_titration.py](scripts/create_titration.py) to setup a titration. For example, the command:

```
create_titration.py -f MD.mdp -c NPT.pdb -p topol.top -n index.ndx -pH 1:10:1 -nr 2
```

creates directories corresponding to pH 1 to 9, with each subdirectory containing two replicas (each containing the appropriate input files for `gmx mdrun`).

---

<b>Performing parameterizations</b>

The following section describes a procedure for parameterizing (new) *two-state* ligands for CpHMD simulations. For convenience, we will use the word ligand to refer to any new lambdagrouptype. In this workflow, we will consider parameterizing ATP as an example. Note that performing parameterizations correctly is relatively complicated, and the reader is advised to check [Scalable Constant pH Molecular Dynamics in GROMACS](https://pubs.acs.org/doi/10.1021/acs.jctc.2c00516) for more information on parameterization in CpHMD.

1. Prepare the residuetype topology. [gmx pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) cannot create topologies for non-standard residuetypes, so you most likely have to provide your own .itp file. This .itp has a number of specific requirements, outlined below:

    the [ moleculetype ] name should be equal to the name of the *titratable version* of the ligand. So in the case of ATP, we would have
    ```
    [ moleculetype ]
    ; name	nrexcl
    ATPT	     3
    ```
    in the ATP.itp file. Furthermore, the names in the [ atoms ] section should also reflect this:
    ```
    [ atoms ]
    ; nr	type	resnr	residu	atom	cgnr	charge	mass
    1        CN7      1      ATPT    C4'      1      0.160    12.0110   ; qtot  0.160
    etc.
    ```
    Finally, at the end of the .itp file there should be a position restraining section
    ```
    #ifdef POSRES
    #include "posresATPT.itp"
    #endif
    ```
    Where the ifdef must be POSRES. The posresATPT.itp file should contain the specific atom(s) of the ligand you want to restrain, e.g. posresATPT.itp:
    ```
    [ position_restraints ]
    ;  i funct       fcx        fcy        fcz
    40  1  1000  1000  1000
    ```
    will position restrain atom 40 (of moleculetype ATPT). Position restraining during the calibration is required to avoid (strong) interactions between the titratable atoms and the neutralizing buffer particle. If the ligand and buffer particle accidentally get close to each other in some of the calibration runs, the resulting **dvdl** coefficients will be significantly affected. It is also important to remember while selecting atoms for which positions are restrained, that we want to keep the distance between the titratable group and the buffer particle large, but at the same time we want to sample as much orientational configurations as possible. Thus, in the case of ATPT we only fix the phosphorus of $\gamma$-phosphate. The suitable selection of atoms to restrain is system-dependent and therefore the responsibility of the user.

2. In addition to providing and setting up your .itp file, you will likely have to make an addition to the CpHMD force field. This possibly includes new atom, bond, pair, angle, and dihedral types not present in the standard force field, but present in ligand topology. phbuilder does not take care of this and modifying the force field is the responsibility of the user.

3. Place a (copy of the default) lambdagrouptypes.dat file in your working directory. Doing so will override the default file, and this will allow you to define custom lambdagrouptypes. To the lambdagrouptypes.dat (in your working dir) add the parameters corresponding to the new type. For ATP:

    ```
    [ ATPT ]
    incl   = ATP
    atoms  = O3B PG O1G O2G H2G O3G
    qqA    = -0.98 1.50 -0.82 -0.68 0.34 -0.82
    pKa_1  = 4.5
    qqB_1  = -0.86 1.10 -0.90 -0.90 0.00 -0.90
    dvdl_1 = 0
    ```

    dvdl_1 is initially set to zero as this is the parameter we are going to obtain during the calibration.

4. Perform basic workflow steps 1 to 7 to obtain a solvate structure. Next, perform the neutralization step setting `-nbufs = 1`:

    ```
    phbuilder neutralize -f solvated.pdb -nbufs 1
    ```

5. Generate .mdp files for EM/EQ/MD in calibration mode by setting the additional `-cal` flag for genparams:

    ```
    phbuilder genparams -f phneutral.pdb -ph 4.0 -cal
    ```

    Setting the `-cal` flag will modify the resulting .mdp files. It will not only set 
    ```
    lambda-dynamics-calibration = yes
    ```
    but also add
    ```
    define = -DPOSRES -DPOSRES_BUF
    ```
    position restraints. Here, `DPOSRES` corresponds to the ligand atom(s) as described in step 1, and `DPOSRES_BUF` corresponds to the buffer. Finally, setting the `-cal` flag modifies the range and initial lambda for the buffer.

6. Perform basic workflow steps 11 and 12 (check generic .mdp files and perform EM+EQ).

7. Use the included [create_parameterization.py](scripts/create_parameterization.py) to setup the parameterization runs. For example, the command:

    ```
    create_parameterization.py -f MD.mdp -c NPT.pdb -r NPT.pdb -p topol.top -n index.ndx
    ```

    creates directories corresponding to different $\lambda$-values, each containing a `.tpr` run input file for `gmx mdrun`.

8. Extract cphmd-coordiante ... and use fitting script to obtain final dvdl coefficients...
