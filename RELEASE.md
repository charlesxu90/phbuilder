# 1.2.5

### Features

* Added support for pH-dependent pKa values in `phbuilder` and `create_titration.py`.

### Improvements

* Added `-titr` option to `phbuilder genparams` for preparing parameter files for titrations.

### Fixes

* Fixed dependencies.

# 1.2.4 (February 22nd, 2024)

### Improvements

* Added [README.md](phbuilder/ffield/README.md) for the CpHMD force field. #7
* Improved documentation for the [create_titration.py](scripts/create_titration.py) script. #9

### Fixes

* Changed labeling of the $\lambda_2$ and $\lambda_3$ states of HSPT. To clarify: $\lambda_2$ corresponds to HSE/HIE and $\lambda_3$ HSD/HID. #7

# 1.2.2 (January 25th, 2024)

### Fixes

* Preprint link updated to published JCIM article.
* Simulation input files link updated to Zenodo.

# 1.2.1 (October 16th, 2023)

### Fixes

* The CpHMD charges specified in lambdagrouptypes.dat are now no-longer sanitized (previously only charges in the [-1, 1] interval were allowed).

# 1.2 (October 5th, 2023)

### Features
* Out of the box, phbuilder now also includes single-site parameterizations for arginine and lysine in CHARMM36m. The force field modifications and parameterization for lysine were [previously published](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00517) but not yet included in the phbuilder package, and the single-site parameterization for arginine was performed as part of the phbuilder paper.

### Improvements
* phbuilder no-longer checks/complains about which GROMACS version was sourced, we now only error if nothing was set in the environment at all.
* Improved the [EQ_smart.py](scripts/EQ_smart.py) script.
* Clarified some things in the parameterization section of the README.

### Fixes
* Removed "CHARMM is calibrated for 0.14." .mdp comment. This fourierspacing value is fine, however CHARMM technically wasn't 'calibrated' for it.
* Fixed bug where gentopol would crash if the structure does not contain a TITLE and no residues are recognized.
* Fixed bug where genparams would write the lambda-dynamics-simulation-ph only up to one digit after the decimal point.
* Fixed bug where genparams would write dV/dl coefficients and charge values only up to three digits after the decimal point.

# 1.1 (August 26th, 2023)
* Changed the required Python version from 3.11 to 3.8 to make the installation more compatible.
* Environment variables are now passed to the GROMACS subprocess.
* phbuilder now uses the GROMACS path set in `lambdagrouptypes.dat` to check whether the correct GROMACS install was loaded in the environment phbuilder was called from, rather than using this path to execute GROMACS and disregarding any environment variables.
* groupnames in lambdagrouptypes.dat are now sanitized properly (they should be maximum 4 chars in order to adhere to the `.pdb` format).
* Add `-ignw` option for `neutralize`. When set, it passes `-maxwarn -1` to the encapsulated `grompp` calls.
* Significantly improved user and tool description messages.
* Changed default number of buffers from N_sites/2q_max to 2N_sites + 1. This is safer for inexperienced users.
* Output a trajectory frame every 50 ps instead of every 10 ps as a default value in the generic MD.mdp.
* Change the default simulation time from 100 ps to 100 ns, and to 1 ns when `-cal` is set.

# 1.0.2 (August 11th, 2023)
* Initial release (not available on PyPi).
