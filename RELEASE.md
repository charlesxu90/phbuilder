## 1.2 (DATE)
* TODO : add two-state ARG, LYS, and TYR parameters to lambdagrouptypes.dat and update the force field and residuetypes.dat.
* Improved EQ_smart.py.
* phbuilder no-longer checks/complains about which GROMACS version was sourced, we now only error if nothing was set in the environment at all.
* Removed fourier spacing comment "CHARMM is calibrated for 0.14." to prevent confusion. A value of `0.14` is perfectly fine, however CHARMM technically wasn't 'calibrated' for this value.
* Clarified some things in the parameterization section in the README.
* Fixed bug where gentopol would crash if the structure does not contain a TITLE and no residues are recognized.

## 1.1 (August 26th, 2023)
* Changed the required Python version from 3.11 to 3.8 to make the installation more compatible.
* Environment variables are now passed to the GROMACS subprocess.
* phbuilder now uses the GROMACS path set in `lambdagrouptypes.dat` to check whether the correct GROMACS install was loaded in the environment phbuilder was called from, rather than using this path to execute GROMACS and disregarding any environment variables.
* groupnames in lambdagrouptypes.dat are now sanitized properly (they should be maximum 4 chars in order to adhere to the `.pdb` format).
* Add `-ignw` option for `neutralize`. When set, it passes `-maxwarn -1` to the encapsulated `grompp` calls.
* Significantly improved user and tool description messages.
* Changed default number of buffers from N_sites/2q_max to 2N_sites + 1. This is safer for inexperienced users.
* Output a trajectory frame every 50 ps instead of every 10 ps as a default value in the generic MD.mdp.
* Change the default simulation time from 100 ps to 100 ns, and to 1 ns when `-cal` is set.

## 1.0.2 (August 11th, 2023)
* Initial release.
