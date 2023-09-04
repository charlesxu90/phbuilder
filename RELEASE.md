## 1.0.2
* Initial release.

## 1.1
* Changed the required Python version from 3.11 to 3.8 to make the installation more compatible.
* Environment variables are now passed to the GROMACS subprocess.
* phbuilder now uses the GROMACS path set in `lambdagrouptypes.dat` to check whether the correct GROMACS install was loaded in the environment phbuilder was called from, rather than using this path to execute GROMACS and disregarding any environment variables.
* groupnames in lambdagrouptypes.dat are now sanitized properly (they should be maximum 4 chars in order to adhere to the `.pdb` format).
* Add `-ignw` option for `neutralize`. When set, it passes `-maxwarn -1` to the encapsulated `grompp` calls.
* Significantly improved user and tool description messages.
* Changed default number of buffers from N_sites/2q_max to 2N_sites + 1. This is safer for inexperienced users.
* Output a trajectory frame every 50 ps instead of every 10 ps as a default value in the generic MD.mdp.
* Change the default simulation time from 100 ps to 100 ns, and to 1 ns when `-cal` is set.

## 1.2
* Change errors about LD_LIBRARY_PATH and PATH to warnings as some may install the GROMACS library in a non-default location.
