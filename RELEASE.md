## 1.0.2
* Initial release.

## 1.1
* Changed required Python version from 3.11 to 3.8 to make it more compatible.
* Environment variables are now passed to the GROMACS subprocess.
* phbuilder now uses the GROMACS path set in `lambdagrouptypes.dat` to check whether the correct GROMACS install was loaded in the environment phbuilder was called from, rather than using this path to execute GROMACS and disregarding any environment variables.
* groupnames in lambdagrouptypes.dat are now sanitized properly; they should be maximum 4 chars in order to adhere to the `.pdb` format.
* Add `-ignw` option for `neutralize`. When set, it passes `-maxwarn -1` to the encapsulated `grompp` calls.
