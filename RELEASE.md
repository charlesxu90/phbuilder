## 1.0.2
* Initial release.

## 1.1
* Environment variables are now being passed to the subprocess running GROMACS and checked/modified.
* groupnames in lambdagrouptypes.dat are now sanitized properly; they should be maximum 4 chars in order to adhere to the .pdb format.
* Add an `-ignw` option for `neutralize` that when set passes `-maxwarn -1` to `grompp`.
