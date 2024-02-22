# Description of the force field folder

The force field folder includes:
- The CHARMM36 force field with modified torsions for ASP, GLU, HIS, and C-terminal. Those changes were required due to the issues, described in [this paper](https://doi.org/10.1021/acs.jctc.2c00517).
- The corresponding `lambdagrouptypes.dat` file, which is used by `phbuilder` to prepare topologies for GROMACS constant pH runs. Currently, the following groups are included by default:
    - ASP
    - GLU
    - HSP
    - ARG
    - LYS
    - BUF

**NOTE** For histidine, currently $\lambda_1$ corresponds to protonated histidine (HSP), $\lambda_2$ - HSE/HIE, $\lambda_3$ - HSD/HID.

