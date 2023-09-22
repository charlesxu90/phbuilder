firstLine = True  # For formatting of title.

def gen_mdp(Type: str, nsteps: int, nstxout: int, posRes=False):
    """Write a generic .mdp file.

    Args:
        Type (str): the .mdp file type. Choose EM, NVT, NPT, MD.
        nsteps (int): number of simulation steps.
        nstxout (int): output frequency (write once every nstxout simulation steps).
        posRes (bool, optional): enable position restraints (for in calibration mode). Defaults to False.
    """

    assert Type in ['EM', 'NVT', 'NPT', 'MD'], "type should be be EM, NVT, NPT, MD"

    file = open("{0}.mdp".format(Type), 'w')

    def addTitle(title: str):
        """Formatting function for Title.

        Args:
            title (str): title.
        """

        global firstLine
        if firstLine:
            file.write("; {0}\n".format(title.upper()))
            firstLine = False
        else:
            file.write("\n; {0}\n".format(title.upper()))

    def addParam(name: str, value, comment=''):
        """Formatting function for adding a parameter.

        Args:
            name (str): parameter name.
            value (any): parameter value.
            comment (str, optional): parameter inline comment. Defaults to ''.
        """

        if comment == '':
            file.write("{:20s} = {:13s}\n".format(name, str(value)))
        else:
            file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    # POSITION RESTRAINTS DURING CALIBRATION

    if posRes:
        addTitle('Position restraints')
        addParam('define', "-DPOSRES -DPOSRES_BUF", 'Position restrain if calibrating.')

    # RUN CONTROL

    addTitle("Run control")

    if Type in ['EM']:
        dt = 0.01
        addParam('integrator', 'steep', 'Use steep for EM.')
        addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
        addParam('emstep', dt, 'Time step (ps).')

    if Type in ['NVT', 'NPT', 'MD']:
        dt = 0.002
        addParam('integrator', 'md')
        addParam('dt', dt, 'Time step (ps).')

    addParam('nsteps', nsteps, "{:d} ps.".format(int(dt * nsteps)))

    # OUTPUT

    addTitle("Output control")
    addParam('nstxout-compressed', nstxout, 'Write .xtc frame every {:d} ps.'.format(int(dt * nstxout)))

    # NEIGHBOUR SEARCHING

    addTitle("Neighbour searching")
    addParam('cutoff-scheme', 'Verlet', 'Related params are inferred by GROMACS.')

    # BONDED

    if Type in ['NVT', 'NPT', 'MD']:
        addTitle("Bond parameters")
        addParam('constraints', 'h-bonds', 'Constrain H-bond vibrations.')
        addParam('constraint_algorithm', 'lincs', 'Holonomic constraints.')
        addParam('lincs_iter', 1, 'Related to accuracy of LINCS.')
        addParam('lincs_order', 4, 'Related to accuracy of LINCS.')

    # ELECTROSTATICS

    addTitle("Electrostatics")
    addParam('coulombtype', 'PME', 'Particle Mesh Ewald electrostatics.')
    addParam('rcoulomb', 1.2, 'CHARMM is calibrated for 1.2 nm.')
    addParam('fourierspacing', 0.14)

    # VAN DER WAALS

    addTitle("Van der Waals")
    addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')
    addParam('rvdw', 1.2, 'CHARMM is calibrated for 1.2 nm.')
    addParam('vdw-modifier', 'force-switch', 'Specific for CHARMM.')
    addParam('rvdw-switch', 1.0, 'Specific for CHARMM.')

    # TEMPERATURE COUPLING

    if Type in ['NVT', 'NPT', 'MD']:
        addTitle("Temperature coupling")
        addParam('tcoupl', 'v-rescale')
        addParam('tc-grps', 'SYSTEM')
        addParam('tau-t', 0.5, 'Coupling time (ps).')
        addParam('ref-t', 300, 'Reference temperature (K).')

    # PRESSURE COUPLING

    if Type in ['NPT', 'MD']:
        addTitle('Pressure coupling')
        addParam('pcoupl', 'C-rescale', 'Use C-rescale barostat.')
        addParam('pcoupltype', 'isotropic', 'Uniform scaling of box.')
        addParam('tau_p', 5.0, 'Coupling time (ps).')
        addParam('ref_p', 1.0, 'Reference pressure (bar).')
        addParam('compressibility', 4.5e-05, 'Isothermal compressbility of water.')

        if Type == 'NPT' or posRes:
            addParam('refcoord_scaling', 'com', 'Required with position restraints.')

    # PERIODIC BOUNDARY CONDITIONS

    addTitle("Periodic boundary condition")
    addParam('pbc', 'xyz', 'Apply periodic boundary conditions.')

    # WRAP UP

    file.close()
    global firstLine
    firstLine = True
