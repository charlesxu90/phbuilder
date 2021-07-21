firstLine = True # For formatting of title

# Write a default .mdp file
def gen_mdp(Type, nsteps, nstxout, restrainCharge):

    # Sanitize input for Type
    if (Type not in ['EM', 'NVT', 'NPT', 'MD']):
        raise Exception("Unknown .mdp Type specified. Types are: EM, NVT, NPT, MD.")

    # Open file
    file = open("{0}.mdp".format(Type), 'w')

    # Formatting function for titles
    def addTitle(title):
        global firstLine
        if firstLine:
            file.write("; {0}\n".format(title.upper()))
            firstLine = False
        else:
            file.write("\n; {0}\n".format(title.upper()))

    # Formatting function for parameters
    def addParam(name, value, comment=''):
        if (comment == ''):
            file.write("{:20s} = {:13s}\n".format(name, str(value)))
        else:
            file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    # POSITION RESTRAINTS
    
    if Type in ['EM', 'MD'] and restrainCharge:
        addTitle('Position restraints')
        addParam('define', '-DPOSRES_BUF', 'Position restrain buffer(s).')

    if Type in ['NVT', 'NPT']:
        addTitle('Position restraints')
        if restrainCharge:
            addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restrain protein and buffer(s).')
        else:
            addParam('define', '-DPOSRES', 'Position restrain protein.')

    # RUN CONTROL

    addTitle("Run control")

    if (Type in ['EM']):
        dt = 0.01
        addParam('integrator', 'steep', 'Use steep for EM.')
        addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
        addParam('emstep', dt, 'Time step (ps).')

    if (Type in ['NVT', 'NPT', 'MD']):
        dt = 0.002
        addParam('integrator', 'md')
        addParam('dt', dt, 'Time step (ps).')

    addParam('nsteps', nsteps, "{:.1f} ns.".format((dt * nsteps)/1000.0))

    # CENTER OF MASS MOTION REMOVAL (USE WHEN WE HAVE BUFFERS)

    if (Type in ['MD'] and restrainCharge):
        addParam('comm-mode', 'Linear-acceleration-correction', 'Remove COM velocity and correct COM positions.')
        addParam('comm-grps', 'Protein Non-Protein')
        addParam("nstcomm", 10, 'Frequence for COM motion removal.')
        addParam("nstcalcenergy", 10, 'Required as nstcomm is otherwise set to nstcalcenergy.')

    # OUTPUT

    addTitle("Output control")
    addParam('nstxout-compressed', nstxout, 'Write .xtc frame every {:.1f} ps.'.format(dt * nstxout))

    # NEIGHBOUR SEARCHING

    addTitle("Neighbour searching")
    addParam('cutoff-scheme', 'Verlet', 'Related params are inferred by GROMACS.')

    # BONDED

    if (Type in ['NVT', 'NPT', 'MD']):
        addTitle("Bond parameters")
        addParam('constraints', 'h-bonds', 'Constrain H-bond vibrations.')
        addParam('constraint_algorithm', 'lincs', 'Holonomic constraints.')
        addParam('lincs_iter', 1, 'Related to accuracy of LINCS.')
        addParam('lincs_order', 4, 'Related to accuracy of LINCS.')

    # ELECTROSTATICS

    addTitle("Electrostatics")
    addParam('coulombtype', 'PME', 'Particle Mesh Ewald electrostatics.')
    addParam('rcoulomb', 1.2, 'CHARMM is calibrated for 1.2 nm.')
    addParam('fourierspacing', 0.14, 'CHARMM is calibrated for 0.14.')

    # VAN DER WAALS

    addTitle("Van der Waals")
    addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')
    addParam('rvdw', 1.2, 'CHARMM is calibrated for 1.2 nm.')
    addParam('vdw-modifier', 'force-switch', 'Specific for CHARMM.')
    addParam('rvdw-switch', 1.0, 'Specific for CHARMM.')

    # TEMPERATURE COUPLING

    if (Type in ['NVT', 'NPT', 'MD']):
        addTitle("Temperature coupling")
        addParam('tcoupl', 'v-rescale')
        addParam('tc-grps', 'SYSTEM')
        addParam('tau-t', 0.5, 'Coupling strength.')
        addParam('ref-t', 300, 'Reference temperature (K).')

    # PRESSURE COUPLING
    
    if (Type in ['NPT', 'MD']):
        addTitle('Pressure coupling')
        
        if (Type in ['NPT']):
            addParam('pcoupl', 'Berendsen', 'Use Berendsen for NPT.')
        else:
            addParam('pcoupl', 'Parrinello-Rahman')
        
        addParam('pcoupltype', 'isotropic', 'Uniform scaling of box.')
        addParam('tau_p', 5.0, 'Coupling strength.')
        addParam('ref_p', 1.0, 'Reference pressure (bar).')
        addParam('compressibility', 4.5e-5, 'Isothermal compressbility of water.')
        addParam('refcoord_scaling', 'all', 'Required with position restraints.')

    # PERIODIC BOUNDARY CONDITIONS
    
    addTitle("Periodic boundary condition")
    addParam('pbc', 'xyz', 'Apply periodic boundary conditions.')

    # WRAP UP

    file.close()
