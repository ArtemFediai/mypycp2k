def set_global(cp2k_input,
               project_name='my project',
               print_level='LOW',
               run_type='ENERGY',
               extended_fft_length=False):
    GLOBAL = cp2k_input.GLOBAL
    GLOBAL.Print_level = print_level
    GLOBAL.Project_name = project_name
    GLOBAL.Extended_fft_lengths = extended_fft_length
    GLOBAL.Run_type = run_type


def set_geo_opt(
        cp2k_input,
        optimizer='BFGS',
        max_geo_opt_iter=200,
        max_force=4.5e-4,  # turbomole value
        rms_force=3.0e-4,  # turbomole value
        max_dr=3.0e-3
):
    MOTION = cp2k_input.MOTION
    GEO_OPT = MOTION.GEO_OPT
    GEO_OPT.Type = 'MINIMIZATION'
    GEO_OPT.Optimizer = optimizer
    GEO_OPT.Max_iter = max_geo_opt_iter

    # Setting convergence criteria
    GEO_OPT.Max_force = max_force
    GEO_OPT.Rms_force = rms_force
    GEO_OPT.Max_dr = max_dr
