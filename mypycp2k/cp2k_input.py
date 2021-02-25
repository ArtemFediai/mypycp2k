def set_global(cp2k_input,
               project_name='my project',
               print_level='MEDIUM',
               run_type='ENERGY'):
    GLOBAL = cp2k_input.GLOBAL
    GLOBAL.Print_level = print_level
    GLOBAL.Project_name = project_name
    GLOBAL.Extended_fft_lengths = True
    GLOBAL.Run_type = run_type


# def set_force_eval(FORCE_EVAL):
#     FORCE_EVAL.Method = 'QUICKSTEP'
#
#
# def set_global_and_force_eval(cp2k_input, project_name='my project'):
#     set_global(cp2k_input, project_name=project_name)
#     set_force_eval(cp2k_input)
