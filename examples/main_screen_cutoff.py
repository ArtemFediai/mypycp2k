from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, print_mo, set_qs
from mypycp2k.xc import set_pbe, set_pbe0, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_ot_never_fail
from mypycp2k.subsys import add_elements, set_unperiodic_cell, set_topology, center_coordinates
from check_convergence.check_convergence \
    import general_convergence, change_calc_abc, change_calc_cutoff, change_calc_rel_cutoff
import copy

# from check_convergence.check_convergence import check_abc_convergence, check_cutoff_convergence

"""
check, which values to select for: abc, cutoff, rel_cutoff
set_XXX means imperatively set some section
add_XXX means optionally set some section. Optionally == if activate_XXX is set to True
activate_XXX means activate "add_XXX"
"""
def main():
    # My input_from_yaml #
    ## base settings ##
    basis_set_base_path = '/home/artem/soft/cp2k/cp2k-7.1/data/'
    my_abc = '16.0 16.0 16.0'
    # my_basis_set_file_name = basis_set_base_path + 'BASIS_RI_cc-TZ'G
    my_basis_set_file_name = basis_set_base_path + 'BASIS_def2_QZVP_RI_ALL'
    my_vdw_parameters_file = basis_set_base_path + 'dftd3.dat'
    my_basis_set = 'def2-QZVP'
    my_potential_file_name = basis_set_base_path + 'POTENTIAL'
    my_potential = 'ALL'
    my_project_name = "methane_GW_PBE"
    my_xyz_file_name = 'methane.xyz'
    my_ri_aux_basis_set = 'RI-5Z'  #
    # organic_elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'B', 'I']
    # elements = ['H', 'O']  # for test
    # my_element = ['H']  # for test
    my_elements = ['H', 'C']
    activate_vdw = False

    ############################################# pycp2k ###############################################################
    calc = CP2K()
    calc.working_directory = './'
    calc.project_name = 'artem_gw_project'
    calc.mpi_n_processes = 1

    # pycp2k objects
    CP2K_INPUT = calc.CP2K_INPUT
    FORCE_EVAL = CP2K_INPUT.FORCE_EVAL_add()
    FORCE_EVAL.Method = 'QUICKSTEP'
    SUBSYS = FORCE_EVAL.SUBSYS
    DFT = FORCE_EVAL.DFT
    XC = DFT.XC
    SCF = DFT.SCF
    ################################################ mypycp2k ##########################################################

    # GLOBAL #
    # FORCE EVAL #
    set_global(CP2K_INPUT, project_name=my_project_name)
    #set_force_eval(FORCE_EVAL)

    ## SUBSYS ##
    set_unperiodic_cell(SUBSYS, abc=my_abc)
    set_nonperiodic_poisson(DFT)
    set_topology(SUBSYS, xyz_file_name=my_xyz_file_name)
    add_elements(SUBSYS,
                 elements=my_elements,
                 basis=my_basis_set,
                 aux_basis=my_ri_aux_basis_set,
                 pot=my_potential)
    center_coordinates(SUBSYS)
    ## END SUBSYS ##

    ## DFT ##
    set_dft(DFT,
            potential_file_name=my_potential_file_name,
            basis_set_file_name=my_basis_set_file_name)
    set_cutoff(DFT, cutoff=900, rel_cutoff=100, ngrids=5)
    set_scf(DFT, max_scf=1)
    add_ot(SCF)
    set_pbe(XC)  # alter: set_pb0, etc.
    set_qs(DFT,
           eps_default=1.0E-15,
           eps_pgf_orb=1.0E-200)
    # add_gw_ver_0(XC)  # GW!
    #print_mo_cubes(DFT.PRINT)
    #print_mo(DFT.PRINT)
    if activate_vdw:
        add_vdw(XC, vdw_parameters_file=my_vdw_parameters_file)
    ## END DFT ##

############################################### CHECK CONVERGENCE ######################################################

    my_target_accuracy_eV = 1.E-3
    if True:
        my_abcs = [8, 10, 12, 14, 16]
        my_cutoff = [100, 200, 300, 400, 500, 600, 700, 800, 900]
        my_rel_cutoff = [20, 40, 60, 70, 80, 90, 100]

    if False:
        my_abcs = [6, 8, 10]
        my_cutoff = [100, 200, 300]
        my_rel_cutoff = [40, 50, 60]

    print(f"Reference parameters are: \n "
          f"cutoff = {calc.CP2K_INPUT.FORCE_EVAL_list[0].DFT.MGRID.Cutoff},\n "
          f"rel_cut = {calc.CP2K_INPUT.FORCE_EVAL_list[0].DFT.MGRID.Rel_cutoff},\n "
          f"abc = {calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc}\n")

    # check abc (cell size) convergence
    general_convergence(calc=copy.deepcopy(calc),
                        params=my_abcs,
                        func_to_accept_param=change_calc_abc,
                        target_accuracy_eV=my_target_accuracy_eV,
                        param_name='abc')

    # check cutoff (absolute) convergence
    general_convergence(calc=copy.deepcopy(calc),
                        params=my_cutoff,
                        func_to_accept_param=change_calc_cutoff,
                        target_accuracy_eV=my_target_accuracy_eV,
                        param_name='cutoff')

    # check cutoff (relative) convergence
    general_convergence(calc=copy.deepcopy(calc),
                        params=my_rel_cutoff,
                        func_to_accept_param=change_calc_rel_cutoff,
                        target_accuracy_eV=my_target_accuracy_eV,
                        param_name='rel_cutoff')

# rel_cutoff: 40; cutoff: 300; abc = 10

if __name__ == '__main__':
    main()