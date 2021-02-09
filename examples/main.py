from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, print_mo, set_qs
from mypycp2k.xc import set_pbe, set_pbe0, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_ot_never_fail, add_diagonalization, add_mixing, add_smear, add_mos
from mypycp2k.outer_scf import add_outer_scf
from mypycp2k.subsys import add_elements, set_unperiodic_cell, set_topology, center_coordinates
"""
Set cp2k input file GW@DFT for the methane
set_XXX means imperatively set some section
add_XXX means optionally set some section. Optionally == if activate_XXX is set to True
activate_XXX means activate "add_XXX"
"""
def main():
    # My input #

    # rel_cutoff: 40; cutoff: 300; abc = 10
    my_abc = '10.0 10.0 10.0'
    cutoff = 300
    rel_cutoff = 40

    ## base settings ##
    basis_set_base_path = '/home/artem/soft/cp2k/cp2k-7.1/data/'
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
    # my_elements = ['H', 'O']  # for test
    # my_element = ['H']  # for test
    my_elements = ['H', 'C']
    inp_file_name = 'GW_PBE_for_methane.inp'
    activate_vdw = False
    activate_outer_scf = False
    activate_ot = False
    activate_diagonalization = not activate_outer_scf
    wf_corr_num_proc = 4  # 16 in the ref paper; -1 to use all

    ####################################################################################################################
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
    OUTER_SCF = DFT.SCF.OUTER_SCF
    ####################################################################################################################

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
    set_cutoff(DFT, cutoff=cutoff, rel_cutoff=rel_cutoff, ngrids=5)
    set_scf(DFT, eps_scf=1.0E-10)
    if activate_ot:
        add_ot(SCF)
    else:
        add_diagonalization(SCF)
        add_smear(SCF)
        add_mixing(SCF)
        add_mos(SCF)
    #
    add_outer_scf(OUTER_SCF)
    set_pbe(XC)  # alter: set_pb0, etc.
    set_qs(DFT,
           eps_default=1.0E-15,
           eps_pgf_orb=1.0E-200)
    add_gw_ver_0(XC)  # GW!
    print_mo_cubes(DFT.PRINT)
    print_mo(DFT.PRINT)
    if activate_vdw:
        add_vdw(XC, vdw_parameters_file=my_vdw_parameters_file)
    ## END DFT ##

########################################################################################################################
    calc.write_input_file(inp_file_name)
    # calc.run()
    print("I am done")
########################################################################################################################


if __name__ == '__main__':
    main()