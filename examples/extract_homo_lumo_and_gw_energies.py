"""
Set cp2k input_from_yaml file GW@DFT for the methane
set_XXX means imperatively set some section
add_XXX means optionally set some section. Optionally == if activate_XXX is set to True
activate_XXX means activate "add_XXX"
"""
from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, print_mo, set_qs
from mypycp2k.xc import set_pbe, set_pbe0, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_ot_never_fail, add_diagonalization, add_mixing, add_smear, add_mos, remove_ot
from mypycp2k.outer_scf import add_outer_scf
from mypycp2k.subsys import add_elements, set_unperiodic_cell, set_topology, center_coordinates
from extract_functions.extract_from_output import return_homo_lumo, return_gw_energies
import os
from util.units_conversion import eV_to_Hartree
from cp2k_run.cp2k_run import cp2k_run
# import rdkit

def main():
    # My input_from_yaml #

    # rel_cutoff: 40; cutoff: 300; abc = 10
    my_abc = '15.0 15.0 15.0'
    cutoff = 300
    rel_cutoff = 40

    ## base settings ##
    basis_set_base_path = '/home/artem/soft/cp2k/cp2k-7.1/data/'
    # my_basis_set_file_name = basis_set_base_path + 'BASIS_RI_cc-TZ'G
    my_basis_set_file_name = basis_set_base_path + 'BASIS_def2_QZVP_RI_ALL'
    my_basis_set_file_name = basis_set_base_path + 'BASIS_CC_AUG_RI_NEW'
    my_vdw_parameters_file = basis_set_base_path + 'dftd3.dat'
    my_basis_set = 'def2-QZVP'
    my_basis_set = 'cc-pVQZ'
    my_potential_file_name = basis_set_base_path + 'POTENTIAL'
    my_potential = 'ALL'
    my_project_name = "methane_GW_PBE"
    my_xyz_file_name = 'methane.xyz'
    my_ri_aux_basis_set = 'RI-5Z'  #
    my_ri_aux_basis_set = 'cc-pVQZ-RIFIT'  #
    organic_elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'B', 'I']
    # elements = ['H', 'O']  # for test
    # my_element = ['H']  # for test
    my_elements = ['H', 'C']
    my_elements = organic_elements
    inp_file_name = 'GW_PBE_for_methane.inp'
    activate_vdw = False
    activate_outer_scf = False
    wf_corr_num_proc = 4  # 16 in the ref paper; -1 to use all

    ########################################### CREATE TEMPLATE FOR TWO RUNS ###########################################
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
    add_ot(SCF)
    #
    add_outer_scf(OUTER_SCF)
    set_pbe(XC)  # alter: set_pb0, etc.
    #set_pbe0(XC)
    set_qs(DFT,
           eps_default=1.0E-10,
           eps_pgf_orb=1.0E-5)

    # print_mo(DFT.PRINT)
    if activate_vdw:
        add_vdw(XC, vdw_parameters_file=my_vdw_parameters_file)
    ## END DFT ##

######################################## END: CREATE TEMPLATE ##########################################################


######################################## BEGIN: RUN CP2K TWO TIMES #####################################################
    # begin: input_from_yaml
    cp2k_exe_path = '/home/artem/soft/cp2k_test/cp2k-7.1/exe/local/cp2k.popt'
    # cp2k_exe_path = '/usr/bin/cp2k.popt'
    #cp2k_exe_path = 'cp2k.popt'
    run_folder = 'my_run_folder'
    output_file = 'out.out'
    ot_file_name = 'OT_' + inp_file_name
    diag_file_name = 'DIAG_' + inp_file_name
    my_xyz_file_name = 'methane.xyz'
    threads = 8
    my_run_type = 'mpi'
    # end: input_from_yaml

    if not os.path.exists(run_folder):
        os.mkdir(run_folder)

    # OT run to converge quickly
    calc.write_input_file(run_folder + '/' + ot_file_name)
    # first run
    print("Running cp2k with OT ...")
    cp2k_run(input_file=ot_file_name,
             xyz_file=my_xyz_file_name,
             run_type=my_run_type,
             np=threads,
             output_file='out1.out',
             cp2k_executable=cp2k_exe_path,
             execution_directory=run_folder)
    # end: first run
    print("I have finished cp2k with OT")

    # remove the OT method
    remove_ot(SCF)

    # change calculations to a diagonalization
    add_diagonalization(SCF)
    add_smear(SCF)
    add_mixing(SCF)
    add_mos(SCF)
    # plot homo/lumo
    print_mo_cubes(DFT.PRINT, nhomo=10, nlumo=10)  #  all HOMOs are typicall plotted
    # add G0W0!
    add_gw_ver_0(XC, ev_sc_iter=1, size_freq_integ_group=threads)  # GW!

    # DIAGONALIZATION RUN to reliably compute HOMO
    calc.write_input_file(run_folder + '/' + diag_file_name)
    # second run
    print("Running cp2k with DIAG ...")
    my_out_file2 = 'out2.out'
    cp2k_run(input_file=diag_file_name,
             xyz_file=my_xyz_file_name,
             output_file=my_out_file2,
             run_type=my_run_type,
             np=threads,
             cp2k_executable=cp2k_exe_path,
             execution_directory=run_folder)
    print("I have finished cp2k with DIAG")

    # extract homo/lumo and gw homo/lumo from the cp2k output file:
    path_to_out2_file = run_folder + '/' + my_out_file2

    homos, lumos = [], []
    homos, lumos = return_homo_lumo(path_to_out2_file)

    print('homo = ', homos[-1]*eV_to_Hartree(), ' eV')
    print('lumo = ', lumos[0]*eV_to_Hartree(), ' eV')

    gw_occ, gw_vir, homo, lumo = return_gw_energies(path_to_out2_file)


    print('gw homo = ', gw_occ, ' eV')
    print('gw lumo = ', gw_vir, ' eV')

    print("\nI am done")

######################################## END: RUN CP2K TWO TIMES #######################################################


if __name__ == '__main__':
    main()
