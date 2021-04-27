"""
2,3,4,5-cc-pV
or
aug-cc-pV
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
import numpy as np
# import rdkit
from copy import copy
from copy import deepcopy
def main():
    # My input_from_yaml #

    # rel_cutoff: 40; cutoff: 300; abc = 10
    my_abc = '10.0 10.0 10.0'
    cutoff = 300
    rel_cutoff = 40

    # GLOBAL SETTINGS
    threads = 121  # cpus used to compute
    #


    ## base settings ##
    #basis_set_base_path = '/home/artem/soft/cp2k/cp2k-7.1/data/'
    basis_set_base_path = ''
    # my_basis_set_file_name = basis_set_base_path + 'BASIS_RI_cc-TZ'G
    #my_basis_set_file_name = basis_set_base_path + 'BASIS_def2_QZVP_RI_ALL'
    basis_set_file_name = 'BASIS_CC_AUG_RI'  # RI5, 2-5 cc, all aug-cc
    my_vdw_parameters_file = basis_set_base_path + 'dftd3.dat'
    #my_basis_set = 'def2-QZVP'  # not used
    # my_basis_sets = ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z']
    my_basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z']
    # my_potential_file_name = basis_set_base_path + 'POTENTIAL'
    my_potential_file_name = 'POTENTIAL'
    my_potential = 'ALL'
    my_project_name = "this_is_template"
    my_xyz_file_name = 'H2O.xyz'
    my_ri_aux_basis_set = 'RI-5Z'  #
    organic_elements = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'B', 'I']
    my_elements = organic_elements
    inp_file_name = 'test_2345.inp'
    activate_vdw = False
    activate_outer_scf = False
    wf_corr_num_proc = 1  # 16 in the ref paper; -1 to use all

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
    # add_elements(SUBSYS,
    #              elements=elements,
    #              basis=my_basis_set,
    #              aux_basis=my_ri_aux_basis_set,
    #              pot=my_potential)
    center_coordinates(SUBSYS)
    ## END SUBSYS ##

    ## DFT ##
    set_dft(DFT,
            potential_file_name=my_potential_file_name,
            basis_set_file_name=basis_set_file_name)
    set_cutoff(DFT, cutoff=cutoff, rel_cutoff=rel_cutoff, ngrids=5)
    set_scf(DFT, eps_scf=1.0E-8, max_scf=500)
    add_ot(SCF)
    #
    # add_outer_scf(OUTER_SCF)
    set_pbe(XC)  # we start with pbe
    # set_pbe0(XC) no pbe0 in the beginning
    set_qs(DFT,
           eps_default=1.0E-10,
           eps_pgf_orb=np.sqrt(1.0E-10))

    # print_mo(DFT.PRINT)
    if activate_vdw:
        add_vdw(XC, vdw_parameters_file=my_vdw_parameters_file)
    ## END DFT ##

######################################## END: CREATE TEMPLATE ##########################################################


######################################## BEGIN: RUN CP2K TWO TIMES #####################################################
    suffix = ['2', '3', '4', '5']  # of out folde

    # begin: input_from_yaml
    # cp2k_exe_path = '/home/artem/soft/cp2k/cp2k-7.1/exe/local/cp2k.popt'
    cp2k_exe_path = '/home/ws/bh5670/cp2k/cp2k-7.1/exe/local/cp2k.popt' 
    run_folder = '2345_run_folder'
    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
    my_xyz_file_name = 'H2O.xyz'
    my_run_type = 'mpi'

    for i_bs, suffix in enumerate(suffix):
        # bs

        calc_ = deepcopy(calc)

        #
        CP2K_INPUT_ = calc_.CP2K_INPUT
        FORCE_EVAL_ = CP2K_INPUT_.FORCE_EVAL_list[0]
        SUBSYS_ = FORCE_EVAL_.SUBSYS
        DFT_ = FORCE_EVAL_.DFT
        XC_ = DFT_.XC
        SCF_ = DFT_.SCF
        OUTER_SCF_ = DFT_.SCF.OUTER_SCF
        #

        set_global(CP2K_INPUT_, project_name=suffix)
        add_elements(SUBSYS_,
                     elements=my_elements,
                     basis=my_basis_sets[i_bs],
                     aux_basis=my_ri_aux_basis_set,
                     pot=my_potential)
        # bs
        output_file = f'out_{suffix}.out'
        ot_file_name = 'OT_' + f'{suffix}_' + inp_file_name
        diag_file_name = 'DIAG_' + f'{suffix}_' + inp_file_name
        # end: input_from_yaml

        # OT run to converge quickly
        calc_.write_input_file(run_folder + '/' + ot_file_name)
        # first run
        print(f"Running PBE with OT (basis set = {suffix})...")
        if True:
            cp2k_run(input_file=ot_file_name,
                     xyz_file=my_xyz_file_name,
                     run_type=my_run_type,
                     np=threads,
                     output_file=f'out_ot_{suffix}.out',
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=run_folder)
            # end: first run
        print(f"I have finished cp2k with OT (basis set = {suffix})")

        # remove the OT method
        remove_ot(SCF_)

        # change calculations to a diagonalization
        add_diagonalization(SCF_)
        # add_smear(SCF)  # add or not?
        # add_mixing(SCF)  # add or not?
        # add_mos(SCF)  # add or not?
        # plot homo/lumo

        set_pbe0(XC_)  # we want G-W0@PBE0. no pbe0 in the beginning
        print_mo_cubes(DFT_.PRINT, nhomo=10, nlumo=10)  #  all HOMOs are typicall plotted
        set_scf(DFT_, eps_scf=1E-6)
        # add G0W0!
        add_gw_ver_0(XC_,
                     ev_sc_iter=1,
                     wf_corr_num_proc=wf_corr_num_proc,
                     rpa_num_quad_points=100,
                     )  # GW!

        # DIAGONALIZATION RUN to reliably compute HOMO and then GW
        calc_.write_input_file(run_folder + '/' + diag_file_name)
        # second run
        print(f"Running G0W0 with DIAG (basis set = {suffix})...")
        my_out_file2 = f'out_diag_{suffix}.out'
        if True:
            cp2k_run(input_file=diag_file_name,
                     xyz_file=my_xyz_file_name,
                     output_file=my_out_file2,
                     run_type=my_run_type,
                     np=threads,
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=run_folder)
            print(f"I have finished cp2k with DIAG (basis set = {suffix})")

            # extract homo/lumo and gw homo/lumo from the cp2k output file:
            path_to_out2_file = run_folder + '/' + my_out_file2

            try:

                homos, lumos = [], []
                homos, lumos = return_homo_lumo(path_to_out2_file)

                print('basis set = {suffix} ', 'homo = ', homos[-1]*eV_to_Hartree(), ' eV')
                print('basis set = {suffix} ', 'lumo = ', lumos[0]*eV_to_Hartree(), ' eV')

                gw_occ, gw_vir = return_gw_energies(path_to_out2_file)


                print('basis set = {suffix} ', 'gw homo = ', gw_occ, ' eV')
                print('basis set = {suffix} ', 'gw lumo = ', gw_vir, ' eV')
            except:
                print("HOMO/LUMO and GW were not extracted")

        print("\nI am done")
        del calc_


######################################## END: RUN CP2K TWO TIMES #######################################################


if __name__ == '__main__':
    main()
