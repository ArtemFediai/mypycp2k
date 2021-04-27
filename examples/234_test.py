"""
!!! Edit, test, push only at int-nano !!!
2,3,4-cc-pV
or
aug-cc-pV
"""
from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, set_qs
from mypycp2k.xc import set_pbe, set_pbe0, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_diagonalization, add_mixing, add_smear, add_mos, remove_ot
from mypycp2k.subsys import add_elements, set_unperiodic_cell, set_topology, center_coordinates
from extract_functions.extract_from_output import return_homo_lumo, return_gw_energies, \
    extract_number_of_independent_orbital_function
import os
from util.units_conversion import eV_to_Hartree
from cp2k_run.cp2k_run import cp2k_run
import numpy as np
# import rdkit
from copy import deepcopy
import argparse
from util.xyz import XYZ
import yaml
from stand_alone_scripts.db_create.test_yaml import Cp2kOutput
from shutil import rmtree, copytree


def main():

    sim = 'sim'  # simulation folder name
    database_file = 'db'
    bh5670 = 'bh5670'  # the outermost folder in the scratch where all other data are put
    my_offset = 15  # from the molecule to a vacuum. 15-20 is recommended!

    dummy_run = False  # does not invoke cp2k if true

    #  parser begin
    parser = argparse.ArgumentParser(description='rank and num of cpus')
    parser.add_argument('-rank')  # array job number
    parser.add_argument('-num_cpus')  # number of cpus you request for every array job
    args = parser.parse_args()
    #  parser end

    #  parsing input
    threads = args.num_cpus  # cpus used to compute
    rank = '{:0>6}'.format(args.rank)  # transform rank from '1' to '000001' format
    xyz_file_name = f'dsgdb9nsd_{rank}.xyz'
    xyz_file_location = f'dsgdb9nsd/{xyz_file_name}'
    sim_folder_scratch = f'/scratch/{bh5670}/{sim}/{rank}'
    if not os.path.exists(sim_folder_scratch):
        os.mkdir(sim_folder_scratch)
    else:
        rmtree(sim_folder_scratch)  # leftovers from previous simulations will be removed
        os.mkdir(sim_folder_scratch)  # and the new folder will be created


    # xyz object created, normal xyz file is created at scratch
    my_xyz_file_obj = XYZ.from_file(xyz_file_location)  # object created using the file from home
    xyz_at_scratch = sim_folder_scratch + '/' + xyz_file_name  #
    my_xyz_file_obj.write(xyz_at_scratch)  # writes a normal xyz (into scratch)


    # my molecule object is created. It will serve as a DB record
    my_new_mol = Cp2kOutput(rank)

    # rel_cutoff: 40; cutoff: 300; abc = 10
    my_abc = str(my_xyz_file_obj.compute_box_size(offset=my_offset))[1:-2]
    cutoff = 300
    rel_cutoff = 40

    # GLOBAL SETTINGS
    #

    ## base settings ##
    #my_basis_set_file_name = basis_set_base_path + 'BASIS_def2_QZVP_RI_ALL'
    basis_set_file_name = 'BASIS_CC_AUG_RI'  # RI5, 2-5 cc, all aug-cc
    my_vdw_parameters_file = 'dftd3.dat'
    # my_basis_sets = ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z']
    my_basis_sets = ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z']
    # my_potential_file_name = basis_set_base_path + 'POTENTIAL'
    my_potential_file_name = 'POTENTIAL'
    my_potential = 'ALL'
    my_project_name = "this_is_template"
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
    set_topology(SUBSYS, xyz_file_name=xyz_file_name)
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
    set_scf(DFT, eps_scf=1.0E-10, max_scf=500, scf_guess='ATOMIC')
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
    suffix = ['2', '3', '4']  # cardinal numbers of the database

    # begin: input
    # cp2k_exe_path = '/home/artem/soft/cp2k/cp2k-7.1/exe/local/cp2k.popt'
    cp2k_exe_path = '/home/ws/bh5670/cp2k/cp2k-7.1/exe/local/cp2k.popt' 
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
        # end: input

        # OT run to converge quickly
        calc_.write_input_file(sim_folder_scratch + '/' + ot_file_name)
        # first run
        print(f"Running PBE with OT (basis set = {suffix})...")
        if not dummy_run:
            cp2k_run(input_file=ot_file_name,
                     xyz_file=xyz_file_name,
                     run_type=my_run_type,
                     np=threads,
                     output_file=f'out_ot_{suffix}.out',
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=sim_folder_scratch)
            # end: first run
        print(f"I have finished cp2k with OT (basis set = {suffix})")

        # DIAGONALIZATION RUN to reliably compute HOMO and then GW
        # remove the OT method
        remove_ot(SCF_)

        # change calculations to a diagonalization
        add_diagonalization(SCF_)
        add_smear(SCF_)  # add or not?
        add_mixing(SCF_)  # add or not?
        add_mos(SCF_)  # add or not?
        # plot homo/lumo

        set_pbe0(XC_)  # we want G0W0@PBE0. no pbe0 in the beginning
        print_mo_cubes(DFT_.PRINT, nhomo=10, nlumo=10)  # all HOMOs are typicall plotted
        set_scf(DFT_, eps_scf=1E-6, max_scf=200)
        # add G0W0!
        add_gw_ver_0(XC_,
                     ev_sc_iter=1,
                     wf_corr_num_proc=wf_corr_num_proc,
                     rpa_num_quad_points=100,
                     )  # GW!

        calc_.write_input_file(sim_folder_scratch + '/' + diag_file_name)
        # second run
        print(f"Running G0W0 with DIAG (basis set = {suffix})...")
        my_out_file2 = f'out_diag_{suffix}.out'
        if not dummy_run:
            cp2k_run(input_file=diag_file_name,
                     xyz_file=xyz_file_name,
                     output_file=my_out_file2,
                     run_type=my_run_type,
                     np=threads,
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=sim_folder_scratch)
            print(f"I have finished cp2k with DIAG (basis set = {suffix})")

            # extract homo/lumo and gw homo/lumo from the cp2k output file:
            path_to_out2_file = sim_folder_scratch + '/' + my_out_file2

            # extract from the output
            try:
                num_orb = extract_number_of_independent_orbital_function(path_to_out2_file)
                print(f'basis set = {suffix}, number of independent orbital functions: {num_orb}')
            except:
                print('number of orbatals was not extracted')
                num_orb = 'not extracted'

            try:
                homos, lumos = [], []
                homos, lumos = return_homo_lumo(path_to_out2_file)
                print(f'basis set = {suffix} ', 'homo = ', homos[-1]*eV_to_Hartree(), ' eV')
                print(f'basis set = {suffix} ', 'lumo = ', lumos[0]*eV_to_Hartree(), ' eV')
                homo = homos[-1]*eV_to_Hartree()
                lumo = lumos[0]*eV_to_Hartree()
            except:
                print(f'Homo/Lumo were not extracted')
                homo = 'not extracted'
                lumo = 'not extracted'

            try:
                gw_occ, gw_vir, homo_, lumo_ = return_gw_energies(path_to_out2_file)
                if isinstance(homo, str) or isinstance(lumo, str):
                    homo = homo_
                    lumo = lumo_
                print(f'basis set = {suffix} ', 'homo = ', homo, ' eV')
                print(f'basis set = {suffix} ', 'lumo = ', lumo, ' eV')
                print(f'basis set = {suffix} ', 'gw homo = ', gw_occ, ' eV')
                print(f'basis set = {suffix} ', 'gw lumo = ', gw_vir, ' eV')
            except:
                print("GW energies were not extracted")
                # gw_occ = 'not extracted'
                # gw_vir = 'not extracted'

            del calc_

            #  put computed data into the molecule object
            my_new_mol.add_energies(int(suffix), homo, lumo, gw_occ, gw_vir)
            my_new_mol.add_num_orbitals(int(suffix), num_orb)
            my_new_mol.extrapolate_energy()
            db_record = my_new_mol.yield_dict()  # this dict will be written into yaml. it will be a record in the global library
            #

    print("\nI am done\n")
    print('saving to DB...')

    with open(f'{database_file}/DB_{rank}.yaml', 'w') as stream:
        yaml.safe_dump(db_record, stream)

    print(f"saved to {database_file}/DB_{rank}.yaml")
    print('I will remove the content the sim folder')

    # Clean up before leave
    status = my_new_mol.status()
    if status == 'all_extracted':  # all quantities are extracted
        print(f'status: {status} ==> will remove {sim_folder_scratch}')
        try_to_remove_folder(sim_folder_scratch)
    else:
        print(f'status: {status} ==> will copy failed sim folder from scratch')
        sim_folder_home = f'{sim}/{rank}'  # sim folder at home exists. you create {rank} folder
        #if not os.path.exists(sim_folder_home):
        #os.mkdir(sim_folder_home)   # will overwrite if exists
        copytree(sim_folder_scratch, sim_folder_home)  # will rewrite the folder
        print(f"I have copied {sim_folder_scratch} to {sim_folder_home}")
        try_to_remove_folder(sim_folder_scratch)

######################################## END: RUN CP2K TWO TIMES #######################################################

def try_to_remove_folder(folder):
    try:
        rmtree(folder)
        print(f"I have removed {folder}")
    except OSError as e:
        print("Error: %s : %s" % (folder, e.strerror))
        print(f"Folder {folder} could not be deleted")


if __name__ == '__main__':
    main()
