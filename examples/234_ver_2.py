"""
VERSION 2: new element is added to parser: the list of the molecules to be simulated (numbers)

The sh file simply set the number of array jobs. That is it!

here, the output for the TZ is the initial guess for DZ wfn, etc.

!!! Edit, test, push only at int-nano or other cluster!!!

2,3,4-cc-pV
or
aug-cc-pV

"""
from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, set_qs
from mypycp2k.xc import set_pbe, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_diagonalization, add_mixing, add_mos, remove_ot
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
from extract_functions.basis_set_extrapolation import Cp2kOutput
from shutil import rmtree, copy, copytree
from util.general import parse_mixed_list
import csv

def main():

    #  parser begin
    parser = argparse.ArgumentParser(description='rank and num of cpus')
    parser.add_argument('-rank')  # array job number
    parser.add_argument('-num_cpus')  # number of cpus you request for every array job
    parser.add_argument('-i')  # input_from_yaml yaml file
    parser.add_argument('-mol_ids')  # mol_ids to simulate (without prefix and suffix)
    args = parser.parse_args()
    #  parser end

    #  yaml file
    yaml_file_name = args.i
    with open(yaml_file_name) as stream:
        input = yaml.load(stream=stream)
    #  end: yaml file

    #  todo: think over because it is imported twice
    #  end: run-or-check settings
    debug = input['debug']
    dummy_run = input['dummy_run']

    #  if not at cluster: test
    # debug = True
    # dummy_run = True
    # end: if not at cluster

    #  folders names
    sim = input['folder_names']['simulations']
    db = input['folder_names']['database']
    bh5670 = input['folder_names'][
        'scratch']  # the outermost folder in the scratch folder where all other data are put
    prefix_xyz_file_name = input['prefix_xyz_file_name']
    my_offset = input['molecule_vacuum_offset']
    try:
        type_mpi = input['mpi']
    except:
        type_mpi = 'openmpi'


    #  parsing input_from_yaml
    threads = int(args.num_cpus)  # cpus used to compute. I do not subtract 1. This does not help
    # mol_id = parse_mixed_list()
    with open('db/missing_num.csv', 'r') as stream:
        csv_reader = csv.reader(stream)
        all_numbers = csv_reader.__next__()  # only one line in this csv format file, so we do not loop over
        rank = all_numbers[int(args.rank)-1]
    # rank = '{:0>6}'.format(mol_id)  # transform rank from '1' to '000001' format. This is not a general thing
    xyz_file_name = f'{prefix_xyz_file_name}_{rank}.xyz'
    xyz_file_location = f'{prefix_xyz_file_name}/{xyz_file_name}'

    db_record_path = f'{db}/DB_{rank}.yaml'  # file where the results will be saved todo: raeum es alles auf!
    #  check is the output exists
    if os.path.exists(db_record_path):
        print(f'The simulation results of mol. {rank} is already in the folder of reference')
        exit()
    #  end: check if the output exists

    if True: #not dummy_run:
        sim_folder_scratch = f'/scratch/{bh5670}/{sim}/{rank}'
    else:
        sim_folder_scratch = f'scratch/{bh5670}/{sim}/{rank}'

    sim_folder_home = f'{sim}/{rank}'  # sim folder at home exists. you create later {rank} folder
    if not os.path.exists(sim_folder_scratch):
        os.mkdir(sim_folder_scratch)
    else:
        rmtree(sim_folder_scratch)  # leftovers from previous simulations will be removed
        os.mkdir(sim_folder_scratch)  # and the new folder will be created


    #  xyz object created, normal xyz file is created at scratch
    try:
        my_xyz_file_obj = XYZ.from_file(xyz_file_location)  # object created using the file from home
    except:  #  test
        my_xyz_file_obj = XYZ.from_file('H2O.xyz')  # object created using the file from home

    xyz_at_scratch = sim_folder_scratch + '/' + xyz_file_name  #
    my_xyz_file_obj.write(xyz_at_scratch)  # writes a normal xyz (into scratch)

    # my molecule object is created. It will serve as a DB record
    my_new_mol = Cp2kOutput(rank)

    # rel_cutoff: 40; cutoff: 300; abc = 10
    my_abc = str(my_xyz_file_obj.compute_box_size(offset=my_offset))[1:-2]
    input['my_abc'] = my_abc
    input['xyz_file_name'] = xyz_file_name

    # misc
    wf_corr_num_proc = 0  # 16 in the ref paper; -1 to use all
    inp_file_name = 'test_2344.inp'  # base file name

######################################## BEGIN: RUN CP2K TWO TIMES #####################################################
    # suffix = ['2', '3', '4']  # cardinal numbers of the database

    # begin: input_from_yaml
    # cp2k_exe_path = '/home/artem/soft/cp2k/cp2k-7.1/exe/local/cp2k.popt'
    #cp2k_exe_path = '/home/ws/bh5670/cp2k/cp2k-7.1/exe/local/cp2k.popt'
    cp2k_exe_path = input['cp2k_exe_path']
    my_run_type = 'mpi'

    suffix = input['basis_set_suffix']  # todo: fix DZ --> 2, TZ --> 3, QZ --> 4

    InputFactory.set_constants(input_from_yaml=input)

    for i_bs, suffix in enumerate(suffix):

        # start: ot dft
        # I/O
        output_file = f'out_{suffix}.out'
        ot_file_name = 'OT_' + f'{suffix}_' + inp_file_name  # for DFT (OT)
        diag_file_name = 'DIAG_' + f'{suffix}_' + inp_file_name  # for GW (DIAG)
        # end: I/O

        # OT dft simulation to converge quickly: create the simulation object
        dft_ot_simulation = InputFactory.new_dft_ot(i_bs)
        # OT dft: write input file
        dft_ot_simulation.write_input_file(sim_folder_scratch + '/' + ot_file_name)
        # OT dft run below ...

        # ... but before, we copy the RESTART from the previous basis set (it exists unless for the smallest basis set)
        if i_bs != 0:
            try:
                copy(sim_folder_scratch + '/' + f'{int(suffix) - 1}-RESTART.wfn',
                     sim_folder_scratch + '/' + f'{suffix}-RESTART.wfn')
                print('copied restart file 2->3 or 3->4')
            except:
                print('not succesfull copy of the restart file')
        elif i_bs == 0:
            pass
        #

        print(f"Running PBE with OT (basis set = {suffix})...")
        if not dummy_run:
            cp2k_run(input_file=ot_file_name,
                     xyz_file=xyz_file_name,
                     run_type=my_run_type,
                     np=threads,
                     output_file=f'out_ot_{suffix}.out',
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=sim_folder_scratch,
                     type_mpi=type_mpi)
            # end: first run
        print(f"I have finished cp2k with OT (basis set = {suffix})")

        # DIAGONALIZATION RUN to reliably compute HOMO and then GW

        # gw: create the simulation object
        gw_diag_simulations = InputFactory.new_gw(i_bs)
        # gw: write the input file
        gw_diag_simulations.write_input_file(sim_folder_scratch + '/' + diag_file_name)
        # gw run
        print(f"Running G0W0 with DIAG (basis set = {suffix})...")
        my_out_file2 = f'out_diag_{suffix}.out'
        if not dummy_run:
            cp2k_run(input_file=diag_file_name,
                     xyz_file=xyz_file_name,
                     output_file=my_out_file2,
                     run_type=my_run_type,
                     np=threads,
                     cp2k_executable=cp2k_exe_path,
                     execution_directory=sim_folder_scratch,
                     type_mpi=type_mpi)
            print(f"I have finished cp2k with DIAG (basis set = {suffix})")

            # extract homo/lumo and gw homo/lumo from the cp2k output file:
            path_to_out2_file = sim_folder_scratch + '/' + my_out_file2

            # the method to extract?
            # extract from the output
            try:
                num_orb = extract_number_of_independent_orbital_function(path_to_out2_file)
                print(f'basis set = {suffix}, number of independent orbital functions: {num_orb}')
            except:
                print('number of orbitals was not extracted')
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
                if isinstance(homo, str) and isinstance(lumo, str):
                    homo = homo_
                    lumo = lumo_
                print(f'basis set = {suffix} ', 'homo = ', homo, ' eV')
                print(f'basis set = {suffix} ', 'lumo = ', lumo, ' eV')
                print(f'basis set = {suffix} ', 'gw homo = ', gw_occ, ' eV')
                print(f'basis set = {suffix} ', 'gw lumo = ', gw_vir, ' eV')
            except:
                print("GW energies were not extracted")
                gw_occ = 'not extracted'
                gw_vir = 'not extracted'

            del dft_ot_simulation, gw_diag_simulations

            #  put computed data into the molecule object
            my_new_mol.add_energies(int(suffix), homo, lumo, gw_occ, gw_vir)
            my_new_mol.add_num_orbitals(int(suffix), num_orb)
            my_new_mol.extrapolate_energy()  # level up?
            db_record = my_new_mol.yield_dict()  # this dict will be written into yaml. it will be a record in the global library
            #

    ######################################## END: RUN CP2K TWO TIMES #######################################################
    print("\nI am done\n")

    if not dummy_run:
        print('saving to DB...')

        with open(f'{db}/DB_{rank}.yaml', 'w') as stream:
            yaml.safe_dump(db_record, stream)

        print(f"saved to {db}/DB_{rank}.yaml")

    print('I will remove the content of the sim folder')
    # Clean up before leave
    status = my_new_mol.status()
    if status == 'all_extracted':  # all quantities are extracted
        if debug:
            print(f'status: {status}, but debug is on ==> will move {sim_folder_scratch} to {sim_folder_home}')
            copytree(sim_folder_scratch, sim_folder_home, dirs_exist_ok=True)  # will rewrite the folder
        else:
            print(f'status: {status} ==> will remove {sim_folder_scratch}')
            try_to_remove_folder(sim_folder_scratch)
    else:
        print(f'status: {status} ==> will copy failed sim folder from scratch')
        #if not os.path.exists(sim_folder_home):
        #os.mkdir(sim_folder_home)   # will overwrite if exists
        try:
            copytree(sim_folder_scratch, sim_folder_home)  # will rewrite the folder?
            print(f"I have copied {sim_folder_scratch} to {sim_folder_home}")
        except:
            print(f"I could not copy {sim_folder_scratch} to {sim_folder_home}")
        try_to_remove_folder(sim_folder_scratch)


def try_to_remove_folder(folder):
    try:
        rmtree(folder)
        print(f"I have removed {folder}")
    except OSError as e:
        print("Error: %s : %s" % (folder, e.strerror))
        print(f"Folder {folder} could not be deleted")


class InputFactory:
    @classmethod
    def set_constants(cls, input_from_yaml):
        #  folders names
        cls.sim = input_from_yaml['folder_names']['simulations']
        cls.db = input_from_yaml['folder_names']['database']
        cls.bh5670 = input_from_yaml['folder_names'][
            'scratch']  # the outermost folder in the scratch folder where all other data are put
        cls.prefix_xyz_file_name = input_from_yaml['prefix_xyz_file_name']
        # other settings
        cls.offset = input_from_yaml['molecule_vacuum_offset']
        cls.cutoff = input_from_yaml['cutoff']
        cls.rel_cutoff = input_from_yaml['rel_cutoff']
        cls.basis_set_file_name = input_from_yaml['basis_set_file']
        cls.basis_sets = input_from_yaml['basis_set_list']
        cls.ri_basis_sets = input_from_yaml['ri_basis_set_list']
        cls.debug = input_from_yaml['debug']
        cls.dummy_run = input_from_yaml['dummy_run']
        cls.cp2k_exe_path = input_from_yaml['cp2k_exe_path']
        cls.suffix = input_from_yaml['basis_set_suffix']  # DZ --> 2, TZ --> 3, QZ --> 4
        cls.potential_file_name = input_from_yaml['potential_file_name']
        cls.potential = input_from_yaml['potential']  # pseudo. 'ALL' for all-electron
        cls.elements = input_from_yaml['elements']
        cls.my_abc = input_from_yaml['my_abc']  # added afterwards
        cls.xyz_file_name = input_from_yaml['xyz_file_name']  # added afterwards

    ################################## CREATE CONSTANTS FOR ALL TEMPLATES ###############################################

    ########################################### CREATE GW TEMPLATE ######################################################
    @classmethod
    def dft_prototype(cls):
        #  creates CP2K object from scratch. Note that gw prototype takes dft_prototype as a starting point
        calc_dft = CP2K()  # cp2k object
        calc_dft.working_directory = './'

        # pycp2k objects: hierarchy
        CP2K_INPUT = calc_dft.CP2K_INPUT  # CP2K_INPUT is what we need
        FORCE_EVAL = CP2K_INPUT.FORCE_EVAL_add()
        FORCE_EVAL.Method = 'QUICKSTEP'
        SUBSYS = FORCE_EVAL.SUBSYS
        DFT = FORCE_EVAL.DFT
        XC = DFT.XC
        SCF = DFT.SCF
        ####################################################################################################################
        # GLOBAL #
        # FORCE EVAL #
        set_global(CP2K_INPUT)
        ## SUBSYS ##
        set_unperiodic_cell(SUBSYS, abc=cls.my_abc)
        set_nonperiodic_poisson(DFT)
        set_topology(SUBSYS, xyz_file_name=cls.xyz_file_name)
        center_coordinates(SUBSYS)
        ## END SUBSYS ##

        ## DFT ##
        set_dft(DFT,
                potential_file_name=cls.potential_file_name,
                basis_set_file_name=cls.basis_set_file_name)
        set_cutoff(DFT, cutoff=cls.cutoff, rel_cutoff=cls.rel_cutoff, ngrids=5)
        set_scf(DFT, eps_scf=1.0E-8, max_scf=500, scf_guess='RESTART')
        # <-- even if no actual wfn saved, will not collapse, but start with ATOMIC guess
        add_ot(SCF, stepsize=0.05)
        #
        # add_outer_scf(OUTER_SCF)
        set_pbe(XC)  # we start with pbe
        # set_pbe0(XC) no pbe0 in the beginning
        set_qs(DFT,
               eps_default=1.0E-10,
               eps_pgf_orb=np.sqrt(1.0E-10))
        ## END DFT ##
        return calc_dft

    @classmethod
    def gw_prototype(cls):

        calc_gw = deepcopy(InputFactory.dft_prototype())  # we could also make it from scratch,
        # but it is better to inherit it from the dft to decrease the probability of the error

        # pycp2k objects: hierarchy
        CP2K_INPUT = calc_gw.CP2K_INPUT
        FORCE_EVAL = CP2K_INPUT.FORCE_EVAL_list[0]
        FORCE_EVAL.Method = 'QUICKSTEP'
        SUBSYS = FORCE_EVAL.SUBSYS
        DFT = FORCE_EVAL.DFT
        XC = DFT.XC
        SCF = DFT.SCF
        ################################################################################################################

        # remove the OT method
        remove_ot(SCF)
        # change calculations to a diagonalization
        add_diagonalization(SCF)
        # add_smear(SCF_)  # uses final T.
        add_mixing(SCF)  # add or not?
        add_mos(SCF, added_mos=1000)
        # plot homo/lumo

        #set_pbe0(XC_)  # we want G0W0@PBE0. no pbe0 in the beginning
        print_mo_cubes(DFT.PRINT, nhomo=10, nlumo=10)  # all HOMOs are typicall plotted
        set_scf(DFT, eps_scf=1E-6, max_scf=200)
        # add G0W0!
        add_gw_ver_0(XC,
                     ev_sc_iter=1,
                     rpa_num_quad_points=100,
                     max_memory_wf=4000,
                     max_memory_hf=500,
                     corr_occ=1,
                     corr_virt=1)  # GW!
        # it is important to keep WF memory smaller than HF memory, otherwise, it crashes
        return calc_gw

    @classmethod
    def __new_simulation(cls,
                         prototype,
                         i_bs):

        calc_ = deepcopy(prototype)

        #
        CP2K_INPUT_ = calc_.CP2K_INPUT
        FORCE_EVAL_ = CP2K_INPUT_.FORCE_EVAL_list[0]
        SUBSYS_ = FORCE_EVAL_.SUBSYS
        DFT_ = FORCE_EVAL_.DFT
        #XC_ = DFT_.XC
        #SCF_ = DFT_.SCF
        #OUTER_SCF_ = DFT_.SCF.OUTER_SCF
        #

        suffix = cls.suffix[i_bs]
        set_global(CP2K_INPUT_, project_name=suffix)
        add_elements(SUBSYS_,
                     elements=cls.elements,
                     basis=cls.basis_sets[i_bs],
                     aux_basis=cls.ri_basis_sets[i_bs],
                     pot=cls.potential)
        return calc_

    @classmethod
    def new_dft_ot(cls, i_bs):
        return InputFactory.__new_simulation(
            prototype=cls.dft_prototype(),
            i_bs=i_bs)

    @classmethod
    def new_gw(cls, i_bs):
        return InputFactory.__new_simulation(
            prototype=cls.gw_prototype(),
            i_bs=i_bs)


######################################## END: CREATE TEMPLATE ##########################################################


if __name__ == '__main__':
    main()
