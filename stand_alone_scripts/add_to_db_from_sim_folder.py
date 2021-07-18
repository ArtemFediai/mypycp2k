"""
Use case:
you have folders in sim folder with non-converged results.
you converged it there.
now you need to extract gw/homo/lumo energies and put them to db

Very important script
works only with 234_ver_2

"""
import os
from os import path
from extract_functions.basis_set_extrapolation import Cp2kOutput



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
from util.exceptions import SCQPSolutionNotFound, SCFNotConvergedNotPossibleToRunMP2



def main():

    def my_out_file(suf, ot_or_diag):
        return f'out_{ot_or_diag}_{suf}.out'


    SUFFIX = ['2', '3']


    my_path = os.path.abspath('')

    print(my_path)

    all_folders = os.listdir('sim')

    print(f'all folders of interest = {all_folders}')

    path_to_folders = [os.path.abspath('sim' + '/' + folder) for folder in all_folders]

    print(path_to_folders)

    for folder, path in zip(all_folders, path_to_folders):

        db_record_path = f'db/DB_{folder}.yaml'

        if os.path.exists(db_record_path):
            print(f'The simulation results of mol. {folder} is already in the folder of reference')
            exit()

        print(f'\nTHIS IS MOL_NUM {folder}\n')


        my_new_mol = Cp2kOutput(folder)

        for suf in SUFFIX:
            diag_out_file = f"sim/{folder}/{my_out_file(suf=suf, ot_or_diag='diag')}"

    #

            try:
                num_orb = extract_number_of_independent_orbital_function(diag_out_file)
                print(f'basis set = {suf}, number of independent orbital functions: {num_orb}')
            except:
                print('number of orbitals was not extracted')
                num_orb = 'not extracted'
            # extract energies:
            try:
                homos, lumos = [], []
                homos, lumos = return_homo_lumo(diag_out_file)
                print(f'basis set = {suf} ', 'homo = ', homos[-1] * eV_to_Hartree(), ' eV')
                print(f'basis set = {suf} ', 'lumo = ', lumos[0] * eV_to_Hartree(), ' eV')
                homo = homos[-1] * eV_to_Hartree()
                lumo = lumos[0] * eV_to_Hartree()
            except:
                print(f'Homo/Lumo were not extracted')
                homo = 'not extracted'
                lumo = 'not extracted'

            try:
                gw_occ, gw_vir, homo_, lumo_ = return_gw_energies(diag_out_file)
                homo, lumo = redefine_homo_lumo_if_not_extracted_before(homo_, lumo_, homo, lumo)
                print_extracted_energies(suf, homo, lumo, gw_occ, gw_vir)  # on a screen
            except SCFNotConvergedNotPossibleToRunMP2:
                print("GW is not extracted, because SCFNotConvergedNotPossibleToRunMP2. Calling fallback ...")
                print('NOT IMPLEMENTED')
            finally:
                try:
                    gw_occ, gw_vir, homo_, lumo_ = return_gw_energies(diag_out_file)
                    homo, lumo = redefine_homo_lumo_if_not_extracted_before(homo_, lumo_, homo, lumo)
                    print_extracted_energies(suf, homo, lumo, gw_occ, gw_vir)  # on a screen
                # <---
                except:
                    print("GW energies were not extracted even in the fallback")
                    gw_occ = 'not extracted'
                    gw_vir = 'not extracted'


            #  put computed data into the molecule object
            my_new_mol.add_energies(int(suf), homo, lumo, gw_occ, gw_vir)
            my_new_mol.add_num_orbitals(int(suf), num_orb)
        my_new_mol.extrapolate_energy()  # level up?
        db_record = my_new_mol.yield_dict()  # this dict will be written into yaml. it will be a record in the global library
        # <-- EMD: GW run and extraction


        print('saving to DB...')

        with open(f'db/DB_{folder}.yaml', 'w') as stream:
            yaml.safe_dump(db_record, stream)

        print(f"saved to db/DB_{folder}.yaml")

        #



        print('done')


####################################################
#                                                  #
#                     END                          #
#                                                  #
####################################################


def redefine_homo_lumo_if_not_extracted_before(homo_, lumo_, homo, lumo):
    if isinstance(homo, str) and isinstance(lumo, str):
        return homo_, lumo_  #  todo: actually, not only string has it to be but also "not extracted"
    else:  # todo: actually, one has to check if these are floats...
        return homo, lumo


def print_extracted_energies(suf, homo, lumo, gw_occ, gw_vir):
    print(f'basis set = {suf} ', 'homo = ', homo, ' eV')
    print(f'basis set = {suf} ', 'lumo = ', lumo, ' eV')
    print(f'basis set = {suf} ', 'gw homo = ', gw_occ, ' eV')
    print(f'basis set = {suf} ', 'gw lumo = ', gw_vir, ' eV')

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
        try:
            cls.eps_scf_dft = input_from_yaml['eps_scf_dft']
        except KeyError:
            cls.eps_scf_dft = [1.0E-8, 1.0E-8, 1.0E-8]
        print('here')

    ################################## CREATE CONSTANTS FOR ALL TEMPLATES ###############################################

    ########################################### CREATE GW TEMPLATE ######################################################
    @classmethod
    def dft_prototype(cls, i_bs=0):
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
        set_scf(DFT, eps_scf=cls.eps_scf_dft[i_bs], max_scf=500, scf_guess='RESTART')
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
            prototype=cls.dft_prototype(i_bs=i_bs),
            i_bs=i_bs)

    @classmethod
    def new_gw(cls, i_bs):
        return InputFactory.__new_simulation(
            prototype=cls.gw_prototype(),
            i_bs=i_bs)


def try_to_copy_previous_restart_file(i_bs, sim_folder_scratch, suf):
    if i_bs != 0:
        try:
            copy(sim_folder_scratch + '/' + f'{int(suf) - 1}-RESTART.wfn',
                 sim_folder_scratch + '/' + f'{suf}-RESTART.wfn')
            print('copied restart file 2->3 or 3->4')
        except:
            print('not succesfull copy of the restart file')
    elif i_bs == 0:
        pass
######################################## END: CREATE TEMPLATE ##########################################################


if __name__ == '__main__':
    main()
