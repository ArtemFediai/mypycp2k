from pycp2k import CP2K
from mypycp2k.cp2k_input import set_global
from mypycp2k.dft import set_dft, set_scf, set_nonperiodic_poisson, set_cutoff, print_mo_cubes, print_mo, set_qs
from mypycp2k.xc import set_pbe, set_pbe0, add_vdw, add_gw_ver_0
from mypycp2k.scf import add_ot, add_ot_never_fail
from mypycp2k.subsys import add_elements, set_unperiodic_cell, set_topology, center_coordinates
import re
"""
Set cp2k input file GW@DFT for the methane
set_XXX means imperatively set some section
add_XXX means optionally set some section. Optionally == if activate_XXX is set to True
activate_XXX means activate "add_XXX"
"""
def main():
    # My input #
    ## base settings ##
    basis_set_base_path = '/home/artem/soft/cp2k/cp2k-7.1/data/'

    ## specific settings ##
    my_abc = '12.0 12.0 12.0'
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
    inp_file_name = 'main_screen_cutoff.inp'
    activate_vdw = False
    activate_outer_scf = False
    activate_ot = True
    activate_diagonalization = not activate_outer_scf
    wf_corr_num_proc = 1  # 16 in the ref paper; -1 to use all

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
    set_cutoff(DFT, cutoff=900, rel_cutoff=60, ngrids=5)
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

########################################################################################################################
    # # my_cutoffs = [100, 300, 600, 900]
    # my_cutoffs = [100, 200]
    # for my_cutoff in my_cutoffs:
    #     set_cutoff(DFT, cutoff=my_cutoff)
    #     # calc.write_input_file(inp_file_name)
    #     calc.run()
    # print("I am done")

    ######################################## check cutoff ##############################################################
    if True:
        my_cutoffs = [50, 100, 200, 300, 400, 500, 600, 700, 800]
        # my_cutoffs = [50, 100]
        my_rel_cutoff = 400
        my_target_accuracy_eV = 1.E-6
        check_cutoff_convergence(DFT, calc,
                                 cutoffs=my_cutoffs,
                                 rel_cutoff=my_rel_cutoff,
                                 target_accuracy_eV=my_target_accuracy_eV)
    ####################################################################################################################

    # check abc
    if False:
        my_abcs = [5, 7, 9, 11, 13, 15, 17]
        my_target_accuracy_eV = 1.E-6
        check_abc_convergence(SUBSYS, DFT, calc,
                                 abcs=my_abcs,
                                 target_accuracy_eV=my_target_accuracy_eV)


def check_abc_convergence(SUBSYS, DFT, calc,
                          abcs,
                          target_accuracy_eV):
    print("I will check cutoff convergence")
    import os
    energies = []

    all_out_files = list_files1('./', 'out')
    print(f'found {len(all_out_files)} out files. Will delete them')
    for file in all_out_files:
        os.remove(file)

    for i, abc in enumerate(abcs):
        set_unperiodic_cell(SUBSYS, (str(abc)+' ')*3)
        calc.output_path = calc.working_directory + "/" + calc.project_name + str(abc) + ".out"
        calc.write_input_file("test.inp")
        calc.run()
        try:
            os.remove('methane_GW_PBE-RESTART.wfn')
        except:
            pass
        with open(calc.output_path, "r") as fin:
            regex = re.compile(" ENERGY\| Total FORCE_EVAL \( QS \) energy \(a\.u\.\):\s+(.+)\n")
            for line in fin:
                match = regex.match(line)
                if match:
                    energies.append(match.groups()[0])

    print(energies)
    plot_energies(data=energies, cutoffs=abcs, name='ene.png', target_accuracy_ev=target_accuracy_eV)

def check_cutoff_convergence(DFT, calc, cutoffs, rel_cutoff, target_accuracy_eV):
    print("I will check cutoff convergence")
    import os
    energies = []

    set_scf(DFT=DFT, scf_guess='ATOMIC')
    all_out_files = list_files1('./', 'out')
    print(f'found {len(all_out_files)} out files. Will delete them')
    for file in all_out_files:
        os.remove(file)

    for cutoff in cutoffs:
        DFT.MGRID.Cutoff = cutoff
        DFT.MGRID.Rel_cutoff = rel_cutoff
        calc.output_path = calc.working_directory + "/" + calc.project_name + str(cutoff) + ".out"
        calc.run()
        try:
            os.remove('methane_GW_PBE-RESTART.wfn')
        except:
            pass
        with open(calc.output_path, "r") as fin:
            regex = re.compile(" ENERGY\| Total FORCE_EVAL \( QS \) energy \(a\.u\.\):\s+(.+)\n")
            for line in fin:
                match = regex.match(line)
                if match:
                    energies.append(match.groups()[0])

    print(energies)
    plot_energies(data=energies, cutoffs=cutoffs, name='ene.png', target_accuracy_ev=target_accuracy_eV)


def plot_energies(data, cutoffs, name, target_accuracy_ev):
    """
    :param data: computed energies
    :param cutoffs: parameter that has been screened in Ry!
    :param name: name of the picture that will be generated
    :param target_accuracy_ev: accuracy that will be checked for: in eV!
    :return:
    save the picture
    return the message about the cutoff to be used
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.constants import physical_constants as c
    eV_to_Hartree = c['hartree-electron volt relationship'][0]
    ene = [float(string) for string in data]
    ene = np.array(ene)
    ene *= eV_to_Hartree
    plt.figure(figsize=[3, 4])
    plt.plot(cutoffs[0:-1], np.abs(ene[0:-1]-ene[-1]))
    plt.xlabel("Cutoff")
    plt.ylabel("energy, eV")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(name)

    ene_diffs = np.abs(ene[0:-1] - ene[-1])
    print(f"All energy differences: {ene_diffs*eV_to_Hartree}")
    for i, ene_dif in enumerate(ene_diffs):
        if ene_dif < target_accuracy_ev:
            use_cutoff = cutoffs[i]
            print(f"use cutoff: {use_cutoff}")
            break
        else:
            print("not (yet) reached")


def list_files1(directory, extension):
    from os import listdir
    return [f for f in listdir(directory) if f.endswith('.' + extension)]

if __name__ == '__main__':
    main()