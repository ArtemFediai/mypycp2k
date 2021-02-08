import re
from mypycp2k.dft import set_scf
from mypycp2k.subsys import set_unperiodic_cell
from pycp2k import CP2K

def change_calc_abc(calc, param):
    calc.CP2K_INPUT.FORCE_EVAL_list[0].SUBSYS.CELL.Abc = (str(param)+' ')*3


def change_calc_cutoff(calc, param):
    calc.CP2K_INPUT.FORCE_EVAL_list[0].DFT.MGRID.Cutoff = param


def change_calc_rel_cutoff(calc, param):
    calc.CP2K_INPUT.FORCE_EVAL_list[0].DFT.MGRID.Rel_cutoff = param


def general_convergence(calc,
                        params,
                        func_to_accept_param=change_calc_abc,
                        target_accuracy_eV=1.0E-3,
                        param_name='abc'):

    remove_files_with_extension('./', 'out')
    energies = []

    for param in params:
        func_to_accept_param(calc=calc, param=param)
        calc.output_path = calc.working_directory + "/" + calc.project_name + str(param) + ".out"
        calc.run()
        remove_files_with_extension('./', "wfn")
        remove_files_with_extension('./', "wfn.bak-1")  # ugly
        extract_total_energy(energies, path_to_file=calc.output_path)
    print(f"Extracted energies: {energies} a.u.")
    plot_energies_general(energies=energies,
                          params=params,
                          name=param_name,
                          target_accuracy_ev=target_accuracy_eV)


def plot_energies_general(energies, params, name, target_accuracy_ev):
    """
    :param energies: computed energies
    :param params: parameter that has been screened in Ry!
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
    ene = [float(string) for string in energies]
    ene = np.array(ene)
    ene *= eV_to_Hartree
    plt.figure(figsize=[3, 4])
    plt.plot(params[0:-1], np.abs(ene[0:-1] - ene[-1]))
    plt.xlabel(f"{name}")
    plt.ylabel("E - E$_0$, eV")
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(name+'.png', dpi=600)

    ene_diffs = np.abs(ene[0:-1] - ene[-1])
    print(f"All energy differences: {ene_diffs*eV_to_Hartree}")
    for i, ene_dif in enumerate(ene_diffs):
        if ene_dif < target_accuracy_ev:
            use_cutoff = params[i]
            print(f"use cutoff: {use_cutoff}")
            break
        else:
            print("not (yet) reached")


# def check_abc_convergence(SUBSYS, calc,
#                           abcs,
#                           target_accuracy_eV):
#
#
#     remove_files_with_extension('./', 'out')
#
#     print("I will check abc (cell) convergence")
#     energies = []
#
#
#     for i, abc in enumerate(abcs):
#         set_unperiodic_cell(SUBSYS, (str(abc)+' ')*3)
#         calc.output_path = calc.working_directory + "/" + calc.project_name + str(abc) + ".out"
#         calc.run()
#         remove_files_with_extension("wfn")
#
#     extract_total_energy(energies, path_to_file=calc.output_path)
#     print(energies)
#     plot_energies(data=energies, cutoffs=abcs, name='ene.png', target_accuracy_ev=target_accuracy_eV)
#
#
# def check_cutoff_convergence(DFT, calc, cutoffs, rel_cutoff, target_accuracy_eV):
#     print("I will check cutoff convergence")
#     energies = []
#     set_scf(DFT=DFT, scf_guess='ATOMIC')
#     remove_files_with_extension('./', 'out')
#
#     for cutoff in cutoffs:
#         DFT.MGRID.Cutoff = cutoff
#         DFT.MGRID.Rel_cutoff = rel_cutoff
#         calc.output_path = calc.working_directory + "/" + calc.project_name + str(cutoff) + ".out"
#         calc.run()
#         remove_files_with_extension('./', 'wfn')
#         extract_total_energy(energies, path_to_file=calc.output_path)
#
#     print(energies)
#     plot_energies(data=energies, cutoffs=cutoffs, name='energy_vs_cutoff.png', target_accuracy_ev=target_accuracy_eV)
#
#
# def plot_energies(data, cutoffs, name, target_accuracy_ev):
#     """
#     :param data: computed energies
#     :param cutoffs: parameter that has been screened in Ry!
#     :param name: name of the picture that will be generated
#     :param target_accuracy_ev: accuracy that will be checked for: in eV!
#     :return:
#     save the picture
#     return the message about the cutoff to be used
#     """
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from scipy.constants import physical_constants as c
#     eV_to_Hartree = c['hartree-electron volt relationship'][0]
#     ene = [float(string) for string in data]
#     ene = np.array(ene)
#     ene *= eV_to_Hartree
#     plt.figure(figsize=[3, 4])
#     plt.plot(cutoffs[0:-1], np.abs(ene[0:-1]-ene[-1]))
#     plt.xlabel("Cutoff, Ry")
#     plt.ylabel("Energy, eV")
#     plt.yscale('log')
#     plt.tight_layout()
#     plt.savefig(name)
#
#     ene_diffs = np.abs(ene[0:-1] - ene[-1])
#     print(f"All energy differences: {ene_diffs*eV_to_Hartree}")
#     for i, ene_dif in enumerate(ene_diffs):
#         if ene_dif < target_accuracy_ev:
#             use_cutoff = cutoffs[i]
#             print(f"use cutoff: {use_cutoff}")
#             break
#         else:
#             print("not (yet) reached")


def extract_total_energy(energies, path_to_file):
    """
    extract total energy from cp2k file
    :param energies: list, possible empty to be finned wth energies
    :param path_to_file: path to output cp2k file
    :return: append energies list with an extracted energy
    """
    with open(path_to_file, "r") as fin:
        regex = re.compile(" ENERGY\| Total FORCE_EVAL \( QS \) energy \(a\.u\.\):\s+(.+)\n")
        for line in fin:
            match = regex.match(line)
            if match:
                energies.append(match.groups()[0])


def list_files1(directory, extension):
    from os import listdir
    return [f for f in listdir(directory) if f.endswith('.' + extension)]


def remove_files_with_extension(directory, extension):
    from os import remove
    files = list_files2(directory, extension)
    for file in files:
        remove(file)

def list_files2(directory, extension):
    import re
    from os import listdir
    regex = re.compile('.*\.'+extension)
    print([f for f in listdir(directory) if regex.match(f)])
    return [f for f in listdir(directory) if regex.match(f)]