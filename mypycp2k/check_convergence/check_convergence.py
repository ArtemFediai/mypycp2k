# in a perfect world this should be a class
import re

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
    plt.figure(figsize=[4, 3])
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
            use_param_value = params[i]
            print(f"use {name}: {use_param_value}")
            break
        else:
            print("not (yet) reached")


def extract_total_energy(energies, path_to_file):
    """
    takes energy as a list and append it with a new energy extracted from file path_to_file
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
    """
    this version catches only files ended with {extension}. Will ignore "wfn.bak-1. Deprecated"
    :param directory:
    :param extension:
    :return:
    """
    from os import remove
    files = list_files2(directory, extension)
    for file in files:
        remove(file)


def list_files2(directory, extension):
    """
    this version catches anything with ".{extension} inside including i.e. *.wfn.bak-1 files"
    :param directory:
    :param extension:
    :return:
    """
    from os import listdir
    regex = re.compile('.*\.'+extension)
    # print([f for f in listdir(directory) if regex.match(f)])
    return [f for f in listdir(directory) if regex.match(f)]