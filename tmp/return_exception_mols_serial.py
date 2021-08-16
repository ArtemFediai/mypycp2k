from extract_functions.extract_from_output import return_gw_energies_advanced
from util.exceptions import *
import csv
import os
import numpy as np
from util.general import timeit

#  make it general for all esception?

PATH = '/shared/user_data/New_GW_FromArtemToPatrick/data/sim'  # path to folders with results, 012345
LEN = 133885
OUTPUT_FILE_NAMES = ('out_diag_2.out', 'out_diag_3.out')
# list_of_sigc_folders, list_of_scqp_sol_folders, list_scf_not_converged, list_nan,  list_iteration_limit = [],[],[]
MAX_LEN = 10
thousend = 1000

NAMES_OF_FOLDERS = ['list_of_sigc_folders_',
                    'list_of_scqp_sol_folders_',
                    'list_scf_not_converged_',
                    'list_nan_',
                    'list_iteration_limit_']

DICT_OF_FOLDERS = {LargeSigc: 'list_of_sigc_folders_',
                   SCQPSolutionNotFound: 'list_of_scqp_sol_folders_',
                   SCFNotConvergedNotPossibleToRunMP2: 'list_scf_not_converged_',
                   NaNInGW: 'list_of_sigc_folders_',
                   IterationLimit: 'list_iteration_limit_'}

for x in DICT_OF_FOLDERS.values():
    exec(f'{x} = []')
# thsnd = 1

def main(thsnd=1):
    dir_content_raw = os.listdir(PATH)
    # print(dir_content_raw)

    dirs = [x for x in dir_content_raw if (x.__len__() == 6) & x.isdecimal()]  # six decimal digits
    dirs = np.sort(dirs)  # this works LOL
    # print(dir_content)
    actual_num_dirs = len(dirs)
    print(len(dirs))
    if actual_num_dirs == LEN:
        print(f"the number of dirs is {LEN} as expected")
    else:
        print(f"the number of dirs is {actual_num_dirs} != {LEN}")
    # dirs = ['000001', ..., '133885']

    # MAIN -->
    large_sigc_folders, list_of_scqp_sol_folders, list_scf_not_converged, list_nan, list_iteration_limit = \
        extract(dirs, thsnd=thsnd)  # wo arg --> LEN
    # <-- MAIN

    # --> write
    print('\n')
    for folder in DICT_OF_FOLDERS.values():
        write_into_csv_file(name_of_folder=folder)
    for folder in DICT_OF_FOLDERS.values():
        print(f"NUMBERS {folder}: {eval(folder)}")
        print(f"LENGTH of {folder}: {len(set(eval(folder)))}")
        print('\n')

    return large_sigc_folders, list_of_scqp_sol_folders, list_scf_not_converged, list_nan, list_iteration_limit 
    # <-- write


#

# @timeit
def extract(dirs,  thsnd=0):
    print(f'make {thsnd} thousands')
    for i, current_dir in enumerate(dirs[thousend*thsnd:thousend*thsnd+thousend]):
        for output_file_name in OUTPUT_FILE_NAMES:
            try:
                return_gw_energies_advanced(path_to_file=f'{PATH}/{current_dir}/{output_file_name}', silent=True)
            except LargeSigc:
                # print(f"SigC folder: {current_dir}/{output_file_name}")
                eval(DICT_OF_FOLDERS[LargeSigc]).append(current_dir)
            except SCQPSolutionNotFound:
                eval(DICT_OF_FOLDERS[SCQPSolutionNotFound]).append(current_dir)
            except SCFNotConvergedNotPossibleToRunMP2:
                eval(DICT_OF_FOLDERS[SCFNotConvergedNotPossibleToRunMP2]).append(current_dir)
            except NaNInGW:
                eval(DICT_OF_FOLDERS[NaNInGW]).append(current_dir)
            except IterationLimit:
                eval(DICT_OF_FOLDERS[IterationLimit]).append(current_dir)
    return eval(DICT_OF_FOLDERS[LargeSigc]), \
           eval(DICT_OF_FOLDERS[SCQPSolutionNotFound]), \
           eval(DICT_OF_FOLDERS[SCFNotConvergedNotPossibleToRunMP2]), \
           eval(DICT_OF_FOLDERS[NaNInGW]), \
           eval(DICT_OF_FOLDERS[IterationLimit])


def write_into_csv_file(name_of_folder=list(DICT_OF_FOLDERS.values())[0]):
    with open(f'{name_of_folder}.csv', 'w') as scv_stream:
        writer = csv.writer(scv_stream)
        writer.writerow(np.sort(list(set(eval(name_of_folder)))))


if __name__ == '__main__':
    main()
