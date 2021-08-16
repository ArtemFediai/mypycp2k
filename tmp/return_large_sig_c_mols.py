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
MAX_LEN = LEN

NAMES_OF_FOLDERS = ['list_of_sigc_folders',
                    'list_of_scqp_sol_folders',
                    'list_scf_not_converged',
                    'list_nan',
                    'list_iteration_limit']

for x in NAMES_OF_FOLDERS:
    exec(f'{x} = []')



def main():

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
        extract(dirs, max_len=MAX_LEN)  # wo arg --> LEN
    # <-- MAIN

    # --> write
    print('\n')
    for folder in NAMES_OF_FOLDERS:
        write_into_csv_file(name_of_folder=folder)
    for folder in NAMES_OF_FOLDERS:
        print(f"NUMBERS {folder}: {eval(folder)}")
        print(f"LENGTH of {folder}: {len(set(eval(folder)))}")
        print('\n')
    # <-- write


#

@timeit
def extract(dirs, max_len=LEN):
    for i, current_dir in enumerate(dirs[0:max_len]):
        if not i % 1000:
            print(f'finished {i // 1000} thousands')
        for output_file_name in OUTPUT_FILE_NAMES:
            try:
                return_gw_energies_advanced(path_to_file=f'{PATH}/{current_dir}/{output_file_name}', silent=True)
            except LargeSigc:
                #print(f"SigC folder: {current_dir}/{output_file_name}")
                list_of_sigc_folders.append(current_dir)
            except SCQPSolutionNotFound:
                list_of_scqp_sol_folders.append(current_dir)
            except SCFNotConvergedNotPossibleToRunMP2:
                list_scf_not_converged.append(current_dir)
            except NaNInGW:
                list_nan.append(current_dir)
            except IterationLimit:
                list_iteration_limit.append(current_dir)
    return list_of_sigc_folders, list_of_scqp_sol_folders, list_scf_not_converged, list_nan, list_iteration_limit


def write_into_csv_file(name_of_folder=NAMES_OF_FOLDERS[0]):
    with open(f'{name_of_folder}.csv', 'w') as scv_stream:
        writer = csv.writer(scv_stream)
        writer.writerow(np.sort(list(set(eval(name_of_folder)))))


if __name__ == '__main__':
    main()

