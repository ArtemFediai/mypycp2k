import csv
import numpy as np
from joblib import Parallel, delayed
from tmp.return_large_sig_c_mols import main
from util.exceptions import *
from util.general import timeit


N_JOBS = 32  # should be == number of processors you have

DICT_OF_FOLDERS = {LargeSigc: 'list_of_sigc_folders_',
                   SCQPSolutionNotFound: 'list_of_scqp_sol_folders_',
                   SCFNotConvergedNotPossibleToRunMP2: 'list_scf_not_converged_',
                   NaNInGW: 'list_of_sigc_folders_',
                   IterationLimit: 'list_iteration_limit_'}

for x in DICT_OF_FOLDERS.values():
    exec(f'{x} = []')


def process(i):
    return main(i)

@timeit
def res():
    return Parallel(n_jobs=N_JOBS)(delayed(process)(i) for i in range(134))

results = res()

# print('this is results: ', results)

# print('len of results: ', len(results))

for resa in results:
    for i, exc_dict in enumerate(list(DICT_OF_FOLDERS.keys())):
        eval(DICT_OF_FOLDERS[exc_dict]).append(np.array(resa[i]))

for i, exc_dict in enumerate(list(DICT_OF_FOLDERS.keys())):
    tmp = eval(DICT_OF_FOLDERS[exc_dict])
    tmp = np.sort(list(set(np.hstack(tmp))))
    with open(f'{i}.csv', 'w') as stream:
        writer = csv.writer(stream)
        writer.writerow(tmp)