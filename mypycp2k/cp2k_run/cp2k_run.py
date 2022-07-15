from util.general import timeit


@timeit
def cp2k_run(cp2k_executable='cp2k.popt ',
             run_type='serial',
             np=1,
             input_file='input_from_yaml.inp',
             xyz_file='coord.xyz',
             output_file='out.out',
             error_file='err.err',
             execution_directory='./',
             type_mpi='openmpi'):
    import os
    import subprocess
    import shlex
    from shutil import copy
    import shutil

    execution_directory = os.path.abspath(execution_directory)
    outstream = open(execution_directory + '/' + output_file, 'w')
    errstream = open(execution_directory + '/' + error_file, 'w')
    environ = os.environ.copy()

    if not execution_directory == './':
        if not os.path.exists(execution_directory + '/' + xyz_file):
            try:
                copy(xyz_file, execution_directory, )
            except shutil.SameFileError:
                "Same file error: it s okay."

    if run_type == 'serial':
        command = f'{cp2k_executable} -i {input_file}'
    elif run_type == 'mpi':
        if type_mpi == 'openmpi':
            command = f'mpirun -np {str(np)} --oversubscribe --bind-to none {cp2k_executable} -i {input_file}'
        elif type_mpi == 'mpich':
            command = f'mpirun -np {str(np)} {cp2k_executable} -i {input_file}'


    process = subprocess.Popen(shlex.split(command),
                               stdin=subprocess.PIPE,
                               stdout=outstream,
                               stderr=errstream,
                               cwd=execution_directory,
                               encoding="utf-8",
                               env=environ)

    process.wait()  # <-- wait until the process is finished
    # return process
