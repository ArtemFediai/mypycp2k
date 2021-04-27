"""
subprocess mpi run of cp2k
1. setulimit so far made no effect (local PC)
2. commented lines and the function may be important for the cluster

run with an argument:
... -num_cpus 4

"""
import os
import subprocess
import shlex
import resource
import argparse
import shutil

def main():


    parser = argparse.ArgumentParser(description='number of cpus')
    parser.add_argument('-num_cpus')  # number of cpus you request for every array job
    parser.add_argument('-num_array_job')  # number of cpus you request for every array job
    args = parser.parse_args()
    num_cpus = args.num_cpus
    num_array_job = args.num_array_job

    # command = f'mpirun -np {num_cpus} /home/artem/soft/cp2k_test/cp2k-7.1/exe/local/cp2k.popt -i input.inp'

    created_dir = str(num_array_job)
    os.mkdir(created_dir)
    shutil.copyfile('input.inp', created_dir + '/' + 'input.inp')
    os.chdir(created_dir)


    command = f'mpirun -np {num_cpus} /home/ws/bh5670/cp2k/cp2k-7.1/exe/local/cp2k.popt -i input.inp'
    outstream = open('cp2k.out', 'w')
    execution_directory = os.path.abspath('.')
    print('execution_directory', execution_directory)

    environ = os.environ.copy()
    # environ["PARNODES"] = str(num_cpus)
    # environ["PARA_ARCH"] = "SMP"
    # environ["TM_PAR_FORK"] = "on"

    # process = subprocess.Popen(shlex.split(command),
    #                            stdin=subprocess.PIPE,
    #                            stdout=outstream,
    #                            stderr=outstream,
    #                            preexec_fn=setulimit,
    #                            cwd=execution_directory,
    #                            encoding="utf-8",
    #                            env=environ)

    process = subprocess.Popen(shlex.split(command),
                               stdin=subprocess.PIPE,
                               stdout=outstream,
                               stderr=outstream,
                               cwd=execution_directory,
                               encoding="utf-8",
                               env=environ)
    process.wait()
    # return process  # <-- this will crash!
    print('done')


def setulimit():
    resource.setrlimit(resource.RLIMIT_STACK, (-1, resource.RLIM_INFINITY))


if __name__ == '__main__':
    main()