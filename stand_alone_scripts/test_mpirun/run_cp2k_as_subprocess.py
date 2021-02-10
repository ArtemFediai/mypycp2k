"""
subprocess mpi run of cp2k
1. setulimit so far made no effect (local PC)
2. commented lines and the function may be important for the cluster
"""
import os
import subprocess
import shlex
import resource


def main():
    threads = 4
    command = f'mpirun -np {threads} /home/artem/soft/cp2k/cp2k-7.1/exe/local/cp2k.popt -i input.inp'
    outstream = open('my_out_file.out', 'w')
    execution_directory = os.path.abspath('.')

    environ = os.environ.copy()
    # environ["PARNODES"] = str(threads)
    # environ["PARA_ARCH"] = "SMP"
    # environ["TM_PAR_FORK"] = "on"

    process = subprocess.Popen(shlex.split(command),
                               stdin=subprocess.PIPE,
                               stdout=outstream,
                               stderr=outstream,
                               preexec_fn=setulimit,
                               cwd=execution_directory,
                               encoding="utf-8",
                               env=environ)

    print('done')


def setulimit():
    resource.setrlimit(resource.RLIMIT_STACK, (-1, resource.RLIM_INFINITY))


if __name__ == '__main__':
    main()