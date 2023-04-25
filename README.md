# mypycp2k

_Note_: this document describes how to use the software before summer 2022. After this point, the folder structure was changed to make it packagable and the relevant part of the repository became 

## Folders

### 1. mypycp2k
High-level scripts. Every file name corresponds to a CP2K section name and takes the respective pycp2k object as the first argument. For instance:

```python

def set_cutoff(DFT, cutoff=900, rel_cutoff=60, ngrids=5):
    MGGRID = DFT.MGRID
    MGGRID.Cutoff = cutoff
    MGGRID.Rel_cutoff = rel_cutoff
    MGGRID.Ngrids = ngrids
```

This is placed in the ```mypycp2k/DFT``` folder in the file ```dft.py```.
It takes the DFT ```pycp2k``` object and changes it by adding ```MGGRID``` information.
It is supposed to be called in custom scripts.

### 2. examples.
Python files and XYZ files that use ```mypycp2k``` functions to create, save, and run CP2K simulation files from Python. For GW in the basis set limit, the most comprehensive file is ```~/mypycp2k/examples/234_ver_3.py```

## Main scripts in ```examples``` folder

### 1. ```234_ver_3.py```: GW@PBE in basis set limit.
Function: set up and run simulation of the specified molecules from some database.
How to use this script.
##### 1.1. Settings file.
Specifies everything except for specific identifiers of the molecules, which are to be simulated.
The example (let us say this file is called ```input.yaml```) is below:
```yaml
---
# example: 2,3,4 zeta for basis set extrapolation GW
basis_set_file: BASIS_CC_AUG_RI_NEW
basis_set_suffix:
  - 2
  - 3
  - 4
basis_set_list:
  - aug-cc-pVDZ
  - aug-cc-pVTZ
  - aug-cc-pVQZ
ri_basis_set_list:
  - aug-cc-pVDZ-RIFIT
  - aug-cc-pVTZ-RIFIT
  - aug-cc-pVQZ-RIFIT
cutoff: 750
rel_cutoff: 75
molecule_vacuum_offset: 10
eps_scf_dft:
  - 1.0E-9
  - 1.0E-9
  - 1.0E-9
folder_names:
  simulations: sim
  database: db
  scratch: bh5670
debug: True
dummy_run: False
prefix_xyz_file_name: dsgdb9nsd  # --> {}_{rank}.xyz
cp2k_exe_path: /home/ws/bh5670/cp2k/cp2k-7.1/exe/local/cp2k.popt
potential_file_name: 'POTENTIAL'
potential: 'ALL'
elements: ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'B', 'I']
mpi: openmpi

```


##### 1.2. Submit file.
It is desgigned to run on a cluster. Hence, the submit script has to be used to submit the job. 
Attention: it is used array jobs
This is an example (let us say this is called ```run.sh```):

```bash
#!/bin/bash 
#SBATCH --time=03:00:00
#SBATCH --output=./out/out_%A_%a.out
#SBATCH --error=./out/out_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
#SBATCH --tasks-per-node=128
#SBATCH --partition=short
#SBATCH --job-name=7
#SBATCH --array=1-7
#SBATCH --exclude=bionano01,bionano02,bionano03,bionano04,bionano05,bionano06

# first non-empty non-comment line ends SLURM options
# cp2k.8 is in your PATH

module load gnu8/8.3.0
module load openmpi3

CP2K_popt=cp2k.popt


cd ${SLURM_SUBMIT_DIR}
echo slurm_submit_dir, $SLURM_SUBMIT_DIR
export SCRATCH=/scratch

mkdir $SCRATCH/bh5670  # if no such folder --> scratch
mkdir $SCRATCH/bh5670/sim  # simulations --> scratch

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=FALSE
export SLURM_CPU_BIND=none
unset SLURM_CPU_BIND

source activate mypycp2k_new

which python

python ~/mypycp2k/examples/234_ver_3.py -rank $SLURM_ARRAY_TASK_ID -num_cpus $SLURM_NTASKS -i input.yaml -mol_ids 7.csv

```

Some important points:
1. Create the folder ```out``` before you run the script in the work directory. Otherwise, you may not save stdout and stderr, and the program crashes. Yes, this is inconvenient.
2. `7.csv` is a file that contains one line of comma-separated numbers of the molecules from the database that are to be simulated, for instance (content of `csv` file): 

```csv
005852,007613,045353,075965,097131,106026,132458,
```
3. You need to manually specify the position of the identifiers in the mentioned csv file this way:
```bash
#SBATCH --array=1-7
```
This, means that the entries in `csv` file at positions 1,2,3,4,5,6,7 are to be simulated. Would we have more comma-separated identifiers in `csv` file, would these ignored. If the list is not continues, the syntax of the SLURM applies. You can set: `1,5-6`, for instance.

#### 1.3 Output structure.
 Three folders will be created or filled with generated data as follows:
 * `out`. stdout and stderr are redirected here. If something crashes, you can see the logs there. Filenames generated there will be something like: `out_2461885_87.out` (first number is the job_id, second number is the position of the identifier in the `csv` file). Yes, this is not convinient, yet preserve from clashes.
 * `db`. The output of the GW simulation are saved in these folders under the name `DB_{mol_id}.yaml`, e.g. `DB_061102.yaml`.
 * `sim`. If the respecting option is opts, here you will have raw data copy backed from scratch. For every molecule, the folder with `mol_id` will be created and all files from the respecting `scratch` directory will be copied there. The content will be something like:

```bash
2-RESTART.wfn         3-RESTART.wfn.bak-1   4-RESTART.wfn.bak-2   err.err               out_diag_3.out        
2-RESTART.wfn.bak-1   3-RESTART.wfn.bak-2   diag_2.inp            ot_2.inp              out_diag_4.out        
2-RESTART.wfn.bak-2   3-RESTART.wfn.bak-3   diag_3.inp            ot_3.inp              out_ot_2.out          
2-RESTART.wfn.bak-3   4-RESTART.wfn         diag_4.inp            ot_4.inp              out_ot_3.out          
3-RESTART.wfn         4-RESTART.wfn.bak-1   dsgdb9nsd_069817.xyz  out_diag_2.out        out_ot_4.out          
```

