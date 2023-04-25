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
