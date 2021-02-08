# mypycp2k

high level scripts based on pycp2k soft
allows to set up cp2k input files in several lines of python code

Folders:

1.
mypycp2k: high level scripts. every file name corresponds to a cp2k section name and takes the respective pycp2k object as a first argumet.
for instance:

def set_cutoff(DFT, cutoff=900, rel_cutoff=60, ngrids=5):\n
    MGGRID = DFT.MGRID
    MGGRID.Cutoff = cutoff
    MGGRID.Rel_cutoff = rel_cutoff
    MGGRID.Ngrids = ngrids

is placed in the 'mypycp2k/DFT' folder in the file dft.py.

It taked DFT pycp2k object and change it by adding mggrid information.

It is supposed to be called in the custom scripts.

2.
examples: python files and xyz files that uses mypycp2k function to create, save and run cp2k simulation files from python.
the input of these files may be to a large extend duplicating.
