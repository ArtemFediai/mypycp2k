"""
the object corresponding to the xyz file: initializes, reads, writes, etc
"""

import numpy as np
from itertools import cycle
from collections import namedtuple


def main():
    """
    this serves as the test
    """
    # import the molecule (xyz)
    xyz_object = XYZ.from_file('test/H2O.xyz')
    # check the box size
    box_size = xyz_object.compute_box_size(offset=8.0)
    print('box sizes: ', box_size)
    xyz_object.write('testfile.xyz')
    print('I am done')
    # set up the input file


class XYZ:
    @classmethod
    def from_file(cls, fin_name):
        with open(fin_name) as fin:
            natoms = int(fin.readline())
            title = fin.readline()[:-1]
            coords = np.zeros([natoms, 3], dtype="float64")
            atom_types = []
            for x in coords:
                line = fin.readline().split()
                atom_types.append(line[0])
                x[0:3] = list(map(float, line[1:4]))

        return cls(coords=coords, atom_types=atom_types, title=title)

    def __init__(self, coords, atom_types, title):
        self.title = title
        self.coords = coords
        self.atom_types = atom_types
        self.n_atoms = np.shape(coords)[0]
        self.box_size = None

        self.namedtuple_xyz = namedtuple("XYZFile", ["coords", "title", "atom_types"]) (coords, title, atom_types)
        # yes, it is duplicated

    def compute_box_size(self, offset=4.0):
        """
        computes the box size where molecule is placed into
        :param offset: offest between the molecule box to the simulation box
        :return: simulation box size
        """
        def di(coords, i):
            dif = coords[:, i].max() - coords[:, i].min()
            return dif

        self.box_size = np.array([di(self.coords, i) for i in [0, 1, 2]]) + offset
        return self.box_size

    def write(self, filename):
        with open(filename, 'w') as stream:
            stream.write("%d\n%s\n" % (self.coords.size / 3, self.title))
            for x, atom_type in zip(self.coords.reshape(-1, 3), cycle(self.atom_types)):
                stream.write("{}   {:.6f}   {:.6f}   {:.6f}\n".format(atom_type, x[0], x[1], x[2]))


if __name__ == '__main__':
    main()
