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
    # set up the input_from_yaml file

    xyz_object_error = XYZ.from_file('test/dsgdb9nsd_104557.xyz')
    print('I have saved to: test/cured_xyzfile.xyz')
    xyz_object_error.write('test/cured_xyzfile.xyz')
    xyz_object_error.identify_atom_types()
    print(f'Unqie atom types: {xyz_object_error.unique_atom_types}')
    print('I am done')

class XYZ:
    @classmethod
    def from_file(cls, fin_name):
        with open(fin_name) as fin:
            natoms = int(fin.readline())
            title = fin.readline()[:-1]
            # noinspection SpellCheckingInspection
            coords = np.zeros([natoms, 3], dtype="float64")
            atom_types = []
            for x in coords:
                line = fin.readline().split()
                atom_types.append(line[0])
                try:
                    x[0:3] = list(map(float, line[1:4]))
                except ValueError:
                    import re
                    print('I have captured an error in xyz file format, but I think I know how to cure it!')
                    xyz_broken = line[1:4]
                    print(xyz_broken)
                    #--> float pattern
                    f_begin = '(^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'  # f = float. expr: any float
                    f_middle = '(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*'
                    f_end = '(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$)'
                    any_float = re.compile(f_begin + f_middle + f_end)
                    #<-- float pattern
                    broken_float = re.compile('^.*\*\^.*$')
                    pattern_with_3_groups = r'(^.*)(\*\^)(.*$)'
                    x_string = []
                    for i in xyz_broken:
                        if any_float.match(i):
                            x_string.append(i)
                        elif broken_float.match(i):
                            i_new = re.sub(pattern_with_3_groups, r'\1E\3', i)
                            x_string.append(i_new)
                            print("I have cured the problem")
                        else:
                            print('I cannot cure the problem. Fix your xyz format! Exiting!')
                            raise Exception
                    x[0:3] = list(map(float, x_string[0:3]))
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

    def identify_atom_types(self):
        """
        Returns unique atom types
        @return:
        set of unique atom types
        """
        self.unique_atom_types = set(self.atom_types)

if __name__ == '__main__':
    main()
