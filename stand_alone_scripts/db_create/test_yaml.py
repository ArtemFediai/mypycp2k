import yaml

"""
Structure of the database:
Mol_number:
  homos:
    cardinal_number_1: its_ value_of_homo
  lumos:
    ...: ...
  occs:
    ...: ...
  virs:
    ...: ..

"""


def main():
    # create dict
    rank = 1

    with open('my_yaml_db.yaml', 'w+') as stream:
        db = yaml.load(stream, Loader=yaml.SafeLoader)
        my_new_molecule = Cp2kOutput(12345)
        db_upd = my_new_molecule.yield_dict()
        cardinal_number, homo, lumo, occ, vir = 2, 1, 2, 3, 4
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        #
        yaml.dump(db_upd, stream, Dumper=yaml.SafeDumper)
        print(yaml.dump(db_upd))


class Cp2kOutput:
    """
    this must produce a library
    """

    def __init__(self, mol_num):
        self.mol_num = mol_num
        self.homos = {}
        self.lumos = {}
        self.occs = {}
        self.virs = {}

    def add_energies(self, cardinal_number, homo, lumo, occ, vir):
        """
        adds h/l/o/v with a specified cardinal number.
        cardinal number of cc-pvDZ = 2; cc-pvTZ = 3, etc.
        """
        assert cardinal_number in (
        1, 2, 3, 4, 5), "Cardinal Number is different from 1,2,3,4, or 5. Did you expect this?"
        self.homos[cardinal_number] = homo
        self.lumos[cardinal_number] = lumo
        self.occs[cardinal_number] = occ
        self.virs[cardinal_number] = vir

    def add_num_orbitals(self, num_orb):
        self.num_orb = num_orb
        pass

    def yield_dict(self):
        return {self.mol_num:
            {
                'homos': self.homos,
                'lumos': self.lumos,
                'occs': self.occs,
                'virs': self.virs
            }
        }
if __name__ == '__main__':
    main()