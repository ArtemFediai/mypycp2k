import yaml

"""
Structure of the database:
Mol_number:
  homos:
    cardinal_number_1: its_ value_of_homo
    cardinal_number_2: its_ value_of_homo
    cardinal_number_3: its_ value_of_homo
  lumos:
    ...: ...
  occs:
    ...: ...
  virs:
    ...: ..
  homo:
    extrapolation_1
    extrapolation_2
  lumo:
    extrapolation_1
    extrapolation_2
  vir:
    extrapolation_1
    extrapolation_2
  occ:
    extrapolation_1
    extrapolation_2
"""


def main():
    # create dict
    rank = 1

    with open('my_yaml_db.yaml', 'w+') as stream:
        db = yaml.load(stream, Loader=yaml.SafeLoader)
        my_new_molecule = Cp2kOutput(12345)
        #
        cardinal_number, homo, lumo, occ, vir = 2, -9.068935085092464, -0.31584310176186037, -11.748, 1.03
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        cardinal_number, homo, lumo, occ, vir = 3, -9.086503028020836, -0.328626466686336, -12.165, 0.88
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        cardinal_number, homo, lumo, occ, vir = 4, -9.091535229637515, -0.3355201992206428, -12.355, 0.827
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        #
        my_new_molecule.add_num_orbitals(cardinal_number=2, num_orb=41)
        my_new_molecule.add_num_orbitals(cardinal_number=3, num_orb=92)
        my_new_molecule.add_num_orbitals(cardinal_number=4, num_orb=172)
        #
        my_new_molecule.extrapolate_energy()
        #
        db_upd = my_new_molecule.yield_dict()
        yaml.dump(db_upd, stream, Dumper=yaml.SafeDumper)
        print(yaml.dump(db_upd))

        my_new_molecule.status()
        print('status: ', my_new_molecule.status())

        #my_new_molecule.plot_it()

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
        self.num_orbs = {}
        self.homo = [None, None]
        self.lumo = [None, None]
        self.occ = [None, None]
        self.vir = [None, None]

    def add_energies(self, cardinal_number, homo, lumo, occ, vir):
        """
        adds h/l/o/v with a specified cardinal number.
        cardinal number of cc-pvDZ = 2; cc-pvTZ = 3, etc.
        """
        assert cardinal_number in (
        1, 2, 3, 4, 5), "Cardinal Number is different from 1,2,3,4, or 5. Did you expect this?"
        # assert is too rude. Warning would suffice
        self.homos[cardinal_number] = homo
        self.lumos[cardinal_number] = lumo
        self.occs[cardinal_number] = occ
        self.virs[cardinal_number] = vir

    def add_num_orbitals(self, cardinal_number, num_orb):
        self.num_orbs[cardinal_number] = num_orb
        pass

    def extrapolate_energy(self):
        """
        method 1: basis_functions. Energy vs. BF**-1
        method 2: cardinal number. Energy vs. CN**-3
        """

        # try:
        # yaml hates numpy ==> float()
        # [0] --> basis function, [1] --> cardinal number

        x = [num_orb**-1.0 for num_orb in list(self.num_orbs.values())]
        self.homo[0] = self.get_intersect(X=x, Y=list(self.homos.values()))
        self.lumo[0] = self.get_intersect(X=x, Y=list(self.lumos.values()))
        self.occ[0] = self.get_intersect(X=x, Y=list(self.occs.values()))
        self.vir[0] = self.get_intersect(X=x, Y=list(self.virs.values()))
        #
        x = [num_orb**-3.0 for num_orb in list(self.num_orbs.keys())]
        self.homo[1] = self.get_intersect(X=x, Y=list(self.homos.values()))
        self.lumo[1] = self.get_intersect(X=x, Y=list(self.lumos.values()))
        self.occ[1] = self.get_intersect(X=x, Y=list(self.occs.values()))
        self.vir[1] = self.get_intersect(X=x, Y=list(self.virs.values()))

        # except:
        #     self.homo = 'not computed'


    @staticmethod
    def get_intersect(X, Y):
        import numpy as np
        if all(isinstance(x, float) for x in X) and all(isinstance(y, float) for y in Y):
            return float(np.polyfit(X, Y, deg=1)[1])  # y = a*x+b. b --> [1]
        else:
            print('could not extrapolate. some are not floating point numbers')
            return None


    def yield_dict(self):
        return {self.mol_num:
            {
                'homos': self.homos,
                'lumos': self.lumos,
                'occs': self.occs,
                'virs': self.virs,
                'num_orb': self.num_orbs,
                'homo': self.homo,
                'lumo': self.lumo,
                'occ': self.occ,
                'vir': self.vir,
            }
        }

    def status(self):
        status = 'all_extracted'
        all_self_attr = self.__dict__
        all_dicts_and_lists = [attr_value for attr_value in all_self_attr.values()
                               if isinstance(attr_value, list) or isinstance(attr_value, dict)]
        for my_dict_or_list in all_dicts_and_lists:
            for item in my_dict_or_list:
                if not (isinstance(item, float) or isinstance(item, int)):
                    status = f'not everything is extracted'
        return status

    def plot_it(self):
        import numpy as np
        import matplotlib.pyplot as plt

        # occ vs num_orb
        plt.plot()
        x = [num_orb**-1.0 for num_orb in list(self.num_orbs.values())]
        y = list(self.occs.values())
        [a, b] = np.polyfit(x, y, deg=1)
        plt.plot(x, y, color='C0')
        plt.plot(0, b, '*', color='C0', label='orb')
        # occ vs card_num
        x = [n**-3.0 for n in list(self.num_orbs.keys())]
        y = list(self.occs.values())
        [a, b] = np.polyfit(x, y, deg=1)
        plt.plot(x, y, color='C1', label='car')
        plt.plot(0, b, '*', color='C1')
        plt.legend()
        plt.savefig('occs.png')
        plt.close()

        # vir vs num_orb
        plt.plot()
        x = [num_orb ** -1.0 for num_orb in list(self.num_orbs.values())]
        y = list(self.virs.values())
        [a, b] = np.polyfit(x, y, deg=1)
        plt.plot(x, y, color='C0')
        plt.plot(0, b, '*', color='C0', label='orb')
        # occ vs card_num
        x = [n ** -3.0 for n in list(self.num_orbs.keys())]
        y = list(self.virs.values())
        [a, b] = np.polyfit(x, y, deg=1)
        plt.plot(x, y, color='C1', label='car')
        plt.plot(0, b, '*', color='C1')
        plt.legend()
        plt.savefig('virs.png')

if __name__ == '__main__':
    main()