import yaml
from sklearn.metrics import r2_score

"""
IDEA:
(BEGIN)

lumos: 
    2: ...,
    3: ..., 
    4: ...
    extrapolated: 
        -: ... (method1)
        -: ... (method2)
    err:
        -: ...
        -: ...
    r_2:
        -: ...
        -: ...
homos:

    ...

occs:

    ...

virs:

    ...

        
(END)    
"""


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
    extrapolation_method1  (energy vs. basis_functions**-2)
    extrapolation_method2 (energy vs. cardinal_number**-3)
  lumo:
    extrapolation_method1
    extrapolation_method2
  vir:
    extrapolation_method1
    extrapolation_method2
  occ:
    extrapolation_method1
    extrapolation_method2
  lumo_r:
    method1
    method2
  homo_r:
  ...
  occ_r:
  ...
  vir_r:
  ...
  homo_err:
  ...
  vir_err:
  ...
"""


def main():
    # create dict

    with open('test/my_yaml_db.yaml', 'w+') as stream:
        db = yaml.load(stream, Loader=yaml.SafeLoader)
        my_new_molecule = Cp2kOutput(12345)
        #
        cardinal_number, homo, lumo, occ, vir = 2, -9.068935085092464, -0.31584310176186037, -11.748, 1.03
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        cardinal_number, homo, lumo, occ, vir = 3, -9.086503028020836, -0.328626466686336, -12.165, 0.88
        my_new_molecule.add_energies(cardinal_number, homo, lumo, occ, vir)
        cardinal_number, homo, lumo, occ, vir = 4, -9.091535229637515, -0.3355201992206428, -12.355, None
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

        my_new_molecule.plot_it()


class EnergyLevels:
    """
    name = homo/lumo/occ/vir
    """
    allowed_names = ['homo', 'lumo', 'occ', 'vir']
    def __init__(self, name,
                 computed=None,
                 extrapolated=None,
                 r_2=None,
                 std_err=None):
        assert any([name == allowed_name for allowed_name in self.allowed_names])
        if std_err is None:
            std_err = [None, None]
        if r_2 is None:
            r_2 = [None, None]
        if extrapolated is None:
            extrapolated = [None, None]
        if computed is None:
            computed = {2: None, 3: None, 4: None}
        self.name = name
        self.computed = computed
        self.extrapolated = extrapolated
        self.r_2 = r_2
        self.std_err = std_err

    def add_computed_level(self, cardinal_number, homo, lumo, occ, vir):
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

class Cp2kOutput:
    """
    class for handling cp2k output, namely HOMO/LUMO energies for various basis set.
    The purpose is to be used in GW for the basis set extrapolation
    this must produce a library
    """


    def __init__(self, mol_num):
        self.mol_num = mol_num
        self.homos = EnergyLevel()
        self.lumos = EnergyLevel()



        self.homos = {}
        self.lumos = {}
        self.occs = {}
        self.virs = {}
        self.num_orbs = {}
        self.homo = [None, None]
        self.lumo = [None, None]
        self.occ = [None, None]
        self.vir = [None, None]
        self.homo_err = [None, None]
        self.lumo_err = [None, None]
        self.occ_err = [None, None]
        self.vir_err = [None, None]
        self.homo_r2 = [None, None]
        self.lumo_r2 = [None, None]
        self.occ_r2 = [None, None]
        self.vir_r2 = [None, None]

    def __init__old(self, mol_num):
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
        self.homo_err = [None, None]
        self.lumo_err = [None, None]
        self.occ_err = [None, None]
        self.vir_err = [None, None]
        self.homo_r2 = [None, None]
        self.lumo_r2 = [None, None]
        self.occ_r2 = [None, None]
        self.vir_r2 = [None, None]

    @classmethod
    def from_class(cls):
        #  todo: make the object from the yaml file
        #  maybe only a part of it ...
        pass

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

    def extrapolate_energy_old(self):
        """
        method 1: basis_functions. Energy vs. BF**-1
        method 2: cardinal number. Energy vs. CN**-3
        """

        # try:
        # yaml hates numpy ==> float()
        # extrapolation method: [0] --> basis function, [1] --> cardinal number
        try:
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
        except:
            # self.homo[0:1], self.lumo[0:1] = 'not_computed', 'not_computed'
            # self.occ[0:1], self.vir[0:1] = 'not_computed', 'not_computed'
            pass


    def extrapolate_energy(self):
        """
        method 1: basis_functions. Energy vs. BF**-1
        method 2: cardinal number. Energy vs. CN**-3
        """

        # try:
        #  extrapolation method: [0] --> basis function, [1] --> cardinal number
        x = []
        x.append([num_orb**-1.0 for num_orb in list(self.num_orbs.values())])  # num orbs
        x.append([num_orb**-3.0 for num_orb in list(self.num_orbs.keys())])  # cardinal number
        for i, xx in enumerate(x):
            self.homo[i], self.homo_r2[i], self.homo_err[i] = self.my_linregress(xx, list(self.homos.values()))
            self.lumo[i], self.lumo_r2[i], self.lumo_err[i] = self.my_linregress(xx, list(self.lumos.values()))
            self.occ[i], self.occ_r2[i], self.occ_err[i] = self.my_linregress(xx, list(self.occs.values()))
            self.vir[i], self.vir_r2[i], self.vir_err[i] = self.my_linregress(xx, list(self.virs.values()))
        # except:
        #     self.homo[0:1], self.lumo[0:1] = [None, None], [None, None]
        #     self.occ[0:1], self.vir[0:1] = [None, None], [None, None]
        #     self.homo_r2[0:1], self.lumo_r2[0:1] = [None, None], [None, None]
        #     self.occ_r2[0:1], self.vir_r2[0:1] = [None, None], [None, None]

    @staticmethod
    def my_linregress(X, Y):
        """
        linregree from sciepy.stats with removed unnecessary outputs and more importantly
        numpy floats converted to normal floats to be able to save to yaml
        returns:
        intersect (e.g. homo/lumo in the basis set limit)
        R**2 (determination coefficient)
        intersect error (e.g. error of the homo/lumo)
        """
        # import numpy as np
        from scipy.stats import linregress
        if all(isinstance(x, float) for x in X) and all(isinstance(y, float) for y in Y):
            # _, energy, r, _, energy_error = linregress(X, Y)
            result = linregress(X, Y)
            return float(result.intercept), float(result.rvalue)**2, float(result.intercept_stderr)
            # return float(energy), float(r)**2.0, float(energy_error)
        else:
            print('could not extrapolate. some entries are not floating point numbers')
            return None, None, None



    @staticmethod
    def get_intersect(X, Y):
        """
        returns the energy extrapolated to the basis set limit (n_basis_func = Inf) == intersect with y axis
        """
        # yaml hates numpy ==> float()
        import numpy as np
        if all(isinstance(x, float) for x in X) and all(isinstance(y, float) for y in Y):
            return float(np.polyfit(X, Y, deg=1)[1])  # y = a*x+b. b --> [1]
        else:
            print('could not extrapolate. some entries are not floating point numbers')
            return None

    def yield_dict(self):
        """
        turn the content of the instance cp2k_output into a dictionary
        """
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
                'homo_err': self.homo_err,
                'lumo_err': self.lumo_err,
                'occ_err': self.occ_err,
                'vir_err': self.vir_err,
                'homo_r2': self.homo_r2,
                'lumo_r2': self.lumo_r2,
                'occ_r2': self.occ_r2,
                'vir_r2': self.vir_r2,
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
        try:
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
        except:
            print("Cannot create extrapolation figure. Probably, some entries are not floats")

if __name__ == '__main__':
    main()