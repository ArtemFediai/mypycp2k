import yaml
from os import listdir

from extract_functions.basis_set_extrapolation import *

print(os.listdir())

with open(file='DB_023772.yaml', mode='r+') as f:
    db = yaml.load(f, Loader=yaml.SafeLoader)

    my = Cp2kOutput.from_yaml('DB_023772.yaml')
    my.mol_num = '023772'
    my.extrapolate_energy()
    my.yield_dict()

    db_upd = my.yield_dict()
    yaml.dump(db_upd, f, Dumper=yaml.SafeDumper)
    print(yaml.dump(db_upd))

    my.status()
    print('status: ', my.status())

print(my)
