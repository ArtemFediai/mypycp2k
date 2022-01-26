import matplotlib.pyplot as plt
import rdkit
import yaml
from rdkit import Chem
from util.xyz import XYZ
import os
from shutil import copyfile
import numpy as np

DATA_FOLDER = '../data/structures_GW100'
ELEMENTS = set(['C', 'H', 'O', 'N', 'F'])
DATA_FOLDER_FOR_RESULTS = '../data/organic_structures_GW100'

files_list = os.listdir(DATA_FOLDER)

organic_number = 0  # counter
organic_molecules_xyz = []

#identify
for xyz_file in files_list:
    print(xyz_file)
    xyz_object = XYZ.from_file(os.path.join(DATA_FOLDER, xyz_file))
    xyz_object.identify_atom_types()
    print(xyz_object)
    uniq_elements = xyz_object.unique_atom_types
    print(type(uniq_elements))
    if uniq_elements.issubset(ELEMENTS) and set(['C']).issubset(uniq_elements):
        organic_number += 1
        print('organic')
        organic_molecules_xyz.append(xyz_file)
    else:
        print('not organic')


print(f'Number of organic elements: {organic_number}')
print(f'List of organic molecules: {organic_molecules_xyz}')
organic_molecules_names = [name.split('.xyz')[0] for name in organic_molecules_xyz]
print(organic_molecules_names)


# save results
os.makedirs(DATA_FOLDER_FOR_RESULTS, exist_ok=True)  # create if not yet there
for mol in organic_molecules_xyz:
    copyfile(src=os.path.join(DATA_FOLDER, mol), dst=os.path.join(DATA_FOLDER_FOR_RESULTS, mol))


# get from db
db = 'db'
DB_PATH = os.path.join(f'../data/{db}')
db_files = os.listdir(DB_PATH)
db_files = [file for file in db_files if file.endswith('.yaml') and file.startswith('DB_')]
print(db_files)

db_names = [file.split('DB_')[1].split('.yaml')[0] for file in db_files]
print(db_names)




organic_molecules_names
db_names

remain = set(db_names) - set(organic_molecules_names)

print(remain)

organic_molecules_names_yaml = []
for mol in organic_molecules_names:
    if mol not in db_names:
        print(f'{mol} is not in db_names')
    else:
        organic_molecules_names_yaml.append(mol)

print(organic_molecules_names_yaml)

# work with my yaml files
all_gw_homo = {}
for mol_name in organic_molecules_names_yaml:
    path_to_yaml = os.path.join(DB_PATH, f'DB_{mol_name}.yaml')
    with open(path_to_yaml) as fid:
        yaml_dict = yaml.load(fid, Loader=yaml.SafeLoader)
        print(yaml_dict)
        gw_homo_qzvp = yaml_dict[mol_name]['occs'][4]
        if gw_homo_qzvp == None:
            print(f'homo of {mol_name} is None')
        all_gw_homo[mol_name] = gw_homo_qzvp

print(all_gw_homo)


# get data from internet
data_link = 'https://raw.githubusercontent.com/setten/GW100/master/data/G0W0%40PBE_HOMO_Cvx_def2-QZVPN4.json'
import requests
response = requests.get(data_link)
data_from_internet = response.json()

# for mol_name in organic_molecules_names_yaml:
#     with open(data_from_internet):



# compare
gw_homo_comparision = {}

comp_dict = {}

this_work = []
this_work_round = []
gw100_paper = []

for mol_name in organic_molecules_names_yaml:
    my_value = all_gw_homo[mol_name]
    # try:
    internet_value = data_from_internet['data'][mol_name]
    # except KeyError:
    #     print(f'{mol_name} is not in the internet')
    print(internet_value, my_value)
    comp_dict[mol_name] = {'this_work': my_value,
                           'gw100_paper': internet_value,
                           'this_work_round': round(my_value, ndigits=2)}
    this_work.append(my_value)
    this_work_round.append(round(my_value, ndigits=2))
    gw100_paper.append(internet_value)

# plot it

plt.figure(figsize=[3,3])
plt.scatter(this_work, gw100_paper)
min = np.min([np.min(this_work), np.min(gw100_paper)])
max = np.max([np.max(this_work), np.max(gw100_paper)])

plt.plot([min, max],[min, max], linestyle=':')
# plt.show()
plt.savefig(fname='gw100_comp.png', dpi=600)
plt.close()

# deviation
delta = np.array(this_work) - np.array(gw100_paper)
delta_round = np.array(this_work_round) - np.array(gw100_paper)
plt.figure(figsize=[3, 3])
hist = np.histogram(delta)

plt.hist(x=delta)
plt.show()

plt.hist(x=delta_round)
plt.show()

mae = np.sum(np.abs(delta)) / len(delta)
mae_round = np.sum(np.abs(delta_round)) / len(delta_round)
print(f'mae is {mae} eV')
print(f'mae round is {mae_round} eV')



#
#
#
#

delta_round = []
print(delta_round)
# yet_to_simulate  =

# for mol in organic_molecules_names:
#     path_to_mol = os.path.join(DATA_FOLDER_FOR_RESULTS, mol)
#     print(path_to_mol)
#     m = Chem.MolFrom (path_to_mol)




# def filter_by_elements(data_folder, elements=['C', 'H', 'O', 'N', 'F']):
#     mols_after_filtering = None
#     return mols_after_filtering
#
