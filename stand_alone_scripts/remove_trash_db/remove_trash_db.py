"""
Sometimes, one needs to identify those DB_abcdefg.yaml that do not contain data (result of unsuccesfull cp2k run).
These are set of functions and a script to deal with that.
Identify DB_ files without data
Save the list (csv) of such data
Save the number of reference of such data
"""

import os
import os.path
import re
import glob
import numpy as np
import yaml
from util.general import timeit
import shutil
import csv
import argparse


def main():
    # Parsing -->
    parser = argparse.ArgumentParser(description='arguments supplied to the remove trash function')
    parser.add_argument('-format',
                        nargs='?')
    args = parser.parse_args()
    format_of_the_xyz_file = args.format

    # <-- Parsing

    # path = '.'  # current dir
    # only_number = r'^[DB_(\d{6}).yaml$]'
    # only_number_2 = r'(^DB_)(\d{6}).yaml$'
    # valid_db_file_name = re.compile(only_number)
    #
    # prefix = 'DB_'
    # num_digits = 6
    # list_of_db = [f for f in glob.glob("DB_*.yaml")]
    # all_files_and_folders = os.listdir(path)
    # valid_db_files = [item for item in all_files_and_folders if valid_db_file_name.match(item)]
    # six_digit_numbers = [int(re.match(only_number_2, item)[2]) for item in all_files_and_folders]

    # return db files -->
    if format_of_the_xyz_file is None:
        valid_db_files, numbers = get_all_db_files()
    else:
        valid_db_files, numbers = get_all_db_files_general()
    print('list of db: ', valid_db_files)
    try:
        print('numbers: ', *numbers)
    except TypeError:
        print('numbers cannot be returned because general format of xyz file name is used. '
              'Only QM9 dataset will return this')
    # <-- return db files

    # identify trash files -->
    trash_files, trash_numbers = identify_trash_db_files(all_db_files=valid_db_files,
                                                         use_cashe_if_exists=False)  # TODO cash is broken
    print("trash_files: ", trash_files)

    try:
        print(
            f"(number of trash files = {len(trash_files)} which is {100 * len(trash_files) / len(valid_db_files)} \% of all files)")
    except ZeroDivisionError:
        print(
            "No valid db files found. Probably you are in a wrong folder. Go to a folder with the simulation results (db/)")
    write_list_of_trash_db(list_of_trash_db=trash_files, list_of_trash_numbers=trash_numbers)

    move_files_to_trash_folder(file_names_to_put_there=trash_files)
    # <-- identify trash files

    min_max_tuple, missing_db, missing_numbers = get_range_and_missing_items(valid_db_files=valid_db_files,
                                                                             numbers=numbers)
    print('missing numbers:', np.sort(np.array(missing_numbers)))

    broken_num = get_broken_items(min_max_tuple)
    print("Numbers of broken xyz file of the original db", broken_num)


@timeit
def move_files_to_trash_folder(file_names_to_put_there,
                               trash_folder_name='trash',
                               copy=False):
    if not os.path.exists(trash_folder_name):
        os.mkdir(path=trash_folder_name)
    for file_name in file_names_to_put_there:
        if copy:
            operation = shutil.copyfile
        else:
            operation = shutil.move
        operation(file_name, trash_folder_name + '/' + file_name)


@timeit
def get_all_db_files(prefix='DB_',
                     num_digits=6,
                     path='.'):
    # only_number = '^DB_[0-9]{6}.yaml$'
    only_number = '^{}[0-9]{}.yaml$'.format(prefix, '{' + str(num_digits) + '}')
    valid_db_file_name = re.compile(only_number)
    all_files_and_folders = os.listdir(path)
    valid_db_files = [item for item in all_files_and_folders if valid_db_file_name.match(item)]
    # only_number_1 = r'(^DB_)(\d{6}).yaml$' this is what is done below
    only_number_2 = r'(^{})([0-9]{}).yaml$'.format(prefix, '{' + str(num_digits) + '}')
    # numbers = [int(re.match(only_number_2, item)[2]) for item in all_files_and_folders if not os.path.isdir(item)]
    numbers = [int(re.match(only_number_2, item)[2]) for item in valid_db_files]
    return valid_db_files, np.sort(numbers)


@timeit
def get_all_db_files_general(prefix='DB_',
                             path='.'):
    """
    General format. In contrast to one
    @param prefix: every record starts here
    @param num_digits: assumed that the prefix is followed by this number of digits
    @param path: self-explaining
    @return: list of records, sorted numbers
    """
    # only_number = '^DB_[0-9]{6}.yaml$'
    only_number = '^{}[a-zA-Z0-9_]*.yaml$'.format(prefix)
    valid_db_file_name = re.compile(only_number)
    all_files_and_folders = os.listdir(path)
    valid_db_files = [item for item in all_files_and_folders if valid_db_file_name.match(item)]
    # only_number_1 = r'(^DB_)(\d{6}).yaml$' this is what is done below
    # only_number_2 = r'(^{})([0-9]{}).yaml$'.format(prefix, '{' + str(num_digits) + '}')
    # numbers = [int(re.match(only_number_2, item)[2]) for item in all_files_and_folders if not os.path.isdir(item)]
    # numbers = [int(re.match(only_number_2, item)[2]) for item in valid_db_files]
    return valid_db_files, None


@timeit
def identify_trash_db_files(all_db_files,
                            debug_mode=False,
                            _create_cashe=True,
                            _cashe_file_name='cache.txt',
                            use_cashe_if_exists=False):
    # time-consuming part. create cache
    if os.path.exists(_cashe_file_name) and use_cashe_if_exists:
        print("I use cache file to identify trash db file. If it is old, info may be incorrect")
        with open(_cashe_file_name) as stream:
            csv_reader = csv.reader(stream)
            return csv_reader.__next__()  # only one line in this csv format file, so we do not loop over
    else:
        trash_db_files = []
        for file in all_db_files:
            with open(file) as stream:
                dict_from_file = yaml.load(stream, Loader=yaml.SafeLoader)
                not_extracted_values_list = get_recursively(dict_from_file, 'not extracted')
                null_values_list = get_recursively(dict_from_file, None)
                if debug_mode:
                    print('not_extracted_keys_list:', not_extracted_values_list)
                    print('null_values_list:', null_values_list)
                if not_extracted_values_list or null_values_list:
                    trash_db_files.append(file)
        with open(_cashe_file_name, 'w') as stream:
            cvs_writer = csv.writer(stream)
            cvs_writer.writerow(trash_db_files)
            print(f'I wrote the list of trash db file into file {_cashe_file_name}')
        try:
            trash_db_numbers = make_nums_from_name(names=trash_db_files)
        except TypeError:
            trash_db_numbers = None
            print('I am in general mode. db numbers cannot be identified')
        return trash_db_files, trash_db_numbers


def write_list_of_trash_db(list_of_trash_db, list_of_trash_numbers, fname_list_db='trash_db.csv',
                           fname_list_numbers='trash_db_numbers.csv'):
    with open(fname_list_db, 'w') as stream:
        csv_writer = csv.writer(stream)
        print(f"Write list of trash db into {fname_list_db}")
        csv_writer.writerow(list_of_trash_db)
    with open(fname_list_numbers, 'w') as stream:
        csv_writer = csv.writer(stream)
        print(f"Write list of trash db into {fname_list_numbers}")
        csv_writer.writerow(list_of_trash_numbers)


def make_name_from_number(number, num_digits=6, with_prefix=True, prefix='DB_'):
    if not with_prefix:
        return '{:0>{}}'.format(number, num_digits)  # example: '1' --> '000001'
    elif with_prefix:
        return '{}{:0>{}}'.format(prefix, number, num_digits)  # example: 1 --> DB_000001  # todo does not work


def make_nums_from_name(names=[], prefix='DB_', num_digits=6):
    """
    list --> list
    """
    # only_number_2 = r'(^DB_)(\d{6}).yaml$'

    only_number_2 = r'(^{})([0-9]{}).yaml$'.format(prefix, '{' + str(num_digits) + '}')
    numbers = [int(re.match(only_number_2, item)[2]) for item in names]
    return numbers


def get_range_and_missing_items(valid_db_files,
                                numbers,
                                fname_missing_db='missing_db.csv',
                                fname_missing_db_numbers='missing_num.csv'):
    smallest_number = min(numbers)
    largest_number = max(numbers)
    wish_list = np.linspace(smallest_number, largest_number, largest_number - smallest_number + 1, dtype=int)
    actual_list = numbers
    missing_numbers = list(set(wish_list) - set(actual_list))
    missing_numbers = np.sort(missing_numbers)
    # print(set(actual_list).__len__())
    # print(set(wish_list).__len__())
    # print(len(missing_numbers))
    # missing_numbers = [number for number in wish_list if number not in actual_list]
    missing_dbs = [make_name_from_number(number) for number in missing_numbers]
    with open(fname_missing_db, 'w') as stream:
        csv_writer = csv.writer(stream)
        csv_writer.writerow(missing_dbs)
    with open(fname_missing_db_numbers, 'w') as stream:
        csv_writer = csv.writer(stream)
        csv_writer.writerow(missing_numbers)

    return (smallest_number, largest_number), missing_dbs, missing_numbers


@timeit
def get_broken_items(small_large_number, db='dsgdb9nsd'):
    """
    small_large_number = (num_smallest, num_largest)
    """

    min_num, max_num = small_large_number
    broken_num = []

    for number in range(min_num, max_num + 1):
        six_digits_number = make_name_from_number(number, with_prefix=False)
        path_to_orig_xyz_file = f'../{db}/{db}_{six_digits_number}.xyz'
        if exception_found(path_to_orig_xyz_file):
            broken_num.append(six_digits_number)

    fname_list_numbers = 'broken_xyz_files.csv'

    with open(fname_list_numbers, 'w') as stream:
        csv_writer = csv.writer(stream)
        print(f"Write list of broken xyz numbers into {fname_list_numbers}")
        csv_writer.writerow(broken_num)
    print(f"I wrote the numbers of broken xyz files into {fname_list_numbers}")

    return broken_num


def exception_found(fin_name):
    """
    folder is normally outside the folder we run the program from
    """
    import re
    broken_float = re.compile('^.*\*\^.*$')
    with open(fin_name) as fin:
        natoms = int(fin.readline())  # name
        _ = fin.readline()[:-1]  # n_atoms
        for _ in range(natoms):
            line = fin.readline().split()
            for i in line[1:4]:
                if broken_float.match(i):
                    print('this mol contains a broken float. I will add it to the list list_of_broken_floats')
                    return True
    return False


def get_recursively(search_dict, field):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():

        if value == field:
            fields_found.append(key)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)
                if item == field:
                    fields_found.append(key)

    return fields_found


if __name__ == '__main__':
    main()
