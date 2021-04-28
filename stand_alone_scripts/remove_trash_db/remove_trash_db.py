import os
import os.path
import re
import glob
import numpy as np
import yaml
from util.general import timeit
import shutil
import csv

def main():


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

    valid_db_files, numbers = get_all_db_files()

    print('list of db: ', valid_db_files)
    print('numbers: ', *numbers)

    trash_files = identify_trash_db_files(all_db_files=valid_db_files)
    print("trash_files: ", trash_files)


    print(f"(number of trash files = {len(trash_files)} which is {100 * len(trash_files) / len(valid_db_files)} \% of all files)")

    write_list_of_trash_db(list_of_trash_db=trash_files, list_of_trash_numbers=numbers)

    move_files_to_trash_folder(file_names_to_put_there=trash_files)


@timeit
def move_files_to_trash_folder(file_names_to_put_there, trash_folder_name='trash', copy=False):
    if not os.path.exists(trash_folder_name):
        os.mkdir(path=trash_folder_name)
    for file_name in file_names_to_put_there:
        if copy:
            operation = shutil.copyfile
        else:
            operation = shutil.move
        operation(file_name, trash_folder_name + '/' + file_name)

@timeit
def get_all_db_files(prefix='DB_', num_digits=6, path='.'):
    only_number = rf'^[{prefix}(\d{num_digits}).yaml$]'
    valid_db_file_name = re.compile(only_number)
    all_files_and_folders = os.listdir(path)
    valid_db_files = [item for item in all_files_and_folders if valid_db_file_name.match(item)]
    # only_number_1 = r'(^DB_)(\d{6}).yaml$' this is what is done below
    only_number_2 = r'(^{})(\d{}).yaml$'.format(prefix, '{' + str(num_digits) + '}')
    # on = ' prefix + chr(92) + 'd{6}).yaml$'
    # numbers = [int(re.match(only_number_2, item)[2]) for item in all_files_and_folders if not os.path.isdir(item)]
    numbers = [int(re.match(only_number_2, item)[2]) for item in valid_db_files]
    return valid_db_files, np.sort(numbers)


@timeit
def identify_trash_db_files(all_db_files, debug_mode=False, _create_cashe=True, _cash_name='cache.txt'):
    # time-consuming part. create cache
    if os.path.exists(_cash_name):
        print("I use cache file to identify trash db file. If it is old, info may be incorrect")
        with open(_cash_name) as stream:
            csv_reader = csv.reader(stream)
            return csv_reader.__next__()
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
        with open(_cash_name, 'w') as stream:
            cvs_writer = csv.writer(stream)
            cvs_writer.writerow(trash_db_files)
            print(f'I wrote the list of trash db file into file {_cash_name}')
        return trash_db_files


def write_list_of_trash_db(list_of_trash_db, list_of_trash_numbers, fname_list_db='trash_db.csv', fname_list_numbers='trash_db_numbers.csv'):
    with open(fname_list_db, 'w') as stream:
        csv_writer = csv.writer(stream)
        print(f"Write list of trash db into {fname_list_db}")
        csv_writer.writerow(list_of_trash_db)
    with open(fname_list_numbers, 'w') as stream:
        csv_writer = csv.writer(stream)
        print(f"Write list of trash db into {fname_list_numbers}")
        csv_writer.writerow(list_of_trash_numbers)

def make_name_from_number(number:int, num_digits=6, with_prefix=True, prefix='DB_'):
    if not with_prefix:
        return '{:0>num_digits}'.format(number)  # example: '1' --> '000001'
    elif with_prefix:
        return '{}{:0>{}'.format(prefix, number, num_digits)  # example: 1 --> DB_000001


def make_nums_from_name(names = [], prefix='DB_', num_digits=6, path='.'):
    """
    list --> list
    """
    all_files_and_folders = os.listdir(path)
    only_number_2 = r'(^DB_)(\d{6}).yaml$'
    numbers = [int(re.match(only_number_2, item)[2]) for item in all_files_and_folders]
    return numbers

def get_range(path):
    smallest_number = None
    largest_number = None
    return smallest_number, largest_number

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