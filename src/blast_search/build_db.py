from utils import get_date, read_config, check_if_exists

import os


def updated_processed_list(ppi_list, processed_list, config):
    print("[ {} ] Reading the list of PDB IDs from {} file ...".format(get_date(), ppi_list))

    processed_i = 0
    unprocessed_i = 0

    original_list = [x.strip('\n') for x in open(ppi_list)]

    with open(processed_list, 'w') as out_pr:
        for ppi in original_list:
            if not check_if_exists(ppi, config):
                unprocessed_i += 1
            else:
                processed_i += 1
                out_pr.write(ppi + '\n')
        print("Number of processed complexes in {}: {}".format(ppi_list, processed_i))
        print("Number of unprocessed complexes in {}: {}".format(ppi_list, unprocessed_i))

    return None

def build_db(args):

    config=read_config(args)
    db_name = config['db_name']

    target_list_file = config['target_list']
    db_list_file = config['db_list']

    if not os.path.exists(target_list_file):
        raise FileNotFoundError("Target list {} do not exists. Please specify PDB IDs of the target protein".format(target_list_file))
    if not os.path.exists(db_list_file):
        raise FileNotFoundError(
            "Database list {} do not exists. Please specify PDB IDs of the proteins you want to include to the database".format(db_list_file))

    print('[ {} ] Building the database "{}"...'.format(get_date(), db_name))

    updated_processed_list(db_list_file, config['db_list_processed'], config)
    updated_processed_list(target_list_file, config['target_list_processed'], config)

    # print("[ {} ] Creating indicies for MaSIF interaction regions...".format(get_date()))
    # cache_patches(config)


