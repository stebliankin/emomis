from datetime import datetime
from importlib.machinery import SourceFileLoader
import os

def get_date():
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S")

def read_config(args):
    # Read config
    if not args.config:
        from config_default import config
    else:
        config_module = SourceFileLoader("config", args.config).load_module()
        config = config_module.config

    print("[ {} ] Configuration parameters:".format(get_date()))
    print(config)
    return config

def read_ppi_list(args):
    if '.txt' in args.ppi:
        print("[ {} ] Reading list {} ...".format(get_date(), args.ppi))
        ppi_list = [x.strip('\n') for x in open(args.ppi, 'r')]
    else:
        ppi_list = [args.ppi]
    return ppi_list

def check_if_exists_masif(ppi, config):
    pid, ch1 = ppi.split('_')[0], ppi.split('_')[1]
    # pdb.set_trace()
    precompute_dir = config['dirs']['patches']

    if (not os.path.exists('{}/{}/{}_{}_iface_labels.npy'.format(precompute_dir, pid, pid, ch1))):
        return False
    else:
        return True

