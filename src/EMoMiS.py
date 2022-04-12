import pymesh
import argparse
import pdb

from utils import read_config, get_date
from data_prepare.data_prepare import extract_ppi_lists, download, build_blast, ab_contact_map, get_processed_patches, prepare_masif
from time import time
from blast_search.blast_search import run_blast
from blast_search.blast_filter import blast_filter, report_hits
from deep_learning.evaluate_binding import evaluate_binding, filter_MaSIF, compute_native_scores
from deep_learning.compute_descriptors import compute_descriptors
from structural_alignment.structural_alignment import structural_alignment, filter_rmsd

parser = argparse.ArgumentParser()

parser.add_argument('--config', help='optional config file')
parser.add_argument('--reverse', default=False, action="store_true", help='If set True, reverse target and database lists.')
parser.add_argument('--skip_blast', default=False, action="store_true", help='If set True, skip "Blast search part".')
parser.add_argument('--skip_sabdab', default=False, action="store_true", help='If set True, do not download SAbDab".')

args = parser.parse_args()
reverse_flag = args.reverse

start = time()
config = read_config(args)

if "min_seq_len_target" in config['input'].keys():
    min_seq_len = config['input']['min_seq_len_target']
else:
    min_seq_len=None

#
if not args.skip_sabdab:
    ppi_list_db, ppi_list_target = extract_ppi_lists(config, reverse_flag)
else:
    ppi_list_db = [x.strip('\n') for x in open(config['db_list']).readlines()]
    ppi_list_target = [x.strip('\n') for x in open(config['target_list']).readlines()]
#
updated_ppi_list_db = download(ppi_list_db, config, to_write=config['db_list'], min_seq_len=None)
updated_ppi_list_target = download(ppi_list_target, config, to_write=config['target_list'], min_seq_len=min_seq_len)

if not args.skip_blast:
    build_blast(updated_ppi_list_db, config)
    ab_contact_map(updated_ppi_list_db+updated_ppi_list_target, config)
    run_blast(updated_ppi_list_target, config)
    blast_filter(updated_ppi_list_target, config)
    print("*** [ {} ] Reporting blast hits...".format(get_date()))
    updated_ppi_list_db = report_hits(updated_ppi_list_target, updated_ppi_list_db, config)
# #
print("*** [ {} ] Computing structural alignment with TM-align...".format(get_date()))
updated_ppi_list_db = structural_alignment(config)
updated_ppi_list_db, updated_ppi_list_target = filter_rmsd(config)
# # # # #
print("*** [ {} ] Start computing MaSIF patches for {} target proteins...".format(get_date(), len(updated_ppi_list_target)))
updated_ppi_list_target = prepare_masif(updated_ppi_list_target, config)

print("*** [ {} ] Start computing MaSIF patches for {} proteins from DB...".format(get_date(), len(updated_ppi_list_db)))
updated_ppi_list_db = prepare_masif(updated_ppi_list_db, config)
#
# # #
updated_ppi_list_db = get_processed_patches(updated_ppi_list_db, config)
print("Number of database PDBs with computed patches: {}".format(len(updated_ppi_list_db)))
# #
updated_ppi_list_target = get_processed_patches(updated_ppi_list_target, config)
print("Number of target PDBs with computed patches: {}".format(len(updated_ppi_list_target)))
# #
print("*** [ {} ] Computing descriptors for DB...".format(get_date()))
compute_descriptors(updated_ppi_list_db, config, list_type="subject")
print("*** [ {} ] Computing descriptors for target...".format(get_date()))
compute_descriptors(updated_ppi_list_target, config, list_type="target")
# #
evaluate_binding(updated_ppi_list_target, config)
filter_MaSIF(updated_ppi_list_db, config)

print("**********************************")
print("**** Output:")
print("**** Final molecular mimicry hits: {}".format(config['out_files']['MaSIF_scores_filtered']))
print("**********************************")
print("Total time: {:.2f} min".format((time()-start)/60))