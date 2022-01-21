import os

config = {}

######################################################################################################################
# Input files
######################################################################################################################
config['input'] = {}
config['input']['target_name'] = 'severe acute respiratory syndrome coronavirus2' #Scientific species name of the target
config['input']['sabdab_summary'] = os.getcwd()+'../demo-12-2021/metadata/sabdab_summary_all.tsv'
# Species to exclude from the database:
config['input']['db_exclude'] = ('severe acute respiratory syndrome', 'sars coronavirus', 'middle east respiratory')

# Minimum sequence length for target protein:
#config['input']['min_seq_len_target'] = 900

# MaSIF path
config['input']['masif_path'] = '/masif/'
######################################################################################################################
# Constants
######################################################################################################################
config['n_threads'] = 8
config['p_val_threshold'] = 1
config['blast_const'] = {}
config["blast_const"]["surface_acc_threshold"] = 0.2
config['blast_const']["contact_threshold"] = 5

# Filtering criteria
#config['min_tm_score'] = 0.17

# RMSD_threshold_high - high confidence threshold
# RMSD_threshold_medium - medium confidence threshold

# config['rmsd_threshold_dict'] = {'Motif_length': [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
#                                                   22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
#                                  'RMSD_threshold_high': [0.01, 0.24, 0.65, 1.08, 1.29, 1.88, 1.88, 2.17, 2.61, 2.93,
#                                                       3.26, 3.58, 3.9, 4.23, 4.55, 4.88, 5.2, 5.52, 5.85, 6.17, 6.5,
#                                                       6.82, 7.14, 7.47, 7.79, 8.12, 8.44, 8.77, 9.09, 9.41],
#                                 'RMSD_threshold_medium': [0.1, 0.69, 1.12, 1.6, 1.87, 2.38, 2.53, 2.81, 3.36, 3.75, 4.13,
#                                                           4.51, 4.9, 5.28, 5.67, 6.05, 6.43, 6.82, 7.2, 7.58, 7.97, 8.35,
#                                                           8.73, 9.12, 9.5, 9.88, 10.27, 10.65, 11.04, 11.42]
#                                  }
#                                'RMSD_threshold_medium': [0.03, 0.5, 0.93, 1.44, 1.59, 2.13, 2.25, 2.55, 3.04, 3.4, 3.76,
#                                                  4.12, 4.48, 4.84, 5.2, 5.55, 5.91, 6.27, 6.63, 6.99, 7.35, 7.71,
#                                                 8.07, 8.42, 8.78, 9.14, 9.5, 9.86, 10.22, 10.58]
#config['max_masif_score'] = 3


######################################################################################################################
# Directories For Intermediate files
######################################################################################################################

################
### DATA PREPARE
config['dirs'] = {}
config['dirs']['data_prepare'] = os.getcwd() + '/data_preparation/'
config['dirs']['lists'] = './lists/'

config['dirs']['raw_pdb'] = config['dirs']['data_prepare'] + '00-raw_pdbs/'
config['dirs']['protonated_pdb'] = config['dirs']['data_prepare'] + '01-protonated_pdb/'
config['dirs']['fasta'] = config['dirs']['data_prepare'] + '02-AG_fasta/' # directory with FASTA files of antigens
config['dirs']['fasta_maps'] = config['dirs']['data_prepare'] + '02-AG_fasta_maps/' # files that maps fasta sequences with residue ID


config['dirs']['blast_db'] = config['dirs']['data_prepare'] + '03-blast_db/'
config['dirs']['ab_contact_map'] = config['dirs']['data_prepare'] + '04-ab_contact_map/'

config['dirs']['chains_pdb'] = config['dirs']['data_prepare'] + '05-chains_pdbs/'
config['dirs']['surface_ply'] = config['dirs']['data_prepare'] + '06-surface_ply/'
config['dirs']['patches'] = config['dirs']['data_prepare'] + '07-patches/'
# config['dirs']['dssp_maps'] = config['dirs']['data_prepare'] + '05-dssp_maps/'
config['dirs']['map_patch'] = config['dirs']['data_prepare'] + "08-patch_maps/"

################
### BLAST HITS
config['dirs']['blast_hits'] = os.getcwd() + '/blast_hits/'
config['dirs']['raw_blast_hits'] = config['dirs']['blast_hits'] + '00-raw_hits/' #dir with xml files
#config['dirs']['raw_blast_hits_txt'] = config['dirs']['blast_hits'] + '00-raw_hits_txt/'
config['dirs']['blast_hits_filtered'] = config['dirs']['blast_hits'] + '01-hits_filtered/'

################
### DEEP LEARNING
config['dirs']['descriptors'] = config['dirs']['data_prepare'] + '09-descriptors/'
config['dirs']['masif_model'] = config['input']['masif_path'] + "data/masif_ppi_search/nn_models/sc05/all_feat/model_data/"

################
### Output
config['dirs']['output'] = os.getcwd() + '/output/variants/'
# config['dirs']['out_individual_hits'] = config['dirs']['output'] + "individual_hits/"
config['dirs']['motifs_pdb'] = os.getcwd()+  '/motif_pdbs/'

config['out_files'] = {}
config['out_files']['raw_blast_hits'] = config['dirs']['output'] + '1-raw_blast_hits.tsv'
config['out_files']['filtered_blast_hits'] = os.getcwd() + '/lists/filtered_blast_hits.tsv'
config['out_files']['TM_align_all'] = config['dirs']['output'] + '3-TM_align.tsv'
config['out_files']['TM_align_filtered'] = config['dirs']['output'] + '3-TM_align_filtered.tsv'
config['out_files']['MaSIF_scores'] = config['dirs']['output'] + '4-MaSIF_scores.tsv'
config['out_files']['MaSIF_scores_filtered'] = config['dirs']['output'] + '4-MaSIF_scores_filtered.tsv'

# config['dirs']['chains_pdb'] = config['dirs']['data_prepare'] + '01-chains_pdbs/'
# config['dirs']['surface_ply'] = config['dirs']['data_prepare'] + '02-surface_ply/'
# config['dirs']['patches'] = config['dirs']['data_prepare'] + '03-patches/'
# config['dirs']['map_patch'] = config['dirs']['data_prepare'] + '04-map_patch_dir/'
#
# config['dirs']['db_masif_indx'] = config['dirs']['db'] + 'masif_indx/'
#
# config['dirs']['blast_db'] = config['dirs']['db'] + 'blast_db/'
# config['dirs']['query_fasta'] = config['dirs']['db'] + 'query_fasta/'
# config['dirs']['descriptors_target'] = config['dirs']['db'] + 'descriptors_target/'


# Temporary files
config['dirs']['tmp'] = os.getcwd() + '/tmp/'

# Create Directories
for dir in config['dirs'].values():
    if not os.path.exists(dir):
        os.makedirs(dir)

######################################################################################################################
# Intermediate files
######################################################################################################################

### Database configuration
config['db_list'] = config['dirs']['lists']+ 'selected_db.txt'
config['target_list'] = config['dirs']['lists'] + 'target.txt'

# Some proteins can not be preprocessed due to large number of heavy atoms of MSMS failure
# The following lists will contain all complexes that were processed with "prepare" module
config['db_list_processed'] = './lists/db_processed_variants.txt'
config['target_list_processed'] = './lists/target_processed_variants.txt'

######################################################################################################################
# Parameters of MaSIF-search
######################################################################################################################
config['masif_search'] = {}
config['masif_search']['cache_min_per_pdb'] = 3 # minimum number of patches to be extracted from a single complex
config['masif_search']['min_sc_fit'] = 0.2

######################################################################################################################
# Parameters of Blast
######################################################################################################################

