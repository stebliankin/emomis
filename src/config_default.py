import os

# Configuration sections [DO NOT MODIFY]
config = {}
config['input'] = {}
config['blast_const'] = {}
config['dirs'] = {}
config['out_files'] = {}

######################################################################################################################
# Input files
######################################################################################################################
config['input']['sabdab_summary'] = os.getcwd()+'/metadata/sabdab_summary_all.tsv'
config['input']['target_name'] = 'severe acute respiratory syndrome coronavirus2' #Scientific species name of the target
# Species to exclude from the database:
config['input']['db_exclude'] = ('severe acute respiratory syndrome', 'sars coronavirus', 'middle east respiratory')

# Optional parameter:
# Minimum sequence length for target protein:
config['input']['min_seq_len_target'] = 900

################
### Output
config['dirs']['output'] = os.getcwd() + '/output/'

config['out_files']['raw_blast_hits'] = config['dirs']['output'] + '1-raw_blast_hits.tsv'
config['out_files']['filtered_blast_hits'] = config['dirs']['output'] + '2-filtered_blast_hits.tsv'
config['out_files']['TM_align_all'] = config['dirs']['output'] + '3-TM_align.tsv'
config['out_files']['TM_align_filtered'] = config['dirs']['output'] + '3-TM_align_filtered.tsv'
config['out_files']['MaSIF_scores'] = config['dirs']['output'] + '4-DL_scores.tsv'
config['out_files']['MaSIF_scores_filtered'] = config['dirs']['output'] + '4-DL_scores_filtered.tsv'


######################################################################################################################
# Constants
######################################################################################################################
config['n_threads'] = 8
config['p_val_threshold'] = 0.1
config["blast_const"]["surface_acc_threshold"] = 0.2
config['blast_const']["contact_threshold"] = 5


######################################################################################################################
# Directories For Intermediate files
######################################################################################################################

config['dirs']['data_prepare'] = os.getcwd() + '/data_preparation/'
config['dirs']['lists'] = './lists/'

config['dirs']['raw_pdb'] = config['dirs']['data_prepare'] + '00-raw_pdbs/'
config['dirs']['protonated_pdb'] = config['dirs']['data_prepare'] + '01-protonated_pdb/'
config['dirs']['fasta'] = config['dirs']['data_prepare'] + '02-AG_fasta/' # directory with FASTA files of antigens
config['dirs']['fasta_maps'] = config['dirs']['data_prepare'] + '02-AG_fasta_maps/' # files that maps fasta sequences with residue ID


config['dirs']['blast_db'] = config['dirs']['data_prepare'] + '03-blast_db/'
config['dirs']['ab_contact_map'] = config['dirs']['data_prepare'] + '04-ab_contact_map/'

config['dirs']['masif_path'] = '/masif/'
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
### STRUCTURAL SIMILARITY
config['dirs']['motifs_pdb'] = os.getcwd()+  '/motif_pdbs/'

################
### DEEP LEARNING
config['dirs']['descriptors'] = config['dirs']['data_prepare'] + '09-descriptors/'
config['dirs']['masif_model'] = config['dirs']['masif_path'] + "data/masif_ppi_search/nn_models/sc05/all_feat/model_data/"

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
config['db_list'] = config['dirs']['lists']+ 'db.txt'
config['target_list'] = config['dirs']['lists'] + 'target.txt'

# Some proteins can not be preprocessed due to large number of heavy atoms of MSMS failure
# The following lists will contain all complexes that were processed with "prepare" module
config['db_list_processed'] = './lists/db_processed.txt'
config['target_list_processed'] = './lists/target_processed.txt'

######################################################################################################################
# Parameters of MaSIF-search
######################################################################################################################
config['masif_search'] = {}
config['masif_search']['cache_min_per_pdb'] = 3 # minimum number of patches to be extracted from a single complex
config['masif_search']['min_sc_fit'] = 0.2

