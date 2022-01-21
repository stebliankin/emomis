from Bio import Entrez,SeqIO
import os
# Here we will output PDB ids of all complexes required for testing sequence identity.
# The output was used to download FASTA sequences from https://www.rcsb.org/downloads/fasta

# # Download SAbDAB test set
# test_pids = [x.strip('\n').split('_')[0] for x in open('lists/test_filtered.txt', 'r').readlines()]
#
# all_pids = ''
# for pid in test_pids:
#     all_pids+=pid+','
# all_pids = all_pids.strip(',')
# print('SAbDab test PDB IDs:')
# print(all_pids)

masif_pid_list = [x.strip('\n').split('_')[0] for x in open('lists/masif_original_train.txt', 'r').readlines()]

masif_non_redundant = []
masif_batch=[]
all_pids = ''
limit=1000 # output in batches as 1000 is the maximum number of load in rcsb.org
print('Masif train PDB IDs:')

for pid in masif_pid_list:
    if pid not in masif_non_redundant:
        masif_non_redundant.append(pid)
        masif_batch.append(pid)
        all_pids+=pid+','
    if len(masif_batch)>=limit:
        print()
        all_pids = all_pids.strip(',')
        print(all_pids)
        all_pids=''
        masif_batch=[]

all_pids = all_pids.strip(',')
print(all_pids)

# # Download SAbDAB test set
in_pdb_dir = '../04-2021-SAbDab/data_preparation/00-raw_pdbs/'
def extract_fasta(PDBs, out_fasta):
    for pid in PDBs:
        with open(in_pdb_dir+pid+'.pdb', 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                out_fasta.write('>' + record.id+'\n')
                out_fasta.write(str(record.seq) + '\n')

# We will use only non-overlapped SAbDab fasta, as we will exclude overlapped PDBs anyway.
test_pid_list = [x.strip('\n') for x in open('lists/unique_sabdab.txt', 'r').readlines()]

if not os.path.exists('fasta/sabdab_raw'):
    os.mkdir('fasta/sabdab_raw')
out_mimicry_fasta = open('fasta/sabdab_raw/raw.fasta', 'w')

extract_fasta(test_pid_list, out_mimicry_fasta)