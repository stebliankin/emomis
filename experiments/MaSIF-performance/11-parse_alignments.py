# Parse pairwise alignment files and extract unique PDB IDs from SAbDab which were homologous to the MASiF training set

import os

emboss_hits_dir = 'emboss_hits/'

homologs_list_file = 'lists/homologs_masif_sabdab.txt'

homologs_list = []

for hit_file in os.listdir(emboss_hits_dir):
    pid = hit_file.split('_')[0]
    if pid not in homologs_list:
        homologs_list.append(pid)

with open(homologs_list_file, 'w') as out:
    for pid in homologs_list:
        out.write(pid+'\n')

