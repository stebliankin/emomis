# Discretize each fasta record in MaSIF sequence

import os
import shutil

masif_fasta_dir = 'fasta/masif_train_fasta/'
out_fasta_dir = 'fasta/masif_single_fasta/'

if os.path.exists(out_fasta_dir):
    shutil.rmtree(out_fasta_dir)
os.mkdir(out_fasta_dir)

for fasta in os.listdir(masif_fasta_dir):
    with open(masif_fasta_dir+fasta, 'r') as f:
        for line in f.readlines():
            if '>' in line:
                new_fasta_name = line.split('|')[0].strip('>') + '.fasta'
            with open(out_fasta_dir+new_fasta_name, 'a') as out:
                out.write(line)


