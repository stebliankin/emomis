import os
import subprocess
import sys
from tqdm import tqdm
from utils import get_date

def blast_one(ppi, config):
    pid, ch1, ch2 = ppi.split('_')

    query_fasta_dir = config['dirs']['fasta']
    query_fasta = query_fasta_dir + '{}_{}.fasta'.format(pid, ch2)

    db_path = config['dirs']['blast_db'] + 'db.fasta'

    process = subprocess.Popen(
        ['blastp', '-query', query_fasta,
         '-db', db_path,
         '-task', 'blastp-short',
         '-evalue', '999999',
         '-gapopen', '32767',
         '-gapextend', '32767',
         '-outfmt', '5'
         ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    #print(stderr)
   # print(stdout)

    out_file = config['dirs']['raw_blast_hits'] + "{}_{}_blast.xml".format(pid, ch2)
    with open(out_file, 'w') as out:
        out.write(stdout.decode("utf-8"))
    return

def run_blast(updated_ppi_list_target, config):
    print("*** [ {} ] Running blast for {} target sequences...".format(get_date(), len(updated_ppi_list_target)))
    for ppi in tqdm(updated_ppi_list_target):
        blast_one(ppi, config)