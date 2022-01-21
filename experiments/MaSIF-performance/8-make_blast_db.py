import os
import subprocess

fasta_dir = 'fasta/masif_train_fasta/'

curr_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(fasta_dir)

for sp_fasta in os.listdir('.'):
    if '.fasta' in sp_fasta:
        process = subprocess.Popen(['makeblastdb', '-in', sp_fasta, '-dbtype', 'prot'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(stdout)
        print(stderr)
