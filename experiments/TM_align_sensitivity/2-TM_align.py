# In this script we will randomly align motifs of various lengths

import os
import pdb
import random
import subprocess

random.seed(69)

motif_lengh_array = list(range(3,33))
out_motif_dir = 'out_motifs/'

if not os.path.exists('out_rmsd'):
    os.mkdir('out_rmsd')

out_rmsd_table = 'out_rmsd/rmsd_distribution.csv'

def parse_TM_align(stdout):
    for line in str(stdout).split('\n'):
        if "RMSD=" in line:
            rmsd = float(line.split("RMSD=")[1].split(',')[0].strip(' ').strip('\t'))
        if 'TM-score=' in line:
            tm_score = float(line.split('TM-score=')[1].split('(if')[0].strip(' ').strip('\t'))
    return rmsd, tm_score


def save_fasta_file(seq_target, seq_subject, fasta_file):
    with open(fasta_file, 'w') as out:
        out.write('>target\n')
        out.write(seq_target+'\n')
        out.write('>subject\n')
        out.write(seq_subject+'\n')

with open(out_rmsd_table, 'w') as out:
    out.write("Len,PDB1,PDB2,RMSD,TM_score\n")
    for len_i in motif_lengh_array:
        motif_dir = out_motif_dir+ '/len' + str(len_i) + '/'
        all_motifs = os.listdir(motif_dir)
        for i, target_rand in enumerate(all_motifs):
            tmp_all_motifs = all_motifs.pop(i)
            subj_rand = random.sample(all_motifs, 1)[0]
            target_seq = target_rand.split('_')[-1].strip('.pdb')
            subj_seq = subj_rand.split('_')[-1].strip('.pdb')

            assert target_rand!=subj_rand

            fasta_file = motif_dir + target_rand.strip('.pdb')+ "_" + subj_rand.strip('.pdb') + '.fasta'
            save_fasta_file(target_seq, subj_seq, fasta_file)

            process = subprocess.Popen(
                ['TMalign', motif_dir + target_rand, motif_dir+subj_rand,
                 '-I', fasta_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            try:
                rmsd, TM_score = parse_TM_align(stdout)
                out.write("{},{},{},{},{}\n".format(len_i, target_rand.strip('.pdb'), subj_rand.strip('.pdb'), rmsd,
                                                    TM_score))
            except:
                print(stderr)

