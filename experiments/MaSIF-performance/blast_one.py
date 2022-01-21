import os
import subprocess
import sys

query_fasta_i = sys.argv[1]
mimicry_fasta_dir = sys.argv[2]
out_dir_hits = sys.argv[3]
query_fasta_dir = sys.argv[4]

IDENTITY_THRESHOLD=95

if not os.path.exists(out_dir_hits):
    os.makedirs(out_dir_hits)

#for query_fasta_i in os.listdir(query_fasta_dir):

outformat = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos'

for mimicry_fasta in os.listdir(mimicry_fasta_dir):
    print(mimicry_fasta)
    if '.fasta.' not in mimicry_fasta:
        process = subprocess.Popen(
            ['blastp', '-query', query_fasta_dir+query_fasta_i,
             '-db', mimicry_fasta_dir+mimicry_fasta,
            #  '-task', 'blastp-short',
            #  '-gapopen', '32767',
            # '-gapextend', '32767',
            # '-max_target_seqs','200000000',
             '-outfmt', outformat
             ],
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(stdout)
        print(stderr)
        # out_file = out_dir_hits + mimicry_fasta.strip('.fasta') + '_' + query_fasta_i.strip('.fasta') + '.txt'
        # with open(out_file, 'w') as out:
        #     out.write(stdout.decode("utf-8"))

        # Parse blast hit
        stdout = stdout.decode("utf-8")
        stdout = stdout.strip('\n')

        for line in stdout.split('\n'):
            if '#' not in line:
                print(line)
                fields = line.split('\t')
                identity = float(fields[2])
                if identity>=IDENTITY_THRESHOLD:
                    out_file = out_dir_hits+fields[0].replace(':','_')+'-'+fields[1].split('|')[0] + '.txt'
                    with open(out_file, 'w') as out:
                        out.write(line+'\n')
