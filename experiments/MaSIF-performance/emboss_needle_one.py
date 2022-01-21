import os
import subprocess
import sys

blast_hit = sys.argv[1]
sabdab_fasta_dir = sys.argv[2]
masif_fasta_dir = sys.argv[3]
outdir = sys.argv[4]
identity_threshold = int(sys.argv[5])

fasta1 = sabdab_fasta_dir +blast_hit.split('-')[0] + '.fasta'
fasta2 = masif_fasta_dir + blast_hit.split('-')[1].strip('.txt') + '.fasta'

outfile = outdir + '/' + fasta1.split('/')[-1].strip('.fasta') + '-' + fasta2.split('/')[-1].strip('.fasta') + '.txt'

process = subprocess.Popen(
    ['needle', fasta1, fasta2,
    '-gapopen', '10',
    '-gapextend', '0.5',
    '-outfile', outfile
     ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

print(stdout)
print(stderr)

# Parse the output file
with open(outfile, 'r') as f:
    for line in f.readlines():
        if '# Identity:' in line:
            print(line)
            identity = float(line.split('(')[1].strip('%)\n'))
            if identity<identity_threshold:
                os.remove(outfile)
