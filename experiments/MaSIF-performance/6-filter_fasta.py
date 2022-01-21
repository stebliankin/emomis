# Filter out antibody chains from fasta sequences

from Bio import SeqIO
import os
import pdb


def getAntibodyName(pid, chains, species_metadata):
    with open(species_metadata, 'r') as f:
        f.readline()
        for line in f.readlines():
            fields = line.strip('\n').split('\t')
            pid_line = fields[0].upper()
            antibody_chain_line = fields[-1].split('_')[1]
            if pid_line==pid:
                # if pid=='1V7M' and chains=='V':
                #     pdb.set_trace()
                for antibody_chain_i in antibody_chain_line:
                    if antibody_chain_i in chains:
                        return True, pid+ ':'+antibody_chain_i
    return False, ''

def parse_fasta(input_fasta_dir, fasta_name, out_fasta_dir, filter_list):
    in_fasta = input_fasta_dir + fasta_name + '.fasta'
    fasta_sequences = SeqIO.parse(open(in_fasta),'fasta')
    species_metadata = '../04-2021-SAbDab/metadata/sabdab_filtered.tsv'


    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        description = fasta.description

        pid = description.split(':')[0]
        if pid not in filter_list:
            #chains = description.split('|')[1].split('Chain')[1].strip('s').strip(' ').split(',')
            chain = description.split(':')[1]
            antibodyFlag, antibodyName = getAntibodyName(pid, chain, species_metadata)
            if antibodyFlag:
                out_fasta = out_fasta_dir + antibodyName.replace(':','_') + '.fasta'
                with open(out_fasta, 'w') as out_file:
                    print(antibodyName)
                    out_file.write('>'+antibodyName + '\n' + sequence + '\n')
            else:
                print(description)

out_fasta_dir = 'fasta/sabdab_antibodies/'

if not os.path.exists(out_fasta_dir):
    os.mkdir(out_fasta_dir)
input_fasta_dir = 'fasta/sabdab_raw/'

for fasta_file in os.listdir(input_fasta_dir):
    fasta_name = fasta_file.split('.fasta')[0]
    parse_fasta(input_fasta_dir, fasta_name, out_fasta_dir, ())

