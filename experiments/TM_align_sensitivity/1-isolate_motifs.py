
import os
import random
import pdb
import pandas as pd


random.seed(69)

pdb_dir = '../demo-09-2021/data_preparation/00-raw_pdbs/'
contact_map_dir = '../demo-09-2021/data_preparation/04-ab_contact_map/'

out_motif_dir = 'out_motifs/'
list_selected = 'lists/selected.txt'
motif_lengh_array = list(range(3,33))
list_pdbs = 'lists/sabdab.txt'
selected_ppis = [x.strip('\n') for x in open(list_pdbs, 'r').readlines()]

if not os.path.exists(out_motif_dir):
    os.mkdir(out_motif_dir)
def convert_threeAA_oneAA(one_code):
    # convert one amino acid code to three
    aa_dict = {'CYS':'C', 'ASP':'D', 'SER': 'S', 'GLN':'Q', 'LYS':'K',
         'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N',
         'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W',
         'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M', 'UNK':'X', 'XAA':'X'}
    # UNK and XAA denotes unknown amino acids
    return aa_dict[one_code]

def random_motif(contact_file, len_i):
     # get the list of consecutive ranges
    contact_df = pd.read_csv(contact_file)
    if len(contact_df)==0:
        return None, None, None
    random_contact = contact_df.sample()

    random_contact = random_contact.reset_index(drop=True)
    contact_res, contact_ch = random_contact["AG_resid"][0], random_contact["AG_chain"][0]

    random_start = contact_res - random.randint(0, len_i)


    return random_start, random_start+len_i-1, contact_ch

def extract_motif(ppi, len_i, pdb_dir, out_motif_dir_current):

    pid, ch_ab, ch_ag = ppi.split('_')

    pdb_file = pdb_dir + 'pdb{}.ent'.format(pid.lower())
    contact_file = contact_map_dir + "{}.csv".format(pid)

    #randomly pick motif location
    res_start, res_end, ch_ag = random_motif(contact_file, len_i)

    sequence_str=""
    prev_resid_i=-1
    out_text = ''
    n_models = 0 # number of models encountered so far
    if res_start is not None:
        with open(pdb_file, 'r') as f:
            for line in f.readlines():
                if line[:5] == "MODEL":
                    n_models+=1
                    if n_models>1:
                        break # exit the loop if PDB has multiple modules
                if line[:4] == 'ATOM':
                    resid_i = int(line[22:26])
                    ch_i = line[21]
                    resid_i = int(line[22:26])
                    resname_i = line[17:20].strip(' ')
                    if ch_i == ch_ag and resid_i >= res_start and resid_i <= res_end:
                        out_text+=line
                        if resid_i != prev_resid_i:
                            res_letter = convert_threeAA_oneAA(resname_i)
                            sequence_str += res_letter
                        prev_resid_i = resid_i

        if len(sequence_str)==len_i:
            out_motif_file = "{}/{}_{}_{}_{}_len{}_{}.pdb".format(out_motif_dir_current, pid, ch_ag, res_start, res_end, len_i, sequence_str)
            open(out_motif_file, 'w').write(out_text)
    return 1

for len_i in motif_lengh_array:
    print("Isolating motifs of length {}".format(len_i))
    out_motif_dir_current = out_motif_dir+'/len' + str(len_i)
    if not os.path.exists(out_motif_dir_current):
        os.mkdir(out_motif_dir_current)
    random.shuffle(selected_ppis)

    for i in range(0, len(selected_ppis)):
        extract_motif(selected_ppis[i], len_i, pdb_dir, out_motif_dir_current)

