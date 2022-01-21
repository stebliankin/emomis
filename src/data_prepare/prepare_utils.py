import numpy as np
import pandas as pd
from Bio.PDB import *
from scipy.spatial import cKDTree
import pdb
from Bio import Entrez,SeqIO, BiopythonWarning
import warnings
import os
# Ignore biopython warnings

warnings.simplefilter('ignore', BiopythonWarning)


def get_start_res(resid, chain_id):
    chain_curr = ''
    start_res = []
    start_res_curr = resid[0]
    for i, res_i in enumerate(resid):
        if chain_id[i]!=chain_curr:
            start_res_curr=res_i
            chain_curr=chain_id[i]
        start_res.append(start_res_curr)
    return np.array(start_res)

def map_patch_to_atom(pid, ch, config):
    pdb_id = pid
    chain_name = ch

    mappings_dir = config['dirs']['map_patch']
    patch_dir = config['dirs']['patches'] + pdb_id + '/'
    pdb_chain_dir = config['dirs']['chains_pdb']

    out_table = mappings_dir + "/" + pdb_id + "_" + chain_name  + "_map.csv"

    # Read coordinates of
    x_coord = np.load(patch_dir+"/{}_{}_X.npy".format(pdb_id, chain_name))
    y_coord = np.load(patch_dir + "/{}_{}_Y.npy".format(pdb_id, chain_name))
    z_coord = np.load(patch_dir + "/{}_{}_Z.npy".format(pdb_id, chain_name))
    patch_coord = np.column_stack((x_coord,y_coord,z_coord))

    # Read interface
    iface_labels = np.load(patch_dir+"/{}_{}_iface_labels.npy".format(pdb_id, chain_name))

    # Read PDB structure
    pdb_path = "{}/{}_{}.pdb".format(pdb_chain_dir, pdb_id, chain_name)

    parser = PDBParser()
    pdb_struct = parser.get_structure('{}_{}'.format(pdb_id, chain_name), pdb_path)

    ## Get heavy atoms
    heavy_atoms=[]
    heavy_orig_map = {}
    k=0
    for i, atom in enumerate(pdb_struct.get_atoms()):
        if atom.element!='H':
            heavy_orig_map[k]=i #map heavy atom index to original pdb index
            heavy_atoms.append(atom)
            k+=1


    atom_coord = np.array([list(atom.get_coord()) for atom in heavy_atoms])
    atom_names = np.array([atom.get_id() for atom in heavy_atoms])
    residue_id = np.array([atom.parent.id[1] for atom in heavy_atoms])
    residue_name = np.array([atom.parent.resname for atom in heavy_atoms])
    chain_id = np.array([atom.get_parent().get_parent().get_id() for atom in heavy_atoms])

    # get start residue
    start_res = get_start_res(residue_id, chain_id)

    #Create KD Tree
   # patch_tree = cKDTree(patch_coord)
    pdb_tree = cKDTree(atom_coord)

    dist, idx = pdb_tree.query(patch_coord) #idx is the index of pdb heavy atoms that close to every patch from [0 to N patches]
    result_pdb_idx=[]
    for i in idx:
        result_pdb_idx.append(heavy_orig_map[i])
    result_pdb_idx = np.array(result_pdb_idx) #index in original pdb
    #Combine everything to a table:
    df = pd.DataFrame({"patch_ind":range(0, len(result_pdb_idx)),
                       "atom_ind":result_pdb_idx,
                       "res_ind": residue_id[idx],
                       "atom_name":atom_names[idx],
                       "residue_name":residue_name[idx],
                       "chain_id":chain_id[idx],
                       "dist": dist,
                       "iface_label":iface_labels,
                       "start_res": start_res[idx]
                       })


    df.to_csv(out_table, index=False)

# def extract_fasta(pid, ch, config):
#     in_pdb_dir = config['dirs']['protonated_pdb']
#     out_fasta_dir = config['dirs']['fasta']
#
#     out_fasta = out_fasta_dir + '{}_{}.fasta'.format(pid, ch)
#
#     isempty = True
#     with open(out_fasta, 'w') as out:
#         with open(in_pdb_dir + pid + '.pdb', 'r') as pdb_file:
#             try:
#                 for record in SeqIO.parse(pdb_file, 'pdb-atom'):
#                     curr_ch = record.id.split(':')[1]
#                     if curr_ch in ch:
#                         out.write('>' + record.id + '\n')
#                         out.write(str(record.seq) + '\n')
#                         isempty=False
#             except KeyError:
#                 print("Can't extract FASTA from {}".format(pdb_file))
#     if isempty:
#         os.remove(out_fasta)

def convert_threeAA_oneAA(one_code):
    # convert one amino acid code to three
    aa_dict = {'CYS':'C', 'ASP':'D', 'SER': 'S', 'GLN':'Q', 'LYS':'K',
         'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N',
         'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W',
         'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M'}
    return aa_dict[one_code]


def extract_fasta(pid, ch, config, min_seq_len):
    out_fasta_dir = config['dirs']['fasta']

    out_fasta = out_fasta_dir + '{}_{}.fasta'.format(pid, ch)
    out_fasta_map_dir = config['dirs']['fasta_maps']

    pdb_path = config['dirs']['raw_pdb'] + 'pdb{}.ent'.format(pid.lower())

    prev_resid_i = 0
    prev_ch = ""
    n_gaps = 0

    current_fasta_sequence = []
    resid_i_list = []
    resname_i_list = []
    insertion_code_list = []
    chains_list = []

    if not os.path.exists(pdb_path):
        return None

    with open(pdb_path, 'r') as f:
        for line in f.readlines():
            if line[0:4]=='ATOM':
                insertion_code = line[26]
                resid_i = int(line[22:26])
                chain_id_i = line[21]
                if prev_ch=="":
                    prev_ch=chain_id_i

                if chain_id_i in ch:
                    resname_i = line[17:20].strip(' ')
                    try:
                        resname_i_one_letter = convert_threeAA_oneAA(resname_i)
                    except:
                        prev_resid_i = resid_i
                        prev_ch = chain_id_i
                        continue

                    if resid_i!=prev_resid_i:

                        if resid_i - prev_resid_i > 1 and prev_ch == chain_id_i:
                            n_gaps = resid_i - prev_resid_i - 1
                            for i in range(0,n_gaps):
                                current_fasta_sequence.append('X')
                                resid_i_list.append('')
                                resname_i_list.append('')
                                insertion_code_list.append('')
                                chains_list.append(chain_id_i)

                        current_fasta_sequence.append(resname_i_one_letter)
                        resid_i_list.append(resid_i)
                        resname_i_list.append(resname_i)
                        insertion_code_list.append(insertion_code)
                        chains_list.append(chain_id_i)

                prev_resid_i = resid_i
                prev_ch = chain_id_i

    # Create pandas dataframe for extracted sequences:
    df = pd.DataFrame({"resname1":current_fasta_sequence, "resid":resid_i_list, "resname3":resname_i_list,
                       "insertion_codes": insertion_code_list, "chain":chains_list})

    if min_seq_len is not None:
        if len(df)<min_seq_len:
            return None

    if not df.empty:
        with open(out_fasta, 'w') as out:
            unique_chains = df['chain'].unique()
            for ch in unique_chains:
                str_fasta = ""
                df_ch = df[df['chain'] == ch].reset_index(drop=True)
                df_ch.to_csv(out_fasta_map_dir + f'{pid}_{ch}.csv', index_label='indx')
                for index, row in df_ch.iterrows():
                    str_fasta+=row['resname1']
                out.write(f'>{pid}:{ch}\n')
                out.write(str_fasta+'\n')
    return None

