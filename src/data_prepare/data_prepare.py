from Bio.PDB import *
from subprocess import Popen, PIPE
from masif.source.input_output.protonate import protonate
import os
import pdb
from scipy.spatial import cKDTree
import numpy as np
import time
from datetime import datetime
from .prepare_utils import map_patch_to_atom, extract_fasta
from utils import read_config, get_date, read_ppi_list, check_if_exists_masif
from .triangulate import triangulate_one
from .compute_patch import compute_patch
import pandas as pd
import os
import subprocess
from tqdm import tqdm
import subprocess

def protonate_pdb(ppi, config):
    """
    downlaod and add hydrogens to PPI
    """
    pid = ppi.split('_')[0]

    # Download pdb
    pdbl = PDBList()
    pdb_filename = pdbl.retrieve_pdb_file(pid, pdir=config['dirs']['raw_pdb'], file_format='pdb')

    # Protonate downloaded file
    protonated_file = config['dirs']['protonated_pdb']+"/"+pid+".pdb"
    protonate(pdb_filename, protonated_file)

def write_list(alist, filename):
    with open(filename, 'w') as f:
        for ppi in alist:
            f.write(ppi+'\n')

def download(ppi_list, config, to_write=None, min_seq_len=None):
    start = time.time()
    print("*** [ {} ] Start Downloading PDBs and extracting FASTA...".format(get_date()))

    processed_ppi = []
    for i in tqdm(range(len(ppi_list))):
        ppi = ppi_list[i]
        if len(ppi.split('_'))==3:
            pid, ch1, ch2 = ppi.split('_')
        else:
            pid, ch2 = ppi.split('_')

        raw_pdb_filename = config['dirs']['protonated_pdb']+"/"+pid+".pdb"
        fasta_file_name = config['dirs']['fasta']+'/{}_{}.fasta'.format(pid, ch2)

        if not os.path.exists(raw_pdb_filename):
            protonate_pdb(ppi, config)
        # else:
        #     print("PDB file {} already exists. Skipping...".format(pid))
        #

        # Extract FASTA
        if not os.path.exists(fasta_file_name):
            extract_fasta(pid, ch2, config, min_seq_len)  # extract fasta only for antigens
        # else:
        #     print("FASTA file for {} already exists. Skipping...".format(pid))

        if os.path.exists(fasta_file_name):
            processed_ppi.append(ppi)

    if to_write is not None:
        with open(to_write, 'w') as out:
            for ppi in processed_ppi:
                out.write(ppi+'\n')

    print("*** [ {} ] Done with downloading PDBs and FASTA extraction.".format(get_date()))
    print("*** [ {} ] Took {:.2f}min.".format(get_date(), (time.time()-start)/60))
    return processed_ppi


def extract_ppi_lists(config, reverse_flag):
    print("*** [ {} ] Extracting list of PPIs...".format(get_date()))

    target_name = config['input']['target_name']
    metadata_file = config['input']['sabdab_summary']
    print("Inputs:")
    print("\tMetadata file: {}".format(metadata_file))
    print("\tTarget name: {}".format(target_name))

    print("Extracted files:")
    out_list_target = config['target_list']
    out_list_db = config['db_list']
    print("\tDatabase PPI list: {}".format(out_list_db))
    print("\tPPI list of the target species: {}".format(out_list_target))

    target_keyword = target_name

    df = pd.read_csv('metadata/sabdab_summary_all.tsv', sep='\t')
    df = df[["pdb", "Hchain", "Lchain", "antigen_chain", "antigen_species"]]

    df["Hchain"] = df["Hchain"].fillna("")
    df["Lchain"] = df["Lchain"].fillna("")

    df = df.dropna().reset_index(drop=True)
    df["antigen_chain"] = df["antigen_chain"].apply(lambda x: x.replace(" | ", ""))

    df["PPI"] = df.apply(lambda x: "{}_{}_{}".format(x.pdb.upper(), x.Hchain + x.Lchain, x.antigen_chain), axis=1)

    df = df[df["Hchain"]!=df["antigen_chain"]]
    df = df[df["Lchain"]!=df["antigen_chain"]]

    df = df.sort_values("antigen_chain", ascending=True)
    # df = df.drop_duplicates('pdb', keep='first')

    print("Number of unique PPIs in SAbDab: {}".format(len(df)))
    print("Number of unique PIDs in SAbDab: {}".format(len(set(df['pdb']))))

    if not reverse_flag:
        df_target = df[df['antigen_species'].str.contains(target_keyword)]
    else:
        df_target = df[~df['antigen_species'].str.contains(target_keyword)]

    df_db = df[~df["pdb"].isin(list(df_target['pdb']))]
    # Exclude unwanted species:
    for keyword in config['input']['db_exclude']:
        print('Excluding "{}" from DB...'.format(keyword))
        if not reverse_flag:
            df_db = df_db[~df_db['antigen_species'].str.contains(keyword)]
        else:
            df_target = df_target[~df_target['antigen_species'].str.contains(keyword)]

    print('Number of unique target complexes: {}'.format(len(df_target['pdb'].unique())))
    unique_ppi_target = df_target['PPI'].unique()

    df_ppi_target = pd.DataFrame({'PPI': unique_ppi_target})
    df_ppi_target.to_csv(out_list_target, header=False, index=False)

    print("Number of unique PDB in database: {}".format(len(df_db['pdb'].unique())))
    unique_ppi_db = df_db['PPI'].unique()
    df_ppi_db = pd.DataFrame({'PPI': unique_ppi_db})
    df_ppi_db.to_csv(out_list_db, header=False, index=False)

    return list(unique_ppi_db), list(unique_ppi_target)

def combine_fasta(ppi_list, config):
    print("\t[ {} ] Combining FASTA files into a single file...".format(get_date()))
    fasta_dir_single = config['dirs']['fasta']
    blast_db_dir = config['dirs']['blast_db']

    # Combine antigen FASTA files into a single large FASTA
    with open(blast_db_dir+'db.fasta', 'w') as out:
        for ppi in ppi_list:
            pid, ch1, ch2 = ppi.split('_')
            with open(fasta_dir_single+'{}_{}.fasta'.format(pid, ch2)) as fasta:
                for line in fasta.readlines():
                    out.write(line)

def build_blast(updated_ppi_list_db,config):
    print("***[ {} ] Building blast database...".format(get_date()))

    # Combine fasta into a single file
    combine_fasta(updated_ppi_list_db, config)

    fasta_dir = config['dirs']['blast_db']
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(fasta_dir)
    process = subprocess.Popen(['makeblastdb', '-in', 'db.fasta', '-dbtype', 'prot'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout)
    print(stderr)
    os.chdir(curr_dir)

def get_heavy_atoms(pdb_struct, ch, resid=None):
    heavy_atoms=[]
    for i, atom in enumerate(pdb_struct.get_atoms()):
        ch_i = atom.get_parent().get_parent().get_id() # get current chain
        if ch_i in ch:
            if resid is not None:
                if atom.parent.get_id()[1] !=resid:
                    continue # skip the residue

            if atom.element!='H' and atom.get_parent().get_full_id()[3][0]==' ':
                # .get_parent().get_full_id()[3][0] - checks if the residue is not a heteroatom (HETATM)
                # reference: https://stackoverflow.com/questions/25718201/remove-heteroatoms-from-pdb
                heavy_atoms.append(atom)

    return heavy_atoms

def get_start_res_fromPDB(pdb_struct, ch):

    start_res_dict = {}
    for one_ch in ch:
        for resid in pdb_struct.get_residues():
            if resid.get_parent().get_id() in one_ch:
                start_res_dict[one_ch] = resid.get_id()[1]
                break
    return start_res_dict


def ab_contact_map_one(ppi, out_contact_map, config):
    # Create a map: ag_ch, ag_resid, ab_ch, ab_resid, dist, ab_start_res, ag_start_res
    #print("***[ {} ] Creating a contact map for antibody-antigens...".format(get_date()))

    pid, ch_ab, ch_ag = ppi.split("_")
    pdb_dir = config['dirs']['raw_pdb']

    contact_threshold = config['blast_const']["contact_threshold"]

    pdb_struct_ag = PDBParser().get_structure('{}_{}'.format(pid, ch_ag), pdb_dir+'pdb{}.ent'.format(pid.lower()))
    pdb_struct_ab = PDBParser().get_structure('{}_{}'.format(pid, ch_ab), pdb_dir+'pdb{}.ent'.format(pid.lower()))

    start_res_dict = get_start_res_fromPDB(pdb_struct_ag, ch_ag+ch_ab)

    heavy_atoms_ag = get_heavy_atoms(pdb_struct_ag, ch_ag)
    heavy_atoms_ab = get_heavy_atoms(pdb_struct_ab, ch_ab)

    ab_tree = cKDTree([x.get_coord() for x in heavy_atoms_ab])

    ag_heavy_atom_coord = [x.get_coord() for x in heavy_atoms_ag]
    if len(ag_heavy_atom_coord)>0:
        dist, idx = ab_tree.query([x.get_coord() for x in heavy_atoms_ag])
    else:
        # In some cases, antigen doesn't have any residues (ex. 6IDG_HL_A)
        dist, idx = [], []


    sorted_indx =np.argsort(dist)

    processed_ag_pais = []
    with open(out_contact_map, 'w') as out:
        out.write("AG_chain,AG_resid,AB_chain,AB_resid,dist,AG_start,AB_start\n")
        for i in range(0, len(dist)):
            i = sorted_indx[i]

            if dist[i]<=contact_threshold:
                ch_ag_i = heavy_atoms_ag[i].get_parent().get_parent().get_id()
                resid_ag_i = heavy_atoms_ag[i].parent.id[1]
                ch_ab_i = heavy_atoms_ab[idx[i]].get_parent().get_parent().get_id()
                resid_ab_i = heavy_atoms_ab[idx[i]].parent.id[1]

                if "{}:{}".format(ch_ag_i, resid_ag_i) not in processed_ag_pais:
                    # Include residus only wth minimum distatnce
                    processed_ag_pais.append("{}:{}".format(ch_ag_i, resid_ag_i))

                    out.write("{},{},{},{},{},{},{}\n".format(ch_ag_i, resid_ag_i, ch_ab_i, resid_ab_i, dist[i], start_res_dict[ch_ag_i], start_res_dict[ch_ab_i]))


def ab_contact_map(ppi_list, config):
    # precompute residue numbers that are in contact with antibody
    print("***[ {} ] Creating a contact map for antibody-antigens...".format(get_date()))

    for ppi in tqdm(ppi_list):
        pid, ch_ab, ch_ag = ppi.split("_")
        out_contact_map = config['dirs']['ab_contact_map'] + pid + ".csv"

        ab_contact_map_one(ppi, out_contact_map, config)
    return


def helper_func_ppi(f, ppi, config):
    ppi_fields = ppi.split('_')
    f(ppi_fields[0], ppi_fields[1], config)
    if len(ppi_fields)==3:
        f(ppi_fields[0], ppi_fields[2], config)

def prepare_masif(updated_ppi_list_db, config):
    start = time.time()

    processed_ppis = []
    #download(args)

    for ppi in tqdm(updated_ppi_list_db):

        if check_if_exists_masif(ppi, config):
            print("{} is already pre-processed. Skipping...".format(ppi))
            continue

        try:
            print(" \t***[ {} ] Triangulating {}...".format(get_date(), ppi))
            helper_func_ppi(triangulate_one, ppi, config)

            print(" \t***[ {} ] Computing patches for {}...".format(get_date(), ppi))
            compute_patch(ppi, config)

            print(" \t***[ {} ] Mapping patch centers to residue IDs {}...".format(get_date(), ppi))
            helper_func_ppi(map_patch_to_atom, ppi, config)

            processed_ppis.append(ppi)
        except:
            print("Couldn't preprocess {}. Skipping...".format(ppi))
            pass
    return processed_ppis

def get_processed_patches(ppi_list, config):
    processed_ppis = []
    patches_dir = config['dirs']['patches']
    for ppi in ppi_list:
        pid, ch1 = ppi.split('_')[0], ppi.split('_')[1]
        if os.path.exists(f"{patches_dir}/{pid}/{pid}_{ch1}_input_feat.npy"):
            processed_ppis.append(ppi)
    return processed_ppis







