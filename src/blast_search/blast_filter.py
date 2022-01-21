from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
import os
import sys
import pdb
import shutil
import pdb
from utils import get_date
from tqdm import tqdm
from masif.source.input_output.extractPDB import extractPDB
from Bio.PDB import *
from scipy.spatial import cKDTree
import multiprocessing
from functools import partial

def blast_empty(blast_hit):
    blast_records = NCBIXML.parse(open(blast_hit))

    # If the generator 'blast_records' is empty, skip the rest of the code:
    try:
        for blast_records in blast_records:
            pass
    except ValueError:
        print("The file {} is empty.".format(blast_hit))
        return True
    return False

def from_list_to_str(alist):
    out_str=""
    for el in alist:
        out_str+=str(el)+','
    return out_str.strip(",")


# def save_hit(pid_target, ppi_subject, hits_filtered_dir, hsps, resid_subject_f3_ag_list, resid_subject_f3_ab_list, resid_query_f3_list, hit_f3_i, filter1_flag, filter2_flag, filter3_flag):
#
#     # filter1_passed, filter2_passed, filter3_passed, filter3_hit_i, filter3_resid_subject, filter3_resid_query
#     def from_list_to_str(alist):
#         out_str=""
#         for el in alist:
#             out_str+=str(el)+','
#         return out_str.strip(",")
#
#     with open(hits_filtered_dir+pid_target+'.txt', 'a') as out:
#         out.write(pid_target+"\t"+ ppi_subject+"\t"+
#                   hsps.query+"\t"+hsps.match+"\t"+hsps.sbjct +"\t"+ str(hsps.expect) + "\t" +
#                   str(filter1_flag) +"\t"+ str(filter2_flag) +"\t"+ str(filter3_flag) +"\t"+
#                  from_list_to_str(hit_f3_i) +"\t"+
#                   from_list_to_str(resid_subject_f3_ag_list) +"\t"  +
#                   from_list_to_str(resid_subject_f3_ab_list) +"\t"+
#                   from_list_to_str(resid_query_f3_list) + "\n"
#                     )


def extract_all_antigen_chains(pid, config):
    summary_df = pd.read_csv(config['input']['sabdab_summary'], sep='\t')
    summary_df = summary_df[summary_df['pdb']==pid.lower()]
    antigen_chains = list(summary_df['antigen_chain'].unique())
    out_ch = ""
    for ch in antigen_chains:
        if ch:
            try:
                out_ch+=ch.replace(' ','').replace('|','') # edge case: 'C | A'
            except AttributeError:
                print("WARNING: Couldn't process chain {} for {}. All chains: {}".format(ch, pid, antigen_chains))
    return out_ch

def compute_dssp(pid, config):
    # Compute DSSP values as described in https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html

    raw_pdb_dir = config['dirs']['raw_pdb']
    tmp_dir = config['dirs']['tmp']

    # extract PDB without antibody
    ch = extract_all_antigen_chains(pid, config)

    extractPDB(raw_pdb_dir+"pdb{}.ent".format(pid.lower()), tmp_dir+'{}_{}.pdb'.format(pid,ch), chain_ids=ch)

    parser = PDBParser(QUIET=1)
    struct = parser.get_structure(pid, tmp_dir+'{}_{}.pdb'.format(pid,ch))

    model = struct[0]
    try:
        dssp = DSSP(model, tmp_dir+'{}_{}.pdb'.format(pid,ch), dssp='mkdssp') # example of a key: ('A', (' ', 1147, ' '))
    except Exception:
        dssp = None
    #os.remove(tmp_dir+'{}_{}.pdb'.format(pid,ch))
    return dssp

def get_ab_chain(pid, ch, list_path):
    all_ppi = [x.strip('\n') for x in open(list_path)]

    for ppi in all_ppi:
        fields = ppi.split('_')
        if pid==fields[0] and ch in fields[2]:
            return fields[1]


#######################################################################################################################
# Blast filters
#######################################################################################################################
def filter1_neigh_exact(res_i, query_str, subject_str):
    # are neighbors of the residue res_i has exact match
    # if res_i==0 or res_i==len(query_str)-1:
    #     return False
    if (query_str[res_i-1]==subject_str[res_i-1]) and \
            (query_str[res_i] == subject_str[res_i]) and \
            (query_str[res_i + 1] == subject_str[res_i + 1]):
        return True
    else:
        return False

def get_ab_contact(pid, ch, resid, config):
    contact_ab_map_dir = config['dirs']['ab_contact_map']
    contact_df = pd.read_csv(contact_ab_map_dir+pid+'.csv')
    if not contact_df.empty:
        contact_df = contact_df[contact_df['AG_chain']==ch]
    if not contact_df.empty:
        contact_df = contact_df[contact_df["AG_resid"]==resid]

    if not contact_df.empty:
        ab_contact = contact_df.iloc[0]["AB_chain"] + ":" + str(contact_df.iloc[0]["AB_resid"])
        return ab_contact
    else:
        return "None"

def convert_threeAA_oneAA(one_code):
    # convert one amino acid code to three
    aa_dict = {'CYS':'C', 'ASP':'D', 'SER': 'S', 'GLN':'Q', 'LYS':'K',
         'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N',
         'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W',
         'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M'}
    return aa_dict[one_code]


# def verify_start_res(pid, ch, curr_res, res_name, config):
#     # Verify that inferred residue ID is accurate
#     # Use for debugging purposes only
#     pdb_path =  config['dirs']['raw_pdb'] + 'pdb{}.ent'.format(pid.lower())
#
#     prev_resid_i = 0
#     prev_ch = ""
#     n_gaps = 0
#     with open(pdb_path, 'r') as f:
#         for line in f.readlines():
#             if line[0:4]=='ATOM':
#
#                 insertion_code = line[26]
#                 resid_i = int(line[22:26])
#                 chain_id_i = line[21]
#
#                 if chain_id_i == ch and resid_i==curr_res:
#                     resname_i =line[17:20].strip(' ')
#                     resname_i_one_letter = convert_threeAA_oneAA(resname_i)
#                     if resname_i_one_letter==res_name:
#                         return True
#                     else:
#                         return False
#     return False


def get_curr_residue(pid, ch, res_fasta,  config):
    df = pd.read_csv(config['dirs']['fasta_maps']+f"{pid}_{ch}.csv")
    df = df[df['indx']==res_fasta-1]
    return int(df.iloc[0]["resid"])

def filter2_surf_acc(res_i, pid, ch, dssp, hsps, config):
    surface_acc_threshold = config["blast_const"]["surface_acc_threshold"]

    query_start = int(hsps.query_start)
    #start_resid_pdb = compute_start_resid_from_dssp(dssp)

    # n_insertions = count_insertions(pid, ch, query_start+res_i, config, hsps, res_i)
    # resid_curr = query_start+res_i + start_resid_pdb - 1 - n_insertions
    try:
        resid_curr = get_curr_residue(pid, ch, query_start+res_i,  config)
    except ValueError:
        return False, None, None

    # verified_flag = verify_start_res(pid, ch, resid_curr, hsps.query[res_i], config)
    # if not verified_flag:
    #     pdb.set_trace()

    access_flag = True
    for i in range(0, 3):
        try:
            dssp_entry = dssp[(ch, (' ', resid_curr+i-1, ' '))]
            surface_acc = dssp_entry[3]
            if surface_acc<surface_acc_threshold:
                access_flag=False
        except KeyError:
            access_flag=False

    if access_flag:
        return True, ch+":"+str(resid_curr), get_ab_contact(pid, ch, resid_curr, config)
    else:
        return False, None, None

# def count_insertions(pid, chain, seq_res_curr_i, config, hsps=None, res_i=0):
#
#     # Sometimes residues with subscript "A" increment the resid count in sequence, but stays the same in PDB file.
#     # For example, presence of resid 169A 6PZZ:C increment the residue number in 6PZZ_C.fasta
#     # In this function we will count number of residues with insertion code and determine the current residue number in the FASTA hit.
#
#     # Missing residues are included in FASTA as "X", but not included in the PDB
#
#     # pid           - pdb ID of the complex of interest
#     # chain         - chain letter
#     # seq_res_curr_i - position of the current residue in FASTA sequence
#
#     # Returns: number of insertions
#
#     n_isertions = 0
#     unique_aa = set()
#     unique_insertions = []
#
#     pdb_path =  config['dirs']['raw_pdb'] + 'pdb{}.ent'.format(pid.lower())
#
#     prev_resid_i = 0
#     prev_ch = ""
#     n_gaps = 0
#     with open(pdb_path, 'r') as f:
#         for line in f.readlines():
#             if line[0:4]=='ATOM':
#
#                 insertion_code = line[26]
#                 resid_i = int(line[22:26])
#                 chain_id_i = line[21]
#
#                 if chain_id_i == chain:
#                     if resid_i - prev_resid_i>1 and prev_ch==chain_id_i:
#                         n_gaps += resid_i - prev_resid_i - 1
#
#                     if insertion_code != ' ':
#                         if str(resid_i)+insertion_code not in unique_insertions:
#                             n_isertions += 1
#                             unique_insertions.append(str(resid_i)+insertion_code)
#                     unique_aa.add(resid_i)
#
#                 prev_resid_i = resid_i
#                 prev_ch = chain_id_i
#
#             if len(unique_aa)+n_gaps>=seq_res_curr_i:
#                 return n_isertions
#     return None

def filter3_ab_proximity(res_i, subject_pid, subject_chain, hsps, config):
    # Return True if residue "res_i" of "subject_pid" is in close proximity with antibody
    contact_ab_map_dir = config['dirs']['ab_contact_map']

    subj_start = int(hsps.sbjct_start)

    contact_df = pd.read_csv(contact_ab_map_dir+subject_pid+'.csv')

    if not contact_df.empty:
        contact_df = contact_df[contact_df['AG_chain']==subject_chain]
    if not contact_df.empty:

        # n_insertions = count_insertions(subject_pid, subject_chain, subj_start+res_i, config, hsps)

        #curr_res_i_real = res_i + subj_start + contact_df.iloc[0]["AG_start"] - 1 - n_insertions
        curr_res_i_real = get_curr_residue(subject_pid, subject_chain, subj_start+res_i,  config)

        # verified_flag = verify_start_res(subject_pid, subject_chain, curr_res_i_real, hsps.sbjct[res_i], config)
        # if not verified_flag:
        #     pdb.set_trace()
        # if n_insertions > 0:
        #     print("{} INSERTIONS FOR {} ch {}; CURR_RES: {}; QUERY: {}; QUERY_i:{}; RES_NAME: {} !!!!!!!!!!!!!!!!!".format(n_insertions, subject_pid, subject_chain, curr_res_i_real, hsps.sbjct, res_i, hsps.sbjct[res_i]))

        contact_df = contact_df[contact_df["AG_resid"]==curr_res_i_real]

    if not contact_df.empty:
        ab_contact = contact_df.iloc[0]["AB_chain"] + ":" + str(contact_df.iloc[0]["AB_resid"])
        ag_contact = subject_chain + ":" + str(curr_res_i_real)
        return True, ag_contact, ab_contact
    else:
        return False, None, None

#######################################################################################################################

# def compute_start_resid_from_dssp(dssp):
#     start_resid = dssp.keys()[0][1][1]
#     return start_resid # example of a dssp key: ('A', (' ', 14, ' '))

def filter_alignment(hsps, target_pid, target_ch, subject_pid, subject_chain, dssp_target, config):

    def update_flag(filter_flag, filter_flag_i):
        if filter_flag_i:
            return True
        else:
            return filter_flag

    hits_filtered_dir = config['dirs']['blast_hits_filtered']


    #save_hit(raw_hits_subjects_dir, pid_target, subject_pid, hsps.sbjct, hsps.sbjct_start, hsps.sbjct_end)
    filter1_flag, filter2_flag, filter3_flag = False, False, False

    resid_subject_f3_ag_list, resid_subject_f3_ab_list, contact_ab_query_f3_list, contact_ag_query_f3_list, hit_f3_list = [], [], [], [], []

    for res_i in range(1, len(hsps.query)-1):
        filter1_flag_i = filter1_neigh_exact(res_i, hsps.query, hsps.sbjct)
        filter1_flag = update_flag(filter1_flag, filter1_flag_i)
        if filter1_flag_i:
            filter2_flag_i, contact_ag_query_i, contact_ab_query_i = filter2_surf_acc(res_i, target_pid, target_ch, dssp_target, hsps, config)
            filter2_flag = update_flag(filter2_flag, filter2_flag_i)
            if filter2_flag_i:
                filter3_flag_i, contact_subject_ag_i, contact_subject_ab_i = filter3_ab_proximity(res_i, subject_pid, subject_chain, hsps, config)
                filter3_flag = update_flag(filter3_flag, filter3_flag_i)
                if filter3_flag_i:
                    resid_subject_f3_ag_list.append(contact_subject_ag_i)
                    resid_subject_f3_ab_list.append(contact_subject_ab_i)
                    contact_ag_query_f3_list.append(contact_ag_query_i)
                    contact_ab_query_f3_list.append(contact_ab_query_i)

                    hit_f3_list.append(res_i+1)


    out_str = target_pid + "\t" + subject_pid + "_" + subject_chain + "\t" + \
                hsps.query + "\t" + hsps.match + "\t" + hsps.sbjct + "\t" + str(hsps.expect) + "\t" + \
                str(filter1_flag) + "\t" + str(filter2_flag) + "\t" + str(filter3_flag) + "\t" + \
                from_list_to_str(hit_f3_list) + "\t" + \
                from_list_to_str(resid_subject_f3_ag_list) + "\t" + \
                from_list_to_str(resid_subject_f3_ab_list) + "\t" + \
                from_list_to_str(contact_ag_query_f3_list) + "\t" + \
                from_list_to_str(contact_ab_query_f3_list) + "\n"
    return out_str



    #
    # save_hit(pid_target=pid_target,
    #          ppi_subject=subject_pid+"_" + subject_chain,
    #          hits_filtered_dir=hits_filtered_dir,
    #          hsps = hsps,
    #          resid_subject_f3_ag_list=resid_subject_f3_ag_list,
    #          resid_subject_f3_ab_list=resid_subject_f3_ab_list,
    #          resid_query_f3_list=resid_query_f3_list,
    #          hit_f3_i=hit_f3_list,
    #          filter1_flag=filter1_flag,
    #          filter2_flag=filter2_flag,
    #          filter3_flag=filter3_flag)


def filter_one(target_ppi, config):
    target_pid, ch1, ch2 = target_ppi.split('_')
    blast_hits_dir = config['dirs']['raw_blast_hits']
    blast_hit = blast_hits_dir +  "{}_{}_blast.xml".format(target_pid, ch2)

    out_hits_file = config['dirs']['blast_hits_filtered'] + target_ppi + '.txt'
    if os.path.exists(out_hits_file):
        return 1


    dssp_target = compute_dssp(target_pid, config)
    if not dssp_target: # if failed to compute DSSP
        return None

    if blast_empty(blast_hit):
        # Exit filtering function if the blast hit doesn't have any entries
        return None
    # Current blast_records reached the end. Read the generator again:
    blast_records = NCBIXML.parse(open(blast_hit))

    out_hits = set()

    # Iterate trough each query in blast hit
    for blast_record in blast_records:
        target_ch = blast_record.query.split(":")[1]

        alignments = blast_record.alignments
        for alignment in alignments:
            subject_pid = alignment.hit_def.split(':')[0]
            subject_chain = alignment.hit_def.split(':')[1]
            for hsps in alignment.hsps:
                one_hit = filter_alignment(hsps, target_pid, target_ch, subject_pid, subject_chain, dssp_target, config)
                out_hits.add(one_hit)

    with open(out_hits_file, 'w') as out:
        out.write("PDB_target" + "\t" + "PDB_subject" + "\t" + \
                "query_target" + "\t" + "match" + "\t" + "subject" + "\t" + "eval" + "\t" + \
                "filter1_flag" + "\t" + "filter2_flag" + "\t" + "filter3_flag" + "\t" + \
                "hit_i" + "\t" + \
                "contact_subject_ag" + "\t" + \
                  "contact_subject_ab" + "\t" + \
                  "target_resid_ag" + "\t" + \
                  "target_resid_ab" + "\n")
        for hit in out_hits:
            out.write(hit)
    return 1



def blast_filter(updated_ppi_list_target, config):
    print("***[ {} ] Filtering blast hits for {} target sequences...".format(get_date(), len(updated_ppi_list_target)))
    pool = multiprocessing.Pool(processes=1 if not config['n_threads'] else config['n_threads'])
    pool.map(partial(filter_one, config=config), updated_ppi_list_target)
    pool.close()
    pool.join()
    # for ppi in tqdm(updated_ppi_list_target):
    #     pid, ch1, ch2 = ppi.split('_')
    #     filter_one(ppi, config)

def report_hits(updated_ppi_list_target, ppi_list_db, config):
    df_all = pd.DataFrame([])
    for ppi in updated_ppi_list_target:
        filtered_hit_file = config['dirs']['blast_hits_filtered'] + ppi + '.txt'
        if not os.path.exists(filtered_hit_file):
            continue
        tmp_df = pd.read_csv(filtered_hit_file, sep='\t')
        if len(df_all)==0:
            df_all=tmp_df
        else:
            df_all = df_all.append(tmp_df)
    df_all.to_csv(config['out_files']['raw_blast_hits'], sep="\t", index=False)

    print("Number of unique raw hits: {}".format(len(df_all[["PDB_subject", "subject"]].dropna().drop_duplicates())))

    for i in range(1,4):
        filter_df = df_all[df_all["filter{}_flag".format(i)]==True]
        print("Number of unique hits passed filter {}: {}".format(i, len(filter_df[["PDB_subject", "subject"]].dropna().drop_duplicates())))

    def check_target_ab(ab_contact):
        all_contacts = ab_contact.split(',')
        for contact in all_contacts:
            if contact!="None":
                return True
        return False

    filter_df["cross_ab_avail"] = filter_df['target_resid_ab'].apply(lambda x: check_target_ab(x))
    filter_df.to_csv(config['out_files']['filtered_blast_hits'], index=False, sep="\t")

    subject_pids = filter_df['PDB_subject'].unique()
    updated_ppi_db = []
    for pid_ch in subject_pids:
        pid, ch = pid_ch.split('_')
        for ppi in ppi_list_db:
            fields = ppi.split('_')
            if pid==fields[0] and ch==fields[2]:
                updated_ppi_db.append(ppi)
    print("Number of unique subject PDBs: {}".format(len(updated_ppi_db)))

    filter_df["hit_length"] = filter_df["query_target"].apply(lambda x: len(x))
    filter_df["PDB_subject"] = filter_df['PDB_subject'].apply(lambda x: x.split('_')[0].lower())

    metadata_df = pd.read_csv(config['input']['sabdab_summary'], sep='\t')

    tmp_df = filter_df.merge(metadata_df, how="left", left_on='PDB_subject', right_on='pdb')

    tmp_df = tmp_df.drop_duplicates(["subject", "antigen_name", "antigen_species"])
    tmp_df.to_csv(config['dirs']['blast_hits_filtered']+"tmp.txt", index=False, sep="\t")
    return updated_ppi_db









