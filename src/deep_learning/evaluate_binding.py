import pandas as pd
import pdb
import numpy as np
import os
from tqdm import tqdm
from blast_search.blast_filter import compute_dssp
from scipy.stats import zscore, norm
import pymesh
from scipy.spatial import cKDTree
from deep_learning.compute_descriptors import compute_descriptors_from_patch_indx, read_learning_object

def eval_score_one(desc1_str, desc2_flip):
    best_score = float("inf")
    for i in range(0, len(desc1_str)):
        for j in range(0, len(desc2_flip)):
            desc_dist = np.sum((desc1_str[i] - desc2_flip[j]) ** 2)
            if desc_dist < best_score:
                best_score = desc_dist
    return best_score

def get_binding_score_one(pid1, pid2, contact1, contact2, config):

    desc1_str_file = config['dirs']['descriptors'] + '{}/{}_{}_str.npy'.format(pid1, pid1, contact1)
    desc2_flip_file = config['dirs']['descriptors'] + '{}/{}_{}_flip.npy'.format(pid2, pid2, contact2)

    if contact1=="None" or contact2=="None" or not os.path.exists(desc1_str_file) or not os.path.exists(desc2_flip_file):
        return float("inf")

    desc1_str = np.load(desc1_str_file, allow_pickle=True)
    desc2_flip = np.load(desc2_flip_file, allow_pickle=True)
    best_score = eval_score_one(desc1_str, desc2_flip)

    return best_score


def evaluate_binding(ppi_list, config):

    print("Evaluating binding with MaSIF...")

    df = pd.read_csv(config['out_files']['TM_align_filtered'], sep='\t')

    result_contact_scores = {"subject-subject": [], "target-subject": [], "target-target": [],
                                 "subject-target": []}

    for i, row in df.iterrows():
        pid_subject = row['PDB_subject'].split('_')[0]
        pid_target = row['PDB_target'].split('_')[0]
        contact_subject_ag_list = row['contact_subject_ag'].split(',')
        contact_subject_ab_list = row['contact_subject_ab'].split(',')
        contact_target_ag_list = row['target_resid_ag'].split(',')
        contact_target_ab_list = row['target_resid_ab'].split(',')

        result_contact_scores_tmp = {"subject-subject": [], "target-subject": [], "target-target": [],
                                 "subject-target": []}

        for i in range(0, len(contact_subject_ag_list)):
            binding_score = get_binding_score_one(pid_subject, pid_subject, contact_subject_ag_list[i], contact_subject_ab_list[i], config)
            result_contact_scores_tmp["subject-subject"].append(binding_score)

            binding_score = get_binding_score_one(pid_target, pid_subject, contact_target_ag_list[i],
                                                  contact_subject_ab_list[i], config)
            result_contact_scores_tmp["target-subject"].append(binding_score)


            binding_score = get_binding_score_one(pid_target, pid_target, contact_target_ag_list[i],
                                                  contact_target_ab_list[i], config)


            result_contact_scores_tmp["target-target"].append(binding_score)

            binding_score = get_binding_score_one(pid_subject, pid_target, contact_subject_ag_list[i],
                                                  contact_target_ab_list[i], config)

            result_contact_scores_tmp["subject-target"].append(binding_score)
        result_contact_scores["subject-subject"].append(np.min(result_contact_scores_tmp["subject-subject"]))
        result_contact_scores["target-subject"].append(np.min(result_contact_scores_tmp["target-subject"]))
       # result_contact_scores["target-target"].append(np.min(result_contact_scores_tmp["target-target"]))
        #result_contact_scores["subject-target"].append(np.min(result_contact_scores_tmp["subject-target"]))

    df["score_subject_subject"] = result_contact_scores["subject-subject"]
    df["score_target_subject"] = result_contact_scores["target-subject"]
 #   df["score_target_target"] = result_contact_scores["target-target"]
 #   df["score_subject_target"] = result_contact_scores["subject-target"]

    df.to_csv(config['out_files']['MaSIF_scores'], sep='\t', index=False)

    return

def check_target_ab(ab_contact):
    all_contacts = ab_contact.split(',')
    for contact in all_contacts:
        if contact!="None":
            return True
    return False

def compute_secondary_structure_from_contact(pid_list, contact_list, config):

    all_sec_str = []

    for i, pid in enumerate(tqdm(pid_list)):
        dssp = compute_dssp(pid, config)
        sec_str_i = ''
        for contact in contact_list[i].split(','):
            ch = contact.split(':')[0]
            resid = int(contact.split(':')[1])
            dssp_entry = dssp[(ch, (' ', resid, ' '))]

            sec_str_i+=dssp_entry[2] + ',' # add information on secondary structure

        sec_str_i=sec_str_i.strip(',')
        all_sec_str.append(sec_str_i)
    return all_sec_str

def compute_secondary_structure_motif(pid_list, contact_list, hit_i_list, query_list, config):
    all_sec_str = []

    for i, pid in enumerate(tqdm(pid_list)):
        dssp = compute_dssp(pid, config)

        hit_i = int(hit_i_list[i].split(',')[0])
        contact_res_i = int(contact_list[i].split(',')[0].split(':')[1])
        ch = contact_list[i].split(',')[0].split(':')[0]
        start_i = contact_res_i - hit_i + 1
        hit_len = len(query_list[i])
        sec_str = ''

        for j in range(0, hit_len):
            resid = start_i + j
            try:
                dssp_entry = dssp[(ch, (' ', resid, ' '))]
                sec_str += dssp_entry[2]
            except KeyError:
                sec_str += '-'
        all_sec_str.append(sec_str)
    return all_sec_str

def compute_sim(str1, str2):
    assert len(str1)==len(str2)
    n_match = 0
    for i in range(0, len(str1)):
        if str1[i] == str2[i]:
            n_match+=1
    return '{}'.format(int((n_match/len(str1))*100))

def get_hit_interval(match, hit_i, contact_i):
    hit_i = int(hit_i.split(',')[0]) - 1
    contact_res_i = int(contact_i.split(',')[0].split(':')[1])
    ch = contact_i.split(',')[0].split(':')[0]
    start_i = contact_res_i - hit_i + 1

    # go left
    curr_match_res = match[hit_i]
    curr_hit_i = hit_i
    while curr_match_res!=' ' and curr_hit_i>0:
        curr_hit_i-=1
        curr_match_res = match[curr_hit_i]

    left_range=start_i+curr_hit_i-1
    if curr_match_res==' ':
        left_range+=1

    # go_right
    curr_match_res = match[hit_i]
    curr_hit_i = hit_i
    while curr_match_res != ' ' and curr_hit_i < len(match)-1:
        curr_hit_i+=1
        curr_match_res = match[curr_hit_i]
    right_range=start_i+curr_hit_i-1
    if curr_match_res==' ':
        right_range-=1

    return "{}:{}-{}".format(ch, left_range, right_range)


def filter_MaSIF(updated_ppi_list_db, config):
    # Filter results by the MaSIF score

    # Read the distribution of negative scores
    curr_path_split = os.path.realpath(__file__).split('/')
    curr_path = '/'.join(curr_path_split[0:len(curr_path_split) - 1])

    masif_scores_native = compute_native_scores(updated_ppi_list_db, config)

    all_neg_scores = [float(x) for x in open(curr_path + '/MaSIF_neg_scores.csv').readlines()]

    def compute_pval(masif_score):
        if masif_score==float('inf'):
            return 1

        tmp_neg_scores = all_neg_scores.copy()
        #tmp_neg_scores = tmp_neg_scores.sort()
        tmp_neg_scores.append(masif_score)
        #tmp_neg_scores = [-x for x in tmp_neg_scores]
#        score_indx = tmp_neg_scores.index(masif_score)

        # compute z-scores
        z_scores = zscore(tmp_neg_scores)
        pvals = norm.sf(z_scores) # one tailed
        return 1-pvals[-1]

    def compute_zscore(masif_score):
        if masif_score == float('inf'):
            return None
        tmp_neg_scores = all_neg_scores.copy()
        tmp_neg_scores.append(masif_score)
        z_scores = zscore(tmp_neg_scores)
        return z_scores[-1]


    mimicry_df = pd.read_csv(config['out_files']['MaSIF_scores'], sep='\t')

    mimicry_df['score_subject_subject'] = mimicry_df['PPI_subject'].apply(lambda x: float('inf') if x not in masif_scores_native else masif_scores_native[x])
    mimicry_df['pval_target_subject'] = mimicry_df['score_target_subject'].apply(lambda row: compute_pval(row))
    mimicry_df['zscore_target_subject'] = mimicry_df['score_target_subject'].apply(lambda row: compute_zscore(row))

    mimicry_df['pval_subject_subject'] = mimicry_df['score_subject_subject'].apply(lambda row: compute_pval(row))

    mimicry_df.to_csv(config['out_files']['MaSIF_scores'], sep='\t', index=False)

    mimicry_df = mimicry_df[(mimicry_df['pval_target_subject']<=config['p_val_threshold']) | (mimicry_df['pval_target_subject']<mimicry_df['pval_subject_subject'])]
    mimicry_df.to_csv(config['out_files']['MaSIF_scores_filtered'], sep='\t', index=False)

  #  mimicry_df[]

    # min_tm_score = config['min_tm_score']
    # max_masif_score = config['max_masif_score']
    #
    # all_mimicry_file = config['dirs']['output'] + 'all_mimicry_results.tsv'
    # df = pd.read_csv(all_mimicry_file, sep='\t')
    # df = df[df['TM_score']>min_tm_score]
    # df = df[df['masif_score']<max_masif_score]
    #
    # df['hit_interval_subject_ag'] = df.apply(lambda row: get_hit_interval(row['match'], row['hit_i'], row['contact_subject_ag']), axis=1)
    # df['hit_interval_target_ag'] = df.apply(lambda row: get_hit_interval(row['match'], row['hit_i'], row['target_resid_ag']), axis=1)
    #
    # print('Computing secondary structure for subjects ...')
    # secondary_structure_subject_ag = compute_secondary_structure_motif(list(df['PDB_subject']), list(df['contact_subject_ag']),
    #                                                                    list(df['hit_i']),  list(df['subject']), config)
    #
    # print('Computing secondary structure for targets ...')
    # secondary_structure_target_ag = compute_secondary_structure_motif(list(df['PDB_target']),
    #                                                                           list(df['target_resid_ag']),
    #                                                                   list(df['hit_i']), list(df['query_target']), config)
    #
    # df['subject_sec_str'] = secondary_structure_subject_ag
    # df['target_sec_str'] = secondary_structure_target_ag
    #
    # df['second_str_sim_percent'] = df.apply(lambda row: compute_sim(row['subject_sec_str'], row['target_sec_str']), axis=1)
    #
    # df.to_csv(config['dirs']['output']+'filtered_mimicry_results.tsv', sep='\t', index=False)
    #
    # if os.path.exists(config['dirs']['output']+'tmp_mimicry_results.tsv'):
    #     os.remove(config['dirs']['output']+'tmp_mimicry_results.tsv')


def compute_native_scores(updated_ppi_list_db, config):

    masif_scores_native = {}
    learning_obj = read_learning_object(config) # read deep learning model

    for ppi in tqdm(updated_ppi_list_db):
        pid, ch1, ch2 = ppi.split('_')
        labels = np.load(config['dirs']['patches']+'{}/{}_{}_sc_labels.npy'.format(pid,pid,ch1), allow_pickle=True)
        # Take the median of the percentile 25 shape complementarity.
        mylabels = labels[0]
        labels = np.median(mylabels, axis=1)

        ply_fn1 = config['dirs']['surface_ply'] + '{}_{}.ply'.format(pid, ch1)
        ply_fn2 = config['dirs']['surface_ply'] + '{}_{}.ply'.format(pid, ch2)
        pos_labels = np.where(labels > config['masif_search']['min_sc_fit'])[0]  # contact points in p1
        K = len(pos_labels)
        if K < config['masif_search']['cache_min_per_pdb']:
            #pdb.set_trace()
           # print(ppi)
            K = config['masif_search']['cache_min_per_pdb']
            pos_labels = (-labels).argsort()[:K]

        l = np.arange(len(pos_labels))  # l is the index of pos label
        np.random.shuffle(l)
        l = l[:K]
        l = pos_labels[l]  # contact points in p1

        v1 = pymesh.load_mesh(ply_fn1).vertices[l]  # coordinates of contact points in p1
        v2 = pymesh.load_mesh(ply_fn2).vertices  # all coordinates in p2

        # For each point in v1, find the closest point in v2.
        kdt = cKDTree(v2)
        dist, r = kdt.query(v1)  # r is the indicies in p2 that are close to v1 (shape of p1)
        # Contact points: those within a cutoff distance.
        sorted_indx = np.argsort(dist)
        dist, r = dist[sorted_indx], r[sorted_indx]

        contact_points = np.where(dist <= config['blast_const']["contact_threshold"])[0]  # which r pass the interface cutoff (ex. [1,2]) (p2 indicies)

        if len(contact_points)<config['masif_search']['cache_min_per_pdb']:
            contact_points = list(range(0,config['masif_search']['cache_min_per_pdb'])) # take the closest points into consideration

        k1 = l[contact_points]  # contact points in p1 (patch indicies)
        k2 = r[contact_points]  # contact points in p2 (patch indicies)
        assert len(k1) == len(k2)
        n_pos = len(k1)

        if n_pos == 0:
            raise ValueError("There is no positive points")

        desc1_str, desc1_flip = compute_descriptors_from_patch_indx(pid, ch1, k1, learning_obj, config)
        desc2_str, desc2_flip = compute_descriptors_from_patch_indx(pid, ch2, k2, learning_obj, config)

        masif_score = eval_score_one(desc1_str, desc2_flip)
        masif_scores_native[ppi] = masif_score
    return masif_scores_native
