import pandas as pd
import pdb
from masif.source.masif_modules.MaSIF_ppi_search import MaSIF_ppi_search
from masif.source.masif_modules.train_ppi_search import compute_val_test_desc
import numpy as np
from tqdm import tqdm
import os

def extract_contact(ppi, list_type, config):

    def extract_contacts_from_df(df, contact_col):
        contact_list = []
        for contact_i in list(df[contact_col]):
            contact_list_i = contact_i.split(',')
            for contact_j in contact_list_i:
                if contact_j not in contact_list:
                    contact_list.append(contact_j)
        return contact_list

    if list_type=="target":
        contact_ag_col = "target_resid_ag"
        contact_ab_col = "target_resid_ab"
    else:
        contact_ag_col = "contact_subject_ag"
        contact_ab_col = "contact_subject_ab"

    blast_hits_file = config['out_files']['filtered_blast_hits']
    blast_hits_df = pd.read_csv(blast_hits_file, sep="\t")
    blast_hits_df["PDB_{}".format(list_type)] = blast_hits_df["PDB_{}".format(list_type)].apply(lambda x: x.split("_")[0])
    blast_hits_df = blast_hits_df[blast_hits_df["PDB_{}".format(list_type)]==ppi.split('_')[0]]

    contacts_ag = extract_contacts_from_df(blast_hits_df, contact_ag_col)
    contacts_ab = extract_contacts_from_df(blast_hits_df, contact_ab_col)

    return contacts_ag, contacts_ab


def get_patch_index_from_contact(pid, ch, contact_i, config):
    ch_i, resid_i = contact_i.split(':')

    contact_map_file = config['dirs']['map_patch'] + "{}_{}_map.csv".format(pid, ch)

    df = pd.read_csv(contact_map_file)
    df = df[df["chain_id"]==ch_i]
    df = df[df["res_ind"]==int(resid_i)]
    patch_ind = df["patch_ind"].unique()
    return list(patch_ind)


def read_learning_object(config):
    #### Read learning object
    learning_obj = MaSIF_ppi_search(
        12.0,  # radius of neural network
        n_thetas=16,
        n_rhos=5,
        n_rotations=16,
        idx_gpu="/gpu:0",
        feat_mask=[1.0] * 5,
    )
    learning_obj.saver.restore(learning_obj.session, config['dirs']['masif_model'] + "model")
    return learning_obj

def compute_descriptors_from_patch_indx(pid, ch, all_patches_indx, learning_obj, config):


    # Read input features
    in_dir = config['dirs']['patches'] + pid + '/'

    rho_wrt_center = np.load(in_dir + f"{pid}_{ch}" + "_rho_wrt_center.npy")[all_patches_indx]
    theta_wrt_center = np.load(in_dir + f"{pid}_{ch}" + "_theta_wrt_center.npy")[all_patches_indx]
    input_feat = np.load(in_dir + f"{pid}_{ch}" + "_input_feat.npy")[all_patches_indx]
    mask = np.load(in_dir + f"{pid}_{ch}" + "_mask.npy")[all_patches_indx]
    indx_tmp = np.array(range(len(rho_wrt_center)))

    desc1_str = compute_val_test_desc(learning_obj, indx_tmp, rho_wrt_center, theta_wrt_center, input_feat, mask, batch_size=1000, flip=False)
    desc1_flip = compute_val_test_desc(
        learning_obj,
        indx_tmp,
        rho_wrt_center,
        theta_wrt_center,
        input_feat,
        mask,
        batch_size=1000,
        flip=True
    )



    return desc1_str, desc1_flip

def save_descriptors(pid, contact_i, desc_str, desc_flip, all_patches_indx, config):
    descriptors_dir = config['dirs']['descriptors'] + pid + '/'

    if not os.path.exists(descriptors_dir):
        os.mkdir(descriptors_dir)

    np.save(descriptors_dir + f"{pid}_{contact_i}_str.npy", desc_str)
    np.save(descriptors_dir + f"{pid}_{contact_i}_flip.npy", desc_flip)
    np.save(descriptors_dir + f"{pid}_{contact_i}_indx.npy", np.array(all_patches_indx))
    return None

def compute_descriptor_from_contact(pid, ch, contacts_list, learning_obj, config):
    descriptors_dir = config['dirs']['descriptors'] + pid + '/'

    for contact_i in contacts_list:
        if os.path.exists(descriptors_dir + f"{pid}_{contact_i}_str.npy") or contact_i=="None":
            continue
        all_patches_indx = get_patch_index_from_contact(pid, ch, contact_i, config)
        if len(all_patches_indx)>0:
            desc_str, desc_flip = compute_descriptors_from_patch_indx(pid, ch, all_patches_indx, learning_obj, config)
            save_descriptors(pid, contact_i, desc_str, desc_flip, all_patches_indx, config)

    return None

def compute_descriptors(ppi_list, config, list_type):
    if list_type not in ("target", "subject"):
        raise ValueError('The list type should be either "target" or "subject"')



    for ppi in tqdm(ppi_list):
        learning_obj = read_learning_object(config)
        contacts_ag, contacts_ab = extract_contact(ppi, list_type, config)
        if len(ppi.split('_'))==2:
            pid, ch1 = ppi.split('_')
            compute_descriptor_from_contact(pid, ch1, contacts_ag, learning_obj, config)
        else:
            pid, ch1, ch2 = ppi.split('_')
            compute_descriptor_from_contact(pid, ch2, contacts_ag, learning_obj, config)
            compute_descriptor_from_contact(pid, ch1, contacts_ab, learning_obj, config)

    return