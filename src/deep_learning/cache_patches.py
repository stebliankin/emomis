# Header variables and parameters.
import sys
import pymesh
import os
import numpy as np
from IPython.core.debugger import set_trace
import importlib

from scipy.spatial import cKDTree

from masif.source.default_config.masif_opts import masif_opts
import pdb
#from tqdm import tqdm

"""
Source: https://github.com/LPDI-EPFL/masif/blob/master/source/masif_ppi_search/masif_ppi_search_cache_training_data.py
"""

def cache_patches(config):
    params = masif_opts['ppi_search']
    min_per_complex = config['masif_search']['cache_min_per_pdb']
    parent_in_dir = config['dirs']['patches']
    cache_dir = config['dirs']['db_masif_indx']
    ppi_list_db = [x.strip('\n') for x in open(config['db_list_processed'])]
    ply_dir = config['dirs']['surface_ply']
    cached_list_file = config['db_list_processed']

    # Read the positive first

    cache_log = cache_dir + "/cache_log.txt"

    p1target_rho_wrt_center = []
    p1target_theta_wrt_center = []
    p1target_input_feat = []
    p1target_mask = []
    p1target_iface_labels = []

    p2pos_rho_wrt_center = []
    p2pos_theta_wrt_center = []
    p2pos_input_feat = []
    p2pos_mask = []
    p2pos_iface_labels = []

    np.random.seed(0)
    p1target_names = []  # target names
    p2pos_names = []

    idx_count = 0
    processed_complexes_list = []

    for ppi_pair_id in (ppi_list_db):
        pid, ch1, ch2 = ppi_pair_id.split('_')
        in_dir = parent_in_dir + pid + '/'

        try:
            labels = np.load('{}/{}_{}_sc_labels.npy'.format(in_dir, pid, ch1), allow_pickle=True)
            # Take the median of the percentile 25 shape complementarity.
            mylabels = labels[0]
            labels = np.median(mylabels, axis=1)

        except Exception as e:
            print("Couldn't open {}/{}_{}_sc_labels.npy : ".format(in_dir, pid, ch1) + str(e))
            continue

        # Read the corresponding ply files.
        ply_fn1 = ply_dir + "{}_{}.ply".format(pid, ch1)
        ply_fn2 = ply_dir + "{}_{}.ply".format(pid, ch2)

        # pos_labels: points > max_sc_filt and >  min_sc_filt.
        pos_labels = np.where((labels > config['masif_search']['min_sc_fit']))[
            0]  # contact points in p1
        K = int(params['pos_surf_accept_probability'] * len(pos_labels))

        if K < 1 and min_per_complex == 0:
            print("For {} 0 patches found".format(ppi_pair_id))
            # If there are no patch within a threshold of sc, select a patch with the maximum sc
            continue
            # pos_labels = np.array([np.argmax(labels)])
            # K=1
            # pdb.set_trace()
        elif K < min_per_complex:
            print("No complementary patches found for {} with sc<{}.".format(pid, config['masif_search']['min_sc_fit']))
            print("Selecting {} most complementary patches".format(min_per_complex))
            K = min_per_complex
            pos_labels = (-labels).argsort()[:K]

        l = np.arange(len(pos_labels))  # l is the index of pos label
        np.random.shuffle(l)
        l = l[:K]
        l = pos_labels[l]  # contact points in p1

        v1 = pymesh.load_mesh(ply_fn1).vertices[l]  # coordinates of contact points in p1
        v2 = pymesh.load_mesh(ply_fn2).vertices  # all coordinates in p2

        # For each point in v1, find the closest point in v2.
        kdt = cKDTree(v2)
        d, r = kdt.query(v1)  # r is the indicies in p2 that are close to v1 (shape of p1)
        # Contact points: those within a cutoff distance.
        contact_points = np.where(d < params['pos_interface_cutoff'])[
            0]  # which r pass the interface cutoff (ex. [1,2]) (p2 indicies)
        try:
            k1 = l[contact_points]  # contact points in p1
        except:
            set_trace()
        k2 = r[contact_points]  # contact points in p2

        assert len(k1) == len(k2)
        n_pos = len(k1)
        if n_pos == 0:
            continue

        for ii in k1:
            p1target_names.append('{}_{}'.format(ppi_pair_id, ii))
        # for ii in p1target_names:
        #     k1_ind.append(ii)

        for ii in k2:
            p2pos_names.append('{}_{}'.format(ppi_pair_id, ii))

        processed_complexes_list.append(ppi_pair_id)
        rho_wrt_center = np.load(in_dir + '{}_{}_rho_wrt_center.npy'.format(pid, ch1), allow_pickle=True)
        theta_wrt_center = np.load(in_dir + '{}_{}_theta_wrt_center.npy'.format(pid, ch1), allow_pickle=True)
        input_feat = np.load(in_dir + '{}_{}_input_feat.npy'.format(pid, ch1),allow_pickle=True)
        mask = np.load(in_dir + '{}_{}_mask.npy'.format(pid, ch1),allow_pickle=True)
        list_indices = np.load(in_dir  + '{}_{}_list_indices.npy'.format(pid, ch1),allow_pickle=True)
        iface_labels = np.load(in_dir  + '{}_{}_iface_labels.npy'.format(pid, ch1),allow_pickle=True)

        p1target_rho_wrt_center.append(rho_wrt_center[k1])
        p1target_theta_wrt_center.append(theta_wrt_center[k1])
        p1target_input_feat.append(input_feat[k1])
        p1target_mask.append(mask[k1])

        iface_tmp = np.zeros((input_feat[k1].shape[0], input_feat[k1].shape[1]))
        for i in range(0, iface_tmp.shape[0]):
            for neigh_i in range(0, iface_tmp.shape[1]):
                try:
                    iface_tmp[i][neigh_i] = iface_labels[list_indices[k1][i][neigh_i]]
                except IndexError:
                    iface_tmp[i][neigh_i] = 0
        p1target_iface_labels.append(iface_tmp)
        #print(list_indices[k1].shape)

        # Read as positives those points.
        rho_wrt_center = np.load(in_dir + '{}_{}_rho_wrt_center.npy'.format(pid, ch2), allow_pickle=True)
        theta_wrt_center = np.load(in_dir + '{}_{}_theta_wrt_center.npy'.format(pid, ch2), allow_pickle=True)
        input_feat = np.load(in_dir + '{}_{}_input_feat.npy'.format(pid, ch2),allow_pickle=True)
        mask = np.load(in_dir + '{}_{}_mask.npy'.format(pid, ch2),allow_pickle=True)
        list_indices = np.load(in_dir  + '{}_{}_list_indices.npy'.format(pid, ch2),allow_pickle=True)
        iface_labels = np.load(in_dir  + '{}_{}_iface_labels.npy'.format(pid, ch2),allow_pickle=True)

        p2pos_rho_wrt_center.append(rho_wrt_center[k2])
        p2pos_theta_wrt_center.append(theta_wrt_center[k2])
        p2pos_input_feat.append(input_feat[k2])
        p2pos_mask.append(mask[k2])
        #p2pos_list_indices.append(list_indices[k2])

        iface_tmp = np.zeros((input_feat[k2].shape[0], input_feat[k2].shape[1]))
        for i in range(0, iface_tmp.shape[0]):
            for neigh_i in range(0, len(iface_labels[list_indices[k2][i]])):
                iface_tmp[i][neigh_i] = iface_labels[list_indices[k2][i][neigh_i]]
        p2pos_iface_labels.append(iface_tmp)
        idx_count += n_pos

    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    np.save(cache_dir + '/p1target_names.npy', p1target_names)

    p1target_rho_wrt_center = np.concatenate(p1target_rho_wrt_center, axis=0)
    np.save(cache_dir + '/p1target_rho_wrt_center.npy', p1target_rho_wrt_center)
    p1target_rho_wrt_center = None

    p1target_theta_wrt_center = np.concatenate(p1target_theta_wrt_center, axis=0)
    np.save(cache_dir + '/p1target_theta_wrt_center.npy', p1target_theta_wrt_center)
    p1target_theta_wrt_center = None

    p1target_input_feat = np.concatenate(p1target_input_feat, axis=0)
    np.save(cache_dir + '/p1target_input_feat.npy', p1target_input_feat)
    p1target_input_feat = None

    p1target_mask = np.concatenate(p1target_mask, axis=0)
    np.save(cache_dir + '/p1target_mask.npy', p1target_mask)
    p1target_mask = None

    p1target_iface_labels = np.concatenate(p1target_iface_labels, axis=0)
    np.save(cache_dir + '/p1target_iface_labels.npy',p1target_iface_labels)
    p1target_iface_labels = None

    p2pos_rho_wrt_center = np.concatenate(p2pos_rho_wrt_center, axis=0)
    np.save(cache_dir + '/p2pos_rho_wrt_center.npy', p2pos_rho_wrt_center)
    print("Read {} positive shapes".format(len(p2pos_rho_wrt_center)))
    p2pos_rho_wrt_center = None


    p2pos_theta_wrt_center = np.concatenate(p2pos_theta_wrt_center, axis=0)
    np.save(cache_dir + '/p2pos_theta_wrt_center.npy', p2pos_theta_wrt_center)
    p2pos_theta_wrt_center = None

    p2pos_input_feat = np.concatenate(p2pos_input_feat, axis=0)
    np.save(cache_dir + '/p2pos_input_feat.npy', p2pos_input_feat)
    p2pos_input_feat = None

    p2pos_mask = np.concatenate(p2pos_mask, axis=0)
    np.save(cache_dir + '/p2pos_mask.npy', p2pos_mask)
    p2pos_mask = None

    p2pos_iface_labels = np.concatenate(p2pos_iface_labels, axis=0)
    np.save(cache_dir+'/p2pos_iface_labels.npy',p2pos_iface_labels)
    p2pos_iface_labels = None

    ###
    # Added by Vitalii
    np.save(cache_dir + '/p2pos_names.npy', p2pos_names)

    with open(cached_list_file, 'w') as out:
        print("Extracted positive patch pairs from {} out of {} complexes.".format(len(processed_complexes_list), len(ppi_list_db)))
        for ppi in processed_complexes_list:
            out.write(ppi + "\n")

    print("Total number of extracted patches: {}".format(len(p1target_names)))



