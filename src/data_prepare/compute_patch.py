import pymesh
import pdb
import time
import numpy as np
from sklearn.decomposition import PCA

from masif.source.masif_modules.read_data_from_surface import read_data_from_surface, compute_shape_complementarity
from masif.source.default_config.masif_opts import masif_opts

from sklearn.neighbors import KDTree
import os


def get_iface_verticies(mesh):
    iface = mesh.get_attribute('vertex_iface')
    vertices = mesh.vertices
    iface_indx = np.where(iface>0)
    return vertices[iface_indx]


def save_precompute(pid, ch, config, input_feat, rho, theta, mask, neigh_indices, iface_labels, verts):

    out_patch_dir = config['dirs']['patches']
    my_precomp_dir = out_patch_dir + pid + '/'
    if not os.path.exists(my_precomp_dir):
        os.mkdir(my_precomp_dir)

    np.save(my_precomp_dir + pid + '_' + ch + '_rho_wrt_center', rho)
    np.save(my_precomp_dir + pid + '_' + ch + '_theta_wrt_center', theta)
    np.save(my_precomp_dir + pid + '_' + ch + '_input_feat', input_feat)
    np.save(my_precomp_dir + pid + '_' + ch + '_mask', mask)
    np.save(my_precomp_dir + pid + '_' + ch + '_list_indices', neigh_indices)
    np.save(my_precomp_dir + pid + '_' + ch + '_iface_labels', iface_labels)
    # Save x, y, z
    np.save(my_precomp_dir + pid + '_' + ch + '_X.npy', verts[:, 0])
    np.save(my_precomp_dir + pid + '_' + ch + '_Y.npy', verts[:, 1])
    np.save(my_precomp_dir + pid + '_' + ch + '_Z.npy', verts[:, 2])


def compute_patch(ppi, config):

    params = masif_opts['ppi_search']
    ply_dir = config['dirs']['surface_ply']

    if len(ppi.split('_'))==2:
        pid, ch1 = ppi.split('_')
    else:
        pid, ch1, ch2 = ppi.split('_')

    ply_fn1 = ply_dir + '{}_{}.ply'.format(pid, ch1)
    input_feat1, rho1, theta1, mask1, neigh_indices1, iface_labels1, verts1 = read_data_from_surface(ply_fn1, params)

    out_patch_dir = config['dirs']['patches']
    my_precomp_dir = out_patch_dir + pid + '/'
    save_precompute(pid, ch1, config, input_feat1, rho1, theta1, mask1, neigh_indices1, iface_labels1, verts1)

    if len(ppi.split('_')) == 3:
        pid, ch1, ch2 = ppi.split('_')



        ply_fn2 = ply_dir + '{}_{}.ply'.format(pid, ch2)
        input_feat2, rho2, theta2, mask2, neigh_indices2, iface_labels2, verts2 = read_data_from_surface(ply_fn2, params)

        p1_sc_labels, p2_sc_labels = compute_shape_complementarity(ply_fn1, ply_fn2, neigh_indices1,
                                                                   neigh_indices2, rho1, rho2, mask1,
                                                                   mask2, params)
        save_precompute(pid, ch2, config, input_feat2, rho2, theta2, mask2, neigh_indices2, iface_labels2, verts2)
        np.save(my_precomp_dir + pid + '_' + ch1 + '_sc_labels.npy', p1_sc_labels)
        np.save(my_precomp_dir + pid + '_' + ch2 + '_sc_labels.npy', p2_sc_labels)



