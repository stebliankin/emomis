import numpy as np
import os
import pdb

from shutil import copyfile, rmtree

# Local includes
import pymesh

from masif.source.default_config.masif_opts import masif_opts
from masif.source.triangulation.computeMSMS import computeMSMS
from masif.source.triangulation.fixmesh import fix_mesh
from masif.source.input_output.extractPDB import extractPDB
from masif.source.input_output.save_ply import save_ply
from masif.source.triangulation.computeHydrophobicity import computeHydrophobicity
from masif.source.triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from masif.source.triangulation.computeAPBS import computeAPBS
from masif.source.triangulation.compute_normal import compute_normal
from sklearn.neighbors import KDTree


def triangulate_one(pid, ch, config):
    """
    triangulate one chain
    """
    chains_pdb_dir = config['dirs']['chains_pdb']

    tmp_pdb_dir = chains_pdb_dir + pid + '_' + ch + '/'
    if not os.path.exists(tmp_pdb_dir):
        os.mkdir(tmp_pdb_dir)

    pdb_filename = config['dirs']['protonated_pdb'] + pid + '.pdb'

    # Extract chains for each interacting protein
    out_filename1 = tmp_pdb_dir + pid + '_' + ch
    extractPDB(pdb_filename, out_filename1+ '.pdb', ch)

    vertices1, faces1, normals1, names1, areas1 = computeMSMS(out_filename1+ '.pdb', protonate=True)

    #exit()
    vertex_hbond = computeCharges(out_filename1, vertices1, names1)
    vertex_hphobicity = computeHydrophobicity(names1)


    vertices2 = vertices1
    faces2 = faces1

    # Fix the mesh.
    mesh = pymesh.form_mesh(vertices2, faces2)
    regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])


    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    # hbonds
    vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices1,\
                                          vertex_hbond, masif_opts)

    vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices1, \
                                               vertex_hphobicity, masif_opts)


    vertex_charges = computeAPBS(regular_mesh.vertices, out_filename1 + ".pdb", out_filename1)

    iface = np.zeros(len(regular_mesh.vertices))

    v3, f3, _, _, _ = computeMSMS(pdb_filename, protonate=True)
    # Regularize the mesh
    mesh = pymesh.form_mesh(v3, f3)
    # I believe It is not necessary to regularize the full mesh. This can speed up things by a lot.
    full_regular_mesh = mesh
    # Find the vertices that are in the iface.
    v3 = full_regular_mesh.vertices
    # Find the distance between every vertex in regular_mesh.vertices and those in the full complex.
    kdt = KDTree(v3)
    d, r = kdt.query(regular_mesh.vertices)
    d = np.square(d)  # Square d, because this is how it was in the pyflann version.
    assert (len(d) == len(regular_mesh.vertices))
    iface_v = np.where(d >= 2.0)[0]
    iface[iface_v] = 1.0
    # Convert to ply and save.

    outply = config['dirs']['surface_ply'] + pid + '_' + ch

    copyfile(out_filename1+".pdb", chains_pdb_dir + '{}_{}.pdb'.format(pid, ch))

    save_ply(outply + ".ply", regular_mesh.vertices, \
             regular_mesh.faces, normals=vertex_normal, charges=vertex_charges, \
             normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity, \
             iface=iface)

    rmtree(tmp_pdb_dir)

    return