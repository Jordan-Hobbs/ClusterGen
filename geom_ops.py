import numpy as np
from numpy.linalg import norm
from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Geometry import Point3D

def orientate_rod(molecule):
    conf = molecule.GetConformer()
    xyz_pos = conf.GetPositions()
    mean = xyz_pos.mean(axis=0)
    centered = xyz_pos - mean

    pca = PCA(n_components=3)
    pca.fit(centered)
    rod_axis = pca.components_[0]
    z_axis = np.array([0, 0, 1])

    a, b = (rod_axis / norm(rod_axis)).reshape(3), (z_axis / norm(z_axis)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    
    for i in range(molecule.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        vec = np.array([pos.x, pos.y, pos.z])
        rotated_vec = rotation_matrix @ vec
        conf.SetAtomPosition(i, rotated_vec)
    return molecule

def translate_mol(molecule, dx=0, dy=0, dz=0):
    molecule = Chem.Mol(molecule)
    conf = molecule.GetConformer()
    for i in range(molecule.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(pos.x + dx, pos.y + dy, pos.z + dz))
    return molecule

def rotate_mol(molecule, angle_deg):
    molecule = Chem.Mol(molecule)
    conf = molecule.GetConformer()
    angle_rad = np.radians(angle_deg)
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)
    
    for i in range(molecule.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x_new = cos_angle * pos.x - sin_angle * pos.y
        y_new = sin_angle * pos.x + cos_angle * pos.y
        conf.SetAtomPosition(i, Point3D(x_new, y_new, pos.z))
    return molecule

def tilt_mol(molecule, tilt_angle):
    conf = molecule.GetConformer()
    xyz_pos = conf.GetPositions()
    mean = xyz_pos.mean(axis=0)

    tilt_angle_rad = np.radians(tilt_angle)
    azimuth = np.radians(np.random.uniform(0, 360))

    tilt_axis = np.array([np.cos(azimuth), np.sin(azimuth), 0.0])

    v = tilt_axis / norm(tilt_axis)
    s = np.sin(tilt_angle_rad)
    c = np.cos(tilt_angle_rad)
    kmat = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    rotation_matrix = np.eye(3) + s * kmat + (1 - c) * kmat @ kmat

    for i in range(molecule.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i)) - mean
        rotated = rotation_matrix @ pos + mean
        conf.SetAtomPosition(i, rotated)
    return molecule