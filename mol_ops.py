import numpy as np

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import geom_ops

def find_min_conformer(smiles, num_conf: int = 100, max_opt_iters: int = 1000):
    print("\n----------------------------------------------------------------")
    print("Generating initial CREST input structure:\n")

    molecule = Chem.MolFromSmiles(smiles)

    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    print(
        f"Generating {num_conf} conformers for initial sorting and "
        "optimization using RDKit."
    )
    rdDistGeom.EmbedMultipleConfs(
        molecule_h,
        num_conf,
        params=rdDistGeom.ETKDGv3()
    )
    conf_energy = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
        molecule_h,
        maxIters=max_opt_iters,
        ignoreInterfragInteractions=False
    )

    if all(conf_set[0] == 0 for conf_set in conf_energy):
        print("All conformers converged.")
    else:
        print("WARNING! Not all conformers converged.")

    min_energy = float("inf")
    min_index = 0
    for index, (_, energy) in enumerate(conf_energy):
        if energy < min_energy:
            min_energy = energy
            min_index = index
    mol_min = Chem.Mol(molecule_h, False, min_index)

    return mol_min

def place_ring(mol, mols, existing_coords, radius, n_mols,
               min_dist=2.5, max_attempts=2000, flip_indices=None):

    if flip_indices is None:
        flip_indices = []

    for i in range(n_mols):
        for attempt in range(max_attempts):
            angle_deg = (i * 360 / n_mols)
            angle_rad = np.radians(angle_deg)

            dtheta = np.random.uniform(0, 360)
            tilt = np.random.uniform(0, 10)
            dx = radius * np.cos(angle_rad) + np.random.uniform(-1, 1)
            dy = radius * np.sin(angle_rad) + np.random.uniform(-1, 1)
            dz = np.random.uniform(-3, 3)

            mol_copy = Chem.Mol(mol)
            if i in flip_indices:
                geom_ops.flip_mol(mol_copy)
            geom_ops.rotate_mol(mol_copy, dtheta)
            geom_ops.tilt_mol(mol_copy, tilt)
            geom_ops.translate_mol(mol_copy, dx=dx, dy=dy, dz=dz)

            new_coords = np.array(mol_copy.GetConformer().GetPositions())
            existing_coords_np = np.array(existing_coords)

            if existing_coords_np.shape[0] > 0:
                dists = np.linalg.norm(existing_coords_np[:, None] - new_coords, axis=2)
                too_close = np.any(dists < min_dist)
            else:
                too_close = False

            if not too_close:
                mols.append(mol_copy)
                existing_coords.extend(new_coords)
                break
        else:
            return False

    return True

def build_cluster_hex_rings(
    mol, num_rings=1, r_start=6.0, r_step=6.0, min_dist=3, max_attempts=1000, num_ap=0
    ):
    mols = [mol]
    existing_coords = list(mol.GetConformer().GetPositions())
    total_mols = 3*num_rings**2+3*num_rings+1
    flip_indices = np.sort(
    np.random.choice(np.arange(1, total_mols), size=min(num_ap, total_mols - 1), replace=False)
    )


    current_index = 1
    for ring in range(1, num_rings + 1):
        n_mols = 6 * ring
        radius = r_start + (ring - 1) * r_step

        ring_flip_indices = [i - current_index for i in flip_indices if current_index <= i < current_index + n_mols]

        if not place_ring(
            mol, mols, existing_coords, radius, n_mols,
            min_dist=min_dist, max_attempts=max_attempts,
            flip_indices=ring_flip_indices
        ):
            return None

        current_index += n_mols

    combined = mols[0]
    for m in mols[1:]:
        combined = Chem.CombineMols(combined, m)

    return combined