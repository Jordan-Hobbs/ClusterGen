import random
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import geom_ops
import writers

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

def build_cluster(mol, n_mols, radius=6.0, min_dist=3, max_attempts=1000):

    mols = [mol]
    existing_coords = list(mol.GetConformer().GetPositions())

    for i in range(n_mols):
        for attempt in range(max_attempts):
            angle_deg = i * 360/n_mols
            angle_rad = np.radians(angle_deg)
            dx = radius * np.cos(angle_rad) + random.uniform(-2, 2)
            dy = radius * np.sin(angle_rad) + random.uniform(-2, 2)
            dz = random.uniform(-3, 3)
            tilt = random.uniform(0, 10)
            
            dtheta = random.uniform(0, 360)
            rotated = geom_ops.rotate_mol(mol, dtheta)
            tilted = geom_ops.tilt_mol(rotated, tilt)
            translated = geom_ops.translate_mol(tilted, dx=dx, dy=dy, dz=dz)

            new_coords = list(translated.GetConformer().GetPositions())

            # Check for overlaps
            too_close = any(
                np.linalg.norm(pos1 - pos2) < min_dist
                for pos1 in new_coords
                for pos2 in existing_coords
            )

            if not too_close:
                mols.append(translated)
                existing_coords.extend(new_coords)
                break
        else:
            print(f"Warning: Could not place molecule {i+1} without overlap.")
            return None

    # Combine all molecules into one
    combined = mols[0]
    for m in mols[1:]:
        combined = Chem.CombineMols(combined, m)
    
    return combined

def main():
    smiles = "O=C(OC1=CC(F)=C(C2=CC=C(C#N)C(F)=C2)C(F)=C1)C(C(F)=C3)=CC=C3CCC"
    mol = find_min_conformer(smiles, num_conf=100)
    mol = geom_ops.orientate_rod(mol)
    writers.write_xyz(mol, f"molecule")

    num_clusters = 10
    for i in range(num_clusters):
        print(f"Generating cluster {i+1}...")
        mol_cluster = build_cluster(mol, 6, radius=6.0)
        while mol_cluster is None:
            mol_cluster = build_cluster(mol, 6, radius=6.0)
        writers.write_xyz(mol_cluster, f"cluster_{i+1}")
        writers.write_toml(f"cluster_{i+1}")
    writers.write_sh(f"cluster_calc", num_clusters)

if __name__ == "__main__":
    main()
