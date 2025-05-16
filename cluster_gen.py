import gc
import zipfile
import multiprocessing
import logging

import numpy as np
import scipy
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

import geom_ops
import writers

logger = logging.getLogger(__name__)

def cluster_gen(args):
    total_attempts = 0
    smiles = args.SmilesString
    mol = find_min_conformer(smiles, num_conf=args.RDNumConf)
    mol = geom_ops.orientate_rod(mol)

    logger.info("------------------------------------------------------------")
    logger.info(f"Generating {args.NumClusters} clusters.")

    task_args = [(args, i, mol) for i in range(args.NumClusters)]
    total_attempts = 0
    successful_clusters = 0

    with multiprocessing.Pool(4) as pool, zipfile.ZipFile(args.Output, "w") as zipf, tqdm(total=args.NumClusters, desc="Clusters generated") as pbar:
        for cluster_index, mol_data in pool.imap_unordered(generate_single_cluster, task_args):
            if mol_data is None:
                continue

            mol_cluster, cluster_attempts = mol_data
            total_attempts += cluster_attempts
            successful_clusters += 1
            avg_attempts = total_attempts / successful_clusters

            # Update one clean status line below progress bar
            pbar.set_postfix({
                "Cluster Attempts": total_attempts,
                "Avg/Cluster": f"{avg_attempts:.2f}"
            })

            cluster_name = f"cluster_{cluster_index + 1}"
            writers.write_xyz(mol_cluster, cluster_name, zipf=zipf)
            writers.write_toml(
                file_name=cluster_name,
                zipf=zipf,
                num_cpus=args.NumCPUs,
                optlevel=args.CRESTOptLevel,
                gfn_method=args.CRESTMethod,
            )

            del mol_cluster
            gc.collect()
            pbar.update(1)

        writers.write_sh(
            num_jobs=args.NumClusters,
            zipf=zipf,
            job_name="cluster_calc",
            num_cpus=args.NumCPUs,
            run_time=args.RunTime,
            email=args.Email,
        )
    logger.info("------------------------------------------------------------")
    logger.info("Cluster generation complete")
    logger.info(f"Total number of clusters generated: {successful_clusters}")
    logger.info(f"Total placement attempts: {total_attempts}")
    logger.info(f"Average number of attempts: {avg_attempts}")
    logger.info("------------------------------------------------------------")
    logger.info("Summary of parameters:")
    logger.info(f"  SMILES           : {args.SmilesString}")
    logger.info(f"  Num Clusters     : {args.NumClusters}")
    logger.info(f"  Rings/Cluster    : {args.NumRings}")
    logger.info(f"  Layers           : {args.NumLayers or 1}")
    logger.info(f"  Mol Separation   : {args.MolSep} Å")
    logger.info(f"  Min Distance     : {args.MinSep} Å")
    logger.info(f"  Antiparallel mols: {args.NumAP}")
    logger.info(f"  CREST Method     : {args.CRESTMethod}")
    logger.info(f"  CREST OptLevel   : {args.CRESTOptLevel}")



def generate_single_cluster(args_tuple):
    args, cluster_index, mol = args_tuple

    cluster_attempts = 0
    try:
        while True:
            cluster_attempts += 1
            mol_cluster, _ = build_cluster_rings(
                mol,
                num_rings=args.NumRings,
                min_dist=args.MinSep,
                r_start=args.MolSep,
                r_step=args.MolSep,
                num_ap=args.NumAP,
            )
            if mol_cluster is not None:
                return cluster_index, (mol_cluster, cluster_attempts)
    except Exception:
        return cluster_index, None




def find_min_conformer(smiles, num_conf: int = 100, max_opt_iters: int = 1000):
    logger.info("------------------------------------------------------------")
    logger.info("Generating single molecule structure")

    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)
    Chem.rdCoordGen.AddCoords(molecule_h)

    logger.info(f"Generating {num_conf} conformers using RDKit.")
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
        logger.info("All conformers converged.")
    else:
        logger.warning("Not all conformers converged.")

    min_energy = float("inf")
    min_index = 0
    for index, (_, energy) in enumerate(conf_energy):
        if energy < min_energy:
            min_energy = energy
            min_index = index
    mol_min = Chem.Mol(molecule_h, False, min_index)

    logger.info(f"Selected conformer {min_index} with energy {min_energy:.2f}")
    return mol_min


def place_single_ring(mol, mols, existing_coords, radius, n_mols,
               min_dist=2.5, max_attempts=1000, flip_indices=None):

    total_ring_attempts = 0
    if flip_indices is None:
        flip_indices = []

    for i in range(n_mols):
        for _ in range(max_attempts):
            total_ring_attempts += 1
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
                tree = scipy.spatial.cKDTree(existing_coords_np)
                distances, _ = tree.query(new_coords, distance_upper_bound=min_dist)
                too_close = np.any(distances < min_dist)
            else:
                too_close = False

            if not too_close:
                mols.append(mol_copy)
                existing_coords.extend(new_coords)
                break
        else:
            return False, total_ring_attempts
        
    return True, total_ring_attempts


def build_cluster_rings(
    mol, num_rings=1, r_start=6.0, r_step=6.0, min_dist=3, max_attempts=1000, num_ap=0
):
    total_mols = 3 * num_rings**2 + 3 * num_rings + 1
    flip_indices = np.sort(
        np.random.choice(np.arange(total_mols), size=min(num_ap, total_mols), replace=False)
    )

    mol0 = Chem.Mol(mol)
    if 0 in flip_indices:
        geom_ops.flip_mol(mol0)
    mols = [mol0]
    existing_coords = list(mol0.GetConformer().GetPositions())

    total_attempts = 0
    current_index = 1
    for ring in range(1, num_rings + 1):
        n_mols = 6 * ring
        radius = r_start + (ring - 1) * r_step

        ring_flip_indices = [i - current_index for i in flip_indices if current_index <= i < current_index + n_mols]

        success, ring_attempts = place_single_ring(
            mol, mols, existing_coords, radius, n_mols,
            min_dist=min_dist, max_attempts=max_attempts,
            flip_indices=ring_flip_indices
        )
        total_attempts += ring_attempts
        if not success:
            return None, total_attempts

        current_index += n_mols

    combined = mols[0]
    for m in mols[1:]:
        combined = Chem.CombineMols(combined, m)

    return combined,total_attempts
