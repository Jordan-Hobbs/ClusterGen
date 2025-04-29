import zipfile
import sys

import writers
import mol_ops
import geom_ops

def main():
    num_clusters = 5
    num_rings = 2
    num_ap = 1
    mol_sep = 6
    min_sep = 2.5
    total_attempts = 0

    smiles = "O=C(OC1=CC(F)=C(C2=CC=C(C#N)C(F)=C2)C(F)=C1)C(C(F)=C3)=CC=C3CCC"
    mol = mol_ops.find_min_conformer(smiles, num_conf=100)
    mol = geom_ops.orientate_rod(mol)

    with zipfile.ZipFile("test.zip", "w") as zipf:
        writers.write_xyz(mol, "molecule", zipf=zipf)

        for i in range(num_clusters):
            print(f"Generating cluster {i+1}...")

            attempts = 0
            while True:
                mol_cluster = mol_ops.build_cluster_hex_rings(
                    mol, num_rings=num_rings, min_dist=min_sep, 
                    r_start=mol_sep, r_step=mol_sep, num_ap=num_ap
                )
                attempts += 1
                total_attempts += 1

                if mol_cluster is not None:
                    break

                sys.stdout.write(f"\rCluster {i+1}: {attempts} attempt(s) so far...")
                sys.stdout.flush()

            print(f"Cluster {i+1} generated after {attempts} attempt(s).")

            cluster_name = f"cluster_{i+1}"
            writers.write_xyz(mol_cluster, cluster_name, zipf=zipf)
            writers.write_toml(cluster_name, zipf=zipf)

        writers.write_sh(num_jobs=num_clusters, zipf=zipf, job_name="cluster_calc")

    print(f"All clusters generated. Total placement attempts: {total_attempts}")

if __name__ == "__main__":
    main()
