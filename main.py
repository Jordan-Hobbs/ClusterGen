import writers
import mol_ops
import geom_ops

def main():
    smiles = "O=C(OC1=CC(F)=C(C2=CC=C(C#N)C(F)=C2)C(F)=C1)C(C(F)=C3)=CC=C3CCC"
    mol = mol_ops.find_min_conformer(smiles, num_conf=100)
    mol = geom_ops.orientate_rod(mol)
    writers.write_xyz(mol, f"molecule")

    num_clusters = 10
    for i in range(num_clusters):
        print(f"Generating cluster {i+1}...")
        mol_cluster = mol_ops.build_cluster(mol, 6, radius=6.0)
        while mol_cluster is None:
            mol_cluster = mol_ops.build_cluster(mol, 6, radius=6.0)
        writers.write_xyz(mol_cluster, f"cluster_{i+1}")
        writers.write_toml(f"cluster_{i+1}")
    writers.write_sh(f"cluster_calc", num_clusters)

if __name__ == "__main__":
    main()
