import numpy as np

def write_xyz(molecule, file_name):
    xyz_pos = molecule.GetConformer().GetPositions()
    xyz_pos = np.round(xyz_pos, 5)
    no_atoms = molecule.GetConformer().GetNumAtoms()

    with open(f"{file_name}.xyz", "w", newline="\n") as file:
        file.write(f"{no_atoms}\n\n")
        for index, atom in enumerate(molecule.GetAtoms()):
            atom_symbol = atom.GetSymbol()
            x, y, z = xyz_pos[index]
            line = f"{atom_symbol} {x} {y} {z}"
            if index < no_atoms - 1:
                file.write(line + "\n")
            else:
                file.write(line)

    print(f"{file_name}.xyz file written successfully")

def write_sh(file_name, job_num, run_time="24:00:00", num_cpus=4, email="j.l.hobbs@leeds.ac.uk"):
    sh_text =(
        "#!/bin/bash\n"
        f"#SBATCH --job-name={file_name}\n"
        f"#SBATCH --time={run_time}\n"
        "#SBATCH --ntasks=1\n"
        "#SBATCH --mem-per-cpu=1G\n"
        f"#SBATCH --cpus-per-task={num_cpus}\n"
        "#SBATCH --mail-type=BEGIN,END,FAIL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --array=1-{job_num}\n"
        "\n"
        "module load crest\n"
        "crest cluster_$SLURM_ARRAY_TASK_ID.toml > cluster_$SLURM_ARRAY_TASK_ID.out"
    )

    with open(f"{file_name}.sh", "w", newline="\n") as file:
        file.write(sh_text)

    print(f"{file_name}.sh file written successfully")

def write_toml(file_name, num_cpus=4, optlevel="extreme", gfn_method="gfnff"):
    toml_text = ( 
        "# CREST 3 input file\n"
        "input = \"molecule.xyz\"\n"
        f"ensemble_input=\"{file_name}.xyz\"\n"
        "runtype = \"ancopt_ensemble\"\n"
        f"threads = {num_cpus}\n"
        "\n"
        "[calculation]\n"
        f"optlev = \"{optlevel}\"\n"
        "\n"
        "[[calculation.level]]\n"
        f"method = \"{gfn_method}\"\n"
        f"calcspace = \"{file_name}\""
    )

    toml_name = f"{file_name}.toml"
    with open(toml_name, "w", newline="\n") as file:
        file.write(toml_text)

    print(f"{toml_name} file written successfully")