import numpy as np

def write_xyz(molecule, file_name, zipf):
    xyz_pos = molecule.GetConformer().GetPositions()
    xyz_pos = np.round(xyz_pos, 5)
    no_atoms = molecule.GetConformer().GetNumAtoms()

    lines = [f"{no_atoms}", ""]
    for index, atom in enumerate(molecule.GetAtoms()):
        atom_symbol = atom.GetSymbol()
        x, y, z = xyz_pos[index]
        lines.append(f"{atom_symbol} {x} {y} {z}")

    file_content = "\n".join(lines)
    zipf.writestr(f"{file_name}.xyz", file_content)

def write_toml(file_name, zipf, num_cpus, optlevel, gfn_method):
    toml_text = (
        "# CREST 3 input file\n"
        f"input = \"{file_name}.xyz\"\n"
        "runtype = \"ancopt\"\n"
        f"threads = {num_cpus}\n"
        "\n"
        "[calculation]\n"
        f"optlev = \"{optlevel}\"\n"
        "\n"
        "[[calculation.level]]\n"
        f"method = \"{gfn_method}\""
    )
    zipf.writestr(f"{file_name}_stage1.toml", toml_text)

def write_toml2(file_name, zipf, num_cpus, optlevel, gfn_method):
    toml_text = (
        "# CREST 3 input file\n"
        f"input = \"crestopt.xyz\"\n"
        "runtype = \"ancopt\"\n"
        f"threads = {num_cpus}\n"
        "\n"
        "[calculation]\n"
        f"optlev = \"{optlevel}\"\n"
        "\n"
        "[[calculation.level]]\n"
        f"method = \"{gfn_method}\""
    )
    zipf.writestr(f"{file_name}_stage2.toml", toml_text)

def write_sh(num_jobs, zipf, job_name, run_time, num_cpus, email):
    sh_text = (
        "#!/bin/bash\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --time={run_time}\n"
        "#SBATCH --ntasks=1\n"
        "#SBATCH --mem-per-cpu=1G\n"
        f"#SBATCH --cpus-per-task={num_cpus}\n"
        "#SBATCH --mail-type=BEGIN,END,FAIL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --array=1-{num_jobs}\n"
        "\n"
        "module load crest\n"
        "\n"
        "# Make directories\n"
        "mkdir -p opt_clusters\n"
        "mkdir -p cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "cd cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "\n"
        "# Copy stage 1 input files\n"
        "cp \"../cluster_${SLURM_ARRAY_TASK_ID}_stage1.toml\" \"../cluster_${SLURM_ARRAY_TASK_ID}.xyz\" .\n"
        "\n"
        "# Run stage 1\n"
        "crest cluster_${SLURM_ARRAY_TASK_ID}_stage1.toml > stage1.out\n"
        "wait\n"
        "\n"
        "# Save output\n"
        "if [ -f \"crestopt.xyz\" ]; then\n"
        "    cp \"crestopt.xyz\" \"../opt_clusters/opt_cluster_${SLURM_ARRAY_TASK_ID}_stage1.xyz\"\n"
        "fi\n"
        "cp stage1.out \"../opt_clusters/opt_cluster_${SLURM_ARRAY_TASK_ID}_stage1.out\"\n"
        "\n"
        "# Check and run stage 2\n"
        "if [ -f \"../cluster_${SLURM_ARRAY_TASK_ID}_stage2.toml\" ]; then\n"
        "    cp \"../cluster_${SLURM_ARRAY_TASK_ID}_stage2.toml\" .\n"
        "    cp crestopt.xyz input.xyz\n"
        "    crest cluster_${SLURM_ARRAY_TASK_ID}_stage2.toml > stage2.out\n"
        "    wait\n"
        "    if [ -f \"crestopt.xyz\" ]; then\n"
        "        cp \"crestopt.xyz\" \"../opt_clusters/opt_cluster_${SLURM_ARRAY_TASK_ID}_stage2.xyz\"\n"
        "    fi\n"
        "    cp stage2.out \"../opt_clusters/opt_cluster_${SLURM_ARRAY_TASK_ID}_stage2.out\"\n"
        "fi\n"
        "\n"
        "cd ..\n"
        "rm -rf cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "rm -f cluster_${SLURM_ARRAY_TASK_ID}.xyz cluster_${SLURM_ARRAY_TASK_ID}_stage1.toml cluster_${SLURM_ARRAY_TASK_ID}_stage2.toml"
    )
    zipf.writestr(f"{job_name}.sh", sh_text)