import os
import zipfile

def write_inp(FileName, zipf, NumCPUs, MemPerCPU, xyz_data):
    inp_text = (
        "! r2SCAN-3c def2-TZVP\n"
        f"%MAXCORE {int(MemPerCPU*1000*0.75)}\n"
        "%PAL\n"
        f" nprocs {NumCPUs}\n"
        "END\n"
        "* xyz 0 1\n"
        f"{xyz_data}*\n"
    )
    zipf.writestr(FileName, inp_text)  # Use exact filename passed here

def write_sh(num_jobs, zipf, job_name, run_time, num_cpus, mem_cpus, email):
    sh_text = (
        "#!/bin/bash\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --time={run_time}\n"
        f"#SBATCH --ntasks={num_cpus}\n"
        f"#SBATCH --mem-per-cpu={mem_cpus}G\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --mail-type=BEGIN,END,FAIL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --array=1-{num_jobs}\n"
        "\n"
        "module load orca\n"
        "module load openmpi\n"
        "export RSH_COMMAND=\"/usr/bin/ssh\"\n"
        "\n"
        "# Make directories\n"
        "mkdir -p opt_clusters\n"
        "mkdir -p cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "cd cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "\n"
        "# Copy stage 1 input files\n"
        "cp \"../opt_cluster_${SLURM_ARRAY_TASK_ID}.inp\" .\n"
        "\n"
        "# Run stage 1\n"
        "/opt/apps/pkg/applications/orca/orca_6_0_1_linux_x86-64_shared_openmpi416/orca opt_cluster_${SLURM_ARRAY_TASK_ID}.inp > opt_cluster_${SLURM_ARRAY_TASK_ID}.out\n"
        "wait\n"
        "\n"
        "# Save output\n"
        "cp opt_cluster_${SLURM_ARRAY_TASK_ID}.out \"../opt_clusters/opt_cluster_${SLURM_ARRAY_TASK_ID}.out\"\n"
        "cp opt_cluster_${SLURM_ARRAY_TASK_ID}.property.txt \"../opt_clusers/opt_cluster_${SLURM_ARRAY_TASK_ID}.property.txt\"\n"
        "\n"
        "cd ..\n"
        "rm -rf cluster_job_${SLURM_ARRAY_TASK_ID}\n"
        "rm -f opt_cluster_${SLURM_ARRAY_TASK_ID}.inp"
    )
    zipf.writestr(f"{job_name}.sh", sh_text)

def cluster_iter(
    xyz_folder,
    zip_name="cluster_iter.zip",
    NumCPUs=8,
    MemPerCPU=10,
    job_name="orca_array_job",
    run_time="24:00:00",
    email="your.email@example.com"
):
    xyz_files = [f for f in os.listdir(xyz_folder) if f.endswith(".xyz")]
    
    zip_path = os.path.join(xyz_folder, zip_name)
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
        for filename in xyz_files:
            xyz_path = os.path.join(xyz_folder, filename)
            with open(xyz_path, 'r') as xyz_file:
                lines = xyz_file.readlines()
                xyz_data = "".join(lines[2:])  # Skip header
            inp_filename = os.path.splitext(filename)[0] + ".inp"  # exact base name + .inp
            write_inp(inp_filename, zipf, NumCPUs, MemPerCPU, xyz_data)

        write_sh((len(xyz_files)+1), zipf, job_name, run_time, NumCPUs, MemPerCPU, email)

# Example usage:
cluster_iter(
    r"C:\Users\phyjlho\OneDrive - University of Leeds\Computational Data\ClusterGen\RM734\3AP_XTB2",
    NumCPUs=8,
    MemPerCPU=8,
    job_name="cluster_iter",
    run_time="48:00:00",
    email="j.l.hobbs@leeds.ac.uk"
)