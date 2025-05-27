import os
import re

folder_path = r"C:\Users\Jordan\OneDrive - University of Leeds\Computational Data\ClusterGen\RM734\3AP_XTB1XTB2"
energy_values = []
total_runtimes_sec = []

for filename in os.listdir(folder_path):
    match = re.match(r"opt_cluster_\d+_stage" + re.escape("1") + r"\.out$", filename)
    if match:
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                match_energy = re.search(r"^\s*TOTAL ENERGY\s+([-+]?[0-9]*\.?[0-9]+)\s+Eh", line)
                if match_energy:
                    energy_values.append(match_energy.group(1))

                match_runtime = re.search(
                    r"CREST runtime \(total\)\s+(\d+) d,\s+(\d+) h,\s+(\d+) min,\s+([0-9.]+) sec", line)
                if match_runtime:
                    d, h, m, s = map(float, match_runtime.groups())
                    total_seconds = d * 86400 + h * 3600 + m * 60 + s
                    total_runtimes_sec.append(total_seconds)


# Save energies
energy_output = os.path.join(folder_path, "energy_values.txt")
with open(energy_output, 'w') as out_file:
    for value in energy_values:
        out_file.write(value + '\n')

# Compute average runtime
if total_runtimes_sec:
    avg_seconds = sum(total_runtimes_sec) / len(total_runtimes_sec)
    avg_d = int(avg_seconds // 86400)
    avg_h = int((avg_seconds % 86400) // 3600)
    avg_m = int((avg_seconds % 3600) // 60)
    avg_s = avg_seconds % 60

    runtime_output = os.path.join(folder_path, "average_runtime.txt")
    with open(runtime_output, 'w') as f:
        f.write(f"Average Runtime: {avg_d} d, {avg_h} h, {avg_m} min, {avg_s:.3f} sec\n")

    print(f"Saved average runtime to {runtime_output}")

print(f"Saved {len(energy_values)} TOTAL ENERGY values to {energy_output}")
