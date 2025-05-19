import os
import re

# Folder containing the .xyz files
folder_path = r"C:\Users\phyjlho\OneDrive - University of Leeds\Computational Data\ClusterGen\RM734\2AP_XTB2"

# This will store the extracted Etot values
etot_values = []

# Loop through all .xyz files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".xyz"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r") as file:
            lines = file.readlines()
            if len(lines) >= 2:
                # Match 'Etot=' followed by a number (x)
                match = re.search(r"Etot=\s*([-+]?[0-9]*\.?[0-9]+)", lines[1])
                if match:
                    etot_values.append(match.group(1))

# Save the Etot values to a text file
output_file = os.path.join(folder_path, "energy_values.txt")
with open(output_file, 'w') as out_file:
    for value in etot_values:
        out_file.write(value + '\n')

print(f"Saved {len(etot_values)} Etot values to {output_file}")
