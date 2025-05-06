import os

# Set the folder containing the .xyz files
folder_path = r"C:\Users\Jordan\OneDrive - University of Leeds\Computational Data\ClusterGen\RM734\2AP"

# This will hold the values from the comment lines
comment_values = []

# Loop through all .xyz files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".xyz"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r") as file:
            lines = file.readlines()
            if len(lines) >= 2:
                comment_line = lines[1].strip()
                comment_values.append(comment_line)

# Save the values to a text file
output_file = os.path.join(folder_path, "energy_values.txt")
with open(output_file, 'w') as out_file:
    for value in comment_values:
        out_file.write(value + '\n')

print(f"Saved {len(comment_values)} values to {output_file}")
