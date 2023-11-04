output_file = "./ASDEG/results/snpeff/snpeff_concatenated"
directory = "./ASDEG/results/snpeff/"

# Clear the output file first, if it exists
if os.path.exists(output_file):
    os.remove(output_file)

# Get all the .ran files in the directory
files = [f for f in os.listdir(directory) if f.endswith('.ran')]

# For each file
for file in files:
    filepath = os.path.join(directory, file)
    
    # Extract tissue name from filename
    tissue = "_".join(file.split('_')[0:2])
    
    with open(filepath, 'r') as infile, open(output_file, 'a') as outfile:
        lines = infile.readlines()[11:]  # Skip the first 11 lines
        for line in lines:
            cols = line.split()
            new_line = f"{cols[0]} {cols[1]} {cols[3]} {cols[4]} {cols[7]} {tissue}\n"
            outfile.write(new_line)

# Display first few lines of the result for verification
with open(output_file, 'r') as f:
    print("".join(f.readlines()[:5]))
