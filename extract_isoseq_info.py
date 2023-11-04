def path_to_isoseq_directory():
    iso_seq_dir = '/group/ctbrowngrp/finnolab_shared/phased_isoseq'
    iso_seq_file = 'phased.nopartial.cleaned.vcf'
    iso_seq_file2 = 'phased.nopartial.cleaned.human_readable.txt'
    return iso_seq_dir, iso_seq_file, iso_seq_file2

def write_to_samtools(line_parts, out):
    out.write(f"{line_parts[0]}\t{int(line_parts[1]) - 1}\t{line_parts[1]}\n")

def write_to_tIDs(line_parts, tissue_name, out, PBID):
    out.write(f"{line_parts[0]}\t{line_parts[1]}\t{line_parts[3]}\t{line_parts[4]}\t{PBID}\t{tissue_name}\n")

def write_to_snpeff(line, out):
    if line.startswith('#'):
        out.write(line)
    else:
        line_parts = line.strip().split()
        line_parts[0] = line_parts[0].replace("chr", "")
        out.write("\t".join(line_parts) + "\n")

def write_to_haps(lines, out, PBID, tissue):
    # Loop over the lines in the file
    for line in lines:
        data = line.split()[0]
        data = str(data)
        for i in data:
            out.write(#i + "\t" +
                  data + "\t" +
                  PBID + "\t" + 
                  tissue + "\n")

def extract_isoseq_info(isoseq_dir, iso_seq_file, iso_seq_file2):
    # Loop over the directories in the path
    for dir_name in os.listdir(isoseq_dir):
        # Get the full path of the directory
        dir_path = os.path.join(isoseq_dir, dir_name)
        # Get the base name of the parent directory (tissue)
        tissue = os.path.basename(dir_path)
        # Correct names in Iso Seq data to match rna seq data 
        tissue_mapping = {
            'AH1_Heart': 'AH1_LeftVentricle',
            'AH1_Muscle': 'AH1_Longissimus',
            'AH3_Skin': 'AH3_Adipose',
            'AH4_Muscle': 'AH4_LongMuscle',
            'AH1_Skin' : 'AH1_AdiposeLoin'
        }
        print(f"Starting content extraction for {tissue}")
        if tissue in tissue_mapping:
            tissue = tissue_mapping[tissue]
        # Check if the directory is actually a directory (and not a file)
        if os.path.isdir(dir_path):
            # Define the name of the file to be ran with samtools
            file_for_samtools = "./ASDEG/results/samtools/" + tissue + "_extracted_lines_for_samtools"
            # Define the name of the file to be ran on snpeff
            file_for_snpeff = "./ASDEG/results/snpeff/" + tissue +  "_extracted_lines_for_snpeff" 
            # Define the name of the file for transcript IDs
            file_for_tIDs = "./ASDEG/results/tIDs/" + tissue + "_extracted_pb_ID"
            # Define the name of the file for Haplotypes
            file_for_Haps = "./ASDEG/results/Haps/" + tissue + "_extracted_Haplotypes"
            # Define the path of the sub dir
            for sub_dir in os.listdir(dir_path):
                subdir_path = os.path.join(dir_path, sub_dir)

                # Open all output files to write
                with open(file_for_samtools, "a") as samtools_out, \
                        open(file_for_tIDs, "a") as tIDs_out, \
                        open(file_for_snpeff, "a") as snpeff_out, \
                        open(file_for_Haps, "a") as haplotypes_out:

                    for file_name in os.listdir(subdir_path):
                        PBID = os.path.basename(subdir_path)
                        if file_name == iso_seq_file:
                            file_path = os.path.join(subdir_path, file_name)

                            with open(file_path, "r") as f:
                                for line in f:
                                    columns = line.split()

                                    if not line.startswith("#"):
                                        write_to_samtools(columns, samtools_out)
                                        write_to_tIDs(columns, tissue, tIDs_out, PBID)
                                    write_to_snpeff(line, snpeff_out)
                                    
                        if file_name == iso_seq_file2:
                            file_path = os.path.join(subdir_path, file_name)
                            with open(file_path, "r") as f:
                                # skip first two lines
                                lines = f.readlines()[2:]
                                write_to_haps(lines, haplotypes_out, PBID, tissue)
        
        print(f"Finished content extraction for {tissue}")                   

# Call path_to_isoseq_directory() to get the inputs
iso_seq_dir, iso_seq_file, iso_seq_file2 = path_to_isoseq_directory()
# /group/ctbrowngrp/finnolab_shared/phased_isoseq/
# phased.nopartial.cleaned.vcf
# phased.nopartial.cleaned.human_readable.txt
# Call the function with the obtained inputs
