module load samtools
# Set the directory containing the files to be processed as the first argument
# Set the directory to search for matching files as the second argument
dir1=$1
dir2=$2

# Loop through all the files in the first directory
for file in $dir1/*
do
  # Extract the base name of the file
  name=$(basename $file)
 
  # Extract the name of the file up to the second underscore
  name=$(echo $name | awk -F'_' '{print $1"_"$2}')
  
  # Loop through all the files in the second directory
  for file2 in $dir2/*
  do
    # Extract the base name of the file
    name2=$(basename $file2)
    name2=$(echo $name2 | awk -F'.' '{print $1}')

    # Check if the base names match
    if [ "$name" = "$name2" ]
    then
      # Run samtools on each matching file 
      samtools mpileup -l "$dir1/$name""_extracted_lines_for_samtools"  -o "$dir1/$name""_samtools_results"  "$dir2/$name2"".markDup.bam"
      # Break out of the loop
      break
    fi
  done
  
done
