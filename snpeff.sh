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
    
    # Echo the names before checking for matching
    echo "Processing file in directory 1: $name"
    echo "Checking against file in directory 2: $name2"

    # Check if the base names match
    if [ "$name" = "$name2" ]
    then
      echo "Match found: $name matches $name2"
      # Run snpeff on each matching file 
      cd snpEff
      echo "Running snpEff on: $dir1/$name""_extracted_lines_for_snpeff"
      snpEff ann EquCab3.0.105  "$dir1/$name""_extracted_lines_for_snpeff" > "$dir1/$name""_extracted_lines_snpeff.ran" -no-intron -no-intergenic
      cd
      # Break out of the loop
      break
    fi
  done
  
done
