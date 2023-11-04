# FIRST RUN - not Sichong gtf
df = pd.read_table('./ASDEG/results/snpeff/snpeff_concatenated', sep = ' ', 
                          names=[ "Chr", "Pos", "ALT", "QUAL", "Effect", "Tissue"])

def parse_snpeff_annotation(ann_string):
    # Check if "ANN=" is in the string
    if "ANN=" not in ann_string:
        return []
    
    # Extracting the actual annotations after ANN=
    ann_string = ann_string.split("ANN=")[1]
    
    # Splitting by comma to get all annotations for the variant
    annotations = ann_string.split(",")
    parsed_data = []
    
    for annotation in annotations:
        # Splitting by pipe to get each field in the annotation
        fields = annotation.split("|")
        
        # Ensure there are enough fields before accessing them
        if len(fields) > 9:
            parsed_data.append({
                'Alternate Allele': fields[0],
                'Variant Type': fields[1],
                'Impact': fields[2],
                'Gene Symbol': fields[3],
                'Gene ID': fields[4],
                'Feature Type': fields[5],
                'Transcript ID': fields[6],
                'HGVS Notation': fields[9]
            })
    
    return parsed_data


# Parsing the 'Effect' column and adding new columns
df['Parsed_Effects'] = df['Effect'].apply(parse_snpeff_annotation)

# Exploding the DataFrame based on the parsed effects to create multiple rows for variants with multiple annotations
df_exploded = df.explode('Parsed_Effects')

# Splitting the dictionary in the 'Parsed_Effects' column into separate columns
df_exploded = pd.concat([df_exploded.drop('Parsed_Effects', axis=1), df_exploded['Parsed_Effects'].apply(pd.Series)], axis=1)

# Drop the 'Effect' column as it's no longer required
df_exploded.drop('Effect', axis=1, inplace=True)
df_exploded.to_csv('./ASDEG/results/snpeff/snpeff_concatenated_parsed')
