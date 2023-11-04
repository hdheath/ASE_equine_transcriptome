def extract_row_data(df, index):
    """Extract data from the DataFrame row at the given index."""
    row = df.iloc[index]
    return {
        
        'chrom': row['Chr'],
        'pos': row['Pos'],
        'pbid': row['PBID'],
        'tissue': row['Tissue'],
        'gene1': row['Gene Symbol'],
        'gene2': row['Gene Symbol.1'],
        'gene3': row['Gene ID'],
        'hap1': row['Haplotype'],
        'a1_exp': row['Expression'],
        'variant': row['Variant Type'],
        'impact': row['Impact'],
        'feature': row['Feature Type'],
        'tID': row['Transcript ID'],
        'HGVS': row['HGVS Notation']

    }

def append_new_data(current, next_row, new_data, index):
    """Check if genes and tissues match and then append data to new_data list."""
    if current['gene1'] == next_row['gene1'] and current['tissue'] == next_row['tissue']:
        new_data.append([current['chrom'], current['pos'], current['pbid'], current['tissue'],
                         current['gene1'], current['gene2'], current['gene3'],
                         current['hap1'], current['a1_exp'], next_row['hap1'],
                         next_row['a1_exp'], current['variant'], current['impact'],
                        current['feature'], current['tID'], current['HGVS']])
        return True  # Indicate that appending was successful
    
    #print(index)
    #print(f"Gene 1 : {current['gene']} Tissue 1 : {current['tissue']}")
    #print(f"Gene 2 : {next_row['gene']} Tissue 2 : {next_row['tissue']}")
    return False  # Indicate that appending was not successful

def process_dataframe(df):
    new_data = []
    index = 0
    while index < len(df) - 1:
        current_row_data = extract_row_data(df, index)
        next_row_data = extract_row_data(df, index + 1)

        # If append is successful, increase index by 2. Else, delete current row, reset index and increase by 1
        if append_new_data(current_row_data, next_row_data, new_data, index):
            index += 2
        else:
            df.drop(index, inplace=True)
            df.reset_index(drop=True, inplace=True)

    return new_data

# Read the CSV into a DataFrame
input_path = './ASDEG/results/merged/third_merge'
df = pd.read_csv(input_path)

columns_to_drop = [col for col in df.columns if col.startswith('Unnamed')]
df = df.drop(columns=columns_to_drop)

# Columns for the new DataFrame
new_columns = [ 'Chr', 'Pos', 'PBID', 'Tissue', 'Gene Symbol', 'Gene Symbol.1', 'Gene ID', 'Hap1',
               'A1 Exp', 'Hap2', 'A2 Exp', 'Variant Type', 'Impact', 'Feature Type', 'Transcript ID', 'HGVS Notation']

# Process the DataFrame and get new data
new_data = process_dataframe(df)

# Create the new DataFrame from the list of lists
df2 = pd.DataFrame(new_data, columns=new_columns)
df2.to_csv('./ASDEG/results/ASGs/All_ASGs')
print(f"Total Number of unique ASGs analyzed {len(df2)}")
