def process_group(group):

    ordered_group = group.sort_values(by='Haplotype')
    tissue = group['Tissue'].iloc[0]
    pbid = group['PBID'].iloc[0]
    gene = group['Gene Symbol'].iloc[0]

    haplotypes_data = []

    unique_haplotypes = ordered_group['Haplotype'].unique()
    for haplotype in unique_haplotypes:
        try:
            processed_data = process_haplotype(ordered_group, haplotype)
            haplotypes_data.append((pbid, tissue, gene, *processed_data))  
        except Exception as e:
            print(f"Error processing haplotype {haplotype} in group ({pbid}, {tissue}): {str(e)}")
            continue

    return haplotypes_data

def process_haplotype(ordered_group, haplotype):
    
    haplotype_rows = ordered_group[ordered_group['Haplotype'] == haplotype]
    
    values_by_position = [[] for _ in haplotype]

    for position, character in enumerate(haplotype):
        if character in haplotype_rows:
            values = haplotype_rows[character].values
            values_by_position[position].extend(values)

    express_value = sum(sum(values_by_position, []))
    express_value = express_value/(len(haplotype))

    effect_cols = ['Variant Type', 'Impact', 'Gene Symbol', 'Gene ID', 'Feature Type', 'Transcript ID', 'HGVS Notation']
    effect_data = {}
    for col in effect_cols:
        values = haplotype_rows[col].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else str(x)).values
        effect_data[col] = '   '.join(values)

    return haplotype_rows['Chr'].iloc[0], ' '.join(map(str, haplotype_rows['Pos'].values)), haplotype, express_value, *effect_data.values()

def main():
    df = pd.read_csv('./ASDEG/results/merged/second_merge')
    columns_to_drop = [col for col in df.columns if col.startswith('Unnamed')]
    df = df.drop(columns=columns_to_drop)

    grouped_df = df.groupby(['PBID', 'Tissue'], as_index=False)
    all_data = []

    for _, group in grouped_df:
        all_data.extend(process_group(group))

    result_columns = ['PBID', 'Tissue', 'Gene Symbol', 'Chr', 'Pos', 'Haplotype', 'Expression', 
                      'Variant Type', 'Impact', 'Gene Symbol', 'Gene ID', 'Feature Type', 'Transcript ID', 'HGVS Notation']
    
    result_df = pd.DataFrame(all_data, columns=result_columns)
    
    desired_order = ['Chr', 'Pos', 'PBID', 'Tissue', 'Gene Symbol', 'Gene ID', 'Haplotype', 'Expression','Variant Type', 'Impact', 'Feature Type', 'Transcript ID', 'HGVS Notation']
    result_df = result_df[desired_order]
    
    return result_df

result_df = main()
