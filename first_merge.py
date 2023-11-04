expression_est_df = pd.read_csv('./ASDEG/results/samtools/Nuc_countup')

tIDs_df = pd.read_csv('./ASDEG/results/tIDs/merged_results.csv', sep='\t', 
                    names = ['Chr', 'Pos','ALT','QUAL', 'PBID','Tissue'])

Haps_df = pd.read_csv('./ASDEG/results/Haps/merged_results.csv', sep='\t', 
                    names = ['Haplotype', 'PBID','Tissue'])

merged_df = (tIDs_df.merge(Haps_df, how = 'inner', on=['PBID', 'Tissue'])
               .merge(expression_est_df, how = 'inner', on=['Chr','Pos','Tissue']))

# take out 'chr' from chr column 
merged_df['Chr'] = merged_df['Chr'].str.replace('chr', '')

merged_df = merged_df.drop_duplicates()

print(f"Total Number of Positions Analyzed : {((len(merged_df))+1)/2}")

# Get a list of columns to drop
columns_to_drop = [col for col in merged_df.columns if col.startswith('Unnamed')]
# Drop the columns
merged_df = merged_df.drop(columns=columns_to_drop)

merged_df.to_csv('./ASDEG/results/merged/first_merge')
