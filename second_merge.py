merged_df = pd.read_csv('./ASDEG/results/merged/first_merge')
snpeff_df = pd.read_csv('./ASDEG/results/snpeff/snpeff_concatenated_parsed')

# Merging with snpeff_df
merged_df['Pos'] = merged_df['Pos'].astype(str)
snpeff_df['Pos'] = snpeff_df['Pos'].astype(str)
merged_df['Chr'] = merged_df['Chr'].astype(str)
snpeff_df['Chr'] = snpeff_df['Chr'].astype(str)

merged_df2 = merged_df.merge(snpeff_df, how='inner', on=['Chr', 'Pos','ALT','QUAL','Tissue'])
# Get a list of columns to drop
columns_to_drop = [col for col in merged_df2.columns if col.startswith('Unnamed')]
merged_df2 = merged_df2.drop(columns=columns_to_drop)
merged_df2.to_csv('./ASDEG/results/merged/second_merge')

# Find out avg depth per loci 
import pandas as pd 

df = pd.read_csv('./ASDEG/results/merged/second_merge')
average_read_depth = df['Read Depth'].mean()

print(f"Average Read Depth Per Loci: {average_read_depth}")
