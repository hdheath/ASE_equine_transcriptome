# Read depth plot 
import pandas as pd 
import matplotlib.pyplot as plt 

# Read the data
df = pd.read_csv('./ASDEG/results/merged/second_merge')

# Calculate average read depth
average_read_depth = df['Read Depth'].mean()
print(f"Average Read Depth Per Loci: {average_read_depth}")
# Set the x-axis limits
plt.xlim(result_df['Haplotype Length'].min(), average_haplotype_length * 500)
# Plot histogram for Read Depth distribution
plt.hist(df['Read Depth'], bins=50, edgecolor='black') # Adjust the bins as needed
plt.axvline(average_read_depth, color='red', linestyle='dashed', linewidth=1) # Draw a vertical line to represent the average read depth
plt.title("Distribution of Read Depth")
plt.xlabel("Read Depth")
plt.ylabel("Frequency")
plt.legend(["Average Read Depth", "Read Depth"])
plt.show()

# Distribution of haplotype lengths 

# Read the data
result_df = pd.read_csv('./ASDEG/results/ASGs/All_ASGs')

# Calculate the length of each haplotype
result_df['Haplotype Length'] = result_df['Hap1'].apply(len)

# Count how many haplotypes have a length of 1
count_length_1 = len(result_df[result_df['Haplotype Length'] == 1])

# Print the count
print(f"Haplotypes with Length 1: {count_length_1}")

# Calculate average haplotype length
average_haplotype_length = result_df['Haplotype Length'].mean()
print(f"Average Haplotype Length: {average_haplotype_length:.2f}")

# Set the x-axis limits
plt.xlim(result_df['Haplotype Length'].min(), average_haplotype_length * 10)

# Plot histogram for Haplotype Length distribution
plt.hist(result_df['Haplotype Length'], bins=150, edgecolor='black')  # Adjust the bins as needed
plt.axvline(average_haplotype_length, color='red', linestyle='dashed', linewidth=1)  # Draw a vertical line to represent the average haplotype length
plt.title("Distribution of Haplotype Lengths")
plt.xlabel("Haplotype Length")
plt.ylabel("Frequency")
plt

# Volcano Plot 
import pandas as pd
import numpy as np
import plotly.express as px

def calculate_log_adjusted_p_value(df):
    """Calculate log(adjusted_p-value) and handle infinite values."""
    return -np.log10(df['adjusted_p-value'].replace([np.inf, -np.inf], 100))

def assign_color_based_on_threshold(df, positive_threshold, negative_threshold):
    """Assign colors based on log2fold_change thresholds."""
    conditions = [
        df['log2fold_change'] >= positive_threshold,
        df['log2fold_change'] <= negative_threshold
    ]
    
    choices = ['Up-Regulated', 'Down Regulated']
    
    return np.select(conditions, choices, default='Insignificant')

def filter_data_based_on_thresholds(df, positive_threshold, negative_threshold):
    """Filter data based on positive and negative thresholds."""
    positive_df = df[df['log2fold_change'] >= positive_threshold]
    negative_df = df[df['log2fold_change'] <= negative_threshold]
    gray_df = df[(df['log2fold_change'] > negative_threshold) & (df['log2fold_change'] < positive_threshold)]
    
    return positive_df, negative_df, gray_df

def plot_interactive_scatter(df, positive_threshold, negative_threshold):
    """Generate and display the interactive scatter plot."""
    fig = px.scatter(df, x='log2fold_change', y='log_adjusted_p_value', color='Color', text='rank_by_log2fold')
    
    positive_df, negative_df, gray_df = filter_data_based_on_thresholds(df, positive_threshold, negative_threshold)
    
    for data_frame in [positive_df, negative_df, gray_df]:
        fig.add_trace(px.scatter(data_frame, x='log2fold_change', y='log_adjusted_p_value', 
                                 text='rank_by_log2fold', opacity=0.5, custom_data=['rank_by_log2fold']).data[0])
    
    fig.update_layout(
        title='Interactive Scatter Plot of log2fold_change vs. log(adjusted_p-value)',
        xaxis_title='log2fold_change',
        yaxis_title='log(adjusted_p-value)',
        hovermode='closest'
    )
    
    fig.update_traces(
        hovertemplate='<b>Rank by log2fold:</b> %{customdata[0]}<br><b>log2fold_change:</b> %{x}<br><b>log(adjusted_p-value):</b> %{y}',
        selector=dict(type='scatter')
    )
    
    fig.show()


df = pd.read_csv('./ASDEG/results/ASGs/ranked_ASGs')

df = df.head(50)

# Assuming your DataFrame is named 'df'
df['log_adjusted_p_value'] = calculate_log_adjusted_p_value(df)

# Define thresholds
positive_threshold = 2.0
negative_threshold = -2.0

df['Color'] = assign_color_based_on_threshold(df, positive_threshold, negative_threshold)

plot_interactive_scatter(df, positive_threshold, negative_threshold)

# Distribution of AE analysis 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('./ASDEG/results/ASGs/ranked_ASGs')

# Extract the value after '_' from the 'Tissue' column, abbreviate if longer than 7 characters, 
# and create a new column 'Generic Tissue'
def extract_tissue(x):
    tissue = x.split('_')[1] if '_' in x else x
    if tissue == "LongMuscle":
        return "Longiss"
    if tissue == 'LeftVentricle':
        return "Heart"
    if tissue == 'ParietalCortex':
        return "PC"
    return tissue[:7]

df['Generic Tissue'] = df['Tissue'].apply(extract_tissue)

# Create a mask for the conditions
mask = (
    (df['log2fold_change'].abs() >= 2) &
    (df['A1 Exp'] + df['A2 Exp'] >= 50) &
    (df['adjusted_p-value'] <= .05)
)

# Get the counts of 'Generic Tissue' for the original and filtered dataframes
total_counts = df['Generic Tissue'].value_counts()
filtered_counts = df[mask]['Generic Tissue'].value_counts()

# Plotting
plt.figure(figsize=(10,6))
bars = sns.barplot(x=total_counts.index, y=total_counts.values, color='lightblue', label='Analayzed Genes')
bottom_bars = sns.barplot(x=filtered_counts.index, y=filtered_counts.values, color='blue', label='ASDEG')

# Adding legends and labels
plt.legend()
plt.ylabel('Count')
plt.xlabel('Tissue')
plt.title('Distribution of Analyzed Genes')
plt.show()

# Dsitirbution of ASDEGs
# Extract the value after '_' from the 'Tissue' column, abbreviate if longer than 7 characters, 
# and create a new column 'Generic Tissue'
def extract_tissue(x):
    tissue = x.split('_')[1] if '_' in x else x
    if tissue == "LongMuscle":
        return "Longiss"
    if tissue == 'LeftVentricle':
        return "Heart"
    if tissue == 'ParietalCortex':
        return "PC"
    
    return tissue[:7]

df['Generic Tissue'] = df['Tissue'].apply(extract_tissue)

# Get the counts of 'Generic Tissue' for the original and filtered dataframes
total_counts = df['Generic Tissue'].value_counts()

# Plotting
plt.figure(figsize=(10,6))
bars = sns.barplot(x=total_counts.index, y=total_counts.values, color='lightblue', label='ASDEGs')


# Adding legends and labels
plt.legend()
plt.ylabel('Count')
plt.xlabel('Tissue')
plt.title('Distribution of ASDEGs')
plt.show()

# Distribution of Variants in ASDEGs
# Split and flatten the lists
variant_types = df['Variant Type'].explode().str.split('&').explode()

# Remove unwanted characters
cleaned_variant_types = variant_types.str.replace(r"[\[\]\,']", "", regex=True)

# Split values by the space between annotations
split_variant_types = cleaned_variant_types.str.split().explode()

# Replace underscores with spaces AFTER exploding
split_variant_types = split_variant_types.str.replace("_", " ")

# Filter values that have counts greater than 10
final_filtered_variant_types = split_variant_types[split_variant_types.isin(split_variant_types.value_counts()[split_variant_types.value_counts() > 5].index)]

value_counts = final_filtered_variant_types.value_counts()
percentages = (value_counts / len(final_filtered_variant_types)) * 100

# Plotting
sns.countplot(x=final_filtered_variant_types, order=final_filtered_variant_types.value_counts().index)
plt.title('Variant Types Within ASDEGs')
plt.xlabel('Variant Type')
plt.ylabel('Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Variant Distribution 
# Split and flatten the lists
variant_types = df['Variant Type'].explode().str.split('&').explode()

# Remove unwanted characters
cleaned_variant_types = variant_types.str.replace(r"[\[\]']", "", regex=True)

# Split values by the space between annotations
split_variant_types = cleaned_variant_types.str.split().explode()

# Replace underscores with spaces AFTER exploding
split_variant_types = split_variant_types.str.replace("_", " ")

# Filter values that have counts greater than 10
final_filtered_variant_types = split_variant_types[split_variant_types.isin(split_variant_types.value_counts()[split_variant_types.value_counts() > 1].index)]

value_counts = final_filtered_variant_types.value_counts()
percentages = (value_counts / len(final_filtered_variant_types)) * 100

# Plotting
sns.countplot(x=final_filtered_variant_types, order=final_filtered_variant_types.value_counts().index)
plt.title('Variant Types Across All Genes')
plt.xlabel('Variant Type')
plt.ylabel('Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Calculate average log2 fold change
# Use a nicer plot style
plt.style.use('seaborn-darkgrid')

# Increase the size of the figure
plt.figure(figsize=(10,6))
df['log2fold_change'] = df['log2fold_change'].abs()
average_log2_foldchange = df['log2fold_change'].mean()
print(f"Average Log2 Fold Change: {average_log2_foldchange}")

# Set the x-axis limits (I assume you have a similar variable as average_haplotype_length, adjust as necessary)
plt.xlim(df['log2fold_change'].min(), average_log2_foldchange * 8)  # Adjust the multiplier if needed

# Plot histogram for log2 fold change distribution
plt.hist(df['log2fold_change'], bins=15, edgecolor='black')  # Adjust the bins as needed
plt.axvline(average_log2_foldchange, color='red', linestyle='dashed', linewidth=1)  # Draw a vertical line to represent the average log2 fold change
plt.title("Distribution of Allelic Expression Fold Change")
plt.xlabel("aeFC")
plt.ylabel("Frequency")
plt.legend(["Average Allelic Expression Fold Change", "Allelic Expression Fold Change"])
plt.show()

# Tissue Heat Map 
# Split "Gene Symbol" and explode it to multiple rows
df = df.assign(**{'Gene Symbol': df['Gene Symbol'].str.split()}).explode('Gene Symbol')
df['Gene Symbol'] = df['Gene Symbol'].str.replace(',', '')
df = df[df['Gene Symbol'] != 'nan']

# Update values
df.loc[df['Tissue'] == 'LongMuscle', 'Tissue'] = 'Longissimus'
df.loc[df['Tissue'] == 'LeftVentricle', 'Tissue'] = 'Heart'
df.loc[df['Tissue'] == 'AdiposeLoin', 'Tissue'] = 'Adipose'

# Create a binary matrix
binary_matrix = pd.crosstab(df['Tissue'], df['Gene Symbol'])
binary_matrix = binary_matrix.applymap(lambda x: 1 if x > 0 else 0)  # Convert to 1s and 0s based on gene presence/absence

# Compute Jaccard similarity
jaccard_similarity = pdist(binary_matrix.values, metric='jaccard')
jaccard_similarity = squareform(jaccard_similarity)
jaccard_similarity = 1 - jaccard_similarity  # Convert distance to similarity

# Plot heatmap using seaborn
plt.figure(figsize=(10, 6))
sns.heatmap(jaccard_similarity, annot=True, cmap='viridis', cbar=True, xticklabels=binary_matrix.index, yticklabels=binary_matrix.index)
plt.title("ASDEG Expression Across Tissues")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Binary Tissue Matrix 
import pandas as pd
df = df.dropna(subset=['Gene Symbol'])
# Split "Gene Symbol" and explode it to multiple rows
df = df.assign(**{'Gene Symbol': df['Gene Symbol'].str.split()}).explode('Gene Symbol')
df = df[df['Gene Symbol'] != 'nan']
# Create a binary matrix
binary_matrix = pd.crosstab(df['Gene Symbol'], df['Generic Tissue'])

# Count the number of tissues each gene is expressed in
gene_counts = binary_matrix.sum(axis=1)

# Get top 10 genes that are expressed in the most tissues
top_genes = gene_counts.nlargest(30).index

# Extract the rows corresponding to the top genes from the binary matrix
top_genes_matrix = binary_matrix.loc[top_genes]

print(top_genes_matrix)
