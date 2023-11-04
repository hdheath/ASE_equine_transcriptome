def filter_data(df):
    """Filter rows where both A1 Exp and A2 Exp are not zero."""
    mask = (df['A1 Exp'] != 0) | (df['A2 Exp'] != 0)
    return df[mask].copy()

def adjust_values(df):
    """Adjust A1 Exp and A2 Exp values."""
    df.loc[df['A1 Exp'] == 0, 'A1 Exp'] = 1
    df.loc[df['A2 Exp'] == 0, 'A2 Exp'] = 1
    return df

def compute_log2_fold_change(df):
    """Calculate log2 fold change."""
    df['log2fold_change'] = np.log2(df['A1 Exp'] / df['A2 Exp'])
    return df

def apply_binomial_test(df):
    """Apply the binomial test and calculate related statistics."""
    results = df.apply(lambda x: binomtest(k=x['A1 Exp'], n=x['A1 Exp'] + x['A2 Exp']), axis=1)
    df['p-value'] = results.apply(lambda x: x.pvalue)
    df['conf_low'] = results.apply(lambda x: x.proportion_ci(confidence_level=0.95)[0])
    df['conf_high'] = results.apply(lambda x: x.proportion_ci(confidence_level=0.95)[1])
    return df

def adjust_p_values(df):
    """Adjust p-values using Benjamini-Hochberg procedure."""
    m = df.shape[0]  # total number of tests
    q = 0.05  # desired false discovery rate
    df['rank_by_p-value'] = np.arange(1, m + 1)
    df['adjusted_p-value'] = (df['rank_by_p-value'] / m) * q
    df['adjusted_p-value'] = np.minimum.accumulate(df['adjusted_p-value'][::-1])[::-1]
    return df

def process_dataframe(df):
    """Process the dataframe through the defined steps."""
    df = filter_data(df)
    df = adjust_values(df)
    df = compute_log2_fold_change(df)
    df = apply_binomial_test(df)
    df = adjust_p_values(df)
    
    # Rank by absolute value of log2 fold change
    df['rank_by_log2fold'] = df['log2fold_change'].abs().rank(ascending=False)
    
    # Sort DataFrame by absolute log2 fold change instead of p-value
    df.sort_values(by='log2fold_change', key=abs, ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    df = df.drop(columns='p-value')
    return df


# Read the CSV into a DataFrame
input_path = './ASDEG/results/ASGs/All_ASGs'
df = pd.read_csv(input_path)
df['A1 Exp'] = df['A1 Exp'].astype(int)
df['A2 Exp'] = df['A2 Exp'].astype(int)

# Process the DataFrame
subset_df = process_dataframe(df)

# Get a list of columns to drop
columns_to_drop = [col for col in subset_df.columns if col.startswith('Unnamed')]
# Drop the columns
subset_df = subset_df.drop(columns=columns_to_drop)
