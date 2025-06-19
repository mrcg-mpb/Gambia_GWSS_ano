import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


ibd_file = "3R_phased_data.ibd.gz"

meta_file = "sample_data.csv"
output_dir ="./IBD_results"
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']

ibd_df = pd.read_csv(ibd_file, sep='\t', names=columns, compression='gzip')


# 2. Calculate segment length
ibd_df['segment_length'] = ibd_df['END'] - ibd_df['START']

# 3. Define total genome length (adjust based on your organism and SNP coverage)
total_genome_length =  53_200_684

# 4. Group by individual pairs and sum total IBD shared length
ibd_summary = (
    ibd_df.groupby(['ID1', 'ID2'])['segment_length']
    .sum()
    .reset_index()
    .rename(columns={'segment_length': 'total_ibd_length'})
)

# 5. Compute proportion of genome shared IBD
ibd_summary['proportion_ibd'] = ibd_summary['total_ibd_length'] / total_genome_length

# 6. Show or save result
print(ibd_summary.head())
ibd_summary.to_csv("pairwise_ibd_proportion.csv", index=False)



# Load metadata
meta_df = pd.read_csv(meta_file)  # should contain columns like 'ID' and 'pop'

# Merge metadata to get population info for ID1 and ID2
ibd_summary = pd.read_csv("pairwise_ibd_proportion.csv")

# Merge for ID1
ibd_merged = ibd_summary.merge(meta_df[['ID', 'pop']], left_on='ID1', right_on='ID', how='left')
ibd_merged = ibd_merged.rename(columns={'pop': 'pop1'}).drop(columns=['ID'])

# Merge for ID2
ibd_merged = ibd_merged.merge(meta_df[['ID', 'pop']], left_on='ID2', right_on='ID', how='left')
ibd_merged = ibd_merged.rename(columns={'pop': 'pop2'}).drop(columns=['ID'])

# Now group by population pairs and compute mean proportion IBD
pop_ibd = (
    ibd_merged.groupby(['pop1', 'pop2'])['proportion_ibd']
    .mean()
    .reset_index()
)

# Pivot to matrix form for heatmap
ibd_matrix = pop_ibd.pivot(index='pop1', columns='pop2', values='proportion_ibd')

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(ibd_matrix, annot=True, fmt=".3f", cmap="viridis")
plt.title("Average Proportion IBD Between Populations")
plt.tight_layout()
plt.show()
