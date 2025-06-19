import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import seaborn as sns

# --- Load IBD segment data ---
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']
ibd_df = pd.read_csv("3R_phased_data.ibd.gz", sep='\t', names=columns, compression='gzip')

# Compute total IBD per pair (summing cM length)
pairwise_ibd = ibd_df.groupby(['ID1', 'ID2'])['CM_LEN'].sum().reset_index()

# Create a square matrix from pairwise IBD
individuals = sorted(set(pairwise_ibd['ID1']).union(pairwise_ibd['ID2']))
ibd_matrix = pd.DataFrame(0, index=individuals, columns=individuals)

for _, row in pairwise_ibd.iterrows():
    i, j, val = row['ID1'], row['ID2'], row['CM_LEN']
    ibd_matrix.at[i, j] = val
    ibd_matrix.at[j, i] = val  # symmetric

# Fill diagonal (self-self IBD)
np.fill_diagonal(ibd_matrix.values, ibd_matrix.max().max())

# Load population info
pop_df = pd.read_csv("sample_data.csv")  # columns: pop, taxon
pop_df = pop_df.set_index("pop")

# Match population labels to matrix order
pop_labels = pop_df.loc[ibd_matrix.index, "taxon"]

# Normalize: distance = 1 - (IBD / max)
distance_matrix = 1 - (ibd_matrix / ibd_matrix.max().max())

# Apply MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
coords = mds.fit_transform(distance_matrix)

# Assign colors
palette = sns.color_palette("tab10", len(set(pop_labels)))
color_map = dict(zip(set(pop_labels), palette))
colors = [color_map[label] for label in pop_labels]

# Plot
plt.figure(figsize=(10, 7))
for pop in set(pop_labels):
    idx = [i for i, p in enumerate(pop_labels) if p == pop]
    plt.scatter(coords[idx, 0], coords[idx, 1], label=pop, alpha=0.8, s=50)

plt.xlabel("MDS Dimension 1")
plt.ylabel("MDS Dimension 2")
plt.title("PCA-like MDS from IBD")
plt.legend(title="Population")
plt.grid(True)
plt.tight_layout()
# Save the plot
plt.savefig("mds_ibd_plot.png", dpi=300)






#---------------------------------------------------------------------# PCA of IBD Sharing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# --- Step 1: Load IBD data ---
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']
ibd_df = pd.read_csv("3R_phased_data.ibd.gz", sep='\t', names=columns, compression='gzip')

# --- Step 2: Compute pairwise IBD sum ---
pairwise_ibd = ibd_df.groupby(['ID1', 'ID2'])['CM_LEN'].sum().reset_index()

# --- Step 3: Create square IBD matrix ---
individuals = sorted(set(pairwise_ibd['ID1']) | set(pairwise_ibd['ID2']))
ibd_matrix = pd.DataFrame(0.0, index=individuals, columns=individuals)

for _, row in pairwise_ibd.iterrows():
    i, j, val = row['ID1'], row['ID2'], row['CM_LEN']
    ibd_matrix.at[i, j] = val
    ibd_matrix.at[j, i] = val  # Symmetric

# Optional: Fill diagonal (self-IBD) with max observed IBD
np.fill_diagonal(ibd_matrix.values, ibd_matrix.max().max())

# --- Step 4: Normalize and apply PCA ---
X = ibd_matrix / ibd_matrix.max().max()  # Optional normalization
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)

# --- Step 5: Prepare DataFrame for plotting ---
pca_df = pd.DataFrame({
    'PC1': principal_components[:, 0],
    'PC2': principal_components[:, 1],
    'Sample': ibd_matrix.index
})

# --- Step 6: Merge population labels ---
pop_df = pd.read_csv("sample_data.csv")  # Columns: pop, taxon
pop_map = pop_df.set_index("pop")["taxon"].to_dict()
pca_df['Population'] = pca_df['Sample'].map(pop_map)

# --- Step 7: Plot ---
plt.figure(figsize=(10, 7))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Population', palette='tab10', s=60)
plt.title('PCA of IBD Sharing')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
plt.legend(title='Population')
plt.grid(True)
plt.tight_layout()
# Save the plot
plt.savefig("ibd_pca_plot.png", dpi=300)
plt.close()