#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, pdist
from sklearn.decomposition import PCA
import matplotlib.patches as mpatches


# ------------------------
# 1. Load IBD and Metadata
# ------------------------
chr = "3R"
len_chr = {
    "2R": 61545105,
    "2L": 49364325,
    "3R": 53200684,
    "3L": 41963435,
    "X": 24393108,
} 
ibd_file = f"{chr}_phased_data.ibd.gz"

meta_file = "sample_data.csv"
output_dir ="./IBD_results"
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']

ibd = pd.read_csv(ibd_file, sep='\t', names=columns, compression='gzip')
meta = pd.read_csv(meta_file).set_index("pop")
# Proportion of IBD segments shared between populations

# Convert total shared IBD cM to proportion shared

# Convert to centiMorgans assuming 1 cM = 1 Mb (1,000,000 bp)
len_chr_cM = {chr: bp / 1_000_000 for chr, bp in len_chr.items()}

ibd['proportion'] = ibd['CM_LEN'] / len_chr_cM[chr]  # Proportion of genome shared

print(ibd)
# Filter strong segments
#ibd = ibd[(ibd['LOD'] > 3) & (ibd['CM_LEN'] > 0.5)].copy()

# Normalize pair names
ibd["pair"] = ibd.apply(lambda row: tuple(sorted([row["ID1"], row["ID2"]])), axis=1)

# Sum total IBD per unique pair
pair_ibd = ibd.groupby("pair")["CM_LEN"].sum().reset_index()  # Try after uing the mean IBD per unique population
pair_ibd[["ID1", "ID2"]] = pd.DataFrame(pair_ibd["pair"].tolist(), index=pair_ibd.index)
pair_ibd.drop(columns="pair", inplace=True)

#print("Heatmap: ",pair_ibd)

# ------------------------
# 2. Create Symmetric IBD Matrix (individual level)
# ------------------------
all_ids = sorted(set(pair_ibd["ID1"]).union(pair_ibd["ID2"]))
ibd_matrix = pd.DataFrame(0.0, index=all_ids, columns=all_ids)
for _, row in pair_ibd.iterrows():
    ibd_matrix.loc[row["ID1"], row["ID2"]] = row["CM_LEN"]
    ibd_matrix.loc[row["ID2"], row["ID1"]] = row["CM_LEN"]

# ------------------------
# 3. Annotate Individuals by Population
# ------------------------
id_to_pop = meta["taxon"].to_dict()
pop_labels = [id_to_pop.get(i, "unknown") for i in ibd_matrix.index]

#print(ibd_matrix)
# ------------------------
# 4. Plot Annotated Heatmap with Clustering
# ------------------------
row_colors = pd.Series(pop_labels, index=ibd_matrix.index).map({
    "gambiae": "blue",
    "coluzzii": "orange",
    "bissau": "purple",
    "arabiensis": "green",
    "melas": "red",
})

# Create a color legend with enhanced formatting
legend_label_map = {
    "gambiae": r"$\boldsymbol{An.\ gambiae\ s.s}$",
    "coluzzii": r"$\boldsymbol{An.\ coluzzii}$", 
    "bissau": r"$\mathbf{Bissau}$",
    "arabiensis": r"$\boldsymbol{An.\ arabiensis}$",
    "melas": r"$\boldsymbol{An.\ melas}$",
}

legend_handles = [
    mpatches.Patch(color=color, label=legend_label_map.get(label, label))
    for label, color in {
        "gambiae": "blue",
        "coluzzii": "orange",
        "bissau": "purple",
        "arabiensis": "green",
        "melas": "red",
    }.items()
]

g = sns.clustermap(
    ibd_matrix,
    row_cluster=True,
    col_cluster=True,
    row_colors=row_colors,
    col_colors=row_colors,
    cmap="YlOrRd",
    figsize=(12, 12),
    xticklabels=False,
    yticklabels=False,
    method="ward",
    vmin=0,
    vmax=3,
)

# Add enhanced legend
legend = g.ax_heatmap.legend(handles=legend_handles, title='Population', 
                            bbox_to_anchor=(1.05, 1), loc='upper left',
                            fontsize=18, title_fontsize=20)
plt.setp(legend.get_title(), fontweight='bold')

g.fig.suptitle(rf"Pairwise IBD Sharing Heatmap of $\boldsymbol{{Anopheles\ gambiae}}$ complex on Chromosome {chr}", 
               y=1.02, fontsize=16, fontweight='bold')


g.fig.savefig(f"{output_dir}/ibd_annotated_heatmap_{chr}.png", dpi=600, bbox_inches='tight')
g.fig.savefig(f"{output_dir}/ibd_annotated_heatmap_{chr}.svg", bbox_inches='tight', format='svg')



# ------------------------
# 7. IBD Segment Length Distribution
# ------------------------
ibd_df = ibd.copy()
# Add segment length (Mb or cM if available)
ibd_df["length"] = ibd_df["END"] - ibd_df["START"]
ibd_df["length_mb"] = ibd_df["length"] / 1e6
# Load metadata
sample_metadata = pd.read_csv("sample_data.csv")
# create a dictionary to map sample IDs to populations
mapping_dict = sample_metadata.set_index('pop')['taxon'].to_dict()

# Annotate each individual with their population
ibd_df["pop1"] = ibd_df["ID1"].map(mapping_dict)
ibd_df["pop2"] = ibd_df["ID2"].map(mapping_dict)

# Define comparison type
ibd_df["comparison"] = ibd_df.apply(
    lambda row: "within" if row["pop1"] == row["pop2"] else "between", axis=1
)

# Total shared IBD per pair (optional aggregation)
ibd_summary = ibd_df.groupby(["ID1", "ID2", "comparison"])["length_mb"].sum().reset_index()

# Filter for between-population IBD
between_ibd = ibd_df[ibd_df["comparison"] == "between"]
between_ibd_df = between_ibd[["ID1", "ID2", "START", "END", "length_mb"]]

# Summary of the CM_LEN column
min_cm = ibd_df["CM_LEN"].min()
max_cm = ibd_df["CM_LEN"].max()
mean_cm = ibd_df["CM_LEN"].mean()
median_cm = ibd_df["CM_LEN"].median()

print(f"IBD Segment Lengths (cM):")
print(f" - Min:    {min_cm:.4f} cM")
print(f" - Max:    {max_cm:.4f} cM")
print(f" - Mean:   {mean_cm:.4f} cM")
print(f" - Median: {median_cm:.4f} cM")

plt.figure(figsize=(10, 5))
sns.histplot(ibd_df["CM_LEN"], bins=100, kde=True)
plt.title(rf"Distribution of IBD Segment Lengths on Chromosome {chr}", 
          fontsize=12, fontweight='bold', pad=20)
plt.xlabel("IBD Segment Length (cM)", fontsize=12, fontweight='bold')
plt.ylabel("Count", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_segment_lengths_distribution_{chr}.png", dpi=300)

# ------------------------
# NEW: Population-level IBD sharing matrix
# ------------------------
# Map taxon/population
ibd_df = ibd.copy()
ibd_df['present'] = ibd_df['CM_LEN'].apply(lambda x: 1 if x > 0 else 0)
ibd_df["pop1"] = ibd_df["ID1"].map(meta["taxon"])
ibd_df["pop2"] = ibd_df["ID2"].map(meta["taxon"])
ibd_df = ibd_df.dropna(subset=["pop1", "pop2"])
ibd_df["comparison"] = ibd_df.apply(lambda row: "within" if row["pop1"] == row["pop2"] else "between", axis=1)

pop_ibd = ibd_df.copy()
pop_ibd["pop_pair"] = pop_ibd.apply(lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1)

pop_sharing = pop_ibd.groupby("pop_pair")["CM_LEN"].sum().reset_index()

# Pivot to square matrix
pop_names = sorted(set(pop_ibd["pop1"]).union(pop_ibd["pop2"]))
pop_matrix = pd.DataFrame(0, index=pop_names, columns=pop_names, dtype=float)

for row in pop_sharing.itertuples(index=False):
    p1, p2, cm = row.pop_pair[0], row.pop_pair[1], row.CM_LEN
    pop_matrix.loc[p1, p2] += cm
    if p1 != p2:
        pop_matrix.loc[p2, p1] += cm

# Plot heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(pop_matrix, annot=True, fmt=".1f", cmap="YlGnBu")
plt.title(rf"Total IBD Sharing (cM) Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=12, fontweight='bold', pad=20)
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.ylabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_population_heatmap_{chr}.png")

pop_matrix_prop = pop_matrix / len_chr_cM[chr]  # Divide by chromosome length
print("len_chr_cM: ", len_chr_cM)

# Plot proportion heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(pop_matrix_prop, annot=True, fmt=".3f", cmap="YlGnBu")
plt.title(rf"Proportion of Genome Shared IBD Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=14, fontweight='bold', pad=20)
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.ylabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_population_proportion_heatmap_{chr}.png")

ibd_df["length"] = ibd_df["END"] - ibd_df["START"]
ibd_df["length_mb"] = ibd_df["length"] / 1e6
# Compute total IBD per sample pair
pairwise_ibd = ibd_df.groupby(["ID1", "ID2", "pop1", "pop2"])["length_mb"].sum().reset_index()
# Create sorted population pair (for symmetry)
pairwise_ibd["pop_pair"] = pairwise_ibd.apply(
    lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1
)

# Split into pop1/pop2 again
pairwise_ibd["pop1_sorted"] = pairwise_ibd["pop_pair"].apply(lambda x: x[0])
pairwise_ibd["pop2_sorted"] = pairwise_ibd["pop_pair"].apply(lambda x: x[1])

# Average IBD sharing between population pairs
heatmap_data = pairwise_ibd.groupby(["pop1_sorted", "pop2_sorted"])["length_mb"].mean().reset_index()

# Pivot to square matrix
ibd_matrix = heatmap_data.pivot(index="pop1_sorted", columns="pop2_sorted", values="length_mb")

# Fill symmetric half (optional)
ibd_matrix = ibd_matrix.combine_first(ibd_matrix.T)
plt.figure(figsize=(10, 8))
sns.heatmap(ibd_matrix, annot=True, fmt=".2f", cmap="YlGnBu", square=True, linewidths=0.5)
plt.title(rf"Mean Total IBD Sharing (Mb) Within and Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=14, fontweight='bold', pad=20)
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.ylabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_heatmap_share_pop_{chr}.png", dpi=300)

# ------------------------
# 8. Gene Annotation Overlap
# ------------------------

# --- Load gene annotations global ---
gff_path = "Anopheles_gambiae.AgamP4.60.chr.gff3.gz"
gene_data = []
with gzip.open(gff_path, 'rt') as gff:
    for line in gff:
        if line.startswith("#"): continue
        f = line.strip().split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        attr = dict(i.split("=") for i in f[8].split(";") if "=" in i)
        gene_data.append({"chrom": f[0], "start": int(f[3]), "end": int(f[4]), "gene": attr.get("Name", attr.get("ID", "unknown"))})
genes_df = pd.DataFrame(gene_data)

# --- Define region of interest ---
top_bin_start = ibd_df["START"].min()
top_bin_end = ibd_df["END"].max()

# --- Filter between population IBD segments ---
ibd_df_copy = ibd_df[ibd_df["comparison"] == "between"]
ibd_df_copy['pair'] = ibd_df_copy.apply(lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1)

print("Population pairs in between IBD segments:")
print(ibd_df_copy['pair'].unique())  

# --- Filter IBD segments in that region ---
highlight_ibd = ibd_df_copy[
    (ibd_df_copy['START'] <= top_bin_end) & 
    (ibd_df_copy['END'] >= top_bin_start)
].sort_values("START").copy()
highlight_ibd['y'] = range(len(highlight_ibd))

# --- Create population pair label ---
highlight_ibd["pair_type"] = highlight_ibd.apply(
    lambda row: f"{row['pop1']}-{row['pop2']}" if row["pop1"] <= row["pop2"] else f"{row['pop2']}-{row['pop1']}",
    axis=1
)

# --- Define color palette for pair types ---
pair_colors = {
    "bissau-gambiae": "red",
    "coluzzii-gambiae": "blue",
    "bissau-coluzzii": "green",
}

highlight_ibd["color"] = highlight_ibd["pair_type"].apply(lambda x: pair_colors.get(x, "gray"))

# --- Plot ---
plt.figure(figsize=(16, 8))
for _, row in highlight_ibd.iterrows():
    plt.plot([row['START'], row['END']], [row['y'], row['y']], color=row['color'], lw=3)

if chr == '3R':
    plt.axvspan(32000000, 32080000, color='yellow', alpha=0.5)
    plt.axvspan(28500000, 28800000, color='yellow', alpha=0.5)

# Create custom label mapping for between-population pairs
between_pair_label_map = {
    "bissau-gambiae": r"$\mathbf{Bissau}$-$\boldsymbol{An.\ gambiae\ s.s}$",
    "coluzzii-gambiae": r"$\boldsymbol{An.\ coluzzii}$-$\boldsymbol{An.\ gambiae\ s.s}$", 
    "bissau-coluzzii": r"$\mathbf{Bissau}$-$\boldsymbol{An.\ coluzzii}$",
}

handles = [plt.Line2D([0], [0], color=color, lw=3, label=between_pair_label_map.get(label, label)) 
           for label, color in pair_colors.items()]
legend = plt.legend(handles=handles, title="Population Pair", bbox_to_anchor=(1.05, 1), 
                   loc='upper left', fontsize=18, title_fontsize=20)
plt.setp(legend.get_title(), fontweight='bold')

plt.title(rf"IBD Segments by Population Pair on Chromosome {chr}", 
          fontsize=20, fontweight='bold', pad=20)
plt.xlabel("Genomic Position (bp)", fontsize=14, fontweight='bold')
plt.ylabel("IBD Segment (sample pair)", fontsize=14, fontweight='bold')
# Make tick labels bold
plt.xticks(fontsize=14, fontweight='bold')
plt.yticks(fontsize=14, fontweight='bold')

# Make axes lines thicker
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(3)
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_top_region_between_pop_{chr}.png", dpi=600)

# List genes in that region
genes_df = genes_df.query('chrom =="3R"')[(genes_df['start'] >= 32000000) & (genes_df['end'] <= 32100000)]

# --- Plot within population segments ---
# --- Filter within population IBD segments ---
ibd_df_copy = ibd_df[ibd_df["comparison"] == "within"]
ibd_df_copy['pair'] = ibd_df_copy.apply(lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1)

# --- Filter IBD segments in that region ---
highlight_ibd = ibd_df_copy[
    (ibd_df_copy['START'] <= top_bin_end) & 
    (ibd_df_copy['END'] >= top_bin_start)
].sort_values("START").copy()
highlight_ibd['y'] = range(len(highlight_ibd))

# --- Create population pair label ---
highlight_ibd["pair_type"] = highlight_ibd.apply(
    lambda row: f"{row['pop1']}-{row['pop2']}" if row["pop1"] <= row["pop2"] else f"{row['pop2']}-{row['pop1']}",
    axis=1
)

# --- Define color palette for pair types ---
pair_colors = {
    "gambiae-gambiae": "blue",
    "coluzzii-coluzzii": "orange",
    "bissau-bissau": "purple",
    "arabiensis-arabiensis": "green",
    "melas-melas": "red",
}

highlight_ibd["color"] = highlight_ibd["pair_type"].apply(lambda x: pair_colors.get(x, "gray"))

# --- Plot ---
plt.figure(figsize=(14, 6))
for _, row in highlight_ibd.iterrows():
    plt.plot([row['START'], row['END']], [row['y'], row['y']], color=row['color'], lw=2)

# Create custom label mapping for within-population pairs
within_pair_label_map = {
    "gambiae-gambiae": r"$\boldsymbol{An.\ gambiae}$ s.s",
    "coluzzii-coluzzii": r"$\boldsymbol{An.\ coluzzii}$",
    "bissau-bissau": r"$\mathbf{Bissau}$",
    "arabiensis-arabiensis": r"$\boldsymbol{An.\ arabiensis}$",
    "melas-melas": r"$\boldsymbol{An.\ melas}$",
}

# Optional: Add legend with enhanced formatting
handles = [plt.Line2D([0], [0], color=color, lw=3, label=within_pair_label_map.get(label, label)) 
           for label, color in pair_colors.items()]
legend = plt.legend(handles=handles, title="Population Pair", bbox_to_anchor=(1.05, 1), 
                   loc='upper left', fontsize=12, title_fontsize=14)
plt.setp(legend.get_title(), fontweight='bold')

plt.title(rf"IBD Segments in Region Colored by Population Pair on Chromosome {chr}", 
          fontsize=16, fontweight='bold', pad=20)
plt.xlabel("Genomic Position (bp)", fontsize=12, fontweight='bold')
plt.ylabel("IBD Segment (sample pair)", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_top_region_within_pop_{chr}.png", dpi=600)
plt.savefig(f"{output_dir}/ibd_top_region_within_pop_{chr}.svg", bbox_inches='tight', format='svg')

# ------------------------
# Comparison plots
# ------------------------

# Load the uploaded IBD file
columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']
ibd_df = pd.read_csv(ibd_file, sep='\t', header=None, names=columns, compression='gzip')

# Add segment length (Mb or cM if available)
ibd_df["length"] = ibd_df["END"] - ibd_df["START"]
ibd_df["length_mb"] = ibd_df["length"] / 1e6

# Load metadata
sample_metadata = pd.read_csv("sample_data.csv")
# create a dictionary to map sample IDs to populations
mapping_dict = sample_metadata.set_index('pop')['taxon'].to_dict()

# Annotate each individual with their population
ibd_df["pop1"] = ibd_df["ID1"].map(mapping_dict)
ibd_df["pop2"] = ibd_df["ID2"].map(mapping_dict)

# Define comparison type
ibd_df["comparison"] = ibd_df.apply(
    lambda row: "within" if row["pop1"] == row["pop2"] else "between", axis=1
)

# Total shared IBD per pair (optional aggregation)
ibd_summary = ibd_df.groupby(["ID1", "ID2", "comparison"])["length_mb"].sum().reset_index()

# Then group again to get mean within/between
comparison_summary = ibd_summary.groupby("comparison")["length_mb"].describe()

plt.figure(figsize=(8, 6))
sns.violinplot(data=ibd_summary, x="comparison", y="length_mb", inner="box", palette="pastel")
plt.title(rf"IBD Sharing: Within vs. Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=16, fontweight='bold', pad=20)
plt.ylabel("Total IBD (Mb) per Pair", fontsize=12, fontweight='bold')
plt.xlabel("Comparison Type", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_sharing_comparison_{chr}.png", dpi=300)

# Filter for within-population IBD
within_ibd = ibd_df[ibd_df["comparison"] == "within"]
within_ibd["pop"] = within_ibd["pop1"]

# Summarize
within_summary = within_ibd.groupby("pop")["length_mb"].sum().reset_index()

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=within_summary, x="pop", y="length_mb")
plt.title(rf"Total Within-Population IBD on Chromosome {chr}", 
          fontsize=16, fontweight='bold', pad=20)
plt.ylabel("Sum of IBD (Mb)", fontsize=12, fontweight='bold')
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/within_population_ibd_{chr}.png", dpi=300)

# Compute total IBD per sample pair
pairwise_ibd = ibd_df.groupby(["ID1", "ID2", "pop1", "pop2"])["length_mb"].sum().reset_index()
# Create sorted population pair (for symmetry)
pairwise_ibd["pop_pair"] = pairwise_ibd.apply(
    lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1
)

# Split into pop1/pop2 again
pairwise_ibd["pop1_sorted"] = pairwise_ibd["pop_pair"].apply(lambda x: x[0])
pairwise_ibd["pop2_sorted"] = pairwise_ibd["pop_pair"].apply(lambda x: x[1])

# Average IBD sharing between population pairs
heatmap_data = pairwise_ibd.groupby(["pop1_sorted", "pop2_sorted"])["length_mb"].mean().reset_index()

# Pivot to square matrix
ibd_matrix = heatmap_data.pivot(index="pop1_sorted", columns="pop2_sorted", values="length_mb")

# Fill symmetric half (optional)
ibd_matrix = ibd_matrix.combine_first(ibd_matrix.T)
plt.figure(figsize=(10, 8))
sns.heatmap(ibd_matrix, annot=True, fmt=".2f", cmap="YlGnBu", square=True, linewidths=0.5)
plt.title(rf"Mean Total IBD Sharing (Mb) Within and Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=14, fontweight='bold', pad=20)
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.ylabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/ibd_heatmap_share_pop_{chr}.png", dpi=300)

# --- Compute proportion of genome shared IBD between populations ---

# Define population pair as sorted tuple
ibd_df["pop_pair"] = ibd_df.apply(
    lambda row: tuple(sorted([row["pop1"], row["pop2"]])), axis=1
)

# Group by population pair and sum CM_LEN
pop_ibd_sum = ibd_df.groupby("pop_pair")["CM_LEN"].sum().reset_index()

# Convert chromosome length to cM
chr_len_cM = len_chr[chr] / 1_000_000

# Compute proportion of genome shared
pop_ibd_sum["proportion_ibd"] = pop_ibd_sum["CM_LEN"] / chr_len_cM

# Create symmetric matrix
pops = sorted(set(sum(pop_ibd_sum["pop_pair"].tolist(), ())))
prop_matrix = pd.DataFrame(0.0, index=pops, columns=pops)

for _, row in pop_ibd_sum.iterrows():
    p1, p2 = row["pop_pair"]
    prop_matrix.loc[p1, p2] = row["proportion_ibd"]
    prop_matrix.loc[p2, p1] = row["proportion_ibd"]

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(prop_matrix, cmap="YlGnBu", annot=True, fmt=".3f")
plt.title(rf"Proportion of Genome Shared IBD Between $\boldsymbol{{Anopheles}}$ Populations on Chromosome {chr}", 
          fontsize=14, fontweight='bold', pad=20)
plt.xlabel("Population", fontsize=12, fontweight='bold')
plt.ylabel("Population", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{output_dir}/proportion_ibd_heatmap_{chr}.png", dpi=300)