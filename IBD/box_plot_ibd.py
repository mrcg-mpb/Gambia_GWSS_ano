#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# 1. Load your IBD + metadata
# ----------------------------
ibd_df = pd.read_csv("phased_data.ibd.gz", sep="\t", compression="gzip", header=None)
ibd_df.columns = ['ID1', 'HAP1', 'ID2', 'HAP2', 'CHR', 'START', 'END', 'LOD', 'CM_LEN']

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Merge population info from metadata
sample_metadata = pd.read_csv("sample_data.csv")  # contains pop -> taxon
mapping_dict = sample_metadata.set_index("pop")["taxon"].to_dict()

# Annotate ibd_df with population info
ibd_df["pop1"] = ibd_df["ID1"].map(mapping_dict)
ibd_df["pop2"] = ibd_df["ID2"].map(mapping_dict)

# Remove any rows with missing population info
ibd_df = ibd_df.dropna(subset=["pop1", "pop2"])

# Create sorted population pair name (e.g., coluzzii-gambiae, gambiae-gambiae)
ibd_df["pop_pair"] = ibd_df.apply(
    lambda row: "-".join(sorted([row["pop1"], row["pop2"]])),
    axis=1
)

# Add a column to classify as "within" or "between"
ibd_df["comparison"] = ibd_df.apply(
    lambda row: "within" if row["pop1"] == row["pop2"] else "between", axis=1
)

# --- Boxplot: Within Population Only ---
within_df = ibd_df[ibd_df["comparison"] == "within"]

plt.figure(figsize=(10, 5))
sns.boxplot(data=within_df, x="pop_pair", y="CM_LEN", palette="pastel")
plt.title("IBD Segment Lengths (cM) - Within Population")
plt.ylabel("IBD Segment Length (cM)")
plt.xlabel("Population Pair")
plt.tight_layout()
plt.savefig("ibd_within_population_boxplot.png", dpi=300)
plt.show()

# --- Boxplot: Between Populations Only ---
between_df = ibd_df[ibd_df["comparison"] == "between"]

plt.figure(figsize=(10, 5))
sns.boxplot(data=between_df, x="pop_pair", y="CM_LEN", palette="Set2")
plt.title("IBD Segment Lengths (cM) - Between Populations")
plt.ylabel("IBD Segment Length (cM)")
plt.xlabel("Population Pair")
plt.tight_layout()
plt.savefig("ibd_between_population_boxplot.png", dpi=300)
plt.show()
