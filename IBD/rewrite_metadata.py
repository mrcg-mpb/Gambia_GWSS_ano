#!/usr/bin/env python3
import pandas as pd

metadata = pd.read_csv("sample_metadata.csv")

# define pop columns

metadata["pop"] = metadata["sample_id"] + "_" + metadata["taxon"]

# save to new file
metadata.to_csv("sample_data.csv", index=False)