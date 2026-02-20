import pandas as pd
import numpy as np

# Load table
df = pd.read_csv("summary.tsv", sep="\t")

# Remove % signs
for col in ["Unique%", "Multi%", "TooManyLoci%", "UnmappedTooShort%", "UnmappedOther%", "MismatchRate%"]:
    df[col] = df[col].str.replace("%", "").astype(float)

# Convert numeric columns
df["InputReads"] = df["InputReads"].astype(float)
df["AvgMappedLength"] = df["AvgMappedLength"].astype(float)

# Weighted mean function
def weighted_mean(values, weights):
    return np.sum(values * weights) / np.sum(weights)

metrics = [
    "Unique%",
    "Multi%",
    "TooManyLoci%",
    "UnmappedTooShort%",
    "UnmappedOther%",
    "MismatchRate%",
    "AvgMappedLength"
]

with open("Summary.txt", "w") as out:

    for metric in metrics:
        values = df[metric]

        out.write(f"\n{metric}\n")
        out.write(f"Mean: {values.mean():.2f}\n")
        out.write(f"Median: {values.median():.2f}\n")
        out.write(f"Min: {values.min():.2f}\n")
        out.write(f"Max: {values.max():.2f}\n")
        out.write(f"SD: {values.std(ddof=0):.2f}\n")

        # Weighted mean only for mapping percentages
        if metric in ["Unique%", "Multi%", "TooManyLoci%", "UnmappedTooShort%", "UnmappedOther%"]:
            wmean = weighted_mean(values, df["InputReads"])
            out.write(f"Weighted Mean: {wmean:.2f}\n")