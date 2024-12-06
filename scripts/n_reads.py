#!/usr/bin/env python3

# Imports
import pandas as pd
import numpy as np
from datetime import datetime

# Setting directories and S3 buckets
workflow_results_dir = "../data/results"
delivery_metadata_dir = "../data"
table_dir = "../tables"

basic_stats = pd.read_csv(f"{workflow_results_dir}/qc_basic_stats.tsv", sep="\t")
basic_stats["date"] = pd.to_datetime(
    basic_stats["sample"].str.split("-").str[1:4].str.join("-")
)
basic_stats["sequencing_machine"] = np.where(
    basic_stats["date"] >= datetime(2024, 2, 25),
    "NovaSeq X",
    "NovaSeq 6000",
)

basic_stats = basic_stats[basic_stats["stage"] == "raw_concat"]

print(f"Mean: {basic_stats['n_read_pairs'].mean():.2e}")
print(f"Min: {basic_stats['n_read_pairs'].min():.2e}")
print(f"Max: {basic_stats['n_read_pairs'].max():.2e}")
