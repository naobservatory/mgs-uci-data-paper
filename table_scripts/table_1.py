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

basic_stats_summary = (
    basic_stats.groupby(["sequencing_machine"])
    .agg(
        total_read_pairs=("n_read_pairs", "sum"),
        mean_gc_content=("percent_gc", "mean"),
        total_bases=("n_bases_approx", "sum"),
        n_samples=("sample", "nunique"),
        date_range=("date", lambda x: f"{x.min().date()} to {x.max().date()}"),
    )
    .reset_index()
)

basic_stats_summary["# Bases"] = basic_stats_summary["total_bases"].apply(
    lambda x: f"{x / 1e12:.2f}T"
)
basic_stats_summary["# Reads"] = basic_stats_summary["total_read_pairs"].apply(
    lambda x: f"{x / 1e9:.2f}B"
)
basic_stats_summary["# Samples"] = basic_stats_summary["n_samples"]
basic_stats_summary["GC Content"] = basic_stats_summary["mean_gc_content"].apply(
    lambda x: f"{x:.2f}%"
)
basic_stats_summary["Date Range"] = basic_stats_summary["date_range"]
basic_stats_summary = basic_stats_summary[
    [
        "sequencing_machine",
        "# Samples",
        "Date Range",
        "# Reads",
        "# Bases",
        "GC Content",
    ]
]

basic_stats_summary.to_csv(f"{table_dir}/table_1.tsv", sep="\t", index=False)
exit()
