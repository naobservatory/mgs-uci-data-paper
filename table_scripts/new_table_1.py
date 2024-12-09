#!/usr/bin/env python3

import pandas as pd
import csv
import os
import gzip
from datetime import datetime
from io import StringIO

# Setting directories
workflow_results_dir = "../workflow_results"
delivery_metadata_dir = "../delivery_metadata"
table_dir = "../tables"

# Loading mgs-workflow results
datasets = [
    "JR-2024-03-22-a",
    "JR-2024-03-22-b",
    "JR-2024-04-12",
    "JR-2024-04-15",
    "JR-2024-04-16",
    "JR-2024-08-06",
    "JR-2024-08-27",
]

# Dictionary mapping dataset names to UCI identifiers
dataset_to_uci = {
    "JR-2024-03-22-a": "UCI-2024-03-pilot-a",
    "JR-2024-03-22-b": "UCI-2024-03-pilot-b",
    "JR-2024-04-12": "UCI-2024-04-full-a",
    "JR-2024-04-15": "UCI-2024-04-full-b",
    "JR-2024-08-06": "UCI-2024-08-pilot",
    "JR-2024-08-27": "UCI-2024-08-full",
}


def get_delivery_date_range():
    delivery_dates = {}
    for metadata_file in os.listdir(delivery_metadata_dir):
        if metadata_file.startswith("."):  # Skip hidden files like .DS_Store
            continue
        dataset = metadata_file.split(".")[0]
        uci_name = dataset_to_uci[dataset]

        with open(f"{delivery_metadata_dir}/{metadata_file}", "r") as inf:
            df = pd.read_csv(inf, sep="\t")
            df["date"] = pd.to_datetime(df["date"], errors="coerce")
            df = df[df["date"].dt.year.isin([2023, 2024])]

            min_date = df["date"].min().strftime("%Y-%m-%d")
            max_date = df["date"].max().strftime("%Y-%m-%d")
            delivery_dates[uci_name] = f"{min_date} to {max_date}"

    return delivery_dates


def sample_to_date():
    sample_to_date = {}
    for metadata_file in os.listdir(delivery_metadata_dir):
        if metadata_file.startswith("."):  # Skip hidden files like .DS_Store
            continue
        with open(f"{delivery_metadata_dir}/{metadata_file}", "r") as inf:
            df = pd.read_csv(inf, sep="\t")
            df["date"] = pd.to_datetime(df["date"], errors="coerce")
            df = df[df["date"].dt.year.isin([2023, 2024])]
            sample_to_date.update(
                dict(zip(df["sample"], df["date"].dt.strftime("%Y-%m-%d")))
            )
    return sample_to_date


def concatenate_qc_basic_stats():
    dfs = []
    sample_to_date_dict = sample_to_date()

    for dataset in datasets:
        gzipped_file = (
            f"{workflow_results_dir}/{dataset}/output/results/qc/qc_basic_stats.tsv.gz"
        )
        if os.path.exists(gzipped_file):
            with gzip.open(gzipped_file, "rt") as f:
                df = pd.read_csv(f, sep="\t")

                if dataset == "JR-2024-04-16":
                    df["date"] = df["sample"].map(sample_to_date_dict)
                    df["date"] = pd.to_datetime(df["date"])
                    df["uci_name"] = df["date"].apply(
                        lambda x: (
                            dataset_to_uci["JR-2024-04-12"]
                            if x.date() <= datetime(2024, 1, 7).date()
                            else dataset_to_uci["JR-2024-04-15"]
                        )
                    )
                else:
                    df["uci_name"] = dataset_to_uci[dataset]

                dfs.append(df)
        else:
            print(f"File not found for {dataset}: {gzipped_file}")

    return pd.concat(dfs, ignore_index=True)


def generate_summary_table():
    delivery_dates = get_delivery_date_range()
    df = concatenate_qc_basic_stats()
    df_raw = df[df["stage"] == "raw_concat"]

    df_raw_summary = (
        df_raw.groupby("uci_name")
        .agg(
            total_read_pairs=("n_read_pairs", "sum"),
            mean_gc_content=("percent_gc", "mean"),
            total_bases=("n_bases_approx", "sum"),
            total_libraries=("sample", "nunique"),
        )
        .reset_index()
    )

    df_raw_summary["Date Range"] = df_raw_summary["uci_name"].map(delivery_dates)
    df_raw_summary["total_read_pairs"] = df_raw_summary["total_read_pairs"].apply(
        lambda x: f"{x/1e9:,.2f}B"
    )
    df_raw_summary["total_bases"] = df_raw_summary["total_bases"].apply(
        lambda x: f"{x/1e9:,.0f}B"
    )
    df_raw_summary["mean_gc_content"] = df_raw_summary["mean_gc_content"].apply(
        lambda x: f"{x:,.2f}%"
    )

    os.makedirs(table_dir, exist_ok=True)
    df_raw_summary.to_csv(f"{table_dir}/table_1.tsv", sep="\t", index=False)


if __name__ == "__main__":
    generate_summary_table()
