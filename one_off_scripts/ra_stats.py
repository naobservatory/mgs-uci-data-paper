#! /usr/bin/env python3
#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
import numpy as np
from scipy.stats import gmean
from collections import defaultdict

# Setting directories and S3 buckets
workflow_results_dir = "../data/results"
table_dir = "../tables"
os.makedirs(workflow_results_dir, exist_ok=True)


def generate_table():
    data = []

    tax_ids = {
        "Bacteria": 2,
        "Viruses": 10239,
        "Archaea": 2157,
        "Eukaryota": 2759,
        "Unclassified": 0,
        "Classified": 1,
    }

    kraken_path = os.path.join(
        workflow_results_dir,
        "kraken_reports_merged.tsv",
    )
    zip_kraken_path = kraken_path + ".gz"
    if os.path.exists(zip_kraken_path):
        subprocess.run(
            ["gunzip", "-c", zip_kraken_path],
            stdout=open(kraken_path, "w"),
            check=True,
        )

    kraken_df = pd.read_csv(kraken_path, sep="\t")
    for sample in kraken_df["sample"].unique():
        sample_data = kraken_df[kraken_df["sample"] == sample]
        sample_data = sample_data.groupby("taxid").agg({"n_reads_clade": "sum"})
        sample_data = sample_data.reset_index()
        # Get read counts for all taxonomic groups in one step
        reads_by_group = {
            group: sample_data[sample_data["taxid"] == tax_id]["n_reads_clade"].values[
                0
            ]
            for group, tax_id in tax_ids.items()
        }

        n_reads_bacteria = reads_by_group["Bacteria"]
        n_reads_virus = reads_by_group["Viruses"]
        n_reads_archea = reads_by_group["Archaea"]
        n_reads_eukaryota = reads_by_group["Eukaryota"]
        n_reads_unclassified = reads_by_group["Unclassified"]
        n_reads_classified = reads_by_group["Classified"]
        total_reads = n_reads_unclassified + n_reads_classified

        bacteria_ra = n_reads_bacteria / total_reads
        virus_ra = n_reads_virus / total_reads
        archea_ra = n_reads_archea / total_reads
        eukaryota_ra = n_reads_eukaryota / total_reads
        unclassified_ra = n_reads_unclassified / total_reads

        data.append(
            [sample, bacteria_ra, virus_ra, archea_ra, eukaryota_ra, unclassified_ra]
        )

    df = pd.DataFrame(
        data,
        columns=[
            "sample",
            "bacteria_ra",
            "virus_ra",
            "archea_ra",
            "eukaryota_ra",
            "unclassified_ra",
        ],
    )
    # Calculate statistics for each group
    columns = df.columns[1:]
    stats = {}
    for group in columns:
        median = df[group].median()
        min_val = df[group].min()
        max_val = df[group].max()
        print(f"{group}: {median:.2%} ({min_val:.2%}, {max_val:.2%})")

    return df


def start():

    generate_table()


if __name__ == "__main__":
    start()
