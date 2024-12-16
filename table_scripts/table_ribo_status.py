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
s3_base_url = "s3://nao-mgs-wb/"
os.makedirs(workflow_results_dir, exist_ok=True)


def generate_table():
    with open(f"{table_dir}/table_ribo_status.tsv", "w") as outf:
        outf.write(
            "Sample\tRibosomal\tBacteria\tVirus\tArchaea\tEukaryota\tUnclassified\tTotal Reads\n"
        )

        tax_ids = {
            "Bacteria": 2,
            "Viruses": 10239,
            "Archaea": 2157,
            "Eukaryota": 2759,
            "Unclassified": 0,
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
        for ribo_status in kraken_df["ribosomal"].unique():
            for sample in kraken_df["sample"].unique():
                sample_data = kraken_df[
                    (kraken_df["sample"] == sample)
                    & (kraken_df["ribosomal"] == ribo_status)
                ]
                sample_data = sample_data.groupby(["taxid", "ribosomal"]).agg(
                    {"n_reads_clade": "sum"}
                )
                sample_data = sample_data.reset_index()
                # Get read counts for all taxonomic groups in one step
                reads_by_group = {
                    group: (
                        sample_data[sample_data["taxid"] == tax_id][
                            "n_reads_clade"
                        ].values[0]
                        if len(sample_data[sample_data["taxid"] == tax_id]) > 0
                        else 0
                    )
                    for group, tax_id in tax_ids.items()
                }

                n_reads_bacteria = reads_by_group["Bacteria"]
                n_reads_virus = reads_by_group["Viruses"]
                n_reads_archea = reads_by_group["Archaea"]
                n_reads_eukaryota = reads_by_group["Eukaryota"]
                n_reads_unclassified = reads_by_group["Unclassified"]

                total_reads = (
                    n_reads_bacteria
                    + n_reads_virus
                    + n_reads_archea
                    + n_reads_eukaryota
                    + n_reads_unclassified
                )

                bacteria_ra = f"{n_reads_bacteria / total_reads:.2%}"
                virus_ra = f"{n_reads_virus / total_reads:.2%}"
                archea_ra = f"{n_reads_archea / total_reads:.2%}"
                eukaryota_ra = f"{n_reads_eukaryota / total_reads:.2%}"
                unclassified_ra = f"{n_reads_unclassified / total_reads:.2%}"

                outf.write(
                    f"{sample}\t{ribo_status}\t{bacteria_ra}\t{virus_ra}\t{archea_ra}\t{eukaryota_ra}\t{unclassified_ra}\t{total_reads}\n"
                )


def start():

    generate_table()


if __name__ == "__main__":
    start()