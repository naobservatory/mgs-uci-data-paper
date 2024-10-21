#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
import numpy as np
from scipy.stats import gmean
from collections import defaultdict

# Setting directories and S3 buckets
workflow_results_dir = "../workflow_results"
table_dir = "../tables"
s3_base_url = "s3://nao-mgs-wb/"
os.makedirs(workflow_results_dir, exist_ok=True)


datasets = [
    "JR-2024-03-22-a",
    "JR-2024-03-22-b",
    "JR-2024-04-12",
    "JR-2024-04-15",
    "JR-2024-04-16",
    "JR-2024-08-06",
    "JR-2024-08-27",
]


def get_data():
    for dataset in datasets:
        local_folder_path = os.path.join(workflow_results_dir, dataset)

        if not os.path.exists(local_folder_path):
            print(f"Folder {dataset} not found. Downloading from S3...")

            s3_folder_url = f"{s3_base_url}{dataset}/"

            try:
                subprocess.run(
                    [
                        "aws",
                        "s3",
                        "cp",
                        s3_folder_url,
                        local_folder_path,
                        "--recursive",
                    ],
                    check=True,
                )
                print(f"Successfully downloaded {dataset}")
            except subprocess.CalledProcessError as e:
                print(f"Error downloading {dataset}: {e}")


def generate_table():
    with open(f"{table_dir}/table_3.tsv", "w") as outf:
        outf.write("sample\tbacteria\tvirus\tarchaea\teukaryota\tunclassified\n")

        tax_ids = {
            "Bacteria": 2,
            "Viruses": 10239,
            "Archaea": 2157,
            "Eukaryota": 2759,
            "Unclassified": 0,
        }

        for delivery in os.listdir(workflow_results_dir):
            kraken_path = os.path.join(
                workflow_results_dir,
                delivery,
                "output",
                "results",
                "taxonomy",
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

            for sample, sample_data in kraken_df.groupby("sample"):
                n_reads_bacteria = sample_data[
                    sample_data["taxid"] == tax_ids["Bacteria"]
                ]["n_reads_clade"].sum()
                n_reads_virus = sample_data[sample_data["taxid"] == tax_ids["Viruses"]][
                    "n_reads_clade"
                ].sum()
                n_reads_archea = sample_data[
                    sample_data["taxid"] == tax_ids["Archaea"]
                ]["n_reads_clade"].sum()
                n_reads_eukaryota = sample_data[
                    sample_data["taxid"] == tax_ids["Eukaryota"]
                ]["n_reads_clade"].sum()
                n_reads_unclassified = sample_data[
                    sample_data["taxid"] == tax_ids["Unclassified"]
                ]["n_reads_clade"].sum()

                total_reads = (
                    n_reads_bacteria
                    + n_reads_virus
                    + n_reads_archea
                    + n_reads_eukaryota
                    + n_reads_unclassified
                )

                bacteria_ra = f"{n_reads_bacteria / total_reads:.2e}"
                virus_ra = f"{n_reads_virus / total_reads:.2e}"
                archea_ra = f"{n_reads_archea / total_reads:.2e}"
                eukaryota_ra = f"{n_reads_eukaryota / total_reads:.2e}"
                unclassified_ra = f"{n_reads_unclassified / total_reads:.2e}"

                outf.write(
                    f"{sample}\t{bacteria_ra}\t{virus_ra}\t{archea_ra}\t{eukaryota_ra}\t{unclassified_ra}\n"
                )


def start():
    get_data()
    generate_table()


if __name__ == "__main__":
    start()
