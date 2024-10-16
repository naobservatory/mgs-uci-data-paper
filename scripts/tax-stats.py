#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
import numpy as np
from scipy.stats import gmean

# Setting directories and S3 buckets
workflow_results_dir = "../workflow_results"
work_dir = "../work"
s3_base_url = "s3://nao-mgs-wb/"

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

os.makedirs(workflow_results_dir, exist_ok=True)
for dataset in datasets:
    local_folder_path = os.path.join(workflow_results_dir, dataset)

    # Check if the folder exists locally
    if not os.path.exists(local_folder_path):
        print(f"Folder {dataset} not found. Downloading from S3...")

        # Construct the S3 URL for the folder
        s3_folder_url = f"{s3_base_url}{dataset}/"

        # Use AWS CLI to download the folder
        try:
            subprocess.run(
                ["aws", "s3", "cp", s3_folder_url, local_folder_path, "--recursive"],
                check=True,
            )
            print(f"Successfully downloaded {dataset}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {dataset}: {e}")
    else:
        continue


bacteria_fractions = np.array([])
virus_fractions = np.array([])
archea_fractions = np.array([])
eukaryota_fractions = np.array([])
unclassified_fractions = np.array([])

for delivery in os.listdir(workflow_results_dir):
    bracken_path = os.path.join(
        workflow_results_dir,
        delivery,
        "output",
        "results",
        "taxonomy",
        "bracken_reports_merged.tsv.gz",
    )
    if os.path.exists(bracken_path):
        unzipped_bracken = bracken_path.replace(".gz", "")
        subprocess.run(
            ["gunzip", "-c", bracken_path],
            stdout=open(unzipped_bracken, "w"),
            check=True,
        )
        bracken_df = pd.read_csv(unzipped_bracken, sep="\t")
        # print(bracken_df)
        ribo_bracken = bracken_df[bracken_df["ribosomal"] == True]
        bacteria_shares = ribo_bracken[ribo_bracken["name"] == "Bacteria"][
            "fraction_total_reads"
        ]
        virus_shares = ribo_bracken[ribo_bracken["name"] == "Viruses"][
            "fraction_total_reads"
        ]
        archea_shares = ribo_bracken[ribo_bracken["name"] == "Archaea"][
            "fraction_total_reads"
        ]
        eukaryota_shares = ribo_bracken[ribo_bracken["name"] == "Eukaryota"][
            "fraction_total_reads"
        ]
        bacteria_fractions = np.append(bacteria_fractions, bacteria_shares)
        virus_fractions = np.append(virus_fractions, virus_shares)
        archea_fractions = np.append(archea_fractions, archea_shares)
        eukaryota_fractions = np.append(eukaryota_fractions, eukaryota_shares)
    kraken_path = os.path.join(
        workflow_results_dir,
        delivery,
        "output",
        "results",
        "taxonomy",
        "kraken_reports_merged.tsv.gz",
    )
    if os.path.exists(kraken_path):
        unzipped_kraken = kraken_path.replace(".gz", "")
        subprocess.run(
            ["gunzip", "-c", kraken_path],
            stdout=open(unzipped_kraken, "w"),
            check=True,
        )
        kraken_df = pd.read_csv(unzipped_kraken, sep="\t")
        ribo_kraken = kraken_df[kraken_df["ribosomal"] == True]
        # print(ribo_kraken)
        unclassified_kraken = ribo_kraken[ribo_kraken["taxid"] == 0]
        unclassified_share = unclassified_kraken["pc_reads_total"] / 100
        unclassified_fractions = np.append(unclassified_fractions, unclassified_share)


print("Bacteria: ", np.mean(bacteria_fractions))
print("Virus: ", np.mean(virus_fractions))
print("Archaea: ", np.mean(archea_fractions))
print("Eukaryota: ", np.mean(eukaryota_fractions))
print("Unclassified: ", np.mean(unclassified_fractions))
