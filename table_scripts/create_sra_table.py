#!/usr/bin/env python3

import pandas as pd
import os
import numpy as np

# Setting directories and S3 buckets
table_dir = "../tables"
metadata_dir = "../delivery_metadata"

dataset_run_type = {
    "JR-2024-03-22-a": "pilot",
    "JR-2024-03-22-b": "pilot",
    "JR-2024-04-12": "full",
    "JR-2024-04-15": "full",
    "JR-2024-04-16": "full",
    "JR-2024-08-06": "pilot",
    "JR-2024-08-27": "full",
}

dataset_instrument = {
    "JR-2024-03-22-a": "Illumina NovaSeq 6000",
    "JR-2024-03-22-b": "Illumina NovaSeq 6000",
    "JR-2024-04-12": "Illumina NovaSeq 6000",
    "JR-2024-04-15": "Illumina NovaSeq 6000",
    "JR-2024-04-16": "Illumina NovaSeq 6000",
    "JR-2024-08-06": "Illumina NovaSeq X",
    "JR-2024-08-27": "Illumina NovaSeq X",
}


def create_sra_table():
    df = pd.DataFrame()
    for dataset in dataset_run_type.keys():
        metadata_file = os.path.join(metadata_dir, f"{dataset}.metadata.tsv")
        metadata_df = pd.read_csv(metadata_file, sep="\t")

        metadata_df["run_type"] = dataset_run_type[dataset]
        metadata_df["instrument_model"] = dataset_instrument[dataset]

        precise_dates = metadata_df[
            (metadata_df["date"] != "2023") & (metadata_df["date"] != "2024")
        ]["date"]
        print(precise_dates)
        precise_dates = [np.datetime64(date) for date in precise_dates]
        min_date = np.min(precise_dates)
        max_date = np.max(precise_dates)
        metadata_df["date_with_dash"] = [
            f"{min_date}--{max_date}" if date in ("2023", "2024") else date
            for date in metadata_df["date"]
        ]
        metadata_df["date"] = [
            f"{min_date}/{max_date}" if date in ("2023", "2024") else date
            for date in metadata_df["date"]
        ]
        metadata_df["sample_name"] = "HTP-" + metadata_df["date_with_dash"]
        df = pd.concat([df, metadata_df])

    df.reset_index(inplace=True)

    df["library"] = df["sample_name"] + "_" + df["run_type"]

    df["library_count"] = df.groupby("library").cumcount() + 1
    df["duplicate_flag"] = df.groupby("library").transform("size") > 1
    df["library_ID"] = np.where(
        df["duplicate_flag"],
        df["library"] + "_" + df["library_count"].astype(str),
        df["library"],
    )
    df = df.drop("library_count", axis=1)

    # Add additional columns
    df["title"] = "Metatranscriptomic Sequencing of Los Angeles Influent Sewage"
    df["library_strategy"] = "RNA-Seq"
    df["library_source"] = "Metagenomic"
    df["library_selection"] = "cDNA"
    df["library_layout"] = "paired"
    df["platform"] = "ILLUMINA"
    df["design_description"] = (
        "Samples (60-ml) were filtered through 0.22-μm vacuum filters, then ultracentrifugated with 10-kDa Amicon filters until volumes were reduced to <500 μl. These concentrates were stored at −80°C before RNA extraction. Subsequently an Invitrogen PureLink RNA minikit with DNase (Invitrogen, Waltham, MA cite) was used to extract RNA following the manufacturer's protocol. Library preparation was performed by UC GRT Hub using the Illumina RNA prep with enrichment kit."
    )
    df["filetype"] = "fastq"
    df["filename"] = df["sample"] + "_1.fastq.gz"
    df["filename2"] = df["sample"] + "_2.fastq.gz"

    # Set demultiplexed to False for sample with false metadata.
    df.loc[df["sample"] == "JR-2024-08-27-xR007-PrNotRecog", "demultiplexed"] = False

    # Define the columns for the SRA table
    sra_columns = [
        "sample_name",
        "library_ID",
        "title",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "platform",
        "instrument_model",
        "design_description",
        "filetype",
        "filename",
        "filename2",
        "demultiplexed",
    ]

    # Create the SRA table with selected columns
    sra_table = df[sra_columns]

    # Save the SRA table as a TSV file
    sra_table.to_csv(os.path.join(table_dir, "sra_table.tsv"), sep="\t", index=False)
    print("SRA table has been saved as 'sra_table.tsv'")


if __name__ == "__main__":
    create_sra_table()
