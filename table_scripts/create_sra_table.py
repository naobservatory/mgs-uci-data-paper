#!/usr/bin/env python3

import pandas as pd
import os
import numpy as np
from datetime import datetime

# Setting directories and S3 buckets
table_dir = "../tables"
results_dir = "../data/results"


def create_sra_table():
    metadata = pd.read_csv(os.path.join(results_dir, "qc_basic_stats.tsv"), sep="\t")
    metadata = metadata[
        metadata["stage"] == "cleaned"
    ]  # Only want one entry per sample, not two!

    df = pd.DataFrame()
    df["sample_name"] = metadata["sample"]
    df["library_ID"] = metadata["sample"]
    # Add additional columns
    df["title"] = "Metatranscriptomic Sequencing of Los Angeles Influent Sewage"
    df["library_strategy"] = "RNA-Seq"
    df["library_source"] = "Metagenomic"
    df["library_selection"] = "cDNA"
    df["library_layout"] = "paired"
    df["platform"] = "ILLUMINA"
    df["date"] = (
        df["sample_name"].str.split("-").str[1:4].str.join("-").astype("datetime64[s]")
    )
    df.sort_values(by="date", inplace=True)
    df["instrument_model"] = np.where(
        df["date"] >= datetime(2024, 2, 25),
        "Illumina NovaSeq X",
        "Illumina NovaSeq 6000",
    )
    df["design_description"] = (
        "Samples (60-ml) were filtered through 0.22-μm vacuum filters, then ultracentrifugated with 10-kDa Amicon filters until volumes were reduced to <500 μl. These concentrates were stored at −80°C before RNA extraction. Subsequently an Invitrogen PureLink RNA minikit with DNase (Invitrogen, Waltham, MA cite) was used to extract RNA following the manufacturer's protocol. Library preparation was performed by UC GRT Hub using the Illumina RNA prep with enrichment kit."
    )
    df["filetype"] = "fastq"
    df["filename"] = df["sample_name"] + "-part-1_1.fastq.gz"
    df["filename2"] = df["sample_name"] + "-part-1_2.fastq.gz"
    df["filename3"] = df["sample_name"] + "-part-2_1.fastq.gz"
    df["filename4"] = df["sample_name"] + "-part-2_2.fastq.gz"
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
        "filename3",
        "filename4",
    ]

    # Create the SRA table with selected columns
    sra_table = df[sra_columns]

    # Save the SRA table as a TSV file
    sra_table.to_csv(os.path.join(table_dir, "sra_table.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    create_sra_table()
