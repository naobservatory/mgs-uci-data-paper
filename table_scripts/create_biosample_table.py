#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd


def create_biosample_table():
    # Setting directories and S3 buckets
    table_dir = "../tables"
    results_dir = "../data/results"

    metadata = pd.read_csv(os.path.join(results_dir, "qc_basic_stats.tsv"), sep="\t")
    metadata["date"] = metadata["sample"].str.split("-").str[1:4].str.join("-")
    # Create dictionary mapping sample to date
    sample_to_date = dict(zip(metadata["sample"], metadata["date"]))

    with open(f"{table_dir}/bio_sample_table.tsv", "w") as outf:
        outf.write(
            "\t".join(
                [
                    "Sample Name",
                    "Sample Title",
                    "BioProject accession",
                    "Organism",
                    "collection date",
                    "broad-scale environmental context",
                    "local-scale environmental context",
                    "environmental medium",
                    "geographic location",
                    "latitude and longitude",
                ]
            )
            + "\n"
        )
        for name, date in sample_to_date.items():
            outf.write(
                "\t".join(
                    [
                        name,
                        "Influent wastewater from Hyperion Treatment Plant (LA, USA)",
                        "",
                        "wastewater metagenome",
                        date,
                        "wastewater",
                        "influent wastewater",
                        "Composite wastewater",
                        "USA: California",
                        "33.924223 N 118.431516 W",
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    create_biosample_table()
