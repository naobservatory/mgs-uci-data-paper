#!/usr/bin/env python3

# Imports

import os
import pandas as pd

# Setting directories and S3 buckets
data_dir = "../data"
results_dir = f"{data_dir}/results"
table_dir = "../tables"

os.makedirs(table_dir, exist_ok=True)


def start():

    metadata = pd.read_csv(os.path.join(results_dir, "qc_basic_stats.tsv"), sep="\t")
    metadata["date"] = metadata["sample"].str.split("-").str[1:4].str.join("-")
    # Create dictionary mapping sample to date
    sample_to_date = dict(zip(metadata["sample"], metadata["date"]))
    sample_to_n_read_pairs = dict(zip(metadata["sample"], metadata["n_read_pairs"]))

    with open(f"{table_dir}/table_s1.tsv", "w") as outf:
        outf.write(
            "\t".join(
                [
                    "Sample",
                    "Country",
                    "State",
                    "County",
                    "City",
                    "Treatment Plant",
                    "Date",
                    "Reads",
                ]
            )
            + "\n"
        )
        for name, date in sample_to_date.items():
            n_read_pairs = sample_to_n_read_pairs[name]
            if n_read_pairs > 1e9:
                n_read_pairs = f"{n_read_pairs / 1e9:.2f}B"
            else:
                n_read_pairs = f"{n_read_pairs / 1e6:.2f}M"
            outf.write(
                "\t".join(
                    [
                        name,
                        "United States",
                        "California",
                        "Los Angeles County",
                        "Los Angeles",
                        "Hyperion Treatment Plant",
                        date,
                        n_read_pairs,
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    start()
