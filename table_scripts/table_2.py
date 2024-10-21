#!/usr/bin/env python3

# Imports

import os


# Setting directories and S3 buckets
delivery_metadata_dir = "../delivery_metadata"
table_dir = "../tables"

os.makedirs(table_dir, exist_ok=True)

# Loading sample metadata (delivery level metadata, before being touched by mgs-workflow)
datasets = [
    "JR-2024-03-22-a",
    "JR-2024-03-22-b",
    "JR-2024-04-12",
    "JR-2024-04-15",
    "JR-2024-04-16",
    "JR-2024-08-06",
    "JR-2024-08-27",
]

dataset_to_uci = {
    "JR-2024-03-22-a": "UCI-2024-03-pilot-a",
    "JR-2024-03-22-b": "UCI-2024-03-pilot-b",
    "JR-2024-04-12": "UCI-2024-04-full-a",
    "JR-2024-04-15": "UCI-2024-04-full-b",
    "JR-2024-04-16": "UCI-2024-04-full-NA",  # This might need further investigation
    "JR-2024-08-06": "UCI-2024-08-pilot",
    "JR-2024-08-27": "UCI-2024-08-full",
}


def start():
    with open(f"{table_dir}/table_2.tsv", "w") as outf:

        metadata_files = os.listdir(delivery_metadata_dir)
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
                    "Demultiplexed",
                    "Enrichment",
                    "Reads",
                    "NA Type",
                    "Dataset",
                ]
            )
            + "\n"
        )
        for metadata_file in metadata_files:
            dataset = metadata_file.split(".")[0]
            if dataset == "":
                continue
            with open(f"{delivery_metadata_dir}/{metadata_file}", "r") as inf:
                # Skip the header line
                next(inf)
                for line in inf:
                    (
                        sample,
                        country,
                        county,
                        date,
                        demultiplexed,
                        enrichment,
                        fine_location,
                        location,
                        na_type,
                        reads,
                        state,
                    ) = line.strip().split("\t")

                    treatment_plant = (
                        "Hyperion Treatment Plant" if "HTP" == fine_location else "NA"
                    )

                    city = location

                    if int(reads) > 1e9:
                        reads = f"{int(reads) / 1e9:.2f}B"
                    else:
                        reads = f"{int(reads) / 1e6:.0f}M"
                    uci_name = dataset_to_uci[dataset]
                    outf.write(
                        "\t".join(
                            [
                                sample,
                                country,
                                state,
                                county,
                                city,
                                treatment_plant,
                                date,
                                demultiplexed,
                                enrichment,
                                reads,
                                na_type,
                                uci_name,
                            ]
                        )
                        + "\n"
                    )


if __name__ == "__main__":
    start()
