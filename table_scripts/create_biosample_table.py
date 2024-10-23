#!/usr/bin/env python3

import os
import numpy as np

delivery_metadata_dir = "../delivery_metadata"
table_dir = "../tables"

os.makedirs(table_dir, exist_ok=True)


def create_table():
    all_dates = set()
    metadata_files = os.listdir(delivery_metadata_dir)
    for metadata_file in metadata_files:
        dataset = metadata_file.split(".")[0]
        if dataset == "":
            continue

        precise_dates = []
        dates = []

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
                date = str(date).strip()
                if date != "2023" and date != "2024":
                    precise_dates.append(date)
                dates.append(date)

            precise_dates = [np.datetime64(date) for date in precise_dates]
            min_date = np.min(precise_dates)
            max_date = np.max(precise_dates)
            print(min_date, max_date)
            dates = [
                f"{min_date}/{max_date}" if date in ("2023", "2024") else date
                for date in dates
            ]
            # print(dates)
            all_dates.update(dates)

    # print(f"{dataset}: {min_date} - {max_date}")
    sample_names = {}
    for date in sorted(all_dates):
        if "/" in date:
            min_date, max_date = date.split("/")
            sample_names[date] = f"HTP-{min_date}--{max_date}"
        else:
            sample_names[date] = f"HTP-{date}"

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
        for date, name in sample_names.items():
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
    create_table()
