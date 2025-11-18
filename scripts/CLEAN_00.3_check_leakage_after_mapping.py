#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse

def load_train_accessions(train_file):
    accs = set()
    with open(train_file, "r") as f:
        for line in f:
            parts = line.strip().split(";")
            if not parts:
                continue
            acc = parts[0]
            accs.add(acc)
    return accs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--map", required=True)
    args = ap.parse_args()

    train_accs = load_train_accessions(args.train)
    map_df = pd.read_csv(args.map, sep="\t")

    matched = []
    for _, row in map_df.iterrows():
        acc = str(row["accession"]).strip()
        if acc in train_accs:
            matched.append(row["dataset_uniprot_id"])

    print("\n====== LEAKAGE CHECK ======\n")
    print(f"Total mapped UniProt IDs in dataset: {len(map_df)}")
    print(f"Found in train.txt: {len(matched)}")
    if len(map_df) > 0:
        print(f"Percentage overlap: {100 * len(matched) / len(map_df):.2f} %\n")

    if matched:
        print("Overlap IDs:")
        for m in matched:
            print("  ", m)

    print("\n===========================\n")


if __name__ == "__main__":
    main()
