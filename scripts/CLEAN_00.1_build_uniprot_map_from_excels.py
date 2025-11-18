#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CLEAN_00a.build_uniprot_map_from_excels.py

Extract unique UniProt IDs from kinase Excel tables and map them
to primary UniProt accessions using the UniProt REST API.

Output:
    TSV with columns:
        dataset_uniprot_id   accession   returned_uniprot_id

Requirements:
    - internet connection
    - Python packages: pandas, requests
"""

import argparse
from pathlib import Path
from typing import Set, Dict, Tuple

import pandas as pd
import requests
from time import sleep


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def collect_uniprot_ids(excel_paths, uniprot_col: str) -> Set[str]:
    ids: Set[str] = set()
    for p in excel_paths:
        print(f"[i] Reading Excel: {p}")
        try:
            df = pd.read_excel(p)
        except Exception as e:
            print(f"[warn] Cannot read {p}: {e}")
            continue

        if uniprot_col not in df.columns:
            print(f"[warn] Column '{uniprot_col}' not found in {p}, skipping.")
            continue

        # dropna + strip
        vals = (
            df[uniprot_col]
            .dropna()
            .astype(str)
            .str.strip()
        )
        ids.update(v for v in vals if v)
    print(f"[i] Collected {len(ids)} unique UniProt IDs from Excels.")
    return ids


def query_uniprot(id_str: str, sleep_sec: float = 0.1) -> Tuple[str, str]:
    """
    Query UniProt for a single ID (e.g. DSOR1_DROME) and return
    (accession, returned_id). If not found, returns ("NA", "NA").
    """
    params = {
        "query": f"id:{id_str}",
        "fields": "accession,id",
        "format": "tsv",
        "size": 1,
    }
    try:
        r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=10)
        r.raise_for_status()
        text = r.text.strip()
        lines = text.splitlines()
        if len(lines) < 2:
            # nothing found
            print(f"[warn] UniProt ID not found: {id_str}")
            return "NA", "NA"
        header = lines[0].split("\t")
        values = lines[1].split("\t")
        row = dict(zip(header, values))
        acc = row.get("Entry", "NA")
        uid = row.get("Entry Name", "NA")
        # be gentle to the server
        if sleep_sec > 0:
            sleep(sleep_sec)
        return acc, uid
    except Exception as e:
        print(f"[warn] UniProt query failed for {id_str}: {e}")
        return "NA", "NA"


def build_mapping(uniprot_ids: Set[str]) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    for i, uid in enumerate(sorted(uniprot_ids), 1):
        print(f"[{i}/{len(uniprot_ids)}] Mapping {uid} ...", end="", flush=True)
        acc, ret_id = query_uniprot(uid)
        mapping[uid] = (acc, ret_id)
        print(f" -> {acc}")
    return mapping


def save_tsv(mapping: Dict[str, Tuple[str, str]], out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("dataset_uniprot_id\taccession\treturned_uniprot_id\n")
        for uid, (acc, ret_id) in sorted(mapping.items()):
            f.write(f"{uid}\t{acc}\t{ret_id}\n")
    print(f"[ok] Written UniProt mapping TSV: {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--excels",
        nargs="+",
        required=True,
        help="List of Excel files (can use shell glob expansion, e.g. Kinase_*/*.xlsx).",
    )
    ap.add_argument(
        "--uniprot-col",
        default="UniprotID",
        help="Name of the UniProt column in Excel tables (default: UniprotID).",
    )
    ap.add_argument(
        "--out",
        default="kinase_uniprot_map.tsv",
        help="Output TSV file (default: kinase_uniprot_map.tsv).",
    )
    args = ap.parse_args()

    excel_paths = [Path(p) for p in args.excels]
    ids = collect_uniprot_ids(excel_paths, args.uniprot_col)
    mapping = build_mapping(ids)
    save_tsv(mapping, Path(args.out))


if __name__ == "__main__":
    main()
