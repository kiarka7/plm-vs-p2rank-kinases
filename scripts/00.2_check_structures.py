#!/usr/bin/env python3
"""
00.2_check_structures.py

Check that all expected structures listed in an Excel file
have been downloaded into a given directory as .cif files.

- Reads an Excel file with a column "PDB" (e.g. "4CZT" or "4CZT A").
- Extracts the 4-character PDB ID (first 4 characters, uppercased).
- In the given structures directory, looks for *.cif files.
- Compares:
    expected PDB IDs (from Excel)
    vs.
    present PDB IDs (from filenames "<pdb>.cif").
- Prints a small summary and lists missing/extra IDs.
"""

import os
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--excel",
        required=True,
        help="Path to Excel file (e.g. Kinase_Type_I/Kinase_Ligands_Type I.xlsx)",
    )
    p.add_argument(
        "--sheet",
        help="Sheet name or index (default: first sheet)",
    )
    p.add_argument(
        "--pdb-column",
        default="PDB",
        help='Column name with PDB IDs (default: "PDB")',
    )
    p.add_argument(
        "--struct-dir",
        required=True,
        help="Directory where .cif files are stored (e.g. Kinase_Type_I/structures)",
    )
    return p.parse_args()


def load_expected_pdb_ids(excel_path: str, sheet, pdb_column: str):
    if sheet is None:
        df = pd.read_excel(excel_path)
    else:
        df = pd.read_excel(excel_path, sheet_name=sheet)

    if pdb_column not in df.columns:
        raise ValueError(
            f"PDB column '{pdb_column}' not found in {excel_path}. "
            f"Available columns: {list(df.columns)}"
        )

    pdb_ids = set()
    for val in df[pdb_column]:
        if pd.isna(val):
            continue
        s = str(val).strip()
        if not s:
            continue
        pdb = s[:4].upper()
        if len(pdb) == 4 and pdb.isalnum():
            pdb_ids.add(pdb)

    return pdb_ids


def load_present_pdb_ids(struct_dir: str):
    present = set()
    if not os.path.isdir(struct_dir):
        return present

    for fname in os.listdir(struct_dir):
        if not fname.lower().endswith(".cif"):
            continue
        pdb = fname[:4].upper()
        if len(pdb) == 4 and pdb.isalnum():
            present.add(pdb)
    return present


def main():
    args = parse_args()

    expected = load_expected_pdb_ids(args.excel, args.sheet, args.pdb_column)
    present = load_present_pdb_ids(args.struct_dir)

    missing = sorted(expected - present)
    extra   = sorted(present - expected)

    print(f"[i] Expected structures (from Excel): {len(expected)}")
    print(f"[i] Present  .cif files in {args.struct_dir}: {len(present)}")
    print(f"[i] Missing: {len(missing)}")
    print(f"[i] Extra:   {len(extra)}")

    if missing:
        print("\n[missing PDB IDs]")
        print(" ".join(missing))

    if extra:
        print("\n[extra PDB IDs in struct-dir (not in Excel)]")
        print(" ".join(extra))


if __name__ == "__main__":
    main()
