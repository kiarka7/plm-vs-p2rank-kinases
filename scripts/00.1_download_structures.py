#!/usr/bin/env python3
"""
00.1_download_structures.py

Download mmCIF structures listed in an Excel sheet into a given folder.

- Reads an Excel file with a column "PDB" (e.g. "4CZT" or "4CZT A").
- Extracts the 4-character PDB ID (first 4 characters, uppercased).
- Downloads mmCIF files from RCSB:
    https://files.rcsb.org/download/<PDB>.cif
- Saves them as <pdb>.cif (lowercase) into the target directory.
- Skips files that already exist unless --force is used.
- Retries failed downloads a few times before giving up.
- Prints a small summary at the end.

Requires:
    pandas, requests, openpyxl (for .xlsx)
"""

import os
import argparse
import time
import requests
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
        "--out-dir",
        required=True,
        help="Directory where .cif files will be saved (e.g. Kinase_Type_I/structures)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist",
    )
    p.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of retries for each PDB on failure (default: 3)",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="HTTP timeout in seconds (default: 30)",
    )
    return p.parse_args()


def load_pdb_ids(excel_path: str, sheet, pdb_column: str):
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

    return sorted(pdb_ids)


def download_cif_once(pdb: str, out_path: str, timeout: int = 30) -> bool:
    """
    Single attempt to download <PDB>.cif from RCSB.
    Returns True if downloaded successfully, False otherwise.
    """
    url = f"https://files.rcsb.org/download/{pdb}.cif"
    print(f"[dl]   {pdb} -> {out_path}")
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            print(f"[!] HTTP {r.status_code} for {pdb}")
            return False
        content = r.text
        if len(content) < 500:
            print(f"[!] Downloaded file for {pdb} is very small (len={len(content)}), may be invalid")
            with open(out_path, "w") as f:
                f.write(content)
            return False

        with open(out_path, "w") as f:
            f.write(content)
        return True
    except Exception as e:
        print(f"[!] Error downloading {pdb}: {e}")
        return False


def download_cif_with_retries(pdb: str, out_path: str, retries: int, timeout: int) -> bool:
    """
    Try downloading a CIF file multiple times before giving up.
    Returns True if successful at least once, otherwise False.
    """
    for attempt in range(1, retries + 1):
        ok = download_cif_once(pdb, out_path, timeout=timeout)
        if ok:
            return True
        if attempt < retries:
            print(f"[i] Retry {attempt}/{retries} failed for {pdb}, waiting a bit before next attempt...")
            time.sleep(3)
    return False


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    pdb_ids = load_pdb_ids(args.excel, args.sheet, args.pdb_column)
    print(f"[i] Found {len(pdb_ids)} unique PDB IDs in Excel")

    n_dl = 0
    n_skip = 0
    n_fail = 0

    for pdb in pdb_ids:
        out_path = os.path.join(args.out_dir, f"{pdb.lower()}.cif")

        if os.path.isfile(out_path) and not args.force:
            print(f"[skip] {pdb.lower()}.cif already exists")
            n_skip += 1
            continue

        ok = download_cif_with_retries(pdb, out_path, retries=args.retries, timeout=args.timeout)
        if ok:
            n_dl += 1
        else:
            n_fail += 1
            if os.path.isfile(out_path) and os.path.getsize(out_path) < 500:
                print(f"[i] Removing suspicious small file {out_path}")
                try:
                    os.remove(out_path)
                except OSError:
                    pass

    print("\n[summary]")
    print(f"  total unique PDBs: {len(pdb_ids)}")
    print(f"  downloaded:        {n_dl}")
    print(f"  skipped existing:  {n_skip}")
    print(f"  failed:            {n_fail}")


if __name__ == "__main__":
    main()
