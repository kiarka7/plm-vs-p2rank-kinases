#!/usr/bin/env python3
"""
CLEAN_00.build_clean_excel_from_mask.py

Purpose
-------
Given a list of training UniProt accessions (train.txt), create a CLEAN
version of the kinase Excel table for a given dataset.

Inputs
------
  - train.txt
      File from your colleague, each line e.g.
          A0A0B4J272;UNK;UNK;Y97 L68 ...;SEQUENCE
      → the first token before ';' is the UniProt ACC used in training.

  - Excel with ligands, e.g.
        Kinase_Type_I/Kinase_Ligands_Type_I.xlsx
      containing a column with UniProt IDs (entry names or accessions),
      e.g. column 'UniprotID'.

  - (optional) UniProt mapping TSV:
        --uniprot-map-tsv uniprot_entry_to_accession.tsv
        --map-entry-col   EntryName
        --map-acc-col     Accession
      This lets you translate entry names like 'DSOR1_DROME'
      to accessions like 'Q24324' before comparing to train.txt.

Outputs
-------
  - CLEAN Excel, e.g.
        Kinase_Type_I/CLEAN/Kinase_Ligands_Type_I_CLEAN.xlsx
    containing only rows whose UniProt ACC was *not* used for training.

  - CSV with summary statistics, e.g.
        Kinase_Type_I/CLEAN/00.clean_summary.csv

Statistics reported
-------------------
  - Total number of rows
  - Number of rows with a valid UniProt ID (after mapping)
  - Number of rows whose UniProt ACC appears in training
  - Number of unique UniProt ACCs in the table
  - Number of those ACCs present in training
  - Number of CLEAN rows (not in training)
  - If a 'PDB' column exists:
      * number of unique PDB codes,
      * number of PDBs that contain at least one training ACC.
"""

import os
import argparse
import pandas as pd
from typing import Set, Dict, List
import csv

# -----------------------------------------------
# Loading training UniProt accessions
# -----------------------------------------------

def load_train_uniprots(train_path: str) -> Set[str]:
    """
    Load train.txt and return a set of UniProt accessions.

    Assumes each non-empty, non-comment line has format:
        UNIPROT_ACC;...

    Returns all accessions uppercased.
    """
    uniprots: Set[str] = set()
    with open(train_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(";")
            if not parts:
                continue
            up = parts[0].strip()
            if up:
                uniprots.add(up.upper())
    return uniprots

# -----------------------------------------------
# Loading UniProt mapping (entry name → accessions)
# -----------------------------------------------

def load_uniprot_map(tsv_path: str,
                     entry_col: str,
                     acc_col: str) -> Dict[str, Set[str]]:
    """
    Load a TSV file that maps UniProt entry names to accessions.

    Parameters
    ----------
    tsv_path : str
        Path to TSV file.
    entry_col : str
        Column name containing entry names, e.g. 'EntryName' or 'From'.
    acc_col : str
        Column name containing UniProt accessions, e.g. 'Accession' or 'To'.

    Returns
    -------
    mapping : dict
        mapping[ENTRY_NAME_UPPER] = set({ACC1_UPPER, ACC2_UPPER, ...})
    """
    mapping: Dict[str, Set[str]] = {}
    with open(tsv_path, "r", newline="") as f:
        # assume TSV; if your file is comma-separated, just change delimiter=","
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Mapping file {tsv_path!r} has no header.")
        if entry_col not in reader.fieldnames or acc_col not in reader.fieldnames:
            raise ValueError(
                f"Mapping file {tsv_path!r} must have columns "
                f"{entry_col!r} and {acc_col!r}. Found: {reader.fieldnames}"
            )
        for row in reader:
            entry_raw = (row.get(entry_col) or "").strip()
            acc_raw   = (row.get(acc_col) or "").strip()
            if not entry_raw or not acc_raw:
                continue
            entry = entry_raw.upper()
            acc   = acc_raw.upper()
            if entry not in mapping:
                mapping[entry] = set()
            mapping[entry].add(acc)
    return mapping

# -----------------------------------------------
# PDB helper
# -----------------------------------------------

def extract_pdb_id(pdb_cell: str) -> str:
    """
    Extract 4-letter PDB code (uppercased) from a cell.

    Examples
    --------
    '4CZT'   → '4CZT'
    '4CZTA'  → '4CZT'
    '4CZT A' → '4CZT'
    """
    if pdb_cell is None:
        return ""
    s = str(pdb_cell).strip()
    if len(s) < 4:
        return ""
    return s[:4].upper()

# -----------------------------------------------
# Main
# -----------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True,
                    help="Path to train.txt (format: UniProtACC;...;...;...;sequence)")
    ap.add_argument("--excel", required=True,
                    help="Input Excel with kinases (e.g. Kinase_Ligands_Type_I.xlsx)")
    ap.add_argument("--out-excel", required=True,
                    help="Output CLEAN Excel (e.g. CLEAN/Kinase_Ligands_Type_I_CLEAN.xlsx)")
    ap.add_argument("--dataset-name", default="UNKNOWN",
                    help="Dataset name for summary (e.g. Kinase_Type_I)")
    ap.add_argument("--uniprot-col", default="UniprotID",
                    help="Name of the column with UniProt IDs in Excel "
                         "(entry names or accessions).")
    ap.add_argument("--sheet", default=None,
                    help="Excel sheet name (if None, use first sheet).")
    ap.add_argument("--summary-csv", default=None,
                    help="Where to write summary CSV "
                         "(default: <folder_of_out_excel>/00.clean_summary.csv).")

    # NEW: UniProt mapping options
    ap.add_argument("--uniprot-map-tsv", default=None,
                    help="Optional TSV file mapping UniProt entry names to accessions.")
    ap.add_argument("--map-entry-col", default=None,
                    help="Column name in mapping TSV with entry names "
                         "(e.g. 'EntryName'). Required if --uniprot-map-tsv is given.")
    ap.add_argument("--map-acc-col", default=None,
                    help="Column name in mapping TSV with accessions "
                         "(e.g. 'Accession'). Required if --uniprot-map-tsv is given.")

    args = ap.parse_args()

    # 1) Load training UniProt accessions
    train_uniprots = load_train_uniprots(args.train)
    print(f"[i] Loaded {len(train_uniprots)} UniProt ACCs from training file: {args.train}")

    # 2) Load optional UniProt mapping
    uniprot_map: Dict[str, Set[str]] = {}
    use_mapping = args.uniprot_map_tsv is not None
    if use_mapping:
        if not (args.map_entry_col and args.map_acc_col):
            raise ValueError(
                "When using --uniprot-map-tsv, you must also provide "
                "--map-entry-col and --map-acc-col."
            )
        uniprot_map = load_uniprot_map(args.uniprot_map_tsv,
                                       args.map_entry_col,
                                       args.map_acc_col)
        print(f"[i] Loaded UniProt mapping from {args.uniprot_map_tsv}")
        print(f"    Unique entry names in mapping: {len(uniprot_map)}")

    # 3) Load Excel
    print(f"[i] Reading Excel: {args.excel}")
    if args.sheet is not None:
        df = pd.read_excel(args.excel, sheet_name=args.sheet)
    else:
        df = pd.read_excel(args.excel)

    if args.uniprot_col not in df.columns:
        raise ValueError(
            f"Column {args.uniprot_col!r} not found in Excel. "
            f"Available columns: {list(df.columns)}"
        )

    # 4) Normalize UniProt IDs with optional mapping
    df["_uniprot_raw"] = df[args.uniprot_col]

    norm_list: List[str] = []      # primary ACC per row (or None)
    acc_list_str: List[str] = []   # ';'-joined list of mapped ACCs
    in_train_flags: List[bool] = []
    n_mapped_rows = 0

    for val in df[args.uniprot_col].tolist():
        raw = str(val).strip() if pd.notna(val) else ""
        if not raw or raw.upper() in ("NAN", "NONE"):
            norm_list.append(None)
            acc_list_str.append("")
            in_train_flags.append(False)
            continue

        key = raw.upper()

        # Step 1: try mapping entry name → accessions
        accs: List[str] = []
        if use_mapping and key in uniprot_map:
            accs = sorted(uniprot_map[key])
            n_mapped_rows += 1
        else:
            # Treat the value as an accession directly
            accs = [key]

        # Use the first accession as "normalized" representative
        primary_acc = accs[0] if accs else None
        norm_list.append(primary_acc)
        acc_list_str.append(";".join(accs))

        # Flag: is ANY of the accessions in training?
        in_train = any(acc in train_uniprots for acc in accs)
        in_train_flags.append(in_train)

    df["_uniprot_norm"] = norm_list
    df["_uniprot_acc_list"] = acc_list_str
    df["_in_train"] = in_train_flags

    if use_mapping:
        print(f"[i] UniProt mapping applied to {n_mapped_rows} rows "
              f"using entry-name → accession TSV.")

    # 5) Build CLEAN subset = rows whose UniProt ACC is NOT in training
    df_clean = df[~df["_in_train"]].copy()

    # 6) Row-level and UniProt-level statistics
    total_rows = len(df)
    rows_with_uniprot = sum(1 for x in norm_list if x is not None)
    rows_in_train = sum(1 for flag in in_train_flags if flag)

    uniq_unip_all = {u for u in norm_list if isinstance(u, str)}
    uniq_unip_train = uniq_unip_all & train_uniprots
    n_uniq_unip_all = len(uniq_unip_all)
    n_uniq_unip_train = len(uniq_unip_train)

    pct_rows_in_train = (rows_in_train / rows_with_uniprot * 100.0) if rows_with_uniprot else 0.0
    pct_unip_in_train = (n_uniq_unip_train / n_uniq_unip_all * 100.0) if n_uniq_unip_all else 0.0

    print()
    print(f"[summary] Dataset: {args.dataset_name}")
    print(f"  Total rows:                         {total_rows}")
    print(f"  Rows with a valid UniProt ID:       {rows_with_uniprot}")
    if use_mapping:
        print(f"    (rows mapped via TSV:             {n_mapped_rows})")
    print(f"  Rows with UniProt in training:      {rows_in_train} "
          f"({pct_rows_in_train:.1f} % of rows with UniProt)")
    print(f"  Unique UniProt ACCs in table:       {n_uniq_unip_all}")
    print(f"  Unique UniProt ACCs in training:    {n_uniq_unip_train} "
          f"({pct_unip_in_train:.1f} % of UniProt in table)")
    print(f"  CLEAN rows (not in training set):   {len(df_clean)}")

    # 7) PDB-level stats (if column 'PDB' exists)
    pdb_stats = {}
    if "PDB" in df.columns:
        df["_pdb_id"] = df["PDB"].apply(extract_pdb_id)
        total_pdb = {p for p in df["_pdb_id"] if p}
        pdb_with_train = {p for (p, flag) in zip(df["_pdb_id"], in_train_flags) if p and flag}

        n_total_pdb = len(total_pdb)
        n_pdb_train = len(pdb_with_train)
        pct_pdb_train = (n_pdb_train / n_total_pdb * 100.0) if n_total_pdb else 0.0

        pdb_stats = {
            "n_total_pdb": n_total_pdb,
            "n_pdb_with_train": n_pdb_train,
            "pct_pdb_with_train": pct_pdb_train,
        }

        print()
        print("  [PDB statistics]")
        print(f"    Unique PDB IDs in table:              {n_total_pdb}")
        print(f"    PDBs with ≥1 training ACC:            {n_pdb_train} "
              f"({pct_pdb_train:.1f} % of PDBs)")

    # 8) Save CLEAN Excel (drop helper columns)
    drop_cols = ["_uniprot_raw", "_uniprot_norm", "_uniprot_acc_list", "_in_train", "_pdb_id"]
    drop_cols = [c for c in drop_cols if c in df_clean.columns]
    df_clean_out = df_clean.drop(columns=drop_cols)

    os.makedirs(os.path.dirname(args.out_excel), exist_ok=True)
    df_clean_out.to_excel(args.out_excel, index=False)
    print()
    print(f"[ok] CLEAN Excel written to: {args.out_excel}")

    # 9) Save summary CSV
    if args.summary_csv is None:
        base_dir = os.path.dirname(args.out_excel) or "."
        args.summary_csv = os.path.join(base_dir, "00.clean_summary.csv")

    summary_rows = [{
        "dataset": args.dataset_name,
        "excel_path": args.excel,
        "train_path": args.train,
        "uniprot_col": args.uniprot_col,
        "uniprot_map_tsv": args.uniprot_map_tsv or "",
        "map_entry_col": args.map_entry_col or "",
        "map_acc_col": args.map_acc_col or "",
        "total_rows": total_rows,
        "rows_with_uniprot": rows_with_uniprot,
        "rows_mapped_via_tsv": n_mapped_rows if use_mapping else "",
        "rows_in_train": rows_in_train,
        "pct_rows_in_train": pct_rows_in_train,
        "n_uniq_unip_all": n_uniq_unip_all,
        "n_uniq_unip_train": n_uniq_unip_train,
        "pct_unip_in_train": pct_unip_in_train,
        "n_clean_rows": len(df_clean),
        "n_total_pdb": pdb_stats.get("n_total_pdb", ""),
        "n_pdb_with_train": pdb_stats.get("n_pdb_with_train", ""),
        "pct_pdb_with_train": pdb_stats.get("pct_pdb_with_train", ""),
    }]

    os.makedirs(os.path.dirname(args.summary_csv), exist_ok=True)
    write_header = not os.path.isfile(args.summary_csv)
    with open(args.summary_csv, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        if write_header:
            w.writeheader()
        for r in summary_rows:
            w.writerow(r)

    print(f"[ok] Summary statistics appended to: {args.summary_csv}")
    print("Done. CLEAN subset can now be used as input for the 00/01/02/03 pipeline.")

if __name__ == "__main__":
    main()
