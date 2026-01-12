#!/usr/bin/env python3
"""
01.2_gold_audit_summary.py — end-to-end diagnostics for GOLD builder

What this script does:
- Reads:
    * Excel table with PDB(+optional chain) and Ligand columns
    * GOLD JSON produced by 01.chain_gold_labels.py
    * Audit CSV produced by 01.chain_gold_labels.py
    * Local structures folder with .cif / .cif.gz
- Produces:
    * Text summary (counts, sanity checks)
    * CSV breakdowns:
        - rows used / skipped / not_found / file_missing
        - per-(pdb,chain) summary (counts of each status)
        - per-component (alias_used) usage
        - PDBs missing in structures, Excel-vs-GOLD overlaps, etc.

Why this is useful:
- Quickly verify that GOLD is consistent with the input Excel and structures you actually have.
- Spot which PDB/chain/ligand tokens failed to resolve and why.

Outputs (written to ./<KINASE_FOLDER>/diag/):
    - 01.audit_summary.txt
    - 01.audit_used.csv
    - 01.audit_not_found.csv
    - 01.audit_file_missing.csv
    - 01.audit_skipped.csv
    - 01.audit_by_pdb_chain_summary.csv
    - 01.audit_by_component.csv
    - 01.audit_excel_vs_structures_vs_gold.csv
"""

import os
import json
import gzip
import pandas as pd
from collections import Counter, defaultdict

# -------- CONFIG --------
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"  # adjust per dataset
EXCEL_FILE    = os.path.join(BASE_DIR, KINASE_FOLDER, "Kinase_Ligands_Type I.xlsx")
STRUCT_DIR    = os.path.join(BASE_DIR, KINASE_FOLDER, "structures")
GOLD_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
AUDIT_CSV     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.ligand_audit.csv")

OUT_DIR       = os.path.join(BASE_DIR, KINASE_FOLDER, "diag")
os.makedirs(OUT_DIR, exist_ok=True)
# ------------------------


def list_struct_pdbs(struct_dir: str):
    """Return set of PDB IDs (upper) for which a .cif or .cif.gz exists."""
    if not os.path.isdir(struct_dir):
        return set()
    pdbs = set()
    for fn in os.listdir(struct_dir):
        f = fn.lower()
        if f.endswith(".cif") or f.endswith(".cif.gz"):
            pdb = f.split(".")[0][:4].upper()
            if len(pdb) == 4:
                pdbs.add(pdb)
    return pdbs


def load_gold(path: str):
    """Load GOLD JSON; return list and some simple views."""
    data = json.load(open(path))
    # unique (pdb, chain)
    keys = {(str(e["pdb_id"]).upper(), str(e["chain_id"]).upper()) for e in data}
    pdbs = {p for p, _ in keys}
    return data, keys, pdbs


def norm_pdb_chain_from_excel_cell(cell: str):
    """
    Excel column 'PDB' may be '3CTQ' or '3CTQA'. Return (pdb, chain|None).
    If chain missing, return None (meaning "A" was assumed by GOLD builder).
    """
    s = str(cell).strip()
    if not s:
        return None, None
    pdb = s[:4].upper()
    chain = s[4:].strip().upper() or None
    return pdb, chain


def main():
    # ----- Load Excel -----
    df_x = pd.read_excel(EXCEL_FILE)
    # keep rows that actually have a PDB code and a non-empty Ligand cell
    df_x = df_x[df_x["PDB"].notna()]
    df_x = df_x[df_x["Ligand"].notna()]
    df_x["PDB_norm"]   = df_x["PDB"].astype(str).str.strip().str[:4].str.upper()
    df_x["CHAIN_hint"] = df_x["PDB"].astype(str).str[4:].str.strip().str.upper().replace({"": None})
    excel_rows = len(df_x)
    excel_pdbs = set(df_x["PDB_norm"].unique())

    # ----- Structures present -----
    struct_pdbs = list_struct_pdbs(STRUCT_DIR)

    # ----- GOLD -----
    gold, gold_keys, gold_pdbs = load_gold(GOLD_JSON)
    gold_rows = len(gold)

    # ----- Audit CSV -----
    df_a = pd.read_csv(AUDIT_CSV)
    # normalize some columns just in case
    for col in ["pdb", "chain_req", "alias_used", "status", "found_mode"]:
        if col in df_a.columns:
            df_a[col] = df_a[col].astype(str).str.strip()
    audit_rows = len(df_a)

    # ----------------- Build summary text -----------------
    lines = []
    lines.append(f"Excel rows (non-empty Ligand): {excel_rows}")
    lines.append(f"Excel unique PDB:             {len(excel_pdbs)}")
    lines.append(f"Structures present (.cif*):   {len(struct_pdbs)}")
    lines.append(f"GOLD JSON entries:            {gold_rows}")
    lines.append(f"GOLD unique (pdb,chain):      {len(gold_keys)}")
    lines.append(f"GOLD unique PDB:              {len(gold_pdbs)}")
    lines.append(f"Audit CSV rows:               {audit_rows}")
    lines.append("")

    # overlaps
    miss_struct_from_excel = sorted(excel_pdbs - struct_pdbs)
    extra_struct_not_in_x  = sorted(struct_pdbs - excel_pdbs)
    excel_not_in_gold      = sorted(excel_pdbs - gold_pdbs)
    gold_not_in_excel      = sorted(gold_pdbs - excel_pdbs)

    lines.append(f"Excel PDB missing structures: {len(miss_struct_from_excel)}")
    if miss_struct_from_excel[:10]:
        lines.append("  e.g.: " + ", ".join(miss_struct_from_excel[:10]))
    lines.append(f"Structures not in Excel:      {len(extra_struct_not_in_x)}")
    if extra_struct_not_in_x[:10]:
        lines.append("  e.g.: " + ", ".join(extra_struct_not_in_x[:10]))
    lines.append(f"Excel PDB not in GOLD:        {len(excel_not_in_gold)}")
    if excel_not_in_gold[:10]:
        lines.append("  e.g.: " + ", ".join(excel_not_in_gold[:10]))
    lines.append(f"GOLD PDB not in Excel:        {len(gold_not_in_excel)}")
    if gold_not_in_excel[:10]:
        lines.append("  e.g.: " + ", ".join(gold_not_in_excel[:10]))
    lines.append("")

    # audit distributions
    status_counts = df_a["status"].value_counts().to_dict() if "status" in df_a.columns else {}
    mode_counts   = df_a["found_mode"].value_counts().to_dict() if "found_mode" in df_a.columns else {}
    cls_counts    = df_a["class"].value_counts().to_dict() if "class" in df_a.columns else {}

    lines.append("[audit] status counts:")
    for k, v in status_counts.items():
        lines.append(f"  - {k:12s}: {v}")
    lines.append("[audit] found_mode counts:")
    for k, v in mode_counts.items():
        lines.append(f"  - {k:12s}: {v}")
    lines.append("[audit] class counts:")
    for k, v in cls_counts.items():
        lines.append(f"  - {k:12s}: {v}")
    lines.append("")

    # Top offenders: (pdb,chain) with most "not_found"
    if "status" in df_a.columns:
        nf = df_a[df_a["status"] == "not_found"].copy()
        if not nf.empty:
            grp = nf.groupby(["pdb", "chain_req"]).size().reset_index(name="n_not_found")
            grp = grp.sort_values("n_not_found", ascending=False)
            lines.append("Top (pdb,chain) with NOT_FOUND tokens:")
            for _, row in grp.head(10).iterrows():
                lines.append(f"  - {row['pdb']}.{row['chain_req']}: {int(row['n_not_found'])}")
            lines.append("")

    # Write text summary
    summary_path = os.path.join(OUT_DIR, "01.audit_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(lines))
    print(f"[ok] Wrote summary: {summary_path}")

    # ----------------- CSV breakdowns -----------------
    # Basic splits
    if "status" in df_a.columns:
        df_used  = df_a[df_a["status"] == "used"].copy()
        df_nf    = df_a[df_a["status"] == "not_found"].copy()
        df_miss  = df_a[df_a["status"] == "file_missing"].copy()
        df_skip  = df_a[df_a["status"] == "skipped"].copy()
    else:
        df_used = df_nf = df_miss = df_skip = pd.DataFrame()

    def save_csv(df, name):
        path = os.path.join(OUT_DIR, name)
        df.to_csv(path, index=False)
        print(f"[ok] {name}: {len(df)} rows")

    save_csv(df_used,  "01.audit_used.csv")
    save_csv(df_nf,    "01.audit_not_found.csv")
    save_csv(df_miss,  "01.audit_file_missing.csv")
    save_csv(df_skip,  "01.audit_skipped.csv")

    # Per-(pdb,chain) summary
    if not df_a.empty:
        cols = ["pdb", "chain_req", "status"]
        have = [c for c in cols if c in df_a.columns]
        g = df_a[have].groupby(["pdb", "chain_req", "status"]).size().unstack(fill_value=0)
        g["total_tokens"] = g.sum(axis=1)
        g = g.reset_index()
        save_csv(g, "01.audit_by_pdb_chain_summary.csv")

    # Per component usage
    if "alias_used" in df_a.columns:
        comp = df_a.groupby(["alias_used", "status"]).size().unstack(fill_value=0).reset_index()
        save_csv(comp, "01.audit_by_component.csv")

    # Excel vs structures vs GOLD (per-PDB)
    # summarize chains present in GOLD per PDB
    chains_per_pdb = defaultdict(set)
    for e in gold:
        chains_per_pdb[str(e["pdb_id"]).upper()].add(str(e["chain_id"]).upper())
    rows = []
    for pdb in sorted(excel_pdbs | struct_pdbs | gold_pdbs):
        rows.append({
            "pdb_id": pdb,
            "in_excel": int(pdb in excel_pdbs),
            "has_structure": int(pdb in struct_pdbs),
            "gold_n_chains": len(chains_per_pdb.get(pdb, set())),
            "gold_chains": " ".join(sorted(chains_per_pdb.get(pdb, set())))
        })
    df_over = pd.DataFrame(rows)
    save_csv(df_over, "01.audit_excel_vs_structures_vs_gold.csv")


if __name__ == "__main__":
    main()
