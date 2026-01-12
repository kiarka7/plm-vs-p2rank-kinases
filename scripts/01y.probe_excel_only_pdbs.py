#!/usr/bin/env python3
"""
01y.probe_excel_only_pdbs.py

Goal:
- Find PDBs that are present in Excel and structures but have 0 GOLD entries.
- For those PDBs, show:
  * raw Excel 'Ligand' cells (all rows for that PDB)
  * tokens we would parse from those cells
  * chains available in the mmCIF and non-water HET components present

Outputs:
- Prints a concise console report
- Writes ./KINASE_FOLDER/diag/01.excel_only_probe.csv for later review
"""

import os, re, json, csv
import pandas as pd
from collections import Counter, defaultdict
from Bio.PDB import MMCIFParser, is_aa

# ==== CONFIG (match your dataset) ====
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"
EXCEL_FILE    = os.path.join(BASE_DIR, KINASE_FOLDER, "Kinase_Ligands_Type I.xlsx")
STRUCT_DIR    = os.path.join(BASE_DIR, KINASE_FOLDER, "structures")
OVERVIEW_CSV  = os.path.join(BASE_DIR, KINASE_FOLDER, "00.overview_per_pdb.csv")
OUT_CSV       = os.path.join(BASE_DIR, KINASE_FOLDER, "diag", "01.excel_only_probe.csv")
WATER_CODES   = {"HOH","WAT","DOD"}
LIGAND_TOKEN  = re.compile(r"^[A-Za-z0-9_]+:\d+$")
# =====================================

def iter_tokens(cell: str):
    if not cell:
        return
    for raw in str(cell).replace(";", " ").replace(",", " ").split():
        tok = raw.strip()
        if LIGAND_TOKEN.match(tok):
            yield tok

def extract_nameonly_codes(cell: str):
    if not cell:
        return []
    out = []
    for raw in str(cell).replace(";", " ").replace(",", " ").split():
        s = raw.strip()
        if not s or s.lower() in ("none","no_ligand","nan"):
            continue
        if LIGAND_TOKEN.match(s):
            continue
        if re.fullmatch(r"[A-Za-z0-9_]{2,5}", s):
            out.append(s.upper())
    return out

def get_cif_hets(pdb_id: str):
    """Return available chains and non-water HET comp_ids (counts) from mmCIF."""
    cif = os.path.join(STRUCT_DIR, f"{pdb_id.lower()}.cif")
    if not os.path.isfile(cif):
        return [], Counter()
    parser = MMCIFParser(QUIET=True)
    try:
        struct = parser.get_structure(pdb_id, cif)[0]
    except Exception:
        return [], Counter()
    chains = [ch.id for ch in struct]
    het_counts = Counter()
    for ch in struct:
        for res in ch:
            if is_aa(res, standard=False):
                continue
            name = res.get_resname().upper()
            if name in WATER_CODES:
                continue
            het_counts[name] += 1
    return chains, het_counts

def main():
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

    df_over = pd.read_csv(OVERVIEW_CSV)
    # Robust pick: excel_rows + gold_n_chains must exist here
    excel_only = df_over[(df_over["excel_rows"] > 0) & (df_over["gold_n_chains"].fillna(0).astype(int) == 0)]
    pdb_list = excel_only["pdb_id"].astype(str).str.upper().tolist()

    if not os.path.isfile(EXCEL_FILE):
        print(f"[!] Excel file missing: {EXCEL_FILE}")
        return
    df_x = pd.read_excel(EXCEL_FILE)
    # Expect a column 'PDB' containing "4CHAR" or "4CHAR + chain"
    df_x = df_x[df_x["PDB"].notna()].copy()
    df_x["PDB"] = df_x["PDB"].astype(str).str.strip()

    rows_out = []

    for pdb in pdb_list:
        rows = df_x[df_x["PDB"].str.upper().str.startswith(pdb)]
        lig_cells = rows["Ligand"].astype(str).tolist() if "Ligand" in rows.columns else []
        tokens = []
        names_only = []
        for s in lig_cells:
            tokens += list(iter_tokens(s))
            names_only += extract_nameonly_codes(s)

        chains, het_counts = get_cif_hets(pdb)
        top_hets = ", ".join([f"{k}:{v}" for k,v in het_counts.most_common(12)])

        print(f"\n=== {pdb} ===")
        print(f" Excel rows: {len(rows)}")
        print(f" Ligand cells (first 3): {lig_cells[:3]}")
        print(f" Parsed tokens: {tokens if tokens else '—'}")
        print(f" Name-only codes: {names_only if names_only else '—'}")
        print(f" CIF chains: {','.join(chains) if chains else '—'}")
        print(f" Non-water HETs: {top_hets if top_hets else '—'}")

        rows_out.append({
            "pdb_id": pdb,
            "excel_rows": len(rows),
            "excel_ligand_cells_first3": " | ".join(lig_cells[:3]),
            "parsed_tokens": " ".join(tokens),
            "nameonly_codes": " ".join(names_only),
            "cif_chains": ",".join(chains),
            "nonwater_hets_top12": top_hets
        })

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\n[ok] Wrote probe CSV → {OUT_CSV}")

if __name__ == "__main__":
    main()
