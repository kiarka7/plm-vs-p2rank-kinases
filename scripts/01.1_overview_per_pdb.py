#!/usr/bin/env python3
"""
Quick sanity audit: Excel vs downloaded structures vs GOLD chain labels.

What it checks
--------------
- Reads Excel and extracts PDB IDs (robust to formats like '1ABC', '1ABC_A', '1ABCA', spaces).
- Lists structures present in a folder (.cif and .cif.gz).
- Reads GOLD JSON (01.chain_gold_labels.json) and collects (PDB, chain) pairs.

Outputs
-------
- Console summary with counts and deltas.
- CSVs in --outdir:
    * 00.missing_structures_from_excel.csv          (PDB present in Excel but no .cif file)
    * 00.extra_structures_not_in_excel.csv          (PDB with .cif but not in Excel)
    * 00.gold_chains_per_pdb.csv                    (per PDB chain list from GOLD)
    * 00.duplicates_in_excel.csv                    (if multiple Excel rows per PDB)
    * 00.overview_per_pdb.csv                       (Excel rows, has_structure, #GOLD chains, chain_ids)

Notes
-----
- Big discrepancies are often benign: Excel counts "evidences" (rows), structures count PDB files,
  GOLD counts (PDB, chain) pairs near ligands.
"""

import os, re, json, argparse, csv
import pandas as pd
from collections import Counter, defaultdict

DEF_BASE = "."
DEF_FOLDER = "Kinase_Type_I"
DEF_EXCEL = os.path.join(DEF_BASE, DEF_FOLDER, "Kinase_Ligands_Type I.xlsx")
DEF_STRUCT = os.path.join(DEF_BASE, DEF_FOLDER, "structures")
DEF_GOLD = os.path.join(DEF_BASE, DEF_FOLDER, "01.chain_gold_labels.json")
DEF_OUTDIR = os.path.join(DEF_BASE, DEF_FOLDER)

PDB_RE = re.compile(r"([0-9A-Za-z]{4})")

def extract_pdb_id(s: str) -> str | None:
    """Return 4-char PDB ID (upper) from a messy cell like '1CKJ', '1CKJ_A', ' 1CKJ B', '1ckja'."""
    if not isinstance(s, str):
        s = str(s)
    s = s.strip()
    m = PDB_RE.search(s)
    return m.group(1).upper() if m else None

def load_excel_pdbs(excel_path: str) -> list[str]:
    df = pd.read_excel(excel_path)
    # heuristically pick first column that looks like PDB
    candidates = [c for c in df.columns if "pdb" in str(c).lower()]
    col = candidates[0] if candidates else df.columns[0]
    ids = []
    for v in df[col].tolist():
        pid = extract_pdb_id(str(v))
        if pid:
            ids.append(pid)
    return ids  # may contain duplicates (multiple rows per PDB)

def list_structure_pdbs(struct_dir: str) -> set[str]:
    out = set()
    if not os.path.isdir(struct_dir):
        return out
    for fn in os.listdir(struct_dir):
        low = fn.lower()
        if low.endswith(".cif") or low.endswith(".cif.gz"):
            # take first 4 alphanum chars
            m = PDB_RE.search(fn)
            if m:
                out.add(m.group(1).upper())
    return out

def load_gold_pairs(gold_json: str) -> list[tuple[str,str]]:
    data = json.load(open(gold_json))
    pairs = []
    for e in data:
        pdb = str(e["pdb_id"]).upper()
        ch  = str(e["chain_id"]).upper()
        pairs.append((pdb, ch))
    return pairs

def write_csv(path: str, rows: list[dict], fieldnames: list[str]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel", default=DEF_EXCEL, help="Excel file with PDB list")
    ap.add_argument("--structures", default=DEF_STRUCT, help="Folder with .cif / .cif.gz")
    ap.add_argument("--gold", default=DEF_GOLD, help="01.chain_gold_labels.json")
    ap.add_argument("--outdir", default=DEF_OUTDIR, help="Where to write CSV reports")
    args = ap.parse_args()

    # Excel
    excel_ids = load_excel_pdbs(args.excel)
    excel_rows = len(excel_ids)
    excel_unique = sorted(set(excel_ids))
    excel_dupes = [p for p,count in Counter(excel_ids).items() if count > 1]

    # Structures
    struct_pdbs = sorted(list_structure_pdbs(args.structures))

    # GOLD
    gold_pairs = load_gold_pairs(args.gold)  # list of (pdb, chain)
    gold_pdbs = sorted(set(p for p,_ in gold_pairs))
    chains_by_pdb = defaultdict(set)
    for p, c in gold_pairs:
        chains_by_pdb[p].add(c)

    # Deltas
    missing_struct = sorted(set(excel_unique) - set(struct_pdbs))
    extra_struct   = sorted(set(struct_pdbs) - set(excel_unique))
    gold_only_pdbs = sorted(set(gold_pdbs) - set(excel_unique))
    excel_only_pdbs= sorted(set(excel_unique) - set(gold_pdbs))

    # Console summary
    print(f"Excel rows: {excel_rows}  | Excel unique PDB: {len(excel_unique)}")
    print(f"Structures: {len(struct_pdbs)} (unique PDB with .cif/.cif.gz)")
    print(f"GOLD chain entries: {len(gold_pairs)}  | GOLD unique PDB: {len(gold_pdbs)}")
    print()
    print(f"- Excel duplicates (same PDB multiple rows): {len(excel_dupes)}")
    if excel_dupes:
        print(f"  e.g.: {', '.join(excel_dupes[:12])}")
    print(f"- Missing structures for Excel PDBs: {len(missing_struct)}")
    if missing_struct:
        print(f"  e.g.: {', '.join(missing_struct[:12])}")
    print(f"- Extra structures not in Excel: {len(extra_struct)}")
    if extra_struct:
        print(f"  e.g.: {', '.join(extra_struct[:12])}")
    print(f"- GOLD-only PDBs (in GOLD but not Excel): {len(gold_only_pdbs)}")
    if gold_only_pdbs:
        print(f"  e.g.: {', '.join(gold_only_pdbs[:12])}")
    print(f"- Excel-only PDBs (in Excel but not GOLD): {len(excel_only_pdbs)}")
    if excel_only_pdbs:
        print(f"  e.g.: {', '.join(excel_only_pdbs[:12])}")

    # CSVs
    write_csv(
        os.path.join(args.outdir, "00.missing_structures_from_excel.csv"),
        [{"pdb_id": p} for p in missing_struct],
        ["pdb_id"]
    )
    write_csv(
        os.path.join(args.outdir, "00.extra_structures_not_in_excel.csv"),
        [{"pdb_id": p} for p in extra_struct],
        ["pdb_id"]
    )
    write_csv(
        os.path.join(args.outdir, "00.duplicates_in_excel.csv"),
        [{"pdb_id": p, "excel_row_count": Counter(excel_ids)[p]} for p in excel_dupes],
        ["pdb_id", "excel_row_count"]
    )
    # GOLD per-PDB chains
    gold_rows = []
    for p in sorted(chains_by_pdb.keys()):
        chs = sorted(chains_by_pdb[p])
        gold_rows.append({"pdb_id": p, "n_chains": len(chs), "chains": " ".join(chs)})
    write_csv(
        os.path.join(args.outdir, "00.gold_chains_per_pdb.csv"),
        gold_rows,
        ["pdb_id", "n_chains", "chains"]
    )
    # Overview per PDB
    overview = []
    excel_counts = Counter(excel_ids)
    excel_all = sorted(set(excel_unique) | set(struct_pdbs) | set(gold_pdbs))
    struct_set = set(struct_pdbs)
    for p in excel_all:
        overview.append({
            "pdb_id": p,
            "excel_rows": int(excel_counts.get(p, 0)),
            "has_structure": int(p in struct_set),
            "gold_n_chains": int(len(chains_by_pdb.get(p, set()))),
            "gold_chains": " ".join(sorted(chains_by_pdb.get(p, set())))
        })
    write_csv(
        os.path.join(args.outdir, "00.overview_per_pdb.csv"),
        overview,
        ["pdb_id", "excel_rows", "has_structure", "gold_n_chains", "gold_chains"]
    )
    print(f"\nReports written to: {args.outdir}")

if __name__ == "__main__":
    main()
