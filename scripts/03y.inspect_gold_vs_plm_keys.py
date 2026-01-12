#!/usr/bin/env python3
"""
03a.inspect_gold_vs_plm_keys.py

Purpose
-------
Quick sanity-check between GOLD chains (from 01.chain_gold_labels.json) and
raw PLM chains (from 03.plm_predictions.json). It answers:
  - Which (PDB,CHAIN) are in GOLD but missing in PLM?
  - Which (PDB,CHAIN) are in PLM but not in GOLD?
  - For GOLD-missing chains, is the PDB present in PLM at all? Which PLM chains exist?
  - Heuristic suggestions for aliasing (e.g., AAA -> A, BBB -> B, lowercase -> UPPERCASE).

Outputs
-------
  * 03.plm_missing_from_plm.csv   – GOLD keys that PLM doesn't contain (with hints)
  * 03.plm_extra_in_plm.csv       – PLM keys that GOLD doesn't contain
  * 03.plm_alias_suggestions.json – per-PDB suggested chain alias map (heuristic)
  * Console summary (counts + first few examples)

Notes
-----
- All comparisons use UPPERCASE for PDB and CHAIN.
- This script does NOT modify your data. It only reports & suggests.
- If suggestions look correct, you can feed the alias map into your PLM alignment
  step (03.chain_plm_labels.py) in the future to try remapping chains.
"""

import os, json, csv
from collections import defaultdict
from typing import Dict, List, Tuple, Set

# ---------- CONFIG ----------
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"
GOLD_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")

PLM_RAW_JSON_CANDIDATES = [
    os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_predictions.json"),
]

MISSING_CSV  = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_missing_from_plm.csv")
EXTRA_CSV    = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_extra_in_plm.csv")
ALIAS_JSON   = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_alias_suggestions.json")

PRINT_PREVIEW = 12  # how many examples to show in console
# ----------------------------

def up_pdb(x: str) -> str:   return str(x).strip()[:4].upper()
def up_chain(x: str) -> str: return (str(x).strip() or "A").upper()

def load_gold_keys(path: str) -> List[Tuple[str,str]]:
    data = json.load(open(path))
    seen = set()
    order = []
    for e in data:
        key = (up_pdb(e["pdb_id"]), up_chain(e["chain_id"]))
        if key not in seen:
            seen.add(key); order.append(key)
    return order

def load_plm_keys(paths: List[str]) -> List[Tuple[str,str]]:
    fn = None
    for cand in paths:
        if os.path.isfile(cand):
            fn = cand
            break
    if fn is None:
        raise FileNotFoundError("Could not find PLM raw JSON in candidates: " + ", ".join(paths))
    raw = json.load(open(fn))
    keys = []
    seen = set()
    for e in raw:
        pdb = e.get("pdb_id") or e.get("pdb") or e.get("pdbid")
        ch  = e.get("chain_id") or e.get("chain") or e.get("chainid")
        if pdb is None or ch is None:
            continue
        key = (up_pdb(pdb), up_chain(ch))
        if key not in seen:
            seen.add(key); keys.append(key)
    return keys

def group_by_pdb(keys: List[Tuple[str,str]]) -> Dict[str, Set[str]]:
    d = defaultdict(set)
    for pdb, ch in keys:
        d[pdb].add(ch)
    return d

def heuristic_chain_alias(chain: str) -> List[str]:
    """
    Very conservative suggestions:
      - AAA -> A, BBB -> B, CCC -> C (if all letters the same)
      - lowercase -> uppercase (already normalized, so no-op here)
      - try first char if chain has repeated same char patterns (e.g., 'AAAA' -> 'A')
    We return a list of *candidates*; user must validate.
    """
    s = chain.upper()
    cands = set()

    # If chain is all the same letter repeated (AAA, BBB, CCCC), suggest its single letter
    if len(s) >= 2 and len(set(s)) == 1:
        cands.add(s[0])

    # If chain is long and starts with a letter repeated, also suggest first letter
    if len(s) >= 3 and s.count(s[0]) >= 2:
        cands.add(s[0])

    # If chain contains underscores like "A_1" (rare), suggest the part before underscore
    if "_" in s:
        cands.add(s.split("_", 1)[0])

    # Return in deterministic order
    return sorted(cands)

def main():
    gold_keys = load_gold_keys(GOLD_JSON)
    plm_keys  = load_plm_keys(PLM_RAW_JSON_CANDIDATES)

    gold_set = set(gold_keys)
    plm_set  = set(plm_keys)

    gold_only = sorted(gold_set - plm_set)   # in GOLD but not in PLM
    plm_only  = sorted(plm_set - gold_set)   # in PLM but not in GOLD

    # Group by PDB for hints
    gold_by_pdb = group_by_pdb(gold_keys)
    plm_by_pdb  = group_by_pdb(plm_keys)

    # Build rows for missing report and alias suggestions (per PDB)
    missing_rows = []
    alias_map_per_pdb: Dict[str, Dict[str, List[str]]] = defaultdict(dict)
    for pdb, ch in gold_only:
        plm_has_pdb = pdb in plm_by_pdb
        plm_chains  = sorted(plm_by_pdb.get(pdb, []))
        # Heuristic alias candidates (e.g., AAA -> A)
        suggested = heuristic_chain_alias(ch)
        # Keep suggestions only if target chain actually exists in PLM for that PDB
        suggested_present = [c for c in suggested if c in plm_chains]

        if suggested_present:
            alias_map_per_pdb[pdb][ch] = suggested_present

        missing_rows.append({
            "pdb_id": pdb,
            "chain_id_gold": ch,
            "plm_pdb_present": plm_has_pdb,
            "plm_chains_available": " ".join(plm_chains),
            "alias_candidates": " ".join(suggested),
            "alias_candidates_present": " ".join(suggested_present)
        })

    # Build rows for extra report
    extra_rows = [{"pdb_id": pdb, "chain_id_plm": ch} for pdb, ch in plm_only]

    # Write CSVs
    with open(MISSING_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "pdb_id","chain_id_gold","plm_pdb_present","plm_chains_available",
            "alias_candidates","alias_candidates_present"
        ])
        w.writeheader()
        for r in missing_rows: w.writerow(r)

    with open(EXTRA_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id","chain_id_plm"])
        w.writeheader()
        for r in extra_rows: w.writerow(r)

    # Write alias suggestions (per PDB)
    with open(ALIAS_JSON, "w") as f:
        json.dump(alias_map_per_pdb, f, indent=2)

    # Console summary
    print(f"GOLD chains: {len(gold_set)} | PLM chains: {len(plm_set)}")
    print(f"Missing in PLM (GOLD-only): {len(gold_only)}  → {MISSING_CSV}")
    if gold_only:
        preview = ", ".join([f"{p}.{c}" for p,c in gold_only[:PRINT_PREVIEW]])
        print("  e.g.:", preview)
    print(f"Extra in PLM (PLM-only):    {len(plm_only)}  → {EXTRA_CSV}")
    if plm_only:
        preview = ", ".join([f"{p}.{c}" for p,c in plm_only[:PRINT_PREVIEW]])
        print("  e.g.:", preview)
    print(f"Alias suggestions (per PDB): → {ALIAS_JSON}")
    if alias_map_per_pdb:
        # show first pdb with suggestions
        first_pdb = sorted(alias_map_per_pdb.keys())[0]
        print(f"  example for {first_pdb}: {alias_map_per_pdb[first_pdb]}")

if __name__ == "__main__":
    main()
