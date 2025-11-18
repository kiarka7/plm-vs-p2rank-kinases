#!/usr/bin/env python3
"""
99.inspect_alignment_for_pdb.py

Goal
----
Quick sanity check that GOLD, P2Rank and PLM are aligned and consistent
for a single PDB across chains.

It checks, for each (PDB, CHAIN) in GOLD:

  * That GOLD / P2Rank / PLM use the same author residue numbers
    (sequence_residue_numbers).
  * Length and position mismatches (if any).
  * Overlap between:
        - GOLD positives (labels==1)
        - P2Rank positives (labels==1, or by prob threshold)
        - PLM positives (prob >= threshold)

It prints a compact per-chain summary and optional CSV output.

Inputs (per kinase folder)
--------------------------
Expected files inside KINASE_FOLDER:

  01.chain_gold_labels.json
  02.chain_p2rank_labels.json
  03.chain_plm_labels.json

Usage
-----
Example:

  python 99.inspect_alignment_for_pdb.py \
      --folder Kinase_Type_I \
      --pdb 4CZT \
      --plm-thr 0.75 \
      --p2-thr 0.5 \
      --out-csv Kinase_Type_I/diag/99.inspect_4CZT.csv

If --out-csv is omitted, only console output is produced.
"""

import os
import json
import argparse
import csv
from typing import Dict, List, Tuple, Optional, Any, Set


# ----------------- helpers -----------------

def up_pdb(x: Any) -> str:
    return str(x).strip()[:4].upper()


def up_chain(x: Any) -> str:
    s = str(x).strip()
    return (s or "A").upper()


def load_gold(path: str, pdb_id: str):
    """
    Load GOLD JSON and return per-chain info only for the given PDB.

    Returns:
        gold_pos:  dict[(pdb,chain)] -> [positions]
        gold_lab:  dict[(pdb,chain)] -> [0/1 labels] (OR-merged over duplicates)
        chains:    sorted list of chain IDs for this PDB
    """
    pdb_up = up_pdb(pdb_id)
    data = json.load(open(path))
    gold_pos: Dict[Tuple[str, str], List[int]] = {}
    gold_lab: Dict[Tuple[str, str], List[int]] = {}

    for e in data:
        p = up_pdb(e["pdb_id"])
        if p != pdb_up:
            continue
        c = up_chain(e["chain_id"])
        key = (p, c)
        pos = [int(x) for x in e.get("sequence_residue_numbers", [])]
        lab = [int(x) for x in e.get("labels", [])]

        if len(pos) != len(lab):
            # skip malformed entry
            continue

        if key not in gold_pos:
            gold_pos[key] = pos
            gold_lab[key] = lab
        else:
            # OR-merge labels over same canonical positions
            # assume positions are the same list as before
            if gold_pos[key] != pos:
                # if positions differ, we keep the first occurrence but warn
                print(f"[!] GOLD positions differ for {p}.{c} among entries; "
                      f"keeping the first one.")
                continue
            merged = []
            for a, b in zip(gold_lab[key], lab):
                merged.append(1 if (a == 1 or b == 1) else 0)
            gold_lab[key] = merged

    chains = sorted({c for (p, c) in gold_pos.keys()})
    return gold_pos, gold_lab, chains


def load_p2rank(path: str, pdb_id: str):
    """
    Load P2Rank per-chain predictions for given PDB.

    Returns:
        p2_pos: dict[(pdb,chain)] -> [positions]
        p2_probs: dict[(pdb,chain)] -> [probabilities]
        p2_labels: dict[(pdb,chain)] -> [0/1 labels] (if present; else empty)
    """
    pdb_up = up_pdb(pdb_id)
    if not os.path.isfile(path):
        print(f"[i] P2Rank file not found: {path}")
        return {}, {}, {}
    data = json.load(open(path))

    p2_pos: Dict[Tuple[str, str], List[int]] = {}
    p2_probs: Dict[Tuple[str, str], List[float]] = {}
    p2_labels: Dict[Tuple[str, str], List[int]] = {}

    for e in data:
        p = up_pdb(e.get("pdb_id", ""))
        if p != pdb_up:
            continue
        c = up_chain(e.get("chain_id", ""))
        key = (p, c)
        pos = [int(x) for x in e.get("sequence_residue_numbers", [])]
        probs = [float(x) for x in e.get("pred_probs", [])]
        labels = [int(x) for x in e.get("labels", [])] if "labels" in e else []

        p2_pos[key] = pos
        p2_probs[key] = probs
        p2_labels[key] = labels

    return p2_pos, p2_probs, p2_labels


def load_plm(path: str, pdb_id: str):
    """
    Load PLM per-chain predictions for given PDB
    (output of 03.chain_plm_labels.py).

    Returns:
        plm_pos: dict[(pdb,chain)] -> [positions]
        plm_probs: dict[(pdb,chain)] -> [probabilities]
    """
    pdb_up = up_pdb(pdb_id)
    if not os.path.isfile(path):
        print(f"[i] PLM file not found: {path}")
        return {}, {}
    data = json.load(open(path))

    plm_pos: Dict[Tuple[str, str], List[int]] = {}
    plm_probs: Dict[Tuple[str, str], List[float]] = {}

    for e in data:
        p = up_pdb(e.get("pdb_id", ""))
        if p != pdb_up:
            continue
        c = up_chain(e.get("chain_id", ""))
        key = (p, c)
        pos = [int(x) for x in e.get("sequence_residue_numbers", [])]
        probs = [float(x) for x in e.get("pred_probs", [])]
        plm_pos[key] = pos
        plm_probs[key] = probs

    return plm_pos, plm_probs


def jaccard(a: Set[int], b: Set[int]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union > 0 else 0.0


# ----------------- main logic -----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True,
                    help="Kinase folder (e.g. Kinase_Type_I, Kinase_ALLO, ...)")
    ap.add_argument("--pdb", required=True,
                    help="4-char PDB ID to inspect (case-insensitive)")
    ap.add_argument("--plm-thr", type=float, default=0.75,
                    help="Probability threshold for PLM positives (default: 0.75)")
    ap.add_argument("--p2-thr", type=float, default=0.5,
                    help="Probability threshold for P2Rank positives if labels are absent (default: 0.5)")
    ap.add_argument("--out-csv", default=None,
                    help="Optional CSV output path for per-chain summary")
    args = ap.parse_args()

    base = "."
    folder = args.folder
    pdb_id = up_pdb(args.pdb)

    gold_path = os.path.join(base, folder, "01.chain_gold_labels.json")
    p2_path   = os.path.join(base, folder, "02.chain_p2rank_labels.json")
    plm_path  = os.path.join(base, folder, "03.chain_plm_labels.json")

    print(f"[i] Inspecting PDB {pdb_id} in folder {folder}")
    print(f"    GOLD:   {gold_path}")
    print(f"    P2Rank: {p2_path}")
    print(f"    PLM:    {plm_path}")
    print()

    if not os.path.isfile(gold_path):
        raise FileNotFoundError(f"GOLD file not found: {gold_path}")

    gold_pos, gold_lab, chains = load_gold(gold_path, pdb_id)
    if not chains:
        print(f"[!] No GOLD entries found for PDB {pdb_id}")
        return

    p2_pos, p2_probs, p2_labels = load_p2rank(p2_path, pdb_id)
    plm_pos, plm_probs = load_plm(plm_path, pdb_id)

    if not plm_pos:
        print(f"[!] No PLM entries found for PDB {pdb_id} (03.chain_plm_labels.json)")
    if not p2_pos:
        print(f"[!] No P2Rank entries found for PDB {pdb_id} (02.chain_p2rank_labels.json)")

    rows_for_csv = []

    print("Chain-by-chain summary")
    print("======================")

    for ch in chains:
        key = (pdb_id, ch)
        g_pos = gold_pos.get(key, [])
        g_lab = gold_lab.get(key, [])
        L = len(g_pos)

        if L == 0:
            print(f"{pdb_id}.{ch}: [!] Empty GOLD positions, skipping")
            continue

        # GOLD positives
        gold_pos_set = {p for p, y in zip(g_pos, g_lab) if y == 1}

        # P2Rank
        p2_p = p2_pos.get(key, [])
        p2_pr = p2_probs.get(key, [])
        p2_lab = p2_labels.get(key, [])
        p2_pos_set = set()

        p2_len_ok = (len(p2_p) == L and p2_p == g_pos)

        if p2_p:
            if p2_len_ok:
                # If labels present and aligned, use them
                if p2_lab and len(p2_lab) == L:
                    p2_pos_set = {p for p, y in zip(p2_p, p2_lab) if y == 1}
                else:
                    # fallback: threshold on probabilities
                    if len(p2_pr) == L:
                        p2_pos_set = {p for p, pr in zip(p2_p, p2_pr) if pr >= args.p2_thr}
                    else:
                        print(f"{pdb_id}.{ch}: [!] P2Rank length mismatch (len(pos)={len(p2_p)}, len(probs)={len(p2_pr)}, L={L})")
            else:
                print(f"{pdb_id}.{ch}: [!] P2Rank positions not aligned with GOLD "
                      f"(len GOLD={L}, len P2={len(p2_p)})")
        # PLM
        pl_p = plm_pos.get(key, [])
        pl_pr = plm_probs.get(key, [])
        pl_pos_set = set()
        pl_len_ok = (len(pl_p) == L and pl_p == g_pos)

        if pl_p:
            if pl_len_ok:
                if len(pl_pr) == L:
                    pl_pos_set = {p for p, pr in zip(pl_p, pl_pr) if pr >= args.plm_thr}
                else:
                    print(f"{pdb_id}.{ch}: [!] PLM length mismatch (len(pos)={len(pl_p)}, len(probs)={len(pl_pr)}, L={L})")
            else:
                print(f"{pdb_id}.{ch}: [!] PLM positions not aligned with GOLD "
                      f"(len GOLD={L}, len PLM={len(pl_p)})")

        # Overlaps
        p2_vs_gold_j = jaccard(p2_pos_set, gold_pos_set) if p2_pos_set else 0.0
        plm_vs_gold_j = jaccard(pl_pos_set, gold_pos_set) if pl_pos_set else 0.0
        p2_vs_plm_j = jaccard(p2_pos_set, pl_pos_set) if (p2_pos_set and pl_pos_set) else 0.0

        print(f"{pdb_id}.{ch}: L={L}, "
              f"GOLD+={len(gold_pos_set)}, "
              f"P2+={len(p2_pos_set)}, "
              f"PLM+={len(pl_pos_set)}, "
              f"J(P2,GOLD)={p2_vs_gold_j:.3f}, "
              f"J(PLM,GOLD)={plm_vs_gold_j:.3f}, "
              f"J(P2,PLM)={p2_vs_plm_j:.3f}")

        rows_for_csv.append({
            "pdb_id": pdb_id,
            "chain_id": ch,
            "L": L,
            "n_gold_pos": len(gold_pos_set),
            "n_p2_pos": len(p2_pos_set),
            "n_plm_pos": len(pl_pos_set),
            "p2_vs_gold_jaccard": f"{p2_vs_gold_j:.6f}",
            "plm_vs_gold_jaccard": f"{plm_vs_gold_j:.6f}",
            "p2_vs_plm_jaccard": f"{p2_vs_plm_j:.6f}",
            "p2_len_ok": int(p2_len_ok),
            "plm_len_ok": int(pl_len_ok),
        })

    # Optional CSV
    if args.out_csv and rows_for_csv:
        os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
        fieldnames = list(rows_for_csv[0].keys())
        with open(args.out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows_for_csv:
                w.writerow(r)
        print()
        print(f"[ok] Per-chain summary written to: {args.out_csv}")


if __name__ == "__main__":
    main()
