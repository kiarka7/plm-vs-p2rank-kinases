#!/usr/bin/env python3
"""
03.chain_plm_labels.py — robust alignment of PLM predictions to GOLD chains
GOLD-ordered output (record order follows first appearance in 01.chain_gold_labels.json).

What this script does
---------------------
Align raw PLM outputs to the canonical residue order taken from
01.chain_gold_labels.json, and write one JSON entry per (pdb, chain) with:
  - sequence_residue_numbers (author numbering; canonical order)
  - pred_probs (float per residue; same length as GOLD sequence)
  - (optional) pred_labels_fixed (0/1 per residue using a fixed threshold, for inspection only)

Noise control:
  - VERBOSE=False prints only problematic chains (absent/failed or with unmapped indices).
  - Full per-chain diagnostics are written into 03.plm_chain_report.csv
  - Problematic chains are also listed in 03.plm_problems.csv

IMPORTANT:
  - Now GOLD keys (PDB, chain) are normalized to UPPERCASE as well,
    so case mismatches should no longer cause 'PLM absent'.
"""

import os, json, csv
from typing import Any, Dict, List, Tuple, Optional

# -------- CONFIG --------
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_ALLO"
GOLD_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
PLM_RAW_JSON  = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_predictions.json")
OUTPUT_JSON   = os.path.join(BASE_DIR, KINASE_FOLDER, "03.chain_plm_labels.json")
REPORT_CSV    = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_chain_report.csv")
PROBLEMS_CSV  = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_problems.csv")

# Console verbosity
VERBOSE = False              # True = print every chain; False = print only problems
PRINT_ONLY_PROBLEMS = True   # If VERBOSE is False, restrict output to absent/failed/unmapped>0

# Diagnostics
TOP_K = 5  # how many top residues to show in CSV diagnostics

# Optional convenience labels inside JSON (does NOT affect evaluation)
EMIT_LABELS     = True      # set False to emit probabilities only
FIXED_THRESHOLD = 0.75      # used only if EMIT_LABELS is True
# ------------------------


# ---- helpers to read & normalize keys ----
def norm_pdb_id(x: Any) -> str:
    return str(x).strip()[:4].upper()

def norm_chain_id(x: Any) -> str:
    s = str(x).strip()
    return (s or "A").upper()

def clamp01(v: float) -> float:
    if v < 0.0: return 0.0
    if v > 1.0: return 1.0
    return float(v)

def pick_first(d: Dict, keys: List[str], default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default
# ------------------------------------------


def load_gold_positions_and_order(path: str):
    """
    Returns:
      gold_idx: {(PDB, CHAIN): {"positions": [author_numbers]}}
      gold_order: list[(PDB, CHAIN)] in the *first-occurrence* order from GOLD
    NOTE: Both PDB and chain are normalized to UPPERCASE to avoid case mismatches.
    """
    gold = json.load(open(path))
    gold_idx = {}
    order = []
    seen = set()
    for e in gold:
        key = (norm_pdb_id(e["pdb_id"]), norm_chain_id(e["chain_id"]))
        if key not in seen:
            seen.add(key)
            order.append(key)
        # keep the first occurrence per (pdb,chain); evaluation merges ligands anyway
        if key not in gold_idx:
            pos = [int(x) for x in e.get("sequence_residue_numbers", [])]
            gold_idx[key] = {"positions": pos}
    return gold_idx, order


# ---------- extractors for different PLM shapes ----------
def try_dense_probs(entry: Dict, L: int) -> Optional[List[float]]:
    arr = pick_first(entry, ["probabilities", "pred_probs", "per_residue_probs", "scores"])
    if isinstance(arr, list) and len(arr) == L:
        return [clamp01(float(x)) for x in arr]
    return None

def detect_index_base(idxs: List[int], L: int) -> Optional[int]:
    if not idxs:
        return None
    mn, mx = min(idxs), max(idxs)
    if mn == 0: return 0
    if mx == L and mn >= 1: return 1
    if mx < L and mn >= 0: return 0
    return None

def try_sparse_index(entry: Dict, L: int) -> Tuple[Optional[List[float]], int, int]:
    preds = pick_first(entry, ["residues_predicted", "residues", "predictions", "sites", "positions"])
    if not isinstance(preds, list):
        return None, 0, 0

    idxs, probs = [], []
    n_bad = 0
    for p in preds:
        if not isinstance(p, dict): continue
        idx = pick_first(p, ["position", "index", "residue_index", "i", "idx"])
        prob = pick_first(p, ["probability", "prob", "score", "p"])
        if idx is None or prob is None: continue
        try:
            idx = int(idx)
            pr  = clamp01(float(prob))
        except Exception:
            n_bad += 1
            continue
        idxs.append(idx)
        probs.append(pr)

    if not idxs:
        return None, 0, n_bad

    base = detect_index_base(idxs, L)
    if base is None:
        base = 0  # safer default

    dense = [0.0] * L
    n_oor = 0
    for idx, pr in zip(idxs, probs):
        i = idx - 1 if base == 1 else idx
        if 0 <= i < L:
            dense[i] = max(dense[i], pr)
        else:
            n_oor += 1

    return dense, base, (n_oor + n_bad)

def try_by_number(entry: Dict, gold_positions: List[int]) -> Tuple[Optional[List[float]], int]:
    preds = pick_first(entry, ["residues", "predictions", "sites"])
    if not isinstance(preds, list):
        return None, 0

    pos_to_idx = {int(n): i for i, n in enumerate(gold_positions)}
    L = len(gold_positions)
    dense = [0.0] * L
    n_unmapped = 0

    for p in preds:
        if not isinstance(p, dict): continue
        num = pick_first(p, ["residue_number", "resnum", "author_seq_id", "author_pos", "auth_seq_id", "seq_id"])
        prob = pick_first(p, ["probability", "prob", "score", "p"])
        if num is None or prob is None: continue
        try:
            num = int(num); pr = clamp01(float(prob))
        except Exception:
            n_unmapped += 1
            continue
        if num in pos_to_idx:
            i = pos_to_idx[num]
            dense[i] = max(dense[i], pr)
        else:
            n_unmapped += 1

    return dense, n_unmapped
# --------------------------------------------------------


def main():
    # Load GOLD (now normalized to UPPERCASE)
    gold_idx, gold_order = load_gold_positions_and_order(GOLD_JSON)
    print(f"Loaded {len(gold_idx)} GOLD chains; GOLD-ordered output enabled (UPPERCASE PDB/chain)")

    # Load PLM raw JSON (list of per-chain dicts) and index by normalized key
    plm_raw = json.load(open(PLM_RAW_JSON))
    plm_by_key: Dict[Tuple[str,str], Dict] = {}
    for e in plm_raw:
        pdb = pick_first(e, ["pdb_id", "pdb", "pdbid"])
        chain = pick_first(e, ["chain_id", "chain", "chainid"])
        if pdb is None or chain is None:  # skip malformed entries
            continue
        key = (norm_pdb_id(pdb), norm_chain_id(chain))
        if key not in plm_by_key:
            plm_by_key[key] = e

    json_out = []
    report_rows = []
    problems_rows = []   # absent/failed/unmapped>0

    used = 0
    absent = 0

    for key in gold_order:
        pdb_id, chain_id = key
        gold_pos = gold_idx[key]["positions"]
        L = len(gold_pos)

        entry = plm_by_key.get(key)
        if entry is None:
            # still emit zeros to keep alignment with GOLD/evaluation
            probs = [0.0] * L
            labels = [0] * L if EMIT_LABELS else None

            rec = {
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "sequence_residue_numbers": gold_pos,
                "pred_probs": probs
            }
            if EMIT_LABELS:
                rec["pred_labels_fixed"] = labels
            json_out.append(rec)

            row = {
                "pdb_id": pdb_id, "chain_id": chain_id,
                "mode": "absent",
                "index_base": "",
                "n_out_of_range_or_unmapped": 0,
                "L": L,
                "n_nonzero": 0,
                "mean_prob": f"{0.0:.6f}",
                "max_prob": f"{0.0:.6f}",
                f"top{TOP_K}_positions": "",
                "fixed_thr": (f"{FIXED_THRESHOLD:.2f}" if EMIT_LABELS else ""),
                "fixed_pred_pos": (0 if EMIT_LABELS else "")
            }
            report_rows.append(row)
            problems_rows.append({**row, "problem": "absent"})

            if VERBOSE or PRINT_ONLY_PROBLEMS:
                print(f"{pdb_id}.{chain_id}: PLM absent → filled zeros (L={L})")
            absent += 1
            continue

        # Try different extraction modes
        probs = try_dense_probs(entry, L)
        mode = "dense"
        index_base = ""
        n_bad = 0

        if probs is None:
            probs, index_base, n_bad = try_sparse_index(entry, L)
            mode = "sparse_index"

        if probs is None:
            probs, n_bad = try_by_number(entry, gold_pos)
            mode = "by_number"
            index_base = ""

        if probs is None:
            probs = [0.0] * L
            mode = "failed"
            n_bad = 0
            if VERBOSE or PRINT_ONLY_PROBLEMS:
                print(f"[!] {pdb_id}.{chain_id}: could not parse PLM entry; filled zeros")

        # Optional fixed-threshold labels (inspection only)
        labels = [1 if float(p) >= FIXED_THRESHOLD else 0 for p in probs] if EMIT_LABELS else None

        # Diagnostics
        nz = sum(1 for x in probs if x > 0)
        mean_p = (sum(probs)/L) if L else 0.0
        max_p  = max(probs) if probs else 0.0
        top = sorted(zip(gold_pos, probs), key=lambda t: t[1], reverse=True)[:TOP_K]
        top_str = "; ".join(f"{pos}:{pr:.3f}" for pos, pr in top)
        base_str = str(index_base) if index_base != "" else ""
        n_bad_str = int(n_bad)
        fixed_pos = sum(labels) if EMIT_LABELS else None

        # Print only problems unless VERBOSE=True
        if VERBOSE or (PRINT_ONLY_PROBLEMS and (mode in ("absent", "failed") or n_bad_str > 0)):
            print(f"{pdb_id}.{chain_id}: mode={mode}{('/'+base_str) if base_str else ''}, "
                  f"L={L}, nonzero={nz}, meanP={mean_p:.3f}, maxP={max_p:.3f}, bad/unmapped={n_bad_str}")

        rec = {
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "sequence_residue_numbers": gold_pos,
            "pred_probs": probs
        }
        if EMIT_LABELS:
            rec["pred_labels_fixed"] = labels
        json_out.append(rec)

        row = {
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "mode": mode,
            "index_base": base_str,
            "n_out_of_range_or_unmapped": n_bad_str,
            "L": L,
            "n_nonzero": nz,
            "mean_prob": f"{mean_p:.6f}",
            "max_prob": f"{max_p:.6f}",
            f"top{TOP_K}_positions": top_str,
            "fixed_thr": (f"{FIXED_THRESHOLD:.2f}" if EMIT_LABELS else ""),
            "fixed_pred_pos": (fixed_pos if EMIT_LABELS else "")
        }
        report_rows.append(row)
        if mode in ("absent", "failed") or n_bad_str > 0:
            problems_rows.append({**row, "problem": "failed" if mode=="failed" else ("unmapped>0" if n_bad_str>0 else "ok")})
        used += 1

    # Write outputs (GOLD-ordered)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(json_out, f, indent=2)
    # Full CSV report
    fieldnames = [
        "pdb_id","chain_id","mode","index_base",
        "n_out_of_range_or_unmapped","L","n_nonzero","mean_prob","max_prob",
        f"top{TOP_K}_positions",
        "fixed_thr","fixed_pred_pos"
    ]
    with open(REPORT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
        for r in report_rows: w.writerow(r)
    # Problems-only CSV
    if problems_rows:
        prob_fields = ["pdb_id","chain_id","problem","mode","index_base","n_out_of_range_or_unmapped","L","n_nonzero"]
        with open(PROBLEMS_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=prob_fields); w.writeheader()
            for r in problems_rows:
                w.writerow({k: r.get(k, "") for k in prob_fields})

    # Console summary
    print(f"Written {OUTPUT_JSON} with {len(json_out)} chains (used={used}, absent={absent}) — GOLD-ordered/UPPERCASE")
    print(f"Written {REPORT_CSV} with {len(report_rows)} rows")
    if problems_rows:
        print(f"Written {PROBLEMS_CSV} with {len(problems_rows)} problem rows")
        # show first few missing keys to help debugging
        missing_preview = [f"{r['pdb_id']}.{r['chain_id']}" for r in problems_rows if r.get('problem')=='absent'][:10]
        if missing_preview:
            print("Examples of PLM-absent chains:", ", ".join(missing_preview))
    print("Note: pred_labels_fixed are for inspection only. Keep using pred_probs for evaluation.")

if __name__ == "__main__":
    main()
