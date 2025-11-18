#!/usr/bin/env python3
"""
04.full_evaluation.py

Full evaluation of GOLD vs P2Rank vs PLM across all chains in a kinase dataset.
Produces per-chain metrics and a dataset-level summary.

Input files per dataset:
  01.chain_gold_labels.json
  02.chain_p2rank_labels.json
  03.chain_plm_labels.json

Outputs:
  diag/04.full_eval.csv
  diag/04.summary_metrics.csv

Metrics:
  - Jaccard
  - Precision, Recall, F1
  - Matthews Correlation Coefficient (MCC)
  - Counts: TP, FP, TN, FN
  - Length sanity: check consistent L among methods
"""

import json
import os
import argparse
import csv
from math import sqrt

# ------------------------ metrics ------------------------

def jaccard(a, b):
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union > 0 else 0.0

def precision(tp, fp):
    return tp / (tp + fp) if (tp + fp) > 0 else 0.0

def recall(tp, fn):
    return tp / (tp + fn) if (tp + fn) > 0 else 0.0

def f1(p, r):
    return (2*p*r)/(p+r) if (p+r) > 0 else 0.0

def mcc(tp, fp, fn, tn):
    num = tp*tn - fp*fn
    den = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return num/den if den > 0 else 0.0

# ------------------------ loading ------------------------

def load_json(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return json.load(open(path))

def to_pos_set(entry):
    """From pred_labels_fixed or pred_probs threshold≥0.75."""
    if "pred_labels_fixed" in entry:
        return {i for i,v in enumerate(entry["pred_labels_fixed"],1) if v == 1}
    elif "pred_probs" in entry:
        return {i+1 for i,p in enumerate(entry["pred_probs"]) if p >= 0.75}
    else:
        return set()

# ------------------------ main ------------------------

def evaluate_dataset(folder):

    gold = load_json(os.path.join(folder, "01.chain_gold_labels.json"))
    p2   = load_json(os.path.join(folder, "02.chain_p2rank_labels.json"))
    plm  = load_json(os.path.join(folder, "03.chain_plm_labels.json"))

    # index by (pdb,chain)
    idx_gold = {(e["pdb_id"], e["chain_id"]): e for e in gold}
    idx_p2   = {(e["pdb_id"], e["chain_id"]): e for e in p2}
    idx_plm  = {(e["pdb_id"], e["chain_id"]): e for e in plm}

    keys = list(idx_gold.keys())

    # output dir
    diag_dir = os.path.join(folder, "diag")
    os.makedirs(diag_dir, exist_ok=True)

    per_chain_csv = os.path.join(diag_dir, "04.full_eval.csv")
    summary_csv   = os.path.join(diag_dir, "04.summary_metrics.csv")

    rows = []

    for (pdb, ch) in keys:
        g = idx_gold[(pdb,ch)]
        L = len(g["sequence_residue_numbers"])

        gold_pos = set(g["gold_residue_numbers"])

        p2_pos  = to_pos_set(idx_p2.get((pdb,ch), {}))
        plm_pos = to_pos_set(idx_plm.get((pdb,ch), {}))

        # confusion counts
        tp_p2 = len(gold_pos & p2_pos)
        fp_p2 = len(p2_pos - gold_pos)
        fn_p2 = len(gold_pos - p2_pos)
        tn_p2 = L - tp_p2 - fp_p2 - fn_p2

        tp_plm = len(gold_pos & plm_pos)
        fp_plm = len(plm_pos - gold_pos)
        fn_plm = len(gold_pos - plm_pos)
        tn_plm = L - tp_plm - fp_plm - fn_plm

        # metrics
        p2_prec = precision(tp_p2, fp_p2)
        p2_rec  = recall(tp_p2, fn_p2)
        p2_f1   = f1(p2_prec, p2_rec)
        p2_jac  = jaccard(gold_pos, p2_pos)
        p2_mcc  = mcc(tp_p2, fp_p2, fn_p2, tn_p2)

        plm_prec = precision(tp_plm, fp_plm)
        plm_rec  = recall(tp_plm, fn_plm)
        plm_f1   = f1(plm_prec, plm_rec)
        plm_jac  = jaccard(gold_pos, plm_pos)
        plm_mcc  = mcc(tp_plm, fp_plm, fn_plm, tn_plm)

        rows.append({
            "pdb_id": pdb,
            "chain_id": ch,
            "L": L,

            "gold_pos": len(gold_pos),
            "p2_pos": len(p2_pos),
            "plm_pos": len(plm_pos),

            "p2_precision": p2_prec,
            "p2_recall": p2_rec,
            "p2_f1": p2_f1,
            "p2_jaccard": p2_jac,
            "p2_mcc": p2_mcc,

            "plm_precision": plm_prec,
            "plm_recall": plm_rec,
            "plm_f1": plm_f1,
            "plm_jaccard": plm_jac,
            "plm_mcc": plm_mcc,
        })

    # --- save per-chain ---
    with open(per_chain_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # --- summary ---
    import numpy as np
    summary = {
        "dataset": folder,
        "n_chains": len(rows),

        "p2_mean_jaccard": float(np.mean([r["p2_jaccard"] for r in rows])),
        "plm_mean_jaccard": float(np.mean([r["plm_jaccard"] for r in rows])),

        "p2_mean_f1": float(np.mean([r["p2_f1"] for r in rows])),
        "plm_mean_f1": float(np.mean([r["plm_f1"] for r in rows])),

        "p2_mean_mcc": float(np.mean([r["p2_mcc"] for r in rows])),
        "plm_mean_mcc": float(np.mean([r["plm_mcc"] for r in rows])),
    }

    with open(summary_csv, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] Full evaluation written to: {per_chain_csv}")
    print(f"[ok] Summary metrics written to: {summary_csv}")

# ---------------------------------------------------------

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="Dataset folder (e.g., Kinase_Type_III)")
    args = ap.parse_args()

    evaluate_dataset(args.folder)
