#!/usr/bin/env python3
# 04.0_evaluate_pipeline_all_in_one.py

"""
End-to-end evaluation of GOLD vs P2Rank vs PLM on kinase datasets.

Overview
--------
This script loads per–chain annotations and predictions, aligns them on a
canonical author residue numbering, and computes both MICRO and MACRO
classification metrics for several variants of the pipeline:

  • P2Rank ANY-pocket      – per-residue scores aggregated across all pockets
  • P2Rank RANK1-pocket    – residues from the top-ranked pocket only
  • P2Rank ORACLE-pocket   – best single pocket per chain, chosen using GOLD
  • PLM fixed threshold    – PLM probabilities thresholded at a fixed value
  • PLM best-MCC           – PLM probabilities thresholded at the global
                             MCC-optimal threshold (sweep over [0, 1])

Metrics
-------
For each variant the script reports:

  • MICRO metrics: computed on all residues concatenated across chains
      - Accuracy, F1 (weighted), Matthews correlation coefficient (MCC)
      - ROC AUC (AUROC), Average Precision (AUPR)

  • MACRO metrics: per-chain metric → mean / median / IQR across chains
      - reported separately for P2Rank ANY, PLM fixed, PLM best-MCC
      - optionally for RANK1 and ORACLE if data are available

Inputs
------
All paths are configured at the top of the script via BASE_DIR and KINASE_FOLDER.

Required JSON files (chain-level):

  • 01.chain_gold_labels.json
      - GOLD binding-site labels per (pdb_id, chain_id)
      - fields: pdb_id, chain_id, sequence_residue_numbers, labels

  • 02.chain_p2rank_labels.json
      - P2Rank per-residue scores/labels mapped to author residue numbers
      - fields: pdb_id, chain_id, sequence_residue_numbers,
                pred_probs (or probabilities), labels (optional)

  • 03.chain_plm_labels.json
      - PLM per-residue probabilities mapped to author residue numbers
      - fields: pdb_id, chain_id, sequence_residue_numbers,
                pred_probs (or probabilities)

Optional inputs:

  • 02.p2rank_rank1_map.json
      - per-chain map of residues belonging to the top-ranked P2Rank pocket

  • p2rank_predictions/
      - folder containing P2Rank *_predictions.csv files
      - used to build the ORACLE-best-single-pocket baseline
      - the script searches this directory recursively for files whose name
        starts with <pdb_id>.lower() and ends with "predictions.csv"

Outputs
-------
All outputs are written into <BASE_DIR>/<KINASE_FOLDER>/:

  • 04.eval_all_in_one.json
      - main summary with:
          config, data_summary, thresholds, micro, macro

  • 04.p2rank_oracle_chain_report.csv
      - per-chain ORACLE diagnostics (compact numeric summary)

  • 04.p2rank_oracle_chain_details.json
      - per-chain ORACLE diagnostics (positions and detailed metrics)

  • 04.plm_chain_diagnostics.csv
      - per-chain PLM statistics and threshold sweeps (fixed, best-MCC, best-F1)

  • 04.plm_chain_diagnostics.json
      - same information as the PLM CSV, in a JSON-friendly structure

Key configuration knobs
-----------------------
  • THRESHOLD_PLM_FIXED
      Fixed threshold used for the "PLM fixed" baseline (default: 0.75).

  • IMPUTE_MISSING_P2
      If True, missing P2Rank residues are imputed as probability=0, label=0.
      This keeps PLM and P2Rank aligned on the same residue set.
      If False, chains/positions without overlap are skipped.

  • SWEEP_START / SWEEP_STOP / SWEEP_STEP
      Range and step for threshold sweeps when searching for best MCC or F1.

Usage
-----
Adjust BASE_DIR and KINASE_FOLDER at the top of the script (or run from the
project root with the default settings), then:

    python 04.0_evaluate_pipeline_all_in_one.py

The script prints:

  • primary MICRO results for P2Rank ANY and PLM (fixed + best-MCC),
  • exploratory sweeps for P2Rank ANY / RANK1 / ORACLE,
  • PRIMARY MACRO (mean across chains) for selected variants,
  • basic information about P2Rank imputation (if enabled).

"""

import os, re, json, csv, math
import numpy as np
from collections import defaultdict
from sklearn.metrics import (
    accuracy_score, f1_score, matthews_corrcoef,
    roc_auc_score, average_precision_score,
    precision_score, recall_score
)

# ---------- CONFIG ----------
BASE_DIR = "."
KINASE_FOLDER = "Kinase_Type_ALLO"

GOLD_JSON = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels_CLEAN.json")
P2_JSON   = os.path.join(BASE_DIR, KINASE_FOLDER, "02.chain_p2rank_labels.json")
PLM_JSON  = os.path.join(BASE_DIR, KINASE_FOLDER, "03.chain_plm_labels.json")

# Optional: needs *_predictions.csv from P2Rank
RANK1_JSON= os.path.join(BASE_DIR, KINASE_FOLDER, "02.p2rank_rank1_map.json")
P2RANK_DIR= os.path.join(BASE_DIR, KINASE_FOLDER, "p2rank_predictions")

OUT_JSON  = os.path.join(BASE_DIR, KINASE_FOLDER, "04.eval_all_in_one_CLEAN.json")

ORACLE_CHAIN_CSV   = os.path.join(BASE_DIR, KINASE_FOLDER, "04.p2rank_oracle_chain_report_CLEAN.csv")
ORACLE_CHAIN_JSON  = os.path.join(BASE_DIR, KINASE_FOLDER, "04.p2rank_oracle_chain_details_CLEAN.json")

PLM_CHAIN_CSV      = os.path.join(BASE_DIR, KINASE_FOLDER, "04.plm_chain_diagnostics_CLEAN.csv")
PLM_CHAIN_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "04.plm_chain_diagnostics_CLEAN.json")

THRESHOLD_PLM_FIXED = 0.75
IMPUTE_MISSING_P2   = True

SWEEP_START = 0.00
SWEEP_STOP  = 1.00
SWEEP_STEP  = 0.01
# ----------------------------

_leading_int = re.compile(r"^\s*([+-]?\d+)")

def parse_pos(x):
    """Parse integer residue index from strings like '56', '56A', ' -12B', etc."""
    if isinstance(x, (int, np.integer)):
        return int(x)
    m = _leading_int.match(str(x))
    if not m:
        raise ValueError(f"Cannot parse residue position from: {x!r}")
    return int(m.group(1))

def to_int_list(xs):
    return [parse_pos(x) for x in xs]


# ---------- GOLD / P2 / PLM loaders ----------

def load_and_merge_gold(path):
    """
    Load GOLD chain labels and OR-merge multiple ligand entries per (PDB,CHAIN).

    Returns dict:
      (PDB,CHAIN) -> {
          "pos_order":    sorted list of author residue numbers,
          "labels_by_pos": {pos -> 0/1}
      }
    """
    raw = json.load(open(path))
    bucket = defaultdict(list)
    for e in raw:
        bucket[(e["pdb_id"].upper(), e["chain_id"].upper())].append(e)

    merged = {}
    for key, entries in bucket.items():
        labels_by_pos = defaultdict(int)
        for e in entries:
            pos = to_int_list(e["sequence_residue_numbers"])
            lbl = [int(x) for x in e["labels"]]
            if len(pos) != len(lbl):
                raise ValueError(
                    f"GOLD length mismatch inside {key}: "
                    f"pos={len(pos)} labels={len(lbl)}"
                )
            for p, y in zip(pos, lbl):
                labels_by_pos[p] = 1 if (labels_by_pos[p] or y) else 0
        merged[key] = {
            "pos_order": sorted(labels_by_pos.keys()),
            "labels_by_pos": dict(labels_by_pos),
        }
    return merged

def load_index_by_key(path, prob_keys=("pred_probs","probabilities"), labels_key="labels"):
    """
    Load a list-of-dicts JSON (P2 or PLM) and index by (PDB,CHAIN).
    Adds helper keys:
      _probs_key   – name of the probability array field
      _labels_key  – name of the label array field (if present)
    """
    data = json.load(open(path))
    idx = {}
    for e in data:
        key = (e["pdb_id"].upper(), e["chain_id"].upper())
        rec = dict(e)
        for k in prob_keys:
            if k in rec:
                rec["_probs_key"] = k
                break
        if labels_key in rec:
            rec["_labels_key"] = labels_key
        idx[key] = rec
    return idx

def load_rank1_map(path):
    """Load precomputed rank-1 pocket map from 02.p2rank_rank1_map.json (if present)."""
    if not os.path.isfile(path):
        return {}
    data = json.load(open(path))
    out = {}
    for e in data:
        key = (e["pdb_id"].upper(), e["chain_id"].upper())
        out[key] = set(int(x) for x in e.get("rank1_positions", []))
    return out


def collapse_p2_by_numeric_position(p2_pos_raw, p2_probs_list, p2_labels_list):
    """
    Collapse P2Rank per-residue predictions by numeric residue index
    (handles insertion codes, takes max prob and OR of labels).
    """
    prob_by_pos, label_by_pos = {}, {}
    for pos_raw, prob, lbl in zip(p2_pos_raw, p2_probs_list, p2_labels_list):
        try:
            pos = parse_pos(pos_raw)
        except ValueError:
            continue
        prob = float(prob)
        lbl = int(lbl)
        if pos in prob_by_pos:
            if prob > prob_by_pos[pos]:
                prob_by_pos[pos] = prob
            label_by_pos[pos] = 1 if (label_by_pos[pos] or lbl) else 0
        else:
            prob_by_pos[pos] = prob
            label_by_pos[pos] = lbl
    return prob_by_pos, label_by_pos


# ---------- Metrics helpers ----------

def score(y_t, y_pred, y_prob):
    """
    Compute standard classification metrics for binary labels.
    Returns: ACC, F1 (weighted), MCC, AUROC, AUPR
    """
    acc = accuracy_score(y_t, y_pred)
    f1  = f1_score(y_t, y_pred, average='weighted', zero_division=0)
    mcc = matthews_corrcoef(y_t, y_pred)

    if len(np.unique(y_t)) == 2:
        try:
            auc = roc_auc_score(y_t, y_prob)
        except Exception:
            auc = float("nan")
        try:
            ap  = average_precision_score(y_t, y_prob)
        except Exception:
            ap = float("nan")
    else:
        auc = ap = float("nan")

    return float(acc), float(f1), float(mcc), float(auc), float(ap)

def best_threshold(y_true, probs, metric="mcc", start=0.0, stop=1.0, step=0.01):
    """
    Sweep thresholds and pick the best one by MCC or F1.
    Returns: (best_threshold, best_value)
    """
    thresholds = np.arange(start, stop + 1e-9, step, dtype=float)
    best_t, best_val = None, -1.0
    for t in thresholds:
        y_pred = (probs >= t).astype(int)
        if metric == "mcc":
            val = matthews_corrcoef(y_true, y_pred)
        else:
            val = f1_score(y_true, y_pred, zero_division=0)
        if not np.isnan(val) and val > best_val:
            best_t, best_val = float(t), float(val)
    return best_t, best_val

def macro_stats(values):
    """Return (mean, median, IQR) for a list of floats, ignoring NaNs."""
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    q25, q50, q75 = np.percentile(arr, [25, 50, 75])
    return (float(np.mean(arr)), float(q50), float(q75 - q25))

def fmt_line(label, acc, f1, mcc, auc, aupr):
    return (
        f"  {label:<22} ACC={acc:.3f}  F1={f1:.3f}  MCC={mcc:.3f}  "
        f"AUC={auc:.3f}  AUPR={aupr:.3f}"
    )

def _mean(xs):
    xs = np.array(xs, dtype=float)
    if xs.size == 0:
        return float("nan")
    return float(np.nanmean(xs))

def _macro_block(acc, f1, mcc, auc, aupr):
    m_acc = macro_stats(acc)
    m_f1  = macro_stats(f1)
    m_mcc = macro_stats(mcc)
    m_auc = macro_stats(auc)
    m_aupr= macro_stats(aupr)
    return {
        "acc":   {"mean": m_acc[0], "median": m_acc[1], "iqr": m_acc[2]},
        "f1":    {"mean": m_f1[0],  "median": m_f1[1],  "iqr": m_f1[2]},
        "mcc":   {"mean": m_mcc[0], "median": m_mcc[1], "iqr": m_mcc[2]},
        "auroc": {"mean": m_auc[0], "median": m_auc[1], "iqr": m_auc[2]},
        "aupr":  {"mean": m_aupr[0],"median": m_aupr[1],"iqr": m_aupr[2]},
    }


# ---------- P2Rank predictions.csv parsing (for ORACLE / RANK1) ----------

def find_predictions_file(pdb_id, p2rank_dir):
    """
    Recursively search for a *_predictions.csv file for a given PDB in p2rank_dir.

    Matches typical P2Rank layouts, e.g.:
      p2rank_predictions/1s9j.cif/1s9j_predictions.csv
      p2rank_predictions/1s9j_predictions.csv
    """
    pdb_lower = pdb_id.lower()
    for root, _, files in os.walk(p2rank_dir):
        for fn in files:
            fnl = fn.lower()
            if fnl.startswith(pdb_lower) and fnl.endswith("predictions.csv"):
                return os.path.join(root, fn)
    return None

def read_predictions_rows(fn):
    """
    Try multiple delimiters and return DictReader rows for *_predictions.csv.

    Required columns (case-insensitive):
      name, rank, score, probability, residue_ids|residues
    """
    if not fn or not os.path.isfile(fn):
        return []
    for delim in ("\t", ",", ";"):
        with open(fn, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter=delim)
            if not reader.fieldnames:
                continue
            cols = [c.strip().lower() for c in reader.fieldnames]
            required = {"name", "rank", "score", "probability", "residue_ids"}
            # tolerate 'residues' instead of 'residue_ids'
            if "residue_ids" not in cols and "residues" in cols:
                required.remove("residue_ids")
                required.add("residues")
            if not required.issubset(set(cols)):
                continue
            rows = list(reader)
            if rows:
                return rows
    return []

def parse_predictions_csv_for_pockets(pdb_id, p2rank_dir):
    """
    Parse P2Rank *_predictions.csv into a structure:

      {chain_id -> [ {rank, name, score, probability, residues=set(int)}, ... ] }

    Only residues of format 'A_172' are used (chain + integer position).
    """
    fn = find_predictions_file(pdb_id, p2rank_dir)
    if not fn:
        return {}

    rows = read_predictions_rows(fn)
    pockets_by_chain = defaultdict(list)
    for row in rows:
        row_norm = {k.strip().lower(): v for k, v in row.items()}
        resid_col = "residue_ids" if "residue_ids" in row_norm else (
            "residues" if "residues" in row_norm else None
        )
        res_tokens = str(row_norm.get(resid_col, "")).split()
        by_chain = defaultdict(set)
        for tok in res_tokens:
            if "_" not in tok:
                continue
            ch, pos = tok.split("_", 1)
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            by_chain[ch.upper()].add(pos_i)

        # meta
        try:
            rk = int(row_norm.get("rank", "") or 0)
        except ValueError:
            rk = 0
        try:
            sc = float(row_norm.get("score", "") or 0.0)
        except ValueError:
            sc = 0.0
        try:
            pr = float(row_norm.get("probability", "") or 0.0)
        except ValueError:
            pr = 0.0
        name = str(row_norm.get("name", "")).strip()

        for ch, s in by_chain.items():
            if s:
                pockets_by_chain[ch].append({
                    "rank": rk,
                    "name": name,
                    "score": sc,
                    "probability": pr,
                    "residues": s,
                })

    # sort pockets per chain by rank
    for ch in pockets_by_chain:
        pockets_by_chain[ch].sort(key=lambda d: d.get("rank", 0))
    return pockets_by_chain


# ---------- Main evaluation ----------

def main():
    gold = load_and_merge_gold(GOLD_JSON)
    p2   = load_index_by_key(P2_JSON)
    plm  = load_index_by_key(PLM_JSON)
    rank1_map = load_rank1_map(RANK1_JSON)

    common = set(gold) & set(plm)
    print(f"Evaluating (micro) on {len(common)} chains (GOLD+PLM; P2 imputed).")

    # ---- MICRO collectors ----
    y_true_all = []

    y_p2_any_labels_all = []
    y_p2_any_probs_all  = []

    y_plm_labels_all    = []
    y_plm_probs_all     = []

    y_true_r1_all       = []
    y_p2_r1_labels_all  = []
    y_p2_r1_probs_all   = []

    y_true_or_all       = []
    y_p2_or_labels_all  = []
    y_p2_or_probs_all   = []

    # ---- MACRO collectors ----
    macro_p2_any_acc = []
    macro_p2_any_f1  = []
    macro_p2_any_mcc = []
    macro_p2_any_auc = []
    macro_p2_any_aupr= []

    macro_plm_fix_acc = []
    macro_plm_fix_f1  = []
    macro_plm_fix_mcc = []
    macro_plm_fix_auc = []
    macro_plm_fix_aupr= []

    macro_plm_bm_acc  = []
    macro_plm_bm_f1   = []
    macro_plm_bm_mcc  = []
    macro_plm_bm_auc  = []
    macro_plm_bm_aupr = []

    macro_r1_acc = []
    macro_r1_f1  = []
    macro_r1_mcc = []
    macro_r1_auc = []
    macro_r1_aupr= []

    macro_or_acc = []
    macro_or_f1  = []
    macro_or_mcc = []
    macro_or_auc = []
    macro_or_aupr= []

    # ---- Per-chain export collectors ----
    oracle_rows_csv  = []
    oracle_rows_json = []

    plm_rows_csv  = []
    plm_rows_json = []

    skipped_no_overlap = 0
    total_p2_missing_positions = 0
    total_positions = 0
    pockets_cache = {}

    # cache for PLM per-chain (used later for macro bestMCC + PLM diagnostics)
    per_chain_cache = []  # [(key, y_true, plm_probs)]

    for key in sorted(common):
        g    = gold[key]
        plme = plm[key]
        pdb, chain = key

        # Canonical order = PLM positions (fallback GOLD)
        plm_pos = to_int_list(plme.get("sequence_residue_numbers", g["pos_order"]))
        L = len(plm_pos)

        labels_by_pos = g["labels_by_pos"]
        y_true = np.array([int(labels_by_pos.get(int(p), 0)) for p in plm_pos], dtype=int)
        gold_pos_set = {int(p) for p, y in zip(plm_pos, y_true) if y == 1}

        pk = plme.get("_probs_key")
        if pk is None:
            # no probabilities → skip chain
            continue
        plm_probs = np.array([float(x) for x in plme[pk]], dtype=float)
        if plm_probs.size != L:
            # inconsistent length, skip chain
            continue
        y_plm_fixed = (plm_probs >= THRESHOLD_PLM_FIXED).astype(int)

        # ---- P2 ANY (imputed or subset) ----
        if key in p2:
            p2e = p2[key]
            p2_probs_key = p2e.get("_probs_key", "pred_probs")
            p2_labels_key = p2e.get("_labels_key")

            p2_pos_raw    = p2e.get("sequence_residue_numbers", [])
            p2_probs_list = p2e.get(p2_probs_key, [])

            if p2_labels_key is None:
                # we do not have explicit labels – treat as all-zero (no pockets)
                if IMPUTE_MISSING_P2:
                    p2_any_probs = np.zeros(L, dtype=float)
                    p2_any_lbls  = np.zeros(L, dtype=int)
                    total_p2_missing_positions += L
                    total_positions += L
                else:
                    skipped_no_overlap += 1
                    continue
            else:
                p2_labels_list = p2e.get(p2_labels_key, [])
                if not (len(p2_pos_raw) == len(p2_probs_list) == len(p2_labels_list)):
                    # malformed P2 entry for this chain – skip
                    continue
                p2_prob_by_pos, p2_label_by_pos = collapse_p2_by_numeric_position(
                    p2_pos_raw, p2_probs_list, p2_labels_list
                )
                if IMPUTE_MISSING_P2:
                    p2_any_probs = np.array(
                        [p2_prob_by_pos.get(int(p), 0.0) for p in plm_pos],
                        dtype=float
                    )
                    p2_any_lbls  = np.array(
                        [p2_label_by_pos.get(int(p), 0) for p in plm_pos],
                        dtype=int
                    )
                    total_p2_missing_positions += sum(
                        1 for p in plm_pos if int(p) not in p2_prob_by_pos
                    )
                    total_positions += L
                else:
                    mask = np.array(
                        [int(p) in p2_prob_by_pos for p in plm_pos], dtype=bool
                    )
                    if not np.any(mask):
                        skipped_no_overlap += 1
                        continue
                    plm_pos      = np.array(plm_pos)[mask].tolist()
                    y_true       = y_true[mask]
                    plm_probs    = plm_probs[mask]
                    y_plm_fixed  = y_plm_fixed[mask]
                    p2_any_probs = np.array(
                        [p2_prob_by_pos[int(p)] for p in plm_pos],
                        dtype=float
                    )
                    p2_any_lbls  = np.array(
                        [p2_label_by_pos[int(p)] for p in plm_pos],
                        dtype=int
                    )
        else:
            # no P2 entry at all
            if IMPUTE_MISSING_P2:
                p2_any_probs = np.zeros(L, dtype=float)
                p2_any_lbls  = np.zeros(L, dtype=int)
                total_p2_missing_positions += L
                total_positions += L
            else:
                skipped_no_overlap += 1
                continue

        # ---- P2 RANK1 ----
        r1_labels = r1_probs = None
        if key in rank1_map:
            r1_pos = rank1_map[key]
            r1_labels = np.array(
                [1 if int(p) in r1_pos else 0 for p in plm_pos],
                dtype=int
            )
            r1_probs = np.array(
                [p2_any_probs[i] if r1_labels[i] == 1 else 0.0 for i in range(L)],
                dtype=float
            )

        # ---- P2 ORACLE: per-chain best single pocket ----
        if pdb not in pockets_cache:
            pockets_cache[pdb] = parse_predictions_csv_for_pockets(pdb, P2RANK_DIR)
        chain_pockets = pockets_cache.get(pdb, {}).get(chain, [])

        or_labels = or_probs = None
        oracle_meta = {
            "rank": None,
            "name": None,
            "score": None,
            "probability": None,
            "size": None,
            "residues": [],
        }

        if chain_pockets:
            best_mcc = -1.0
            best_mask = None
            best_meta = None

            for pkt in chain_pockets:
                pkt_set = pkt["residues"]
                mask = np.array(
                    [1 if int(p) in pkt_set else 0 for p in plm_pos],
                    dtype=int
                )
                mcc = matthews_corrcoef(y_true, mask)
                if not np.isnan(mcc) and mcc > best_mcc:
                    best_mcc = mcc
                    best_mask = mask
                    best_meta = pkt

            if best_mask is not None and best_meta is not None:
                or_labels = best_mask
                or_probs  = np.array(
                    [p2_any_probs[i] if or_labels[i] == 1 else 0.0 for i in range(L)],
                    dtype=float
                )
                residues_sorted = sorted(
                    [int(p) for p, y in zip(plm_pos, or_labels) if y == 1]
                )
                oracle_meta = {
                    "rank": best_meta.get("rank"),
                    "name": best_meta.get("name"),
                    "score": best_meta.get("score"),
                    "probability": best_meta.get("probability"),
                    "size": len(best_meta.get("residues", [])),
                    "residues": residues_sorted,
                }

        # ---- MICRO aggregation ----
        y_true_all.extend(y_true.tolist())
        y_p2_any_labels_all.extend(p2_any_lbls.tolist())
        y_p2_any_probs_all.extend(p2_any_probs.tolist())
        y_plm_labels_all.extend(y_plm_fixed.tolist())
        y_plm_probs_all.extend(plm_probs.tolist())

        if r1_labels is not None:
            y_true_r1_all.extend(y_true.tolist())
            y_p2_r1_labels_all.extend(r1_labels.tolist())
            y_p2_r1_probs_all.extend(r1_probs.tolist())

        if or_labels is not None:
            y_true_or_all.extend(y_true.tolist())
            y_p2_or_labels_all.extend(or_labels.tolist())
            y_p2_or_probs_all.extend(or_probs.tolist())

        # ---- MACRO PRIMARY (ANY + PLM) ----
        a, f, m, auc, ap = score(y_true, p2_any_lbls, p2_any_probs)
        macro_p2_any_acc.append(a)
        macro_p2_any_f1.append(f)
        macro_p2_any_mcc.append(m)
        macro_p2_any_auc.append(auc)
        macro_p2_any_aupr.append(ap)

        a, f, m, auc, ap = score(y_true, y_plm_fixed, plm_probs)
        macro_plm_fix_acc.append(a)
        macro_plm_fix_f1.append(f)
        macro_plm_fix_mcc.append(m)
        macro_plm_fix_auc.append(auc)
        macro_plm_fix_aupr.append(ap)

        per_chain_cache.append((key, y_true.copy(), plm_probs.copy()))

        # RANK1 macro
        if r1_labels is not None:
            a, f, m, auc, ap = score(y_true, r1_labels, r1_probs)
            macro_r1_acc.append(a)
            macro_r1_f1.append(f)
            macro_r1_mcc.append(m)
            macro_r1_auc.append(auc)
            macro_r1_aupr.append(ap)

        # ORACLE macro + per-chain exports
        if or_labels is not None:
            a, f, m, auc, ap = score(y_true, or_labels, or_probs)
            macro_or_acc.append(a)
            macro_or_f1.append(f)
            macro_or_mcc.append(m)
            macro_or_auc.append(auc)
            macro_or_aupr.append(ap)

            # Per-chain ORACLE CSV/JSON rows (now actually filled)
            any_mcc = matthews_corrcoef(y_true, p2_any_lbls)
            any_f1  = f1_score(y_true, p2_any_lbls, zero_division=0)

            if r1_labels is not None:
                r1_mcc = matthews_corrcoef(y_true, r1_labels)
                r1_f1  = f1_score(y_true, r1_labels, zero_division=0)
                r1_pos_count = int(r1_labels.sum())
                r1_available = 1
            else:
                r1_mcc = float("nan")
                r1_f1  = float("nan")
                r1_pos_count = 0
                r1_available = 0

            tp = int(((y_true == 1) & (or_labels == 1)).sum())
            fp = int(((y_true == 0) & (or_labels == 1)).sum())
            fn = int(((y_true == 1) & (or_labels == 0)).sum())
            tn = int(L - tp - fp - fn)

            oracle_precision = precision_score(y_true, or_labels, zero_division=0)
            oracle_recall    = recall_score(y_true, or_labels, zero_division=0)

            # Jaccard vs GOLD positives using oracle_meta["residues"]
            oracle_res_set = set(oracle_meta["residues"])
            if gold_pos_set or oracle_res_set:
                inter = len(gold_pos_set & oracle_res_set)
                union = len(gold_pos_set | oracle_res_set)
                oracle_jaccard = inter / union if union > 0 else 0.0
            else:
                oracle_jaccard = 0.0

            oracle_rows_csv.append({
                "pdb_id": pdb,
                "chain_id": chain,
                "L": int(L),
                "gold_pos": int(y_true.sum()),
                "any_pos": int(p2_any_lbls.sum()),
                "any_mcc": any_mcc,
                "any_f1": any_f1,
                "rank1_available": r1_available,
                "rank1_pos": r1_pos_count,
                "rank1_mcc": r1_mcc,
                "rank1_f1": r1_f1,
                "oracle_rank": oracle_meta["rank"],
                "oracle_name": oracle_meta["name"],
                "oracle_score": oracle_meta["score"],
                "oracle_probability": oracle_meta["probability"],
                "oracle_size": oracle_meta["size"],
                "oracle_pos": int(or_labels.sum()),
                "oracle_mcc": m,
                "oracle_f1": f,
                "oracle_precision": oracle_precision,
                "oracle_recall": oracle_recall,
                "oracle_jaccard": oracle_jaccard,
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
            })

            oracle_rows_json.append({
                "pdb_id": pdb,
                "chain_id": chain,
                "L": int(L),
                "gold_positive_positions": sorted(int(p) for p in gold_pos_set),
                "p2_any_positive_positions": sorted(
                    int(p) for p, y in zip(plm_pos, p2_any_lbls) if y == 1
                ),
                "rank1_positive_positions": (
                    sorted(int(p) for p, y in zip(plm_pos, r1_labels) if y == 1)
                    if r1_labels is not None else []
                ),
                "oracle": {
                    "meta": oracle_meta,
                    "predicted_positions": oracle_meta["residues"],
                    "metrics": {
                        "mcc": m,
                        "f1": f,
                        "precision": oracle_precision,
                        "recall": oracle_recall,
                        "jaccard_vs_gold": oracle_jaccard,
                        "tp": tp,
                        "fp": fp,
                        "tn": tn,
                        "fn": fn,
                    },
                },
            })

    # ---- Convert MICRO collectors to arrays ----
    y_true_all          = np.array(y_true_all, dtype=int)
    y_p2_any_labels_all = np.array(y_p2_any_labels_all, dtype=int)
    y_p2_any_probs_all  = np.array(y_p2_any_probs_all, dtype=float)
    y_plm_labels_all    = np.array(y_plm_labels_all, dtype=int)
    y_plm_probs_all     = np.array(y_plm_probs_all, dtype=float)

    y_true_r1_all       = np.array(y_true_r1_all, dtype=int) if y_true_r1_all else None
    y_p2_r1_labels_all  = np.array(y_p2_r1_labels_all, dtype=int) if y_p2_r1_labels_all else None
    y_p2_r1_probs_all   = np.array(y_p2_r1_probs_all, dtype=float) if y_p2_r1_probs_all else None

    y_true_or_all       = np.array(y_true_or_all, dtype=int) if y_true_or_all else None
    y_p2_or_labels_all  = np.array(y_p2_or_labels_all, dtype=int) if y_p2_or_labels_all else None
    y_p2_or_probs_all   = np.array(y_p2_or_probs_all, dtype=float) if y_p2_or_probs_all else None

    # ---- PRIMARY MICRO (fixed thresholds) ----
    acc_p2_any, f1_p2_any, mcc_p2_any, auc_p2_any, ap_p2_any = score(
        y_true_all, y_p2_any_labels_all, y_p2_any_probs_all
    )
    acc_plm_fix, f1_plm_fix, mcc_plm_fix, auc_plm, ap_plm = score(
        y_true_all, y_plm_labels_all, y_plm_probs_all
    )

    # PLM bestMCC (global sweep)
    t_plm_best, _ = best_threshold(
        y_true_all, y_plm_probs_all, metric="mcc",
        start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
    )
    y_plm_best = (y_plm_probs_all >= t_plm_best).astype(int)
    acc_plm_bm, f1_plm_bm, mcc_plm_bm, _, _ = score(
        y_true_all, y_plm_best, y_plm_probs_all
    )

    # ---- PRIMARY MACRO: recompute PLM bestMCC per chain ----
    for (_key, y_true_c, plm_probs_c) in per_chain_cache:
        y_plm_b = (plm_probs_c >= t_plm_best).astype(int)
        a, f, m, auc, ap = score(y_true_c, y_plm_b, plm_probs_c)
        macro_plm_bm_acc.append(a)
        macro_plm_bm_f1.append(f)
        macro_plm_bm_mcc.append(m)
        macro_plm_bm_auc.append(auc)
        macro_plm_bm_aupr.append(ap)

    # ---- EXPLORATORY MICRO sweeps ----
    # P2 ANY sweep
    t_p2_any_mcc, _ = best_threshold(
        y_true_all, y_p2_any_probs_all, metric="mcc",
        start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
    )
    y_p2_any_best_mcc = (y_p2_any_probs_all >= t_p2_any_mcc).astype(int)
    acc_p2_any_bm, f1_p2_any_bm, mcc_p2_any_bm, _, _ = score(
        y_true_all, y_p2_any_best_mcc, y_p2_any_probs_all
    )

    t_p2_any_f1, _ = best_threshold(
        y_true_all, y_p2_any_probs_all, metric="f1",
        start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
    )
    y_p2_any_best_f1 = (y_p2_any_probs_all >= t_p2_any_f1).astype(int)
    acc_p2_any_bf, f1_p2_any_bf, mcc_p2_any_bf, _, _ = score(
        y_true_all, y_p2_any_best_f1, y_p2_any_probs_all
    )

    # RANK1 sweep (if available)
    if y_p2_r1_probs_all is not None and y_p2_r1_probs_all.size:
        t_r1_mcc, _ = best_threshold(
            y_true_r1_all, y_p2_r1_probs_all, metric="mcc",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_r1_best_mcc = (y_p2_r1_probs_all >= t_r1_mcc).astype(int)
        acc_r1_bm, f1_r1_bm, mcc_r1_bm, _, _ = score(
            y_true_r1_all, y_r1_best_mcc, y_p2_r1_probs_all
        )

        t_r1_f1, _ = best_threshold(
            y_true_r1_all, y_p2_r1_probs_all, metric="f1",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_r1_best_f1 = (y_p2_r1_probs_all >= t_r1_f1).astype(int)
        acc_r1_bf, f1_r1_bf, mcc_r1_bf, _, _ = score(
            y_true_r1_all, y_r1_best_f1, y_p2_r1_probs_all
        )
    else:
        t_r1_mcc = t_r1_f1 = None
        acc_r1_bm = f1_r1_bm = mcc_r1_bm = np.nan
        acc_r1_bf = f1_r1_bf = mcc_r1_bf = np.nan

    # ORACLE sweep (if available) + fixed ORACLE micro
    if y_p2_or_probs_all is not None and y_p2_or_probs_all.size:
        t_or_mcc, _ = best_threshold(
            y_true_or_all, y_p2_or_probs_all, metric="mcc",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_or_best_mcc = (y_p2_or_probs_all >= t_or_mcc).astype(int)
        acc_or_bm, f1_or_bm, mcc_or_bm, _, _ = score(
            y_true_or_all, y_or_best_mcc, y_p2_or_probs_all
        )

        t_or_f1, _ = best_threshold(
            y_true_or_all, y_p2_or_probs_all, metric="f1",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_or_best_f1 = (y_p2_or_probs_all >= t_or_f1).astype(int)
        acc_or_bf, f1_or_bf, mcc_or_bf, _, _ = score(
            y_true_or_all, y_or_best_f1, y_p2_or_probs_all
        )

        acc_or_fix, f1_or_fix, mcc_or_fix, auc_or_fix, ap_or_fix = score(
            y_true_or_all, y_p2_or_labels_all, y_p2_or_probs_all
        )
    else:
        t_or_mcc = t_or_f1 = None
        acc_or_bm = f1_or_bm = mcc_or_bm = np.nan
        acc_or_bf = f1_or_bf = mcc_or_bf = np.nan
        acc_or_fix = f1_or_fix = mcc_or_fix = auc_or_fix = ap_or_fix = np.nan

    # ----- WRITE MAIN JSON with MICRO + MACRO -----
    data_summary = {
        "n_chains": len(common),
        "n_positions_micro": int(y_true_all.size),
        "impute_missing_p2": IMPUTE_MISSING_P2,
        "p2_missing_positions_imputed": int(total_p2_missing_positions),
        "p2_missing_impute_pct": (
            100.0 * total_p2_missing_positions / total_positions
            if (IMPUTE_MISSING_P2 and total_positions) else 0.0
        ),
    }
    thresholds = {
        "plm_best_mcc": t_plm_best,
        "p2_any_best_mcc": t_p2_any_mcc,
        "p2_any_best_f1": t_p2_any_f1,
        "rank1_best_mcc": t_r1_mcc,
        "rank1_best_f1": t_r1_f1,
        "oracle_best_mcc": t_or_mcc,
        "oracle_best_f1": t_or_f1,
    }
    micro = {
        "p2_any_fixed": {
            "acc": acc_p2_any,
            "f1": f1_p2_any,
            "mcc": mcc_p2_any,
            "auroc": auc_p2_any,
            "aupr": ap_p2_any,
        },
        "plm_fixed": {
            "acc": acc_plm_fix,
            "f1": f1_plm_fix,
            "mcc": mcc_plm_fix,
            "auroc": auc_plm,
            "aupr": ap_plm,
        },
        "plm_best_mcc": {
            "threshold": t_plm_best,
            "acc": acc_plm_bm,
            "f1": f1_plm_bm,
            "mcc": mcc_plm_bm,
        },
        "p2_any_best_mcc": {
            "threshold": t_p2_any_mcc,
            "acc": acc_p2_any_bm,
            "f1": f1_p2_any_bm,
            "mcc": mcc_p2_any_bm,
        },
        "p2_any_best_f1": {
            "threshold": t_p2_any_f1,
            "acc": acc_p2_any_bf,
            "f1": f1_p2_any_bf,
            "mcc": mcc_p2_any_bf,
        },
    }
    if y_p2_r1_probs_all is not None and y_p2_r1_probs_all.size:
        acc_r1_fix, f1_r1_fix, mcc_r1_fix, auc_r1_fix, ap_r1_fix = score(
            y_true_r1_all, y_p2_r1_labels_all, y_p2_r1_probs_all
        )
        micro["rank1_fixed"] = {
            "acc": acc_r1_fix,
            "f1": f1_r1_fix,
            "mcc": mcc_r1_fix,
            "auroc": auc_r1_fix,
            "aupr": ap_r1_fix,
        }
        micro["rank1_best_mcc"] = {
            "threshold": t_r1_mcc,
            "acc": acc_r1_bm,
            "f1": f1_r1_bm,
            "mcc": mcc_r1_bm,
        }
        micro["rank1_best_f1"] = {
            "threshold": t_r1_f1,
            "acc": acc_r1_bf,
            "f1": f1_r1_bf,
            "mcc": mcc_r1_bf,
        }
    if y_p2_or_probs_all is not None and y_p2_or_probs_all.size:
        micro["oracle_fixed"] = {
            "acc": acc_or_fix,
            "f1": f1_or_fix,
            "mcc": mcc_or_fix,
            "auroc": auc_or_fix,
            "aupr": ap_or_fix,
        }
        micro["oracle_best_mcc"] = {
            "threshold": t_or_mcc,
            "acc": acc_or_bm,
            "f1": f1_or_bm,
            "mcc": mcc_or_bm,
        }
        micro["oracle_best_f1"] = {
            "threshold": t_or_f1,
            "acc": acc_or_bf,
            "f1": f1_or_bf,
            "mcc": mcc_or_bf,
        }

    macro = {
        "p2_any_fixed": _macro_block(
            macro_p2_any_acc, macro_p2_any_f1,
            macro_p2_any_mcc, macro_p2_any_auc, macro_p2_any_aupr
        ),
        "plm_fixed": _macro_block(
            macro_plm_fix_acc, macro_plm_fix_f1,
            macro_plm_fix_mcc, macro_plm_fix_auc, macro_plm_fix_aupr
        ),
        "plm_best_mcc": _macro_block(
            macro_plm_bm_acc, macro_plm_bm_f1,
            macro_plm_bm_mcc, macro_plm_bm_auc, macro_plm_bm_aupr
        ),
    }
    if macro_r1_acc:
        macro["rank1_fixed"] = _macro_block(
            macro_r1_acc, macro_r1_f1,
            macro_r1_mcc, macro_r1_auc, macro_r1_aupr
        )
    if macro_or_acc:
        macro["oracle_fixed"] = _macro_block(
            macro_or_acc, macro_or_f1,
            macro_or_mcc, macro_or_auc, macro_or_aupr
        )

    out = {
        "config": {
            "plm_threshold_fixed": THRESHOLD_PLM_FIXED,
            "threshold_sweep": {
                "start": SWEEP_START,
                "stop": SWEEP_STOP,
                "step": SWEEP_STEP,
            },
            "impute_missing_p2": IMPUTE_MISSING_P2,
        },
        "data_summary": data_summary,
        "thresholds": thresholds,
        "micro": micro,
        "macro": macro,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {OUT_JSON}")
    if IMPUTE_MISSING_P2 and total_positions:
        pct = 100.0 * total_p2_missing_positions / total_positions
        print(
            f"  P2 missing positions imputed: "
            f"{total_p2_missing_positions}/{total_positions} ({pct:.2f}%)"
        )

    # -----  CONSOLE (MICRO) -----
    print("\n[PRIMARY micro] (global concatenation)")
    print(fmt_line("P2 ANY fixed", acc_p2_any, f1_p2_any, mcc_p2_any, auc_p2_any, ap_p2_any))
    print(fmt_line(
        f"PLM fixed (thr={THRESHOLD_PLM_FIXED:.2f})",
        acc_plm_fix, f1_plm_fix, mcc_plm_fix, auc_plm, ap_plm
    ))
    print(fmt_line(
        f"PLM bestMCC (thr={t_plm_best:.2f})",
        acc_plm_bm, f1_plm_bm, mcc_plm_bm, auc_plm, ap_plm
    ))

    print("\n===== EXPLORATORY (sweeps & pocket variants) =====")
    print("-- P2Rank ANY-pocket (sweep, micro) --")
    print(fmt_line(
        f"best-MCC (thr={t_p2_any_mcc:.2f})",
        acc_p2_any_bm, f1_p2_any_bm, mcc_p2_any_bm, float('nan'), float('nan')
    ))
    print(fmt_line(
        f"best-F1  (thr={t_p2_any_f1:.2f})",
        acc_p2_any_bf, f1_p2_any_bf, mcc_p2_any_bf, float('nan'), float('nan')
    ))

    if y_p2_r1_probs_all is not None and y_p2_r1_probs_all.size:
        acc_r1_fix, f1_r1_fix, mcc_r1_fix, auc_r1_fix, ap_r1_fix = score(
            y_true_r1_all, y_p2_r1_labels_all, y_p2_r1_probs_all
        )
        print("-- P2Rank RANK1-pocket --")
        print(fmt_line("fixed", acc_r1_fix, f1_r1_fix, mcc_r1_fix, auc_r1_fix, ap_r1_fix))
        print(fmt_line(
            f"best-MCC (thr={t_r1_mcc:.2f})",
            acc_r1_bm, f1_r1_bm, mcc_r1_bm, float('nan'), float('nan')
        ))
        print(fmt_line(
            f"best-F1  (thr={t_r1_f1:.2f})",
            acc_r1_bf, f1_r1_bf, mcc_r1_bf, float('nan'), float('nan')
        ))
    else:
        print("-- P2Rank RANK1-pocket -- not available")

    if y_p2_or_probs_all is not None and y_p2_or_probs_all.size:
        print("-- P2Rank ORACLE-best-single-pocket --")
        print("  note: uses GOLD to select per-chain best pocket (optimistic upper bound)")
        print(fmt_line(
            "fixed", acc_or_fix, f1_or_fix, mcc_or_fix, auc_or_fix, ap_or_fix
        ))
        print(fmt_line(
            f"best-MCC (thr={t_or_mcc:.2f})",
            acc_or_bm, f1_or_bm, mcc_or_bm, float('nan'), float('nan')
        ))
        print(fmt_line(
            f"best-F1  (thr={t_or_f1:.2f})",
            acc_or_bf, f1_or_bf, mcc_or_bf, float('nan'), float('nan')
        ))
    else:
        print("-- P2Rank ORACLE-best-single-pocket -- not available")

    # ----- PRIMARY MACRO (print means across chains) -----
    mean_p2_acc  = _mean(macro_p2_any_acc)
    mean_p2_f1   = _mean(macro_p2_any_f1)
    mean_p2_mcc  = _mean(macro_p2_any_mcc)

    mean_plf_acc = _mean(macro_plm_fix_acc)
    mean_plf_f1  = _mean(macro_plm_fix_f1)
    mean_plf_mcc = _mean(macro_plm_fix_mcc)

    mean_plb_acc = _mean(macro_plm_bm_acc)
    mean_plb_f1  = _mean(macro_plm_bm_f1)
    mean_plb_mcc = _mean(macro_plm_bm_mcc)

    print("\n[PRIMARY macro] (mean across chains)")
    print(fmt_line(
        "P2 ANY fixed (macro)",
        mean_p2_acc, mean_p2_f1, mean_p2_mcc, float('nan'), float('nan')
    ))
    print(fmt_line(
        "PLM fixed (macro)",
        mean_plf_acc, mean_plf_f1, mean_plf_mcc, float('nan'), float('nan')
    ))
    print(fmt_line(
        "PLM bestMCC (macro)",
        mean_plb_acc, mean_plb_f1, mean_plb_mcc, float('nan'), float('nan')
    ))

    # ----- Per-chain ORACLE CSV/JSON (already collected during loop) -----
    oracle_fields = [
        "pdb_id","chain_id","L","gold_pos",
        "any_pos","any_mcc","any_f1",
        "rank1_available","rank1_pos","rank1_mcc","rank1_f1",
        "oracle_rank","oracle_name","oracle_score","oracle_probability",
        "oracle_size","oracle_pos","oracle_mcc","oracle_f1",
        "oracle_precision","oracle_recall","oracle_jaccard",
        "tp","fp","tn","fn"
    ]
    with open(ORACLE_CHAIN_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=oracle_fields)
        w.writeheader()
        for row in oracle_rows_csv:
            w.writerow(row)
    with open(ORACLE_CHAIN_JSON, "w") as f:
        json.dump(oracle_rows_json, f, indent=2)
    print(f"\nSaved per-chain ORACLE CSV:  {ORACLE_CHAIN_CSV}")
    print(f"Saved per-chain ORACLE JSON: {ORACLE_CHAIN_JSON}")

    # ----- Per-chain PLM diagnostics (CSV + JSON) -----
    for key, y_true_c, plm_probs_c in per_chain_cache:
        pdb, chain = key
        L = int(len(y_true_c))
        gold_pos = int(np.sum(y_true_c))

        mean_p = float(np.mean(plm_probs_c))
        std_p  = float(np.std(plm_probs_c))
        min_p  = float(np.min(plm_probs_c))
        max_p  = float(np.max(plm_probs_c))
        q25, med, q75 = np.percentile(plm_probs_c, [25, 50, 75])
        iqr_p = float(q75 - q25)

        y_fix = (plm_probs_c >= THRESHOLD_PLM_FIXED).astype(int)
        acc_f, f1_f, mcc_f, auc_f, ap_f = score(y_true_c, y_fix, plm_probs_c)
        npos_fix = int(np.sum(y_fix))

        t_mcc, _ = best_threshold(
            y_true_c, plm_probs_c, metric="mcc",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_bm = (plm_probs_c >= t_mcc).astype(int)
        acc_bm, f1_bm, mcc_bm, _, _ = score(y_true_c, y_bm, plm_probs_c)
        npos_bm = int(np.sum(y_bm))

        t_f1, _ = best_threshold(
            y_true_c, plm_probs_c, metric="f1",
            start=SWEEP_START, stop=SWEEP_STOP, step=SWEEP_STEP
        )
        y_bf = (plm_probs_c >= t_f1).astype(int)
        acc_bf, f1_bf, mcc_bf, _, _ = score(y_true_c, y_bf, plm_probs_c)
        npos_bf = int(np.sum(y_bf))

        plm_rows_csv.append({
            "pdb_id": pdb,
            "chain_id": chain,
            "L": L,
            "gold_pos": gold_pos,
            "mean": mean_p,
            "std": std_p,
            "min": min_p,
            "q25": float(q25),
            "median": float(med),
            "q75": float(q75),
            "max": max_p,
            "iqr": iqr_p,
            "fixed_thr": THRESHOLD_PLM_FIXED,
            "fixed_acc": acc_f,
            "fixed_f1": f1_f,
            "fixed_mcc": mcc_f,
            "fixed_auroc": auc_f,
            "fixed_aupr": ap_f,
            "fixed_pred_pos": npos_fix,
            "best_mcc_thr": t_mcc,
            "best_mcc_acc": acc_bm,
            "best_mcc_f1": f1_bm,
            "best_mcc_mcc": mcc_bm,
            "best_mcc_pred_pos": npos_bm,
            "best_f1_thr": t_f1,
            "best_f1_acc": acc_bf,
            "best_f1_f1": f1_bf,
            "best_f1_mcc": mcc_bf,
            "best_f1_pred_pos": npos_bf,
        })

        plm_rows_json.append({
            "pdb_id": pdb,
            "chain_id": chain,
            "L": L,
            "gold_positive_count": gold_pos,
            "probs_stats": {
                "mean": mean_p,
                "std": std_p,
                "min": min_p,
                "q25": float(q25),
                "median": float(med),
                "q75": float(q75),
                "max": max_p,
                "iqr": iqr_p,
            },
            "fixed": {
                "threshold": THRESHOLD_PLM_FIXED,
                "accuracy": acc_f,
                "f1": f1_f,
                "mcc": mcc_f,
                "auroc": auc_f,
                "aupr": ap_f,
                "predicted_positive": npos_fix,
            },
            "best_mcc": {
                "threshold": t_mcc,
                "accuracy": acc_bm,
                "f1": f1_bm,
                "mcc": mcc_bm,
                "predicted_positive": npos_bm,
            },
            "best_f1": {
                "threshold": t_f1,
                "accuracy": acc_bf,
                "f1": f1_bf,
                "mcc": mcc_bf,
                "predicted_positive": npos_bf,
            },
        })

    plm_fields = [
        "pdb_id","chain_id","L","gold_pos",
        "mean","std","min","q25","median","q75","max","iqr",
        "fixed_thr","fixed_acc","fixed_f1","fixed_mcc",
        "fixed_auroc","fixed_aupr","fixed_pred_pos",
        "best_mcc_thr","best_mcc_acc","best_mcc_f1","best_mcc_mcc","best_mcc_pred_pos",
        "best_f1_thr","best_f1_acc","best_f1_f1","best_f1_mcc","best_f1_pred_pos"
    ]
    with open(PLM_CHAIN_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=plm_fields)
        w.writeheader()
        for row in plm_rows_csv:
            w.writerow(row)
    with open(PLM_CHAIN_JSON, "w") as f:
        json.dump(plm_rows_json, f, indent=2)
    print(f"Saved per-chain PLM CSV:    {PLM_CHAIN_CSV}")
    print(f"Saved per-chain PLM JSON:   {PLM_CHAIN_JSON}")

if __name__ == "__main__":
    main()
