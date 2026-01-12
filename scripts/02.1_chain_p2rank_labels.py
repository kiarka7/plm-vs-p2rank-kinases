#!/usr/bin/env python3
"""
02.chain_p2rank_labels.py

Export P2Rank per-residue predictions with GOLD-ordered output and
an extra "best pocket" selection per (pdb, chain).

Key points (pipeline-friendly):
- Output records are in GOLD order (same (pdb,chain) sequence as GOLD).
- Residue arrays (sequence_residue_numbers, pred_probs, labels, residue_pocket_ids)
  are aligned to GOLD's canonical author numbering for that (pdb,chain).
  -> If P2Rank misses a position, we fill prob=0.0, label=0, pocket_id=-1.
- Chain IDs are normalized to UPPERCASE to match GOLD keys.

Writes:
- 02.chain_p2rank_labels.json  (aligned predictions)
- 02.p2rank_chain_report.csv   (small summary per chain)
"""

import os, re, json, glob, csv
import pandas as pd
from collections import defaultdict

# --- CONFIG ---
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"
GOLD_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
P2RANK_DIR    = os.path.join(BASE_DIR, KINASE_FOLDER, "p2rank_predictions")
OUTPUT_JSON   = os.path.join(BASE_DIR, KINASE_FOLDER, "02.chain_p2rank_labels.json")
REPORT_CSV    = os.path.join(BASE_DIR, KINASE_FOLDER, "02.p2rank_chain_report.csv")
# --------------

_leading_int = re.compile(r"^\s*([+-]?\d+)")

def parse_pos(x):
    """Extract leading integer from labels like '56', '56A', '  -12B' etc."""
    if isinstance(x, int): 
        return int(x)
    s = str(x); m = _leading_int.match(s)
    if not m:
        raise ValueError(f"Cannot parse residue position: {x!r}")
    return int(m.group(1))

def load_gold_labels_and_order(path):
    """
    Load GOLD and build:
      order: GOLD-ordered list of (PDB, CHAIN) keys (UPPERCASE).
      gold_pos: dict key -> list of positions (canonical author numbering, in order).
      gold_pos_set: key -> set(positions)
      gold_pos_positive: key -> set(positions where GOLD label == 1), OR-merged across duplicates.
    """
    gold = json.load(open(path))
    order = []
    seen = set()
    gold_pos = {}
    gold_pos_set = {}
    gold_pos_positive = {}

    for e in gold:
        key = (str(e["pdb_id"]).upper(), str(e["chain_id"]).upper())
        if key not in seen:
            seen.add(key); order.append(key)

        # canonical positions for the chain (take the first time we see this key)
        if key not in gold_pos:
            pos = [parse_pos(x) for x in e.get("sequence_residue_numbers", [])]
            gold_pos[key] = pos
            gold_pos_set[key] = set(pos)

        # OR-merge positives across repeated entries for the same chain
        lbl = [int(x) for x in e.get("labels", [])]
        if len(lbl) == len(gold_pos[key]):
            pos_set = set(p for p, y in zip(gold_pos[key], lbl) if y == 1)
            gold_pos_positive[key] = gold_pos_positive.get(key, set()) | pos_set

    return order, gold_pos, gold_pos_set, gold_pos_positive

def find_residue_csvs(p2dir):
    """Find *_residues.csv both flat and one-level nested."""
    return sorted(set(glob.glob(os.path.join(p2dir, "*_residues.csv")) +
                      glob.glob(os.path.join(p2dir, "*", "*_residues.csv"))))

def read_residue_csv(path):
    """
    Read P2Rank per-residue CSV. Expected columns (case-insensitive):
      chain, residue_label, pocket, probability
    Returns normalized DataFrame with correct dtypes.
    """
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    need = {"chain","residue_label","pocket","probability"}
    if not need.issubset(set(df.columns)):
        raise ValueError(f"Unexpected columns in {os.path.basename(path)}: {list(df.columns)}")
    out = pd.DataFrame({
        "chain": df["chain"].astype(str).str.strip(),  # case will be normalized later
        "residue_label": df["residue_label"].astype(str).str.strip(),
        "pocket": pd.to_numeric(df["pocket"], errors="coerce").fillna(0).astype(int),
        "probability": pd.to_numeric(df["probability"], errors="coerce").fillna(0.0).astype(float),
    })
    return out

def main():
    order, gold_pos, gold_pos_set, gold_pos_positive = load_gold_labels_and_order(GOLD_JSON)
    print(f"✓ Loaded GOLD for {len(order)} chains; GOLD-ordered output enabled")

    # accumulators per (PDB, CHAIN) — CHAIN UPPERCASE
    prob_by_pos   = defaultdict(dict)   # key -> {pos -> prob (max)}
    label_by_pos  = defaultdict(dict)   # key -> {pos -> 0/1 (OR over pockets)}
    pocket_by_pos = defaultdict(dict)   # key -> {pos -> pocket_id that gave the max prob}
    n_csv = 0; n_rows_total = 0; n_rows_after = 0

    for csv_path in find_residue_csvs(P2RANK_DIR):
        n_csv += 1
        df = read_residue_csv(csv_path)

        # infer PDB id from path
        folder = os.path.basename(os.path.dirname(csv_path))
        fname  = os.path.basename(csv_path)
        pdb_guess = folder if folder and folder != "p2rank_predictions" else fname.split("_")[0]
        pdb = pdb_guess[:4].upper()

        # group by chain (normalize to UPPERCASE to match GOLD)
        df["chain"] = df["chain"].astype(str).str.strip().str.upper()

        for chain_up, grp in df.groupby("chain"):
            key = (pdb, chain_up)
            labels_raw = grp["residue_label"].tolist()
            pockets    = grp["pocket"].tolist()
            probs      = grp["probability"].tolist()
            n_rows_total += len(labels_raw)

            # collapse to numeric pos: prob=max, label=OR, keep the pocket giving max prob
            tmp_prob = {}
            tmp_lab  = {}
            tmp_poc  = {}
            for lab, poc, pr in zip(labels_raw, pockets, probs):
                try:
                    pos = parse_pos(lab)
                except ValueError:
                    continue
                pr = float(pr); lb = 1 if int(poc) > 0 else 0
                if (pos not in tmp_prob) or (pr > tmp_prob[pos]):
                    tmp_prob[pos] = pr
                    tmp_poc[pos]  = int(poc)
                tmp_lab[pos] = 1 if (tmp_lab.get(pos,0) or lb) else 0

            n_rows_after += len(tmp_prob)

            # merge into global maps
            for pos, pr in tmp_prob.items():
                if (pos not in prob_by_pos[key]) or pr > prob_by_pos[key][pos]:
                    prob_by_pos[key][pos]   = pr
                    pocket_by_pos[key][pos] = tmp_poc[pos]
                label_by_pos[key][pos] = 1 if (label_by_pos[key].get(pos,0) or tmp_lab[pos]) else 0

    # build JSON in GOLD order + compute best pocket vs GOLD positives
    records = []
    report  = []
    missing_keys = 0

    for key in order:
        pdb, chain = key
        gold_positions = gold_pos.get(key, [])
        gold_positions_set = gold_pos_set.get(key, set())
        gold_positives = gold_pos_positive.get(key, set())

        pos_map  = prob_by_pos.get(key, {})   # may be empty
        lab_map  = label_by_pos.get(key, {})
        poc_map  = pocket_by_pos.get(key, {})

        if not pos_map:
            missing_keys += 1

        # Align strictly to GOLD positions and order (critical for downstream eval)
        positions = list(gold_positions)
        probs  = [float(pos_map.get(p, 0.0)) for p in positions]
        labels = [int(lab_map.get(p, 0))     for p in positions]   # predicted binary (pocket>0)
        pockets= [int(poc_map.get(p, -1))    for p in positions]   # -1 if unknown/no pocket

        # choose best pocket vs GOLD positives (use only positions that exist in GOLD)
        pocket_ids = sorted({pid for pid in pockets if pid >= 0})
        best_pid = None
        best_jac = -1.0
        pocket_stats = []

        # build per-pocket sets over GOLD-indexed positions
        pocket_to_pos = defaultdict(set)
        for p, pid in zip(positions, pockets):
            if pid >= 0:
                pocket_to_pos[pid].add(p)

        for pid in pocket_ids:
            P = pocket_to_pos[pid]
            inter = len(P & gold_positives)
            union = len(P | gold_positives) if gold_positives else 0
            jac = (inter/union) if union>0 else 0.0
            pocket_stats.append({
                "pocket_id": pid,
                "n_residues": len(P),
                "overlap_with_gold": inter,
                "jaccard_vs_gold": jac
            })
            if jac > best_jac:
                best_jac = jac; best_pid = pid

        # summary line
        n_pos = len(positions)
        n_pospos = sum(labels)
        mean_p = (sum(probs)/n_pos) if n_pos else 0.0
        max_p  = max(probs) if probs else 0.0
        tk = sorted(zip(positions, probs), key=lambda t: t[1], reverse=True)[:5]
        tk_str = "; ".join(f"{pp}:{pr:.3f}" for pp,pr in tk)
        cov = (len(set(positions) & gold_positions_set)/len(gold_positions_set)*100.0) if gold_positions_set else 0.0
        print(f"{pdb}.{chain}: N={n_pos}, +=#{n_pospos}, meanP={mean_p:.3f}, maxP={max_p:.3f}, "
              f"cover_vs_GOLD={cov:.1f}%, best_pocket={best_pid}, top5=[{tk_str}]")

        # JSON record
        records.append({
            "pdb_id": pdb,
            "chain_id": chain,
            "sequence_residue_numbers": positions,   # exactly GOLD positions/order
            "labels": labels,                        # predicted binary from P2Rank (pocket>0)
            "pred_probs": probs,                     # per-residue probabilities aligned to GOLD
            "residue_pocket_ids": pockets,           # -1 for unknown
            "best_pocket_id": best_pid,
            "pocket_overlap_stats": pocket_stats
        })

        report.append({
            "pdb_id": pdb, "chain_id": chain,
            "n_positions_p2": n_pos,
            "n_positive_p2": n_pospos,
            "mean_prob": f"{mean_p:.6f}",
            "max_prob": f"{max_p:.6f}",
            "coverage_vs_gold": f"{cov:.3f}",
            "best_pocket_id": (best_pid if best_pid is not None else ""),
            "top5_positions": tk_str
        })

    # write outputs
    with open(OUTPUT_JSON, "w") as f:
        json.dump(records, f, indent=2)
    with open(REPORT_CSV, "w", newline="") as f:
        if report:
            w = csv.DictWriter(f, fieldnames=list(report[0].keys()))
            w.writeheader(); w.writerows(report)

    print(f"Saved {OUTPUT_JSON} with {len(records)} chains")
    print(f"Saved {REPORT_CSV} with {len(report)} rows")
    if n_rows_total > 0:
        collapsed = n_rows_total - n_rows_after
        print(f"Collapsed insertion-code duplicates: {collapsed}/{n_rows_total} ({100*collapsed/n_rows_total:.2f}%)")
    print(f"Parsed P2Rank residue CSV files: {n_csv}")
    if missing_keys:
        print(f"[i] Chains present in GOLD but with no P2Rank CSV rows: {missing_keys} (filled with zeros)")

if __name__ == "__main__":
    main()
