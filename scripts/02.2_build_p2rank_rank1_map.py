#!/usr/bin/env python3
"""
02.2_build_p2rank_rank1_map.py  (more robust)

Build per-(PDB, chain) RANK-1 pocket map from P2Rank outputs.

Fixes vs v3:
- Recursive discovery of *_predictions.csv (os.walk).
- Robust PDB id inference:
    1) from filename (e.g., '4czt_predictions.csv', '4czt.cif_predictions.csv')
    2) from parent folder (e.g., '4czt.cif/')
    3) from CSV 'name' column (often contains '<pdb>.cif')
- Flexible column names: accepts residue_ids | residues | residueIds, etc.
- Does NOT require *_residues.csv; falls back to GOLD position order per chain.
- Better diagnostics: counts, samples, and reasons.

Outputs:
  <KINASE_FOLDER>/02.p2rank_rank1_map.json
  <KINASE_FOLDER>/02.p2rank_rank1_report.csv
"""

import os, re, json, sys, csv
import pandas as pd
from collections import defaultdict

# ---------- CONFIG ----------
BASE_DIR = "."
KINASE_FOLDER = "Kinase_Type_ALLO"
GOLD_JSON = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
P2RANK_DIR = os.path.join(BASE_DIR, KINASE_FOLDER, "p2rank_predictions")

OUT_JSON = os.path.join(BASE_DIR, KINASE_FOLDER, "02.p2rank_rank1_map.json")
OUT_CSV  = os.path.join(BASE_DIR, KINASE_FOLDER, "02.p2rank_rank1_report.csv")
# ----------------------------

# 4-char PDB token like '1CKJ' anywhere in text
PDB_TOKEN = re.compile(r'([0-9][A-Za-z0-9]{3})')

# chain/pos token like "A_172", "B:110", "C 265"
RESID_TOKEN = re.compile(r'^\s*([A-Za-z])[_:\s]?([+-]?\d+)')

def parse_pos(x):
    s = str(x).strip()
    m = re.match(r'^\s*([+-]?\d+)', s)
    if not m: raise ValueError(f"Cannot parse residue position from: {x!r}")
    return int(m.group(1))

def load_gold_index(path):
    """Index GOLD by (PDB, chain) -> {pos: 0/1}, OR-merge duplicates."""
    raw = json.load(open(path))
    idx = {}
    for e in raw:
        key = (e["pdb_id"].upper(), e["chain_id"].upper())
        pos = [parse_pos(x) for x in e["sequence_residue_numbers"]]
        lbl = [int(x) for x in e["labels"]]
        if len(pos) != len(lbl): continue
        if key not in idx: idx[key] = {}
        mp = idx[key]
        for p,y in zip(pos,lbl):
            mp[p] = 1 if (mp.get(p,0) or y) else 0
    return idx

def guess_pdb_from_path(pred_path):
    """
    Try to guess PDB id:
      a) from filename
      b) from parent folder (strip extension like '.cif')
      c) as last resort: read CSV 'name' column and regex there
    """
    fn = os.path.basename(pred_path)
    m = PDB_TOKEN.search(fn)
    if m: return m.group(1).upper()

    parent = os.path.basename(os.path.dirname(pred_path))
    # parent like '4czt.cif' or '4czt'
    m = PDB_TOKEN.search(parent)
    if m: return m.group(1).upper()

    # fallback: peek CSV 'name' column
    try:
        head = pd.read_csv(pred_path, nrows=5, sep=None, engine="python")
        for col in head.columns:
            if str(col).strip().lower() == "name":
                for v in head[col].astype(str).tolist():
                    mm = PDB_TOKEN.search(v)
                    if mm: return mm.group(1).upper()
    except Exception:
        pass
    return None

def list_predictions_recursive(p2dir):
    """Return list of prediction file paths (endswith '_predictions.csv')."""
    out = []
    for root, _, files in os.walk(p2dir):
        for fn in files:
            if fn.lower().endswith("_predictions.csv"):
                out.append(os.path.join(root, fn))
    return sorted(out)

def read_rank1_row(pred_path):
    """
    Return dict(score, probability, residue_ids_str) for rank==1 (or 0), None on failure.
    Accepts varied column names.
    """
    try:
        df = pd.read_csv(pred_path, sep=None, engine="python")
    except Exception as ex:
        print(f"[!] Cannot read {pred_path}: {ex}", file=sys.stderr)
        return None

    # normalize columns
    cols_map = {c.strip().lower(): c for c in df.columns}
    need_rank = next((cols_map[k] for k in cols_map if k == "rank"), None)
    col_score = next((cols_map[k] for k in cols_map if k in ("score","scores")), None)
    col_prob  = next((cols_map[k] for k in cols_map if k in ("probability","prob","probabilities")), None)
    col_resid = next((cols_map[k] for k in cols_map if k in ("residue_ids","residues","residueids")), None)

    if not (need_rank and col_score and col_prob and col_resid):
        print(f"[!] Missing expected columns in {pred_path}", file=sys.stderr)
        return None

    # choose rank==1, fallback to rank==0; if multiple rows with same rank, take the best score
    try:
        ranks = pd.to_numeric(df[need_rank], errors="coerce").fillna(0).astype(int)
    except Exception:
        ranks = pd.Series([0]*len(df))
    r1 = df[ranks == 1]
    if r1.empty:
        r1 = df[ranks == 0]
    if r1.empty:
        # if ranks weird, just take the top by score
        r1 = df
    try:
        r1 = r1.sort_values(col_score, ascending=False).iloc[0]
    except Exception:
        return None

    return {
        "residue_ids": str(r1[col_resid]),
        "score": float(pd.to_numeric(r1[col_score], errors="coerce") or 0.0),
        "probability": float(pd.to_numeric(r1[col_prob],  errors="coerce") or 0.0),
    }

def parse_residue_ids_by_chain(resid_str):
    """Turn 'A_172 B_110 ...' into {'A': {172,...}, 'B': {110,...}}."""
    out = defaultdict(set)
    for tok in str(resid_str).replace(",", " ").split():
        m = RESID_TOKEN.match(tok)
        if not m: continue
        ch  = m.group(1).upper()
        pos = int(m.group(2))
        out[ch].add(pos)
    return out

def main():
    gold_idx = load_gold_index(GOLD_JSON)
    valid_chains = set(gold_idx.keys())
    gold_pdbs = {p for p,_ in valid_chains}
    print(f"[i] GOLD chains: {len(valid_chains)} across {len(gold_pdbs)} PDBs")

    pred_files = list_predictions_recursive(P2RANK_DIR)
    print(f"[i] Found {len(pred_files)} predictions.csv files in {P2RANK_DIR} (recursive)")

    # Collect rank-1 info per PDB
    rank1_by_pdb = {}
    reasons = defaultdict(int)
    sample_no_pdb = []
    sample_no_rank = []

    for pred_path in pred_files:
        pdb = guess_pdb_from_path(pred_path)
        if not pdb:
            reasons["no_pdb_inferred"] += 1
            if len(sample_no_pdb) < 8: sample_no_pdb.append(os.path.basename(pred_path))
            continue

        r1 = read_rank1_row(pred_path)
        if r1 is None:
            reasons["bad_predictions_csv"] += 1
            if len(sample_no_rank) < 8: sample_no_rank.append(os.path.basename(pred_path))
            continue

        r1_by_chain = parse_residue_ids_by_chain(r1["residue_ids"])
        rank1_by_pdb[pdb] = {
            "score": r1["score"],
            "prob":  r1["probability"],
            "by_chain": r1_by_chain
        }

    # Emit per-(PDB,CHAIN) JSON only for chains present in GOLD
    json_records = []
    report_rows  = []
    overlap_pdbs = sorted(gold_pdbs & set(rank1_by_pdb.keys()))
    if not overlap_pdbs:
        print("[!] No PDB overlap between GOLD and P2Rank predictions. Check naming/layout.", file=sys.stderr)

    for (pdb_id, chain_id) in sorted(valid_chains):
        if pdb_id not in rank1_by_pdb:
            continue
        r1info = rank1_by_pdb[pdb_id]
        rank1_pos = sorted(r1info["by_chain"].get(chain_id, set()))
        gold_map  = gold_idx.get((pdb_id, chain_id), {})
        gold_total_pos = sum(gold_map.values()) if gold_map else 0
        rank1_gold_pos = sum(gold_map.get(p, 0) for p in rank1_pos) if gold_map else 0

        N_r1  = len(rank1_pos)
        cov   = (rank1_gold_pos / gold_total_pos) if gold_total_pos else None
        prec  = (rank1_gold_pos / N_r1) if N_r1 else None

        json_records.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "rank1_positions": rank1_pos,
            "rank1_pocket_score": r1info["score"],
            "rank1_pocket_probability": r1info["prob"]
        })

        report_rows.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "N_rank1": N_r1,
            "rank1_AND_GOLD_pos": rank1_gold_pos,
            "GOLD_total_pos": gold_total_pos,
            "coverage_of_GOLD": cov,
            "precision_like": prec,
            "rank1_score": r1info["score"],
            "rank1_probability": r1info["prob"]
        })

    # Save
    with open(OUT_JSON, "w") as f:
        json.dump(json_records, f, indent=2)
    pd.DataFrame(report_rows).to_csv(OUT_CSV, index=False)

    print(f"Saved {OUT_JSON} with {len(json_records)} chain entries")
    print(f"Saved {OUT_CSV} with {len(report_rows)} rows")
    # Diagnostics
    print("[diag] reasons (predictions skipped):", dict(reasons))
    print(f"[diag] PDB overlap GOLD vs P2Rank: {len(overlap_pdbs)}")
    if sample_no_pdb:
        print("[diag] sample files with no PDB inferred:", ", ".join(sample_no_pdb))
    if sample_no_rank:
        print("[diag] sample files with unreadable ranks/columns:", ", ".join(sample_no_rank))
    # show a few non-empty rank1
    shown = 0
    for r in json_records:
        if r["rank1_positions"]:
            print(f"  {r['pdb_id']}.{r['chain_id']}: "
                  f"N_rank1={len(r['rank1_positions'])}, "
                  f"score={r['rank1_pocket_score']:.3f}, prob={r['rank1_pocket_probability']:.3f}, "
                  f"top5={r['rank1_positions'][:5]}")
            shown += 1
            if shown >= 10: break

if __name__ == "__main__":
    main()
