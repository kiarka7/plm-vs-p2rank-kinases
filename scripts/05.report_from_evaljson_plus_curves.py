#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05.report_from_evaljsons_plus_curves.py

Combined evaluation report (plots + markdown).

What this script does
---------------------
1) Reads per-dataset evaluation JSON files:
       04.eval_all_in_one.json
   and generates:
   - MICRO bar plots for ACC, F1 (weighted), MCC
   - MACRO bar plots (mean ± IQR/2) for ACC, F1 (weighted), MCC

2) Reads raw GOLD and PLM prediction files:
       01.chain_gold_labels.json
       03.plm_predictions.json
   and generates:
   - PLM ROC curves (micro, all residues concatenated)
   - PLM PR curves (micro)
   - AUROC/AUPR bar plots from raw scores
   - sanity CSVs with score statistics and join diagnostics

3) Writes a small Markdown report referencing all figures and
   summarising AUROC/AUPR numbers.

Usage
-----
python 05.report_from_evaljsons_plus_curves.py \
  --base "." \
  --out "_combined_final" \
  --include Kinase_Type_ALLO Kinase_Type_I Kinase_Type_I.5 Kinase_Type_II Kinase_Type_III
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

# Headless plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score
)

# ------------------------------------------------------------
# File names
# ------------------------------------------------------------

EVAL_FILE = "04.eval_all_in_one_CLEAN.json"
GOLD_FILE = "01.chain_gold_labels_CLEAN.json"
PRED_FILE = "03.plm_predictions_CLEAN.json"

# ------------------------------------------------------------
# Visual style / aesthetics
# ------------------------------------------------------------

# Color palette for datasets (consistent across all plots)
DATASET_COLORS = {
    "Kinase_Type_ALLO":"#991313",  # red
    "Kinase_Type_I":   "#1f77b4",  # dark blue
    "Kinase_Type_I.5": "#17becf",  # cyan
    "Kinase_Type_II":  "#2ca02c",  # green
    "Kinase_Type_III": "#bdb76b",  # khaki
}

# Methods distinguished by hatching and slight alpha
METHODS = [
    ("P2 ANY fixed", "p2_any_fixed", "",    0.80),
    ("PLM fixed",    "plm_fixed",    "///", 0.95),
    ("PLM bestMCC",  "plm_best_mcc", "xx",  0.95),
]

# Y-limit for metric bar plots (a bit higher to leave room for labels)
YMAX_BAR = 1.10

# Error-bar settings (IQR/2) – semi-transparent
ERROR_KW = dict(
    elinewidth=1.5,
    capsize=5,
    capthick=1.2,
    ecolor=(0, 0, 0, 0.35),
)

# Global typography (aimed at publication-quality figures)
BASE_FONT_SIZE   = 14
AXIS_LABEL_SIZE  = 16
TICK_LABEL_SIZE  = 14
TITLE_SIZE       = 18
LEGEND_FONT_SIZE = 14
ANNOT_FONT_SIZE  = 13

plt.rcParams.update({
    "font.size": BASE_FONT_SIZE,
    "axes.labelsize": AXIS_LABEL_SIZE,
    "axes.titlesize": TITLE_SIZE,
    "xtick.labelsize": TICK_LABEL_SIZE,
    "ytick.labelsize": TICK_LABEL_SIZE,
    "legend.fontsize": LEGEND_FONT_SIZE,
    "axes.linewidth": 1.4,
})

# Optional suffix for plot titles (e.g. " on CLEAN subset")
TITLE_SUFFIX = ""  # set to " on CLEAN subset" if you wish

# ------------------------------------------------------------
# Helper utilities
# ------------------------------------------------------------

def nice_name(ds: str) -> str:
    """Prettier dataset name for legends (e.g. 'Kinase Type I.5')."""
    return ds.replace("_", " ")

def load_json(p: Path):
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None

def need_cols(d: dict, keys: List[str]) -> bool:
    return all(k in d for k in keys)

def ensure_list(obj):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        return list(obj.values())
    return []

def fig_save(path: Path):
    """Tight layout + reasonably high DPI."""
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# F1 (weighted) helpers
# ------------------------------------------------------------

WEIGHTED_NOTES: List[str] = []

def pick_f1_micro(block: dict, ds: str, method_label: str) -> float:
    """
    For MICRO metrics, prefer weighted F1 if present; fall back to 'f1'.
    """
    for k in ("f1_weighted", "f1_w", "f1weighted"):
        if k in block:
            return float(block[k])
    if "f1" in block:
        WEIGHTED_NOTES.append(
            f"[{ds}] {method_label}: weighted F1 not found -> using 'f1' (fallback)."
        )
        return float(block["f1"])
    return float("nan")

def pick_f1_macro(block: dict, ds: str, method_label: str) -> float:
    """
    For MACRO metrics, prefer mean of 'f1_weighted'; otherwise mean of 'f1'.
    Expects structure { 'f1_weighted': {'mean':..,'median':..,'iqr':..}, ... }.
    """
    src = None
    if isinstance(block.get("f1_weighted"), dict) and "mean" in block["f1_weighted"]:
        src = block["f1_weighted"]
    elif isinstance(block.get("f1"), dict) and "mean" in block["f1"]:
        src = block["f1"]
        WEIGHTED_NOTES.append(
            f"[{ds}] {method_label}: weighted F1 (macro) not found -> using 'f1' (fallback)."
        )
    if src is None:
        return float("nan")
    return float(src["mean"])

def pick_f1_macro_iqr(block: dict) -> float:
    if isinstance(block.get("f1_weighted"), dict) and "iqr" in block["f1_weighted"]:
        return float(block["f1_weighted"]["iqr"])
    if isinstance(block.get("f1"), dict) and "iqr" in block["f1"]:
        return float(block["f1"]["iqr"])
    return float("nan")

# ------------------------------------------------------------
# Bar-plot helpers
# ------------------------------------------------------------

def _annotate_bar(ax, bar, value, y_max=YMAX_BAR):
    """Safely place numeric label above a bar without hitting the top frame."""
    if not np.isfinite(value):
        return
    x = bar.get_x() + bar.get_width() / 2.0
    y = min(bar.get_height() + 0.015, y_max - 0.02)
    ax.text(
        x, y, f"{value:.3f}",
        ha="center",
        va="bottom",
        fontsize=ANNOT_FONT_SIZE,
        zorder=5,
    )

def save_micro_bars_from_json(datasets, evals: Dict[str, dict], outdir: Path):
    metrics = [("acc", "Accuracy"), ("f1", "F1 (weighted)"), ("mcc", "MCC")]
    fnames  = ["bars_micro_acc.png", "bars_micro_f1.png", "bars_micro_mcc.png"]

    for (met, title_label), fname in zip(metrics, fnames):
        fig, ax = plt.subplots(figsize=(16, 6))
        x = np.arange(len(datasets))
        width = 0.25

        for i, (label, key, hatch, alpha) in enumerate(METHODS):
            vals = []
            for ds in datasets:
                try:
                    blk = evals[ds]["micro"][key]
                    if met == "f1":
                        vals.append(pick_f1_micro(blk, ds, label))
                    else:
                        vals.append(float(blk[met]))
                except Exception:
                    vals.append(np.nan)

            colors = [DATASET_COLORS.get(ds, "#888888") for ds in datasets]
            offs = x + (i - 1) * width
            bars = ax.bar(
                offs, vals, width,
                color=colors,
                edgecolor="black",
                linewidth=0.4,
                hatch=hatch,
                alpha=alpha,
                label=label,
            )
            for b, v in zip(bars, vals):
                _annotate_bar(ax, b, v, y_max=YMAX_BAR)

        ax.set_xticks(x)
        ax.set_xticklabels(datasets, rotation=20)
        ax.set_ylabel(title_label)
        ax.set_ylim(0, YMAX_BAR)
        ax.set_title(f"{title_label} (micro){TITLE_SUFFIX}")
        ax.legend(
            loc="lower left",
            ncol=3,
            frameon=True,
            framealpha=0.95,
            fancybox=True,
            edgecolor="black",
        )
        ax.tick_params(axis="both", which="both", length=5)
        fig_save(outdir / fname)

def save_macro_bars_from_json(datasets, evals: Dict[str, dict], outdir: Path):
    metrics = [("acc", "Accuracy"), ("f1", "F1 (weighted)"), ("mcc", "MCC")]
    fnames  = ["bars_macro_acc_mean.png", "bars_macro_f1_mean.png", "bars_macro_mcc_mean.png"]

    for (met, title_label), fname in zip(metrics, fnames):
        fig, ax = plt.subplots(figsize=(16, 6))
        x = np.arange(len(datasets))
        width = 0.25

        for i, (label, key, hatch, alpha) in enumerate(METHODS):
            means, iqrs = [], []
            for ds in datasets:
                try:
                    blk = evals[ds]["macro"][key]
                    if met == "f1":
                        means.append(pick_f1_macro(blk, ds, label))
                        iqrs.append(pick_f1_macro_iqr(blk))
                    else:
                        means.append(float(blk[met]["mean"]))
                        iqrs.append(float(blk[met]["iqr"]))
                except Exception:
                    means.append(np.nan)
                    iqrs.append(np.nan)

            errs = np.array(iqrs) / 2.0
            colors = [DATASET_COLORS.get(ds, "#888888") for ds in datasets]
            offs = x + (i - 1) * width
            bars = ax.bar(
                offs, means, width,
                yerr=errs,
                error_kw=ERROR_KW,
                color=colors,
                edgecolor="black",
                linewidth=0.4,
                hatch=hatch,
                alpha=alpha,
                label=label,
            )
            for b, v in zip(bars, means):
                _annotate_bar(ax, b, v, y_max=YMAX_BAR)

        ax.set_xticks(x)
        ax.set_xticklabels(datasets, rotation=20)
        ax.set_ylabel(title_label)
        ax.set_ylim(0, YMAX_BAR)
        ax.set_title(f"{title_label} (macro mean ± IQR/2){TITLE_SUFFIX}")
        ax.legend(
            loc="lower left",
            ncol=3,
            frameon=True,
            framealpha=0.95,
            fancybox=True,
            edgecolor="black",
        )
        ax.tick_params(axis="both", which="both", length=5)
        fig_save(outdir / fname)

# ------------------------------------------------------------
# RAW curve helpers (PLM ROC / PR)
# ------------------------------------------------------------

def collect_pairs_strict(ds_dir: Path) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Load GOLD and PLM predictions for a dataset directory, and
    return concatenated (y_true, scores) plus a small debug dict.
    """
    gold_p = ds_dir / GOLD_FILE
    pred_p = ds_dir / PRED_FILE

    if not gold_p.is_file() or not pred_p.is_file():
        return np.array([]), np.array([]), {
            "used_records": 0,
            "pairs": 0,
            "n_raw": 0,
            "dropped": 0,
            "reason": "missing_files",
        }

    gold = ensure_list(load_json(gold_p))
    preds = ensure_list(load_json(pred_p))

    gold_map = {}
    for rec in gold:
        if not isinstance(rec, dict):
            continue
        if not need_cols(rec, ["pdb_id", "chain_id", "labels"]):
            continue
        key = (rec["pdb_id"], rec["chain_id"])
        labels = list(map(int, rec["labels"]))
        gold_map[key] = {
            "labels": np.asarray(labels, dtype=int),
            "seq": rec.get("sequence_residue_numbers", None),
        }

    y_all, s_all = [], []
    used = 0
    dropped = 0
    join_rows = []

    for rec in preds:
        if not isinstance(rec, dict):
            continue
        if not need_cols(rec, ["pdb_id", "chain_id", "probabilities"]):
            continue
        key = (rec["pdb_id"], rec["chain_id"])
        probs = np.asarray(rec["probabilities"], dtype=float)

        if key not in gold_map:
            dropped += 1
            continue

        gl = gold_map[key]["labels"]
        if len(gl) != len(probs):
            join_rows.append({
                "pdb_id": key[0],
                "chain_id": key[1],
                "len_gold": len(gl),
                "len_probs": len(probs),
                "dataset": ds_dir.name,
            })
            dropped += 1
            continue

        y_all.append(gl)
        s_all.append(probs)
        used += 1
        join_rows.append({
            "pdb_id": key[0],
            "chain_id": key[1],
            "join": len(gl),
            "len_gold": len(gl),
            "len_probs": len(probs),
            "join_ratio_gold": 1.0,
            "join_ratio_probs": 1.0,
            "dataset": ds_dir.name,
        })

    y_all = np.concatenate(y_all) if y_all else np.array([], dtype=int)
    s_all = np.concatenate(s_all) if y_all.size else np.array([], dtype=float)

    return y_all, s_all, {
        "used_records": used,
        "pairs": int(y_all.size),
        "n_raw": len(preds),
        "dropped": dropped,
        "join_rows": join_rows,
        "reason": "ok",
    }

def plot_overlay(xs: Dict[str, np.ndarray], ys: Dict[str, np.ndarray],
                 xlabel: str, ylabel: str, title: str, out_png: Path,
                 dataset_order: List[str]):
    if not xs:
        return
    plt.figure(figsize=(14, 9))

    for name in dataset_order:
        if name not in xs:
            continue
        c = DATASET_COLORS.get(name, None)
        plt.plot(
            xs[name], ys[name],
            label=nice_name(name),
            color=c,
            linewidth=2.5,
        )

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title + TITLE_SUFFIX)
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.grid(alpha=0.25, linestyle="--", linewidth=0.8)

    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),  
        ncol=3,                       
        frameon=True,
        framealpha=0.95,
        fancybox=True,
        edgecolor="black",
    )

    fig_save(out_png)


def plot_bars(df: pd.DataFrame, value_col: str, title: str,
              ylabel: str, out_png: Path):
    if df.empty:
        return
    order = list(df["dataset"])
    vals = df[value_col].values
    plt.figure(figsize=(max(12, 1.6 * len(order)), 6))
    x = np.arange(len(order))
    colors = [DATASET_COLORS.get(ds, "#888888") for ds in order]
    bars = plt.bar(
        x, vals,
        color=colors,
        edgecolor="black",
        linewidth=0.4,
    )
    for b, v in zip(bars, vals):
        if pd.notnull(v):
            y = min(b.get_height() + 0.015, YMAX_BAR - 0.02)
            plt.text(
                b.get_x() + b.get_width() / 2.0,
                y,
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=ANNOT_FONT_SIZE,
                zorder=5,
            )
    plt.xticks(x, order, rotation=20)
    plt.ylabel(ylabel)
    plt.ylim(0, YMAX_BAR)
    plt.title(title + TITLE_SUFFIX)
    fig_save(out_png)

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="Root folder containing dataset subfolders.")
    ap.add_argument("--out",  required=True, help="Output folder for the combined report.")
    ap.add_argument("--include", nargs="*", default=None,
                    help="Whitelist of dataset subfolder names (if given).")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="Blacklist substrings (if name contains any, it is skipped).")
    args = ap.parse_args()

    base = Path(args.base).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- dataset directories ---
    ds_dirs: List[Path] = []
    for p in sorted(base.iterdir()):
        if not p.is_dir():
            continue
        name = p.name
        if args.include and name not in args.include:
            continue
        if any(ex.lower() in name.lower() for ex in args.exclude):
            continue
        ds_dirs.append(p)
    if not ds_dirs:
        raise SystemExit("No dataset folders found.")

    datasets = [d.name for d in ds_dirs]

    # --- load eval JSONs for bar plots ---
    evals: Dict[str, dict] = {}
    missing = []
    for d in ds_dirs:
        ejp = d / EVAL_FILE
        if ejp.is_file():
            ej = load_json(ejp)
            if isinstance(ej, dict):
                evals[d.name] = ej
            else:
                missing.append(d.name)
        else:
            missing.append(d.name)
    if missing:
        raise SystemExit("Missing eval JSON for: " + ", ".join(missing))

    save_micro_bars_from_json(datasets, evals, out_dir)
    save_macro_bars_from_json(datasets, evals, out_dir)

    # --- ROC/PR from raw scores ---
    roc_x, roc_y = {}, {}
    pr_x, pr_y = {}, {}
    auc_rows, aupr_rows = [], []
    debug_lines = []
    join_rows_all = []
    sanity_rows = []
    mismatch_rows = []

    for d in ds_dirs:
        ds = d.name
        y, s, dbg = collect_pairs_strict(d)

        if dbg.get("join_rows"):
            for jr in dbg["join_rows"]:
                if "join" in jr:
                    join_rows_all.append(jr)
                else:
                    mismatch_rows.append(jr)

        if y.size == 0 or s.size == 0:
            debug_lines.append(
                f"[{ds}] SKIP: {dbg.get('reason')} "
                f"(used={dbg['used_records']}, n_raw={dbg['n_raw']}, dropped={dbg['dropped']})"
            )
            continue

        uniq = np.unique(s)
        sanity_rows.append({
            "dataset": ds,
            "n": len(s),
            "min": float(np.min(s)),
            "max": float(np.max(s)),
            "mean": float(np.mean(s)),
            "unique_vals": len(uniq),
            "frac_exact_0": float(np.mean(s == 0.0)),
            "frac_exact_1": float(np.mean(s == 1.0)),
            "looks_binarized": bool(len(uniq) <= 4),
            "n_pos": int((y == 1).sum()),
            "n_neg": int((y == 0).sum()),
        })

        fpr, tpr, _ = roc_curve(y, s, pos_label=1)
        roc_auc = auc(fpr, tpr)
        roc_x[ds] = fpr
        roc_y[ds] = tpr
        auc_rows.append({"dataset": ds, "auroc_from_scores": float(roc_auc)})

        prec, rec, _ = precision_recall_curve(
            y, s,
            pos_label=1,
            drop_intermediate=False,
        )
        ap = average_precision_score(y, s, pos_label=1)
        pr_x[ds] = rec
        pr_y[ds] = prec
        aupr_rows.append({"dataset": ds, "aupr_from_scores": float(ap)})

        debug_lines.append(
            f"[{ds}] OK used_records={dbg['used_records']} "
            f"pairs={dbg['pairs']} dropped={dbg['dropped']}"
        )

        # Small inspection CSVs for score calibration
        order = np.argsort(s)
        low_ix = order[:50]
        mid_ix = order[len(order) // 2 - 25: len(order) // 2 + 25]
        high_ix = order[-50:][::-1]
        for tag, idx in [("LOW", low_ix), ("MID", mid_ix), ("HIGH", high_ix)]:
            df_ins = pd.DataFrame({
                "rank": np.arange(1, len(idx) + 1),
                "score": s[idx],
                "label": y[idx].astype(int),
            })
            df_ins.to_csv(out_dir / f"inspect_{ds}_top50_{tag}.csv", index=False)

    (out_dir / "_plm_curve_debug.txt").write_text(
        "\n".join(debug_lines),
        encoding="utf-8",
    )
    pd.DataFrame(sanity_rows).to_csv(out_dir / "scores_sanity.csv", index=False)
    if join_rows_all:
        pd.DataFrame(join_rows_all).to_csv(out_dir / "join_stats.csv", index=False)
    if mismatch_rows:
        pd.DataFrame(mismatch_rows).to_csv(out_dir / "mismatch_lengths.csv", index=False)

    if roc_x:
        plot_overlay(
            roc_x, roc_y,
            "False positive rate (FPR)",
            "True positive rate (TPR)",
            "PLM ROC curves (micro)",
            out_dir / "plm_roc_overlay.png",
            datasets,
        )
    if pr_x:
        plot_overlay(
            pr_x, pr_y,
            "Recall",
            "Precision",
            "PLM PR curves (micro)",
            out_dir / "plm_pr_overlay.png",
            datasets,
        )

    df_auc = pd.DataFrame(auc_rows).sort_values("dataset")
    df_aupr = pd.DataFrame(aupr_rows).sort_values("dataset")
    if not df_auc.empty:
        df_auc.to_csv(out_dir / "plm_auroc_from_scores.csv", index=False)
        plot_bars(
            df_auc.rename(columns={"auroc_from_scores": "val"}),
            "val",
            "PLM AUROC (micro) across datasets (from raw scores)",
            "AUROC",
            out_dir / "bars_plm_auroc.png",
        )
    if not df_aupr.empty:
        df_aupr.to_csv(out_dir / "plm_aupr_from_scores.csv", index=False)
        plot_bars(
            df_aupr.rename(columns={"aupr_from_scores": "val"}),
            "val",
            "PLM AUPR (micro) across datasets (from raw scores)",
            "AUPR",
            out_dir / "bars_plm_aupr.png",
        )

    # --------------------------------------------------------
    # Markdown report
    # --------------------------------------------------------
    md = []
    md.append("# Combined report")
    md.append("")
    md.append(f"_Generated: {pd.Timestamp.utcnow().isoformat()}Z_")
    md.append("")
    md.append("## Micro bars (from eval JSON)")
    md += [
        "",
        "![acc](bars_micro_acc.png)",
        "",
        "![f1](bars_micro_f1.png)",
        "",
        "![mcc](bars_micro_mcc.png)",
        "",
    ]
    md.append("## Macro bars (mean ± IQR/2, from eval JSON)")
    md += [
        "",
        "![acc](bars_macro_acc_mean.png)",
        "",
        "![f1](bars_macro_f1_mean.png)",
        "",
        "![mcc](bars_macro_mcc_mean.png)",
        "",
    ]
    if roc_x:
        md.append("## PLM ROC curves (micro) — from raw scores")
        md += ["", "![roc](plm_roc_overlay.png)", ""]
    if pr_x:
        md.append("## PLM PR curves (micro) — from raw scores")
        md += ["", "![pr](plm_pr_overlay.png)", ""]
    if not df_auc.empty:
        md.append("## AUROC (from raw scores)")
        md.append("")
        md.append(pd.DataFrame(df_auc).to_markdown(index=False))
        md.append("")
    if not df_aupr.empty:
        md.append("## AUPR (from raw scores)")
        md.append("")
        md.append(pd.DataFrame(df_aupr).to_markdown(index=False))
        md.append("")
    if WEIGHTED_NOTES:
        md.append("## Notes")
        md.append("")
        md.append("```")
        md.extend(WEIGHTED_NOTES)
        md.append("```")

    (out_dir / "combined_report.md").write_text(
        "\n".join(md),
        encoding="utf-8",
    )

    print(f"Done. Outputs written to: {out_dir}")
    if WEIGHTED_NOTES:
        print("\n".join(WEIGHTED_NOTES))

if __name__ == "__main__":
    main()
