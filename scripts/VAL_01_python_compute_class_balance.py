#!/usr/bin/env python3
import json
from pathlib import Path
import pandas as pd

BASE = Path(".")

datasets = [
    "Kinase_Type_ALLO",
    "Kinase_Type_I",
    "Kinase_Type_I.5",
    "Kinase_Type_II",
    "Kinase_Type_III",
]

rows = []

for ds in datasets:
    gold_path = BASE / ds / "01.chain_gold_labels_CLEAN.json"
    with open(gold_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    n_chains = 0
    n_pos = 0
    n_neg = 0

    for rec in data:
        labels = rec["labels"]
        n_chains += 1
        pos = sum(labels)
        neg = len(labels) - pos
        n_pos += pos
        n_neg += neg

    total = n_pos + n_neg
    prevalence = n_pos / total if total > 0 else float("nan")

    rows.append({
        "dataset": ds,
        "chains": n_chains,
        "positives": n_pos,
        "negatives": n_neg,
        "prevalence": prevalence * 100.0,
    })

df = pd.DataFrame(rows)
df.to_csv("class_balance_clean.csv", index=False)
print(df)
