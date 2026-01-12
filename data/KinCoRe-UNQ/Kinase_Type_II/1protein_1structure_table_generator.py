#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

INPUT_XLSX = "Kinase_Ligands_Type II.xlsx"
SHEET_NAME = 0
OUT_XLSX = "Kinase_Ligands_Type II_oneStructPerProtein.xlsx"
OUT_TSV  = "Kinase_Ligands_Type II_oneStructPerProtein.tsv"

df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)

required = ["UniprotID", "PDB"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in input table: {missing}\nColumns present: {list(df.columns)}")

sort_cols = ["UniprotID"]
if "Resolution" in df.columns:
    sort_cols += ["Resolution"]
sort_cols += ["PDB"]

df_sorted = df.sort_values(sort_cols, kind="mergesort")
df_one = df_sorted.drop_duplicates(subset=["UniprotID"], keep="first").copy()

df_one.to_excel(OUT_XLSX, index=False)
df_one.to_csv(OUT_TSV, sep="\t", index=False)

print(f"Input rows: {len(df)}")
print(f"Unique proteins (UniprotID): {df['UniprotID'].nunique()}")
print(f"Output rows (1 per protein): {len(df_one)}")
print(f"Saved: {OUT_XLSX}")
print(f"Saved: {OUT_TSV}")
