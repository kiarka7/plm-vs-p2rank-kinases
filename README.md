# PLM vs P2Rank on Kinase Binding Sites

This repository contains preprocessing, evaluation and plotting scripts used in the
analysis comparing a fine-tuned Protein Language Model (PLM) against P2Rank for
residue-level binding-site prediction across kinase families.

## Repository structure

- `scripts/`  
  End-to-end pipeline for:
  - preprocessing of GOLD labels (ground truth),
  - parsing and normalisation of P2Rank outputs,
  - inference and post-processing of PLM predictions,
  - computation of micro/macro metrics,
  - generation of all figures reported in the manuscript.

- `figures/`  
  Exported bar plots, ROC and PR overlays.

- `env/`  
  Environment specifications (`requirements.txt`).

- `model/`  
  Placeholder for the fine-tuned PLM checkpoint (not included in the public repo).

## Reproducibility

All evaluations were performed per structure–chain key and aligned by numeric
residue indices. Only the **clean** leakage-controlled subset was used for
computing metrics reported in the manuscript.

A full reproducibility section (including instructions for running the entire
pipeline end-to-end) will be added in the final version.
