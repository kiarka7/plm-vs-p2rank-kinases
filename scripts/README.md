# Analysis Pipeline Scripts

This directory contains all scripts used to preprocess data, generate predictions,
evaluate models, and produce figures reported in the manuscript.

The pipeline operates on individual protein structure–chain entries
identified by (PDB ID, chain ID) and aligns all residue-level annotations
using numeric residue indices.

---

## Overview of the pipeline

Scripts are organized by numeric prefixes reflecting the recommended execution order:

00. Structure download and validation  
01. Construction and audit of ground-truth (GOLD) residue labels  
02. P2Rank prediction and post-processing  
03. PLM prediction and post-processing  
04. Evaluation and metric computation  
05. Figure and report generation  
06. Class prevalence (imbalance) analysis  

Each step can be run independently per dataset (e.g. `Kinase_Type_I`).

---

## 00 – Structure download and validation

### `00.1_download_structures.py`
Downloads protein structures listed in the input Excel table
and saves them in CIF format.

### `00.1_download_structures_run.sh`
Wrapper script to run structure download in batch mode.

### `00.2_check_structures.py`
Validates downloaded CIF files, checks basic consistency,
and reports missing or problematic structures.

### `00.2_check_structures_run.sh`
Batch wrapper for structure validation.

---

## 01 – Ground-truth (GOLD) residue labels

### `01.0_chain_gold_labels.py`
Constructs residue-level binary ground-truth labels based on ligand proximity.
Each structure–chain entry is annotated independently.

### `01.1_overview_per_pdb.py`
Generates summary statistics per PDB entry
(number of chains, labeled residues, ligands).

### `01.2_gold_audit_summary.py`
Produces audit summaries of GOLD labels
(e.g. residue counts, positives vs. negatives).

### `01y.probe_excel_only_pdbs.py`
Diagnostic script identifying PDB entries present in the Excel input
but missing from processed datasets.

---

## 02 – P2Rank predictions

### `02.0_run_p2rank_batch.sh`
Runs P2Rank in batch mode on all structures within a dataset.

### `02.1_chain_p2rank_labels.py`
Parses raw P2Rank outputs and converts them into
residue-level labels aligned to numeric residue indices.

### `02.2_build_p2rank_rank1_map.py`
Builds mapping of top-ranked P2Rank pockets (rank 1)
for downstream evaluation and oracle analyses.

---

## 03 – PLM predictions

### `03.0_plm_predict_from_gold.py`
Runs prediction using a fine-tuned Protein Language Model (PLM)
on residues defined in the GOLD annotations.

### `03.1_chain_plm_labels.py`
Post-processes PLM probability outputs and aligns them
to residue indices for evaluation.

### `03y.inspect_gold_vs_plm_keys.py`
Diagnostic script verifying consistency of structure–chain keys
between GOLD labels and PLM predictions.

---

## 04 – Evaluation

### `04.0_evaluate_pipeline_all_in_one.py`
Computes evaluation metrics for all methods, including:
- accuracy,
- weighted F1,
- Matthews correlation coefficient (MCC),
- AUROC and AUPR (from raw probabilities).

Both micro- and macro-averaged metrics are produced.

---

## 05 – Figure and report generation

### `05.report_from_evaljson_plus_curves.py`
Generates all plots and tables reported in the manuscript
from evaluation JSON files, including ROC and PR curves.

### `05.report_from_evaljson_plus_curves_SVG.py`
SVG-based variant of the plotting script producing vector graphics
for publication-quality figures.

### `05.report_from_evaljson_plus_curves_run.sh`
Wrapper script for batch generation of reports and figures.

---

## 06 – Class prevalence analysis

### `06.0_prevalence.py`
Computes class prevalence (class imbalance) for residue-level annotations.
Prevalence is calculated across all protein–ligand complexes
and reflects the fraction of positive residues among all annotated residues.

---

## Notes

- Scripts are designed to be run sequentially but are modular.
- Large intermediate data files are not stored in this repository
  and are provided via external data storage (see main README).
- All figures and metrics reported in the manuscript
  are directly reproducible from the outputs of these scripts.
