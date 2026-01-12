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
Downloads mmCIF structures listed in the input Excel table.

- Reads an Excel file with a column `PDB` (default), takes the first 4 characters (PDB ID; any chain suffix is ignored).
- Downloads mmCIF files from RCSB (`https://files.rcsb.org/download/<PDB>.cif`).
- Saves them as `<pdb>.cif` (lowercase) into the target directory.
- Skips existing files unless `--force` is used.
- Retries failed downloads and removes suspiciously small files.

### `00.1_download_structures_run.sh`
Wrapper script to run structure download in batch mode on a selected dataset.
Edit paths to the Excel input and output directory as needed for each dataset.

### `00.2_check_structures.py`
Validates downloaded CIF files, checks basic consistency,
and reports missing or problematic structures.

- Reads expected PDB IDs from the Excel input table.
- Scans the structures directory for downloaded `.cif` files.
- Compares expected vs. present structures.
- Reports missing and extra PDB IDs.

This script is intended as a sanity check after structure download.

### `00.2_check_structures_run.sh`
Batch wrapper for structure validation.
Edit paths to the Excel input and output directory as needed for each dataset.

---

## 01 – Ground-truth (GOLD) residue labels

### `01.0_chain_gold_labels.py`
Constructs residue-level GOLD (ground-truth) binding-site annotations.

- Reads ligand annotations from the Excel input table.
- Loads mmCIF structures using author residue numbering.
- Resolves ligand components using a robust multi-step strategy:
  - author name + number
  - author number only
  - label name + number
  - label number only
- Applies controlled name-only fallbacks (including optional manual overrides),
  avoiding duplicate or ambiguous assignments.
- Classifies components (ligand, cap, modified amino acid, polymeric, water).
- Labels all amino-acid residues within a fixed distance cutoff (4.0 Å)
  of matched ligand atoms.

Multiple entries per structure–chain may be produced if multiple ligands
are present; downstream evaluation merges them as required.

### `01.1_overview_per_pdb.py`
Quick diagnostic overview: Excel vs downloaded structures vs GOLD chain labels.

- Extracts PDB identifiers from the Excel input table,
- lists downloaded structure files (`.cif` / `.cif.gz`),
- summarizes structure–chain pairs present in `01.chain_gold_labels.json`.

It reports counts and set differences between these sources and produces
CSV summaries highlighting missing or extra entries.

This script is intended for diagnostic and bookkeeping purposes only;
discrepancies between Excel rows, structure files, and GOLD chain entries
are common and often benign due to differences in data granularity
(evidence rows vs. structures vs. chain-level annotations).

### `01.2_gold_audit_summary.py`

Provides an end-to-end diagnostic overview of the GOLD label construction process.

The script cross-checks:
- the input Excel table (PDB identifiers and ligand annotations),
- locally available structure files (`.cif` / `.cif.gz`),
- the generated GOLD labels (`01.chain_gold_labels.json`),
- and the detailed audit log produced during GOLD construction.

It produces a concise text summary and multiple CSV reports describing
how individual ligand tokens were resolved, skipped, or not found,
and summarizes overlaps and discrepancies between Excel entries,
available structures, and GOLD chain-level annotations.

This script is intended for sanity checking and bookkeeping purposes only
and does not modify any data or enforce assumptions about completeness.

### `01y.probe_excel_only_pdbs.py`

Diagnostic helper script for investigating PDB entries that appear in the input Excel table
but do not produce any GOLD chain-level annotations.

For each such PDB, the script:
- lists raw ligand annotation cells from Excel,
- shows how ligand tokens would be parsed,
- reports available chains and non-water HET components present in the corresponding mmCIF file.

The output is intended to provide contextual information explaining why certain Excel entries
do not result in GOLD labels (e.g., missing ligands, non-resolvable tokens, or absence of
ligand-like components in the structure), without making assumptions about data correctness.

---

## 02 – P2Rank predictions

### `02.0_run_p2rank_batch.sh`
Batch script for running P2Rank predictions on all mmCIF structures in a dataset.

- Iterates over all `.cif` files in `<dataset>/structures/`,
- runs P2Rank (v2.5) in prediction mode for each structure,
- stores raw P2Rank outputs in `<dataset>/p2rank_predictions/`.

The P2Rank binary path and dataset folder can be configured via environment
variables (`P2RANK_BIN`, `KINASE_FOLDER`). Visualizations are disabled to allow
headless and batch execution.

### `02.1_chain_p2rank_labels.py`

Parses raw P2Rank per-residue prediction outputs and aligns them to the
canonical GOLD residue order for each (PDB ID, chain).

- Aggregates P2Rank residue-level probabilities across pockets,
- resolves insertion-code duplicates by taking the maximum probability,
- aligns predictions to GOLD residue numbering and ordering,
- fills missing positions with zero probability,
- optionally identifies the P2Rank pocket with the highest overlap to GOLD
  annotations (for reporting and diagnostic purposes only).

### `02.2_build_p2rank_rank1_map.py`

Extracts the top-ranked (rank-1) P2Rank pocket for each structure and maps
its residues to individual protein chains present in the GOLD annotations.

- Robustly identifies PDB IDs from P2Rank prediction outputs,
- extracts the highest-ranked pocket (rank-1) and its residue set,
- intersects rank-1 pocket residues with GOLD residue annotations,
- reports coverage- and precision-like overlap statistics for diagnostic purposes.

---

## 03 – PLM predictions

### `03.0_plm_predict_from_gold.py`
Runs predictions of a fine-tuned protein language model (PLM) on protein chains
defined by the GOLD annotations.

- Uses GOLD annotations only to determine which (PDB, chain) pairs to process,
- extracts amino-acid sequences from mmCIF structures using author residue numbering,
- performs residue-level inference with a fine-tuned ESM2-based model,
- outputs dense per-residue probability scores aligned to GOLD residue indices.

The output is designed for downstream evaluation and comparison with
structure-based methods (e.g., P2Rank). No ground-truth labels are used during inference.

### `03.1_chain_plm_labels.py`

Aligns raw PLM inference outputs to the canonical residue order defined by
GOLD annotations.

- Normalizes PLM outputs of different possible formats
(dense arrays, sparse indices, residue-number–based predictions)
into a single, consistent per-residue representation aligned to GOLD chains.

Optional fixed-threshold labels are provided
for inspection only and are not used in quantitative evaluation.

### `03y.inspect_gold_vs_plm_keys.py`

Diagnostic helper script that compares the set of (PDB, chain) pairs
present in GOLD annotations with those present in raw PLM prediction outputs.

- Reports missing and extra chains, groups discrepancies per PDB,
and provides conservative heuristic suggestions for possible chain ID aliases.

---

## 04 – Evaluation

### `04.0_evaluate_pipeline_all_in_one.py`
Evaluates residue-level binding-site predictions against GOLD annotations
for all implemented methods in a single, consistent framework.

- Aligns GOLD labels, P2Rank predictions, and PLM predictions
on a common author residue numbering and computes standard classification
metrics from per-residue probabilities.
- For P2Rank, the primary evaluated output is the per-residue prediction
aggregated across all predicted pockets
(`02.chain_p2rank_labels.json`, ANY-pocket variant).
- Additional P2Rank variants (top-ranked pocket and oracle-selected pocket)
are derived internally and used for diagnostic and comparative analysis.

Reported metrics include:
- accuracy,
- weighted F1-score,
- Matthews correlation coefficient (MCC),
- AUROC and AUPR.

Both **micro-averaged** (all residues pooled) and **macro-averaged**
(per-chain) metrics are produced.

The following prediction variants are evaluated:
- P2Rank (all predicted pockets combined),
- P2Rank top-ranked pocket (rank-1),
- P2Rank oracle pocket (best single pocket selected using GOLD; diagnostic upper bound),
- PLM predictions with a fixed probability threshold,
- PLM predictions with an automatically optimized threshold (best MCC).

In addition to global metrics, the script also exports per-chain diagnostic
tables for PLM and P2Rank oracle predictions.

---

## 05 – Final report and visualisation

### `05.report_from_evaljson_plus_curves.py`

Generates a combined, publication-ready evaluation report from previously
computed results and raw prediction scores. 
Aggregates evaluation outputs across multiple kinase datasets
and produces a unified set of figures and tables, including:

- bar plots of MICRO- and MACRO-averaged performance metrics
  (accuracy, weighted F1, MCC) for all evaluated methods,
- ROC and precision–recall (PR) curves for PLM predictions
  computed directly from raw residue-level probabilities,
- summary tables with AUROC and AUPR values derived from raw scores.

In addition to figures, the script automatically creates a concise
Markdown report that references all generated plots and tables,
allowing easy inspection of results without rerunning the full pipeline.

### `05.report_from_evaljson_plus_curves_SVG.py`
SVG-based variant of the plotting script producing vector graphics
for publication-quality figures.

### `05.report_from_evaljson_plus_curves_run.sh`
Wrapper script for batch generation of reports and figures.

---

## 06 – Class prevalence (dataset imbalance)

### `06.0_prevalence.py`
Computes residue-level class prevalence (class imbalance) for each kinase dataset.

The prevalence is computed at the level of individual protein–chain–ligand complexes.
Each complex (i.e. one structure, one protein chain, and one bound ligand) is treated
as an independent sample and contributes its full set of residue-level labels.

As a consequence, the same protein chain may appear multiple times in the dataset if
it is associated with multiple ligands, which explains why the number of evaluated
chains can exceed the number of unique protein structures.

This analysis provides essential context for interpreting evaluation metrics, in
particular MCC and AUPR, and highlights the strong class imbalance typical for
residue-level binding-site prediction tasks.

---
## Auxiliary scripts (tables for SI)
---
Script is located in the `si_helpers/` folder: `si_helpers/select_one_structure_per_protein.py`

### `si_helpers/select_one_structure_per_protein.py`
Utility script used to create a **1 protein = 1 structure** subset of the input table
(e.g., for Supporting Information). The script groups entries by `UniprotID` and keeps
a single representative structure per protein (by default preferring lower `Resolution`
when available; ties are broken by `PDB`).

## Diagnostic and auxiliary scripts

Scripts in this section are intended for debugging, sanity checks, and manual inspection
of individual structures. Diagnostic scripts are located in the `diagnostics/` folder:

 `diagnostics/inspect_alignment_for_pdb.py`
 
 `diagnostics/inspect_alignment_for_pdb_run.sh`
 
 `diagnostics/inspect_in_pymol.pml`

### `diagnostics/inspect_alignment_for_pdb.py`
Sanity-check alignment for a single PDB across all chains present in GOLD.  
For each chain it verifies that GOLD / P2Rank / PLM use the same `sequence_residue_numbers` (author residue numbering), reports length/position mismatches, and summarizes overlaps (counts + Jaccard) between GOLD positives, P2Rank positives, and PLM positives (at given thresholds).  

### `diagnostics/inspect_alignment_for_pdb_run.sh`
Convenience wrapper to run `inspect_alignment_for_pdb.py` for a chosen dataset folder and PDB ID with predefined thresholds and an optional CSV output path.

### `diagnostics/inspect_in_pymol.pml`
PyMOL script for manual visual inspection of a single structure.
Loads the mmCIF structure and highlights GOLD, P2Rank, and PLM positive residues
using different colors (GOLD = red, P2Rank = cyan, PLM = yellow).

Run from the project root using:
```console
pymol inspect_in_pymol.pml
```

## Notes

- Scripts are designed to be run sequentially but are modular.

- Large intermediate data files are not stored in this repository
  and are provided via external data storage (see main README).

- All figures and metrics reported in the manuscript
  are directly reproducible from the outputs of these scripts.
