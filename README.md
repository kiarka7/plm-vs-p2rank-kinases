# PLM vs P2Rank on Kinase Binding Sites

This repository and the accompanying Zenodo archive provide data and scripts used in the study:

**Protein Language Models and Structure-Based Machine Learning for Prediction of Allosteric Binding Sites in Protein Kinases: An Explainable AI Framework Guided by Energy Landscape–Derived Frustration**

---

## 1. Overview

This project compares residue-level binding-site predictions produced by: a fine-tuned Protein Language Model (PLM), and the structure-based pocket prediction method **P2Rank** (version 2.5), across multiple kinase datasets (Type I, I.5, II, III, and Allosteric).

The fine-tuned PLM model can be downloaded from this data storage: https://owncloud.cesnet.cz/index.php/s/8RSJqt60D2uWJNa

 The code and datasets related to this PLM are available at: https://github.com/skrhakv/LBS-pLM

All scripts required to reproduce both PLM and P2Rank predictions reported in this study are publicly available in this repository:
https://github.com/kiarka7/plm-vs-p2rank-kinases

The analysis emphasizes reproducibility and fair evaluation across protein–ligand complexes.

----

## 2. Repository Structure
### scripts/

Python scripts implementing the full analysis pipeline:

- `00.*` – download of structures in .cif format from `.xlsx` input tables 
- `01.*` – construction of GOLD (ground-truth) residue labels
- `02.*` – P2Rank predictions and post-processing
- `03.*` – PLM predictions and post-processing
- `04.*` – evaluation (micro/macro metrics, MCC, AUROC, AUPR)
- `05.*` – figures and report generation
- `06.*` – class prevalence (class imbalance) calculations

Scripts operate per `(PDB ID, chain ID)` and align predictions using numeric residue indices.
Additional control and diagnostic scripts are included for validation and auditing purposes.

----
### data/

Each dataset is stored in a separate folder, e.g.:

- `Kinase_Type_I/`
- `Kinase_Type_I.5/`
- `Kinase_Type_II/`
- `Kinase_Type_III/`
- `Kinase_Type_ALLO/`

Expected essential files per dataset include:
- `Kinase_Ligands_*.xlsx` (source tables; original KinCoRe / KinCoRe-UNQ exports)
- `structures/*.cif` (downloaded on demand from RCSB PDB; not part of the archive)
- `01.chain_gold_labels.json` (ground-truth residue labels derived from ligand annotations)
- `p2rank_predictions/*_predictions.csv` (raw P2Rank outputs)
- `02.chain_p2rank_labels.json` (processed P2Rank outputs aligned to residue indices)
- `02.p2rank_rank1_map.json` (rank-1 P2Rank pocket mapping)
- `03.plm_predictions.json` (PLM raw probability outputs)
- `03.chain_plm_labels.json` (PLM residue-level predictions aligned to residue indices)
- `04.eval_all_in_one.json` (main evaluation results)
- `04.p2rank_oracle_chain_details.json` (oracle evaluation details)

Generated reports and figures are stored in:
`_combined_final_CLEAN/`

See `scripts/README.md` for a detailed description.

### Expected directory layout

All scripts assume a fixed directory layout and are designed to be run
from the **repository root directory**.

The expected structure is:

plm-vs-p2rank-kinases/

```
│ ├── 00.1_download_structures.py
│ ├── 01.0_chain_gold_labels.py
│ ├── 02.0_run_p2rank_batch.sh
│ ├── ...
│ └── 06.0_prevalence.py
│ ├── Kinase_Type_I/
│ ├── Kinase_Type_I.5/
│ ├── Kinase_Type_II/
│ ├── Kinase_Type_III/
│ ├── Kinase_Type_ALLO/
│ └── _combined_final_CLEAN/
```

Scripts reference dataset folders using **relative paths** (e.g. `Kinase_Type_I/`)
and therefore require this directory structure to be preserved. Running scripts from a different working directory is not supported unless paths inside the scripts are adjusted accordingly.

### Pipeline (quick run)

The pipeline is executed **per dataset** (e.g. `Kinase_Type_I`). Dataset names are specified directly in individual scripts where required.

Typical execution order:

```console
bash 00.1_download_structures_run.sh 
```

```console
bash 00.2_check_structures_run.sh
```

```console
python 01.0_chain_gold_labels.py
```

```console
python 01.1_overview_per_pdb.py 
```

```console
python 01.2_gold_audit_summary.py 
```

```console
python 01y.probe_excel_only_pdbs.py 
```

```console
bash 02.0_run_p2rank_batch.sh
```

```console
python 02.1_chain_p2rank_labels.py
```

```console
python 02.2_build_p2rank_rank1_map.py
```

```console
python 03.0_plm_predict_from_gold.py
```

```console
python 03.1_chain_plm_labels.py
```

```console
python 03y.inspect_gold_vs_plm_keys.py 
```

```console
python 04.0_evaluate_pipeline_all_in_one.py
```

```console
bash 05.report_from_evaljson_plus_curves_run.sh
```

```console
python 06.0_prevalence.py
```

## 3. Leakage Control and Evaluation Policy

Sequence similarity control was applied prior to evaluation. Structures with >30% sequence identity to the PLM training set were excluded.
All evaluations reported in the manuscript are performed at the structure–chain level using the original KinCoRe annotations. 
For controlled analyses, only one structure per protein was retained
(KinCoRe-UNQ export).

## 4. Reproducibility

To reproduce the reported results:

Use the provided dataset folders with the original input Excel tables.

Run scripts in numerical order (01 → 06).

Metrics and figures reported in the manuscript correspond directly to: 

- 04.eval_all_in_one.json

- figures generated by 05.report_from_evaljsons_plus_curves.py

Structure files are not archived to avoid data duplication;
all required structures can be retrieved programmatically using the provided scripts.

## 5. Notes

AUROC and AUPR are computed from raw probability scores and are independent of threshold selection.

Fixed-threshold and best-MCC operating points are described in the manuscript Methods section.

Class prevalence reflects residue-level annotations aggregated across all protein–ligand complexes.
Prevalence for unique proteins is reported only for controlled analyses.
