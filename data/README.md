# Data Availability

Due to size limitations of GitHub, no outputs used in this study are stored directly in this repository.

---
### Datasets
The datasets, including:

- residue-level GOLD annotations,
- processed P2Rank outputs,
- PLM raw predictions,
- evaluation JSON files,
- and all generated figures, reports

are available via institutional data storage at:

**[[Charles University storage](https://cunicz-my.sharepoint.com/:f:/g/personal/40001559_cuni_cz/IgBFFuT8WMGHQrGYNnq4ncp9ATPYrBBmPGa3gdEVuP57D50?e=QkmkY8)]**

(access provided for review and archival purposes)

### Dataset variants

Two dataset variants are used throughout the analysis:

- **`KinCoRe/`**
  The full dataset, containing all available protein–ligand complexes
  for each kinase type.

- **`KinCoRe-UNQ/`**
  A reduced dataset in which only **one representative structure per protein**
  (UniProt ID) is retained (Supporting Information tables).

Both dataset variants share the same internal structure and differ only
in the set of included structures.

### Scripts and reproducibility

All scripts provided in the `scripts/` directory operate directly on
the dataset folders (`Kinase_Type_I/`, `Kinase_Type_I.5/`,...).

The only external inputs required to regenerate the datasets are the
original Excel tables availabel in each of the `Kinase_*` directories. These tables contain the protein–ligand information.

The directory `_combined_final_CLEAN/`  (in the sharepoint link) contains aggregated figures and
reports generated from the evaluation outputs.

Atomic protein structures (.cif files) are not included in the sharepoint link. 
All structures are publicly available from the RCSB Protein Data Bank and can be automatically downloaded using the provided script:

```console
bash 00.1_download_structures_run.sh
```

Diagnostic and inspection scripts are documented in the GitHub repository. Their outputs are included in the sharepoint link when relevant for evaluation or reproducibility; additional diagnostic or inspection outputs used only for sanity checks are not included.
