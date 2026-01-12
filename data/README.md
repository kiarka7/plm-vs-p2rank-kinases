# Data Availability

Due to size limitations of GitHub, not outputs used in this study are stored directly in this repository.

---
### Datasets
The datasets, including:
- residue-level GOLD annotations,
- processed P2Rank outputs,
- PLM raw predictions,
- evaluation JSON files,
- and all generated figures, reports

are available via institutional data storage at:

**[[UNIVERSITY SHAREPOINT LINK HERE](https://cunicz-my.sharepoint.com/:f:/g/personal/40001559_cuni_cz/IgBFFuT8WMGHQrGYNnq4ncp9ATPYrBBmPGa3gdEVuP57D50?e=QkmkY8)]**

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

The only external input required to regenerate the datasets is the
**original Excel table** listing protein–ligand structures.
All subsequent processing steps (filtering, prediction mapping,
evaluation, and reporting) are fully scripted and reproducible.

