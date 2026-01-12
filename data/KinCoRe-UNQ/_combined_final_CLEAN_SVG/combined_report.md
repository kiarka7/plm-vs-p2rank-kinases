# Combined report

_Generated: 2026-01-08T20:23:51.855974+00:00Z_

## Micro bars (from eval JSON)

![acc](bars_micro_acc.png)

![f1](bars_micro_f1.png)

![mcc](bars_micro_mcc.png)

## Macro bars (mean ± IQR/2, from eval JSON)

![acc](bars_macro_acc_mean.png)

![f1](bars_macro_f1_mean.png)

![mcc](bars_macro_mcc_mean.png)

## PLM ROC curves (micro) — from raw scores

![roc](plm_roc_overlay.png)

## PLM PR curves (micro) — from raw scores

![pr](plm_pr_overlay.png)

## AUROC (from raw scores)

| dataset          |   auroc_from_scores |
|:-----------------|--------------------:|
| Kinase_Type_ALLO |            0.68599  |
| Kinase_Type_I    |            0.96649  |
| Kinase_Type_I.5  |            0.971667 |
| Kinase_Type_II   |            0.955231 |
| Kinase_Type_III  |            0.871374 |

## AUPR (from raw scores)

| dataset          |   aupr_from_scores |
|:-----------------|-------------------:|
| Kinase_Type_ALLO |          0.0721546 |
| Kinase_Type_I    |          0.619043  |
| Kinase_Type_I.5  |          0.726287  |
| Kinase_Type_II   |          0.688647  |
| Kinase_Type_III  |          0.275881  |

## Notes

```
[Kinase_Type_ALLO] P2 ANY fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I] P2 ANY fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I.5] P2 ANY fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_II] P2 ANY fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_III] P2 ANY fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_ALLO] PLM fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I] PLM fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I.5] PLM fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_II] PLM fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_III] PLM fixed: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_ALLO] PLM bestMCC: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I] PLM bestMCC: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_I.5] PLM bestMCC: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_II] PLM bestMCC: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_III] PLM bestMCC: weighted F1 not found -> using 'f1' (fallback).
[Kinase_Type_ALLO] P2 ANY fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I] P2 ANY fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I.5] P2 ANY fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_II] P2 ANY fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_III] P2 ANY fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_ALLO] PLM fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I] PLM fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I.5] PLM fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_II] PLM fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_III] PLM fixed: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_ALLO] PLM bestMCC: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I] PLM bestMCC: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_I.5] PLM bestMCC: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_II] PLM bestMCC: weighted F1 (macro) not found -> using 'f1' (fallback).
[Kinase_Type_III] PLM bestMCC: weighted F1 (macro) not found -> using 'f1' (fallback).
```