#!/bin/bash

python CLEAN_00.2_build_clean_excel_from_mask_V2.py \
  --train vita-model-modificated-finetuning/train.txt \
  --excel "Kinase_Type_ALLO/Allostreric Kinase LIgands Only.xlsx" \
  --out-excel "Kinase_Type_ALLO/CLEAN/Allostreric Kinase LIgands Only_CLEAN.xlsx" \
  --dataset-name Kinase_Type_ALLO \
  --uniprot-col UniprotID \
  --uniprot-map-tsv Kinase_Type_ALLO/kinase_uniprot_map_typeALLO.tsv \
  --map-entry-col dataset_uniprot_id \
  --map-acc-col accession
