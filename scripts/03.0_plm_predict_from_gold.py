#!/usr/bin/env python3
"""
03.0_plm_predict_from_gold.py  (robust, GOLD-driven inference, with checkpoints + resume)

What this does
--------------
- Loads a finetuned PLM and a tokenizer.
- Uses GOLD list (01.chain_gold_labels.json) to know exactly which (PDB, CHAIN) to predict.
- Extracts chain sequences from local .cif files (author residue numbering).
- Runs the model and writes one JSON with BOTH:
    * "probabilities": full per-residue probs (length L; alignment-friendly)
    * "residues_predicted": sparse list above threshold (human-friendly)
- Keeps output order GOLD-consistent.

Robust bits
-----------
- Case-insensitive chain resolution (handles 'A' in GOLD vs 'a' in mmCIF).
- No dependency on Bio.PDB.Polypeptide.three_to_one (uses a safe map with PTMs).
- Handles tokenizer special tokens; trims logits to length L when needed.
- Quiet logging: prints only problems, periodic heartbeat, and a final summary.
- Periodic checkpoints every N chains to OUTPUT_JSON.tmp
- Resume from OUTPUT_JSON.tmp after crash/reboot (skips already processed keys)
- Atomic writes to avoid corrupted checkpoints

Outputs
-------
- ./****/03.plm_predictions.json         (final)
- ./****/03.plm_predictions.json.tmp     (checkpoint)
"""

import os, sys, json
import torch
import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from Bio.PDB import MMCIFParser, is_aa
from Bio.Data import IUPACData
from transformers import AutoTokenizer

# ================== CONFIG (edit these to taste) ==================

BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"

STRUCT_DIR      = os.path.join(BASE_DIR, KINASE_FOLDER, "structures")
GOLD_JSON       = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
OUTPUT_JSON     = os.path.join(BASE_DIR, KINASE_FOLDER, "03.plm_predictions.json")
CHECKPOINT_JSON = OUTPUT_JSON + ".tmp"   # rolling checkpoint
CHECKPOINT_EVERY = 10                    # save checkpoint every N processed chains
HEARTBEAT_EVERY  = 10                    # print a progress line every N chains
RESUME_FROM_CHECKPOINT = True            # resume if .tmp exists

# Finetuning folder (contains model.pt and e.g. finetuning_utils.py)
DEFAULT_FINETUNE_DIR = os.path.join(BASE_DIR, "vita-model-modificated-finetuning")
DEFAULT_MODEL_PATH   = os.path.join(DEFAULT_FINETUNE_DIR, "model.pt")

# HF tokenizer name (ESM2 650M by default)
TOKENIZER_NAME = "facebook/esm2_t33_650M_UR50D"

# Threshold for the auxiliary sparse output "residues_predicted"
SPARSE_THRESHOLD = 0.75

# Device: CUDA if available
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# =================================================================

# ---- Safe 3-letter → 1-letter mapping (handles common PTMs) ----
AA3_TO_1: Dict[str, str] = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}
AA3_TO_1.update({
    "MSE": "M",  # Selenomethionine
    "SEC": "U",  # Selenocysteine
    "PYL": "O",  # Pyrrolysine
    "SEP": "S",  # Phosphoserine
    "TPO": "T",  # Phosphothreonine
    "PTR": "Y",  # Phosphotyrosine
    "HYP": "P",  # Hydroxyproline
    "CSO": "C", "CSD": "C", "CME": "C",  # oxidized/methylated cysteines
    "M3L": "K", "MLY": "K", "MLZ": "K",  # methyl-lysines
})

def three_to_one_safe(resname: str) -> str:
    """Return 1-letter AA for a 3-letter residue name. Unknown → 'X'."""
    return AA3_TO_1.get(str(resname).strip().upper(), "X")


# ------------------------- I/O helpers ---------------------------

def ensure_finetune_imports(finetune_dir: str) -> None:
    """
    If the pickled model needs auxiliary modules (e.g., 'finetuning_utils'),
    add their directory to sys.path BEFORE torch.load.
    """
    if finetune_dir and os.path.isdir(finetune_dir):
        if finetune_dir not in sys.path:
            sys.path.insert(0, finetune_dir)
        print(f"[i] Added to sys.path for pickled deps: {finetune_dir}")

def load_gold_keys(path: str) -> List[Tuple[str, str]]:
    """
    Read GOLD JSON and return a GOLD-ordered list of (PDB, CHAIN) keys.
    (First occurrence per key wins; order is preserved.)
    """
    data = json.load(open(path))
    keys, seen = [], set()
    for e in data:
        key = (str(e["pdb_id"]).upper(), str(e["chain_id"]).upper())
        if key not in seen:
            seen.add(key)
            keys.append(key)
    return keys


# -------- NEW: robust, case-insensitive chain resolution ----------

def resolve_chain_id_case_insensitive(model, requested_id: str) -> Optional[str]:
    """
    Return a chain id present in `model` that matches `requested_id` case-insensitively.

    Preference order:
      1) Exact match (fast path).
      2) Unique case-insensitive match (e.g., GOLD 'A' vs mmCIF 'a').
      3) If single-letter requested id, try swapped case ('A' -> 'a', 'a' -> 'A').

    If ambiguous or not found, return None.
    """
    req = str(requested_id)
    available = [ch.id for ch in model]
    # 1) exact
    if req in available:
        return req
    # 2) unique CI match
    matches = [cid for cid in available if cid.lower() == req.lower()]
    if len(matches) == 1:
        return matches[0]
    # 3) trivial swap for single-letter ids
    if len(req) == 1:
        alt = req.lower() if req.isupper() else req.upper()
        if alt in available:
            return alt
    return None


def extract_chain_residues(cif_path: str, chain_id: str) -> Optional[Tuple[str, List[int]]]:
    """
    Build a sequence from residues of a chain and collect author residue numbers.
    Returns (sequence, author_positions) or None if chain not found / empty.

    Uses resolve_chain_id_case_insensitive() to be robust to mmCIF lower/upper ids.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    model = structure[0]

    chosen = resolve_chain_id_case_insensitive(model, chain_id)
    if chosen is None:
        # chain not present; print info about available chains
        avail = ",".join(sorted([ch.id for ch in model]))
        print(f"[!] Chain '{chain_id}' not in {os.path.basename(cif_path)}; available: {avail}")
        return None

    chain = model[chosen]
    letters, positions = [], []
    for res in chain:
        if not is_aa(res, standard=False):
            continue
        aa = three_to_one_safe(res.get_resname())
        letters.append(aa)
        positions.append(int(res.get_id()[1]))  # author residue number

    if not letters:
        print(f"[!] Empty/unsuitable sequence for {os.path.basename(cif_path)} chain {chosen}")
        return None
    return ("".join(letters), positions)


# ---------------- small helpers ----------------

def strip_special_tokens(probs: np.ndarray, L: int) -> np.ndarray:
    """
    Align token-level probabilities to residue-level length L by trimming/padding.
    Heuristics for HF tokenizers (CLS/EOS):
      len==L      -> use as-is
      len==L+2    -> drop first & last
      len==L+1    -> drop first
      len> L+2    -> take next L after dropping first (fallback)
      len< L      -> pad zeros to L
    """
    T = int(probs.shape[0])
    if T == L:
        return probs
    if T == L + 2:
        return probs[1:-1]
    if T == L + 1:
        return probs[1:]
    if T > L + 2:
        return probs[1:1+L] if (T - 1) >= L else probs[:L]
    out = np.zeros(L, dtype=float)
    out[:min(T, L)] = probs[:min(T, L)]
    return out

def atomic_write_json(path: str, obj) -> None:
    """
    Write JSON atomically: write to path+'.writing' and then os.replace().
    Prevents corrupted files if process is interrupted mid-write.
    """
    tmp = path + ".writing"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)


# ----------------------- Inference routine -----------------------

def main():
    finetune_dir = DEFAULT_FINETUNE_DIR
    model_path   = DEFAULT_MODEL_PATH

    # Make sure pickled aux modules are importable
    ensure_finetune_imports(finetune_dir)

    # Load model (entire object, not just state_dict)
    print("Loading model…")
    if not os.path.isfile(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    model = torch.load(model_path, map_location=DEVICE, weights_only=False)
    model.to(DEVICE)
    model.eval()
    print(f"[i] Model loaded from: {model_path}")

    # Tokenizer
    print("Loading tokenizer…")
    tokenizer = AutoTokenizer.from_pretrained(TOKENIZER_NAME)
    print(f"[i] Tokenizer: {TOKENIZER_NAME}")

    # GOLD keys
    keys = load_gold_keys(GOLD_JSON)
    total = len(keys)
    print(f"Generating PLM predictions for {total} GOLD chains…")

    # Prepare output dir
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)

    # ---- Resume logic: load partial checkpoint if present ----
    out_records: List[Dict] = []
    processed: Set[Tuple[str, str]] = set()
    if RESUME_FROM_CHECKPOINT and os.path.isfile(CHECKPOINT_JSON):
        try:
            partial = json.load(open(CHECKPOINT_JSON))
            if isinstance(partial, list):
                out_records = partial
                for r in out_records:
                    processed.add((str(r.get("pdb_id","")).upper(), str(r.get("chain_id","")).upper()))
                print(f"[i] Resuming from checkpoint: {CHECKPOINT_JSON} ({len(out_records)} chains already)")
        except Exception as e:
            print(f"[!] Failed to read checkpoint ({CHECKPOINT_JSON}): {e} — starting fresh")

    # Counters
    n_ok = 0
    n_missing_cif = 0
    n_missing_chain = 0
    n_empty = 0
    n_truncated = 0

    # Main loop
    for i, (pdb_id, chain_id) in enumerate(keys, 1):
        if (pdb_id, chain_id) in processed:
            # already in checkpoint — skip computing, keep order as-is
            continue

        cif = os.path.join(STRUCT_DIR, f"{pdb_id.lower()}.cif")
        if not os.path.isfile(cif):
            if n_missing_cif < 12:
                print(f"[!] Missing CIF for {pdb_id} (expected {cif})")
            n_missing_cif += 1
            continue

        got = extract_chain_residues(cif, chain_id)
        if got is None:
            n_missing_chain += 1
            continue
        seq, author_pos = got
        L = len(seq)
        if L == 0:
            n_empty += 1
            continue

        # Tokenize; ESM has 1024 limit — allow truncation with warning counter
        enc = tokenizer(
            seq,
            return_tensors="pt",
            add_special_tokens=True,
            truncation=True,
            max_length=1024
        ).to(DEVICE)

        tok_len = int(enc["input_ids"].shape[1])
        if tok_len > (L + 2):
            n_truncated += 1  # indicates some reshape/truncation happened

        # Optional plDDT channel if the finetuned head expects it
        seq_len_tokens = enc["input_ids"].shape[1]
        plddt = torch.ones(seq_len_tokens, dtype=torch.float16, device=DEVICE).unsqueeze(0).unsqueeze(1)
        enc["plDDTs"] = plddt  # keep the key name as in your finetuning

        # Forward pass
        try:
            with torch.no_grad():
                out = model(enc)
        except Exception as e:
            print(f"[!] Forward error for {pdb_id}.{chain_id}: {e} — skipping")
            continue

        # Accept either pure tensor or a tuple/list with first tensor
        if isinstance(out, torch.Tensor):
            logits = out
        elif isinstance(out, (tuple, list)) and len(out) >= 1 and isinstance(out[0], torch.Tensor):
            logits = out[0]
        else:
            print(f"[!] Unexpected model output for {pdb_id}.{chain_id}: {type(out)} — skipping")
            continue

        probs_token = torch.sigmoid(logits).squeeze().detach().float().cpu().numpy()
        probs_res = strip_special_tokens(probs_token, L)  # align to residue count
        if probs_res.shape[0] != L:
            # final safety
            probs_res = probs_res[:L] if probs_res.shape[0] > L else np.pad(probs_res, (0, L - probs_res.shape[0]))

        # Dense per-residue probs (for evaluation)
        probs_list = [float(x) for x in probs_res]

        # Sparse (thresholded) list for quick inspection (1-based positions)
        residues = [
            {"position": int(ii + 1), "probability": float(p)}
            for ii, p in enumerate(probs_res) if p >= SPARSE_THRESHOLD
        ]

        out_records.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "sequence_residue_numbers": author_pos, # same author ID asi in GOLD
            "probabilities": probs_list,         # dense, length L
            "residues_predicted": residues       # sparse, human-friendly
        })
        processed.add((pdb_id, chain_id))
        n_ok += 1

        # Heartbeat
        if (len(processed) % HEARTBEAT_EVERY) == 0:
            print(f"[i] Processed {len(processed)}/{total} chains…", flush=True)

        # Periodic checkpoint (atomic write)
        if (len(processed) % CHECKPOINT_EVERY) == 0:
            try:
                atomic_write_json(CHECKPOINT_JSON, out_records)
                print(f"[ckpt] Wrote checkpoint: {CHECKPOINT_JSON} ({len(processed)}/{total})", flush=True)
            except Exception as e:
                print(f"[!] Failed to write checkpoint: {e}")

    # Final save (atomic)
    try:
        atomic_write_json(OUTPUT_JSON, out_records)
        print(f"Done. Saved {OUTPUT_JSON} with {len(out_records)} chains.")
    finally:
        # Optionally keep the checkpoint; or clean it up:
        # try: os.remove(CHECKPOINT_JSON)
        # except OSError: pass
        pass

    # Summary
    if n_missing_cif:
        print(f"  - Missing CIF files: {n_missing_cif}")
    if n_missing_chain:
        print(f"  - Chains not found / empty: {n_missing_chain}")
    if n_empty:
        print(f"  - Empty sequences: {n_empty}")
    if n_truncated:
        print(f"  - Tokenized length > L+2 (likely trunc/reshape): {n_truncated}")
    print(f"  - OK: {n_ok}")

if __name__ == "__main__":
    main()
