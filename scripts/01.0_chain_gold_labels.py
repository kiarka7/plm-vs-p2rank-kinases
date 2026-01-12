#!/usr/bin/env python3
"""
01.chain_gold_labels.py — robust gold-label extraction with audit report

What this script does:
- Reads an Excel sheet with columns: PDB (e.g. "4CZT" or "4CZT A"), Ligand (tokens like "ANP:60401")
- Loads mmCIF using author numbering (use_author=True).
- For each token NAME:NUMBER, classifies the component (cap / modified AA / polymeric / water / ligand candidate).
- Skips caps + modified AA by default; optionally skip polymeric components.
- Tries to locate the residue(s) by:
    1) NAME + AUTHOR number
    2) AUTHOR number-only (any non-AA, excluding water and optionally modified AA)
    3) NAME + LABEL number   (map LABEL→AUTHOR via _atom_site)
    4) LABEL number-only     (map LABEL→AUTHOR)
- Fallback-only: If Excel uses "code:xx:yy" like "337:01:00", we extract the leading 3-char code (e.g., "337")
                and search **by name only** across the requested chain — but **only if** no classic token
                for the same (pdb, chain, code) was used in that row.
- Numeric fallback: If the fallback token is purely numeric (e.g. "801") and there is exactly one non-AA HET
                at that author number in the chain, accept it (treated as NAME-ONLY numeric fallback).
- Manual overrides: Optional CSV `01.gold_nameonly_overrides.csv` can provide extra name-only codes per (PDB,CHAIN).
                These are merged with Excel fallback codes and audited as `name_only_override`.
- Labels all amino-acid residues in the chain within DIST_CUTOFF Å of any matched residue atom as 1, else 0.
- Writes:
    - 01.chain_gold_labels.json  (one entry per USED token)
    - 01.ligand_audit.csv        (one row per Excel/override token with diagnostics)

Notes:
- You may see multiple JSON entries for the same (pdb,chain) if Excel lists multiple ligands (different codes/numbers).
  Downstream evaluation OR-merges them; here we just avoid duplicate entries for the same code coming from the fallback.
"""

import os, json, re, csv
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set
from Bio.PDB import MMCIFParser, NeighborSearch, is_aa
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ========= CONFIG =========
BASE_DIR      = "."
KINASE_FOLDER = "Kinase_Type_I"  # <- adjust to your folder
EXCEL_FILE    = os.path.join(BASE_DIR, KINASE_FOLDER, "Kinase_Ligands_Type I.xlsx")
STRUCT_DIR    = os.path.join(BASE_DIR, KINASE_FOLDER, "structures")
OUTPUT_JSON   = os.path.join(BASE_DIR, KINASE_FOLDER, "01.chain_gold_labels.json")
AUDIT_CSV     = os.path.join(BASE_DIR, KINASE_FOLDER, "01.ligand_audit.csv")

# Optional: manual name-only overrides for (PDB, CHAIN) -> list of comp codes
OVERRIDES_CSV = os.path.join(BASE_DIR, KINASE_FOLDER, "01.gold_nameonly_overrides.csv")

DIST_CUTOFF   = 4.0  # Å
EXCLUDE_HYDROGENS = True

# Skip these outright as ligands (terminal caps etc.)
CAP_RES_CODES = {"ACE","NME","NH2","FOR","AME","AIB","PCA"}

# Modified amino acids you usually don't want to treat as small-molecule ligands:
MODIFIED_AA_CODES = {
    "YTH","SEP","TPO","PTR","MSE","SEC","CSO","CME","CSD","MLY","M3L","MLZ",
    "HYP","LLP","TYS","KCX","HSK","FME"
}
EXCLUDE_MOD_AA = True

# Water-like components to ignore as ligands
WATER_CODES = {"HOH","WAT","DOD"}

# Optionally skip any component that appears in that chain as POLYMER (group_PDB='ATOM')
EXCLUDE_POLYMERIC_COMPONENTS = False

# Ligand name aliases if mmCIF uses a different component ID than Excel/PDB
LIGAND_ALIASES: Dict[str,str] = {
    # "YTH": "FJ9",
}

# Diagnostics
PRINT_DIAGNOSTICS = True
DIAG_LIST_LIMIT   = 20

# Accept tokens like "ANP:60401"
LIGAND_TOKEN = re.compile(r'^[A-Za-z0-9_]+:\d+$')
# Accept leading 3-char code from messy cells (e.g., "337:01:00" → "337")
LEAD3  = re.compile(r"\b([A-Z0-9]{3})(?=[:;,\s/]|$)", re.I)
TOKEN3 = re.compile(r"^[A-Z0-9]{3}$", re.I)
# ==========================


def iter_ligand_tokens(cell: str):
    """Yield NAME:NUMBER tokens from the Ligand cell (handles commas/semicolons/whitespace)."""
    if not cell:
        return
    for raw in str(cell).replace(";", " ").replace(",", " ").split():
        tok = raw.strip()
        if LIGAND_TOKEN.match(tok):
            yield tok

def extract_nameonly_codes(cell: str) -> List[str]:
    """
    Extract leading 3-char codes like '337' or 'ATP' from messy cells (e.g., '337:01:00').
    Returned codes are unique and preserve order. Skip explicit waters.
    """
    if not cell:
        return []
    s = str(cell)
    seen, out = set(), []
    for chunk in re.split(r"[;/|]+", s):
        m = LEAD3.search(chunk.strip())
        if not m:
            continue
        code = m.group(1).upper()
        if TOKEN3.match(code) and code != "HOH" and code not in seen:
            seen.add(code)
            out.append(code)
    return out

def make_parser_author():
    """MMCIF parser preferring author numbering."""
    try:
        return MMCIFParser(QUIET=True, use_author=True)
    except TypeError:
        return MMCIFParser(QUIET=True)

def canon_name(name: str) -> str:
    """Apply alias mapping (uppercased)."""
    return LIGAND_ALIASES.get(name, name).upper()

def find_chain(struct, chain_id: str):
    """Find chain by exact, upper or lower ID."""
    if chain_id in struct:
        return struct[chain_id]
    up, lo = chain_id.upper(), chain_id.lower()
    if up in struct: return struct[up]
    if lo in struct: return struct[lo]
    return None

def gen_variants(resnum: int):
    """Tolerate long numbers: include last 3–5 digits."""
    s = str(resnum)
    return {resnum} | {int(s[-i:]) for i in (3,4,5) if len(s) >= i}

# ---------- label↔author mapping ----------
def load_label_author_map(cif_path: str):
    """
    Build per-auth-chain mapping using _atom_site.
    Returns: dict auth_chain_id -> mapping dicts (label2auth, auth2label, comp_at, group_at)
    """
    mm = MMCIF2Dict(cif_path)
    label_asym = list(mm.get("_atom_site.label_asym_id", []))
    auth_asym  = list(mm.get("_atom_site.auth_asym_id", []))
    label_seq  = list(mm.get("_atom_site.label_seq_id", []))
    auth_seq   = list(mm.get("_atom_site.auth_seq_id", []))
    comp_id    = list(mm.get("_atom_site.label_comp_id", []))
    group_pdb  = list(mm.get("_atom_site.group_PDB", []))
    icode      = list(mm.get("_atom_site.pdbx_PDB_ins_code", []))

    n = min(len(label_asym), len(auth_asym), len(label_seq), len(auth_seq), len(comp_id), len(group_pdb))
    use_icode = bool(icode) and (len(icode) >= n)

    per_chain = defaultdict(lambda: {
        "label2auth": {},
        "auth2label": {},
        "comp_at": defaultdict(set),
        "group_at": defaultdict(set)
    })

    for i in range(n):
        la = str(label_asym[i])
        aa = str(auth_asym[i])
        ls = label_seq[i]; as_ = auth_seq[i]
        try:
            ls_i = int(ls); as_i = int(as_)
        except Exception:
            continue
        ic = (icode[i].strip() if (use_icode and icode[i] not in (".","?")) else "")
        comp = str(comp_id[i]).upper()
        grp  = str(group_pdb[i]).upper()

        entry = per_chain[aa]
        entry["label2auth"].setdefault(ls_i, (as_i, ic))
        entry["auth2label"].setdefault(as_i, ls_i)
        entry["comp_at"][ls_i].add(comp)
        entry["group_at"][ls_i].add(grp)

    return per_chain

def is_polymeric_component_in_chain(label_map, auth_chain_id: str, comp_name_up: str) -> bool:
    """True if comp_name shows up at some label_seq_id with group_PDB='ATOM' in this chain."""
    entry = label_map.get(auth_chain_id, None)
    if not entry:
        return False
    for ls, comps in entry["comp_at"].items():
        if comp_name_up in comps and "ATOM" in entry["group_at"].get(ls, set()):
            return True
    return False
# -----------------------------------------

# ---------- search helpers ----------
def residues_by_author_name(chain, comp_name_up: str, author_number: int):
    """Find NON-AA residues at author_seq_id matching name."""
    variants = gen_variants(author_number)
    hits = []
    for r in chain:
        if is_aa(r):
            continue
        het, num, ic = r.get_id()
        if int(num) in variants and r.get_resname().upper() == comp_name_up:
            hits.append(r)
    return hits

def residues_by_author_number_only(chain, author_number: int):
    """Find NON-AA residues at author_seq_id, any name (excluding waters and modAA)."""
    variants = gen_variants(author_number)
    hits = []
    for r in chain:
        if is_aa(r):
            continue
        het, num, ic = r.get_id()
        if int(num) in variants:
            rn = r.get_resname().upper()
            if rn in WATER_CODES:
                continue
            if EXCLUDE_MOD_AA and rn in MODIFIED_AA_CODES:
                continue
            hits.append(r)
    return hits

def residues_by_label_number(struct, label_map, auth_chain_id: str, label_number: int,
                             comp_name_up: Optional[str] = None):
    """Using label_map: find residues by label_seq_id (optionally require name)."""
    chain = find_chain(struct, auth_chain_id)
    if chain is None: return []

    entry = label_map.get(auth_chain_id, None)
    if not entry: return []
    if label_number not in entry["label2auth"]: return []

    auth_num, _ = entry["label2auth"][label_number]
    variants = gen_variants(auth_num)
    hits = []
    for r in chain:
        if is_aa(r): continue
        het, num, ic = r.get_id()
        if int(num) in variants:
            rn = r.get_resname().upper()
            if rn in WATER_CODES:
                continue
            if EXCLUDE_MOD_AA and rn in MODIFIED_AA_CODES:
                continue
            if comp_name_up is None or rn == comp_name_up:
                hits.append(r)
    return hits

def residues_by_name_any_number(chain, comp_name_up: str):
    """
    Find NON-AA residues in the chain that match a given HET name, irrespective of number.
    Used by the name-only fallback.
    """
    hits = []
    for r in chain:
        if is_aa(r):
            continue
        rn = r.get_resname().upper()
        if rn == comp_name_up:
            if rn in WATER_CODES:
                continue
            if EXCLUDE_MOD_AA and rn in MODIFIED_AA_CODES:
                continue
            hits.append(r)
    return hits
# ------------------------------------

def classify_component(comp_up: str, label_map, auth_chain_id: str) -> str:
    """Heuristic classification for audit."""
    if comp_up in CAP_RES_CODES:      return "cap"
    if comp_up in MODIFIED_AA_CODES:  return "modified_aa"
    if comp_up in WATER_CODES:        return "water"
    if is_polymeric_component_in_chain(label_map, auth_chain_id, comp_up):
        return "polymeric"
    return "ligand_candidate"

# ---------- overrides loader ----------
def load_nameonly_overrides(path: str):
    """
    Read a CSV with columns: pdb_id, chain_id, code
    Returns: dict[(PDB, CHAIN)] -> set([CODE,...]) ; all uppercased and 4-char PDB.
    """
    out = defaultdict(set)
    if not os.path.isfile(path):
        return out
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            pdb  = str(row.get("pdb_id","")).strip().upper()[:4]
            ch   = (str(row.get("chain_id","")).strip() or "A").upper()
            code = str(row.get("code","")).strip().upper()
            if pdb and ch and code:
                out[(pdb, ch)].add(code)
    return out
# --------------------------------------

def main():
    df = pd.read_excel(EXCEL_FILE)
    parser = make_parser_author()
    out_entries = []
    audit_rows = []

    # Load optional name-only overrides from CSV
    overrides = load_nameonly_overrides(OVERRIDES_CSV)
    if PRINT_DIAGNOSTICS and overrides:
        total_codes = sum(len(s) for s in overrides.values())
        print(f"[i] Loaded name-only overrides: {total_codes} code(s) for {len(overrides)} PDB.chain")

    # --- de-dup state ---
    seen_classic: Set[Tuple[str,str,str,int]] = set()  # (pdb, chain, comp_up, number)
    seen_nameonly: Set[Tuple[str,str,str]]    = set()  # (pdb, chain, comp_up)

    # counters for a short summary
    n_used_classic = n_used_nameonly = 0
    n_skipped_nameonly_due_to_classic = 0

    for _, row in df.iterrows():
        pdb_full = str(row.get("PDB", "")).strip()
        if not pdb_full or pdb_full.lower() == "nan":
            continue

        pdb = pdb_full[:4].lower()

        # Mírnější varianta: pokud není explicitní chain, řádek přeskočíme
        if len(pdb_full) == 4:
            if PRINT_DIAGNOSTICS:
                print(f"[!] Missing chain for {pdb_full}, skipping row.")
            continue

        chain_in = pdb_full[4:].strip()
        if not chain_in:
            if PRINT_DIAGNOSTICS:
                print(f"[!] Empty chain after parsing PDB field '{pdb_full}', skipping row.")
            continue
        chain_in = chain_in.upper()

        lig_str = str(row.get("Ligand", "")).strip()
        if not lig_str or lig_str.lower() in ("no_ligand", "none", "nan"):
            # Still allow overrides to contribute if Excel ligand cell is empty
            if (pdb.upper(), chain_in) not in overrides:
                continue

        cif = os.path.join(STRUCT_DIR, f"{pdb}.cif")
        if not os.path.isfile(cif):
            if PRINT_DIAGNOSTICS:
                print(f"[!] Missing file: {cif}")
            # audit rows (classic + name-only + overrides) as file_missing
            for tok in iter_ligand_tokens(lig_str):
                name, num_str = tok.split(":", 1)
                try: num = int(num_str)
                except: continue
                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": name, "token_number": num,
                    "alias_used": canon_name(name), "class": "unknown",
                    "status": "file_missing", "found_mode": "",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
            for code in extract_nameonly_codes(lig_str):
                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": code, "token_number": "",
                    "alias_used": canon_name(code), "class": "unknown",
                    "status": "file_missing", "found_mode": "",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
            # also audit overrides for this key as file_missing
            for code in overrides.get((pdb.upper(), chain_in), set()):
                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": code, "token_number": "",
                    "alias_used": canon_name(code), "class": "unknown",
                    "status": "file_missing", "found_mode": "name_only_override",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
            continue

        try:
            struct = parser.get_structure(pdb, cif)[0]
        except Exception as e:
            if PRINT_DIAGNOSTICS:
                print(f"[!] Failed to parse {cif}: {e}")
            continue

        ch = find_chain(struct, chain_in)
        if ch is None:
            if PRINT_DIAGNOSTICS:
                print(f"[!] Chain {chain_in} not found in {pdb.upper()}")
            continue

        # Build mapping once per structure
        label_map = load_label_author_map(cif)

        # ---- (A) CLASSIC tokens on this row ----
        classic_codes_on_row: Set[str] = set()  # track comp_up present as classic for fallback gating

        for tok in iter_ligand_tokens(lig_str):
            name, num_str = tok.split(":", 1)
            try:
                num = int(num_str)
            except ValueError:
                continue

            comp_up = canon_name(name)
            classic_codes_on_row.add(comp_up)

            key_classic = (pdb, chain_in, comp_up, num)
            if key_classic in seen_classic:
                continue  # global de-dup across the whole file
            seen_classic.add(key_classic)

            comp_class = classify_component(comp_up, label_map, chain_in)
            if comp_class in ("cap", "modified_aa", "polymeric"):
                if PRINT_DIAGNOSTICS:
                    why = "cap" if comp_class=="cap" else ("modAA" if comp_class=="modified_aa" else "polymeric")
                    print(f"[i] Skipping {why} component {name}:{num} in {pdb.upper()}.{chain_in}")
                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": name, "token_number": num,
                    "alias_used": comp_up, "class": comp_class,
                    "status": "skipped", "found_mode": f"skip_{comp_class}",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
                continue

            # --- classic search pipeline ---
            found_mode = ""
            lig_res = residues_by_author_name(ch, comp_up, num)
            if lig_res:
                found_mode = "author_name_number"
            else:
                lig_res = residues_by_author_number_only(ch, num)
                if lig_res:
                    found_mode = "author_number_only"
                    if PRINT_DIAGNOSTICS:
                        names = sorted({r.get_resname() for r in lig_res})
                        print(f"    ↳ Using AUTHOR number-only for {pdb.upper()}.{chain_in} {name}:{num} "
                              f"(found: {', '.join(names)})")

            if not lig_res:
                lig_res = residues_by_label_number(struct, label_map, chain_in, num, comp_name_up=comp_up)
                if lig_res:
                    found_mode = "label_name_number"
                    auth_num, _ = label_map.get(chain_in, {}).get("label2auth", {}).get(num, (None, ""))
                    if PRINT_DIAGNOSTICS and auth_num is not None:
                        print(f"    ↳ Resolved via LABEL name+number for {pdb.upper()}.{chain_in} "
                              f"{name}:{num} (label {num} → author {auth_num})")

            if not lig_res:
                lig_res = residues_by_label_number(struct, label_map, chain_in, num, comp_name_up=None)
                if lig_res:
                    found_mode = "label_number_only"
                    auth_num, _ = label_map.get(chain_in, {}).get("label2auth", {}).get(num, (None, ""))
                    if PRINT_DIAGNOSTICS and auth_num is not None:
                        names = sorted({r.get_resname() for r in lig_res})
                        print(f"    ↳ Using LABEL number-only for {pdb.upper()}.{chain_in} {name}:{num} "
                              f"(label {num} → author {auth_num}; names: {', '.join(names)})")

            if not lig_res:
                # Not found
                if PRINT_DIAGNOSTICS:
                    hets = []
                    for r in ch:
                        if not is_aa(r):
                            het, n, ic = r.get_id()
                            hets.append(f"{r.get_resname()}:{int(n)}{(ic or '').strip()}")
                    preview = ", ".join(hets[:DIAG_LIST_LIMIT])
                    extra = " ..." if len(hets) > DIAG_LIST_LIMIT else ""
                    print(f"[!] Ligand {name}:{num} not found in {pdb.upper()}.{chain_in}")
                    print(f"    HET in {pdb.upper()}.{chain_in} (first {min(len(hets), DIAG_LIST_LIMIT)}): "
                          f"{preview}{extra}")

                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": name, "token_number": num,
                    "alias_used": comp_up, "class": comp_class,
                    "status": "not_found", "found_mode": "",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
            else:
                # Build labels
                atoms = [a for a in ch.get_atoms() if (not EXCLUDE_HYDROGENS or a.element != 'H')]
                ns = NeighborSearch(atoms)
                close_res = set()
                for lig in lig_res:
                    for atom in lig:
                        if EXCLUDE_HYDROGENS and atom.element == 'H':
                            continue
                        for nb in ns.search(atom.coord, DIST_CUTOFF):
                            parent = nb.get_parent()
                            if is_aa(parent):
                                close_res.add(parent)

                seq_res   = [r for r in ch if is_aa(r)]
                positions = [r.get_id()[1] for r in seq_res]
                labels    = [1 if r in close_res else 0 for r in seq_res]

                out_entries.append({
                    "pdb_id": pdb.upper(),
                    "chain_id": ch.id,
                    "sequence_residue_numbers": positions,
                    "labels": labels,
                    "ligand_name": comp_up,
                    "ligand_number": num
                })

                matched_names = ",".join(sorted({r.get_resname().upper() for r in lig_res}))
                auth_num = ""
                entry = label_map.get(chain_in, {})
                if entry and num in entry.get("label2auth", {}):
                    auth_num = str(entry["label2auth"][num][0])

                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": name, "token_number": num,
                    "alias_used": comp_up, "class": comp_class,
                    "status": "used", "found_mode": found_mode,
                    "matched_chain": ch.id, "matched_auth_number": auth_num,
                    "matched_names": matched_names
                })
                n_used_classic += 1

        # ---- (B) NAME-ONLY FALLBACK on this row (Excel codes + optional overrides) ----
        # gather Excel name-only codes
        nameonly_list = extract_nameonly_codes(lig_str)
        # add manual overrides (if any) for this (PDB,CHAIN)
        nameonly_override_set = overrides.get((pdb.upper(), chain_in), set())
        nameonly_list += list(nameonly_override_set)

        for code in nameonly_list:
            comp_up = canon_name(code)

            # Fallback policy:
            # - skip if this row already had any classic token with the same comp_up
            # - skip if we've already used name-only for this (pdb, chain, comp_up) before
            if comp_up in classic_codes_on_row:
                n_skipped_nameonly_due_to_classic += 1
                continue
            key_no = (pdb, chain_in, comp_up)
            if key_no in seen_nameonly:
                continue
            seen_nameonly.add(key_no)

            comp_class = classify_component(comp_up, label_map, chain_in)
            if comp_class in ("cap", "modified_aa", "polymeric"):
                if PRINT_DIAGNOSTICS:
                    why = "cap" if comp_class == "cap" else ("modAA" if comp_class == "modified_aa" else "polymeric")
                    print(f"[i] Skipping {why} component {code} in {pdb.upper()}.{chain_in} (name-only)")
                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": code, "token_number": "",
                    "alias_used": comp_up, "class": comp_class,
                    "status": "skipped", "found_mode": f"skip_{comp_class}",
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
                continue

            origin = "name_only_override" if code in nameonly_override_set else "name_only"

            # 1) standard name-only search: any residue in this chain with that component name
            lig_res = residues_by_name_any_number(ch, comp_up)
            found_mode = "author_name_any"

            # 2) numeric-only permissive fallback
            used_numeric_fallback = False
            if not lig_res and code.isdigit():
                auth_num = int(code)
                uniq_hits = residues_by_author_number_only(ch, auth_num)  # excludes waters & modAAs
                if len(uniq_hits) == 1:
                    lig_res = uniq_hits
                    found_mode = "name_only_numeric_unique_author"
                    used_numeric_fallback = True

            if not lig_res:
                # Not found (neither name-only nor numeric-only)
                if PRINT_DIAGNOSTICS:
                    hets = []
                    for r in ch:
                        if not is_aa(r):
                            het, n, ic = r.get_id()
                            hets.append(f"{r.get_resname()}:{int(n)}{(ic or '').strip()}")
                    preview = ", ".join(hets[:DIAG_LIST_LIMIT])
                    extra = " ..." if len(hets) > DIAG_LIST_LIMIT else ""
                    print(f"[!] Ligand {code} ({origin}) not found in {pdb.upper()}.{chain_in}")
                    print(f"    HET in {pdb.upper()}.{chain_in} (first {min(len(hets), DIAG_LIST_LIMIT)}): "
                          f"{preview}{extra}")

                audit_rows.append({
                    "pdb": pdb.upper(), "chain_req": chain_in, "token_name": code, "token_number": "",
                    "alias_used": comp_up, "class": comp_class,
                    "status": "not_found", "found_mode": origin,
                    "matched_chain": "", "matched_auth_number": "", "matched_names": ""
                })
                continue  # next name-only token

            # --- Build labels (shared for both modes) ---
            atoms = [a for a in ch.get_atoms() if (not EXCLUDE_HYDROGENS or a.element != 'H')]
            ns = NeighborSearch(atoms)
            close_res = set()
            for lig in lig_res:
                for atom in lig:
                    if EXCLUDE_HYDROGENS and atom.element == 'H':
                        continue
                    for nb in ns.search(atom.coord, DIST_CUTOFF):
                        parent = nb.get_parent()
                        if is_aa(parent):
                            close_res.add(parent)

            seq_res   = [r for r in ch if is_aa(r)]
            positions = [r.get_id()[1] for r in seq_res]
            labels    = [1 if r in close_res else 0 for r in seq_res]

            # For numeric-fallback, we *know* the author number and actual comp name
            if used_numeric_fallback:
                actual_name = lig_res[0].get_resname().upper()
                actual_num  = int(lig_res[0].get_id()[1])
                out_entries.append({
                    "pdb_id": pdb.upper(),
                    "chain_id": ch.id,
                    "sequence_residue_numbers": positions,
                    "labels": labels,
                    "ligand_name": actual_name,
                    "ligand_number": actual_num
                })
                matched_names = actual_name
                matched_auth_number = str(actual_num)
                fm = found_mode  # numeric fallback mode
            else:
                out_entries.append({
                    "pdb_id": pdb.upper(),
                    "chain_id": ch.id,
                    "sequence_residue_numbers": positions,
                    "labels": labels,
                    "ligand_name": comp_up,   # original name-only token (or override)
                    "ligand_number": None
                })
                matched_names = ",".join(sorted({r.get_resname().upper() for r in lig_res}))
                matched_auth_number = ""
                fm = origin  # name_only or name_only_override

            audit_rows.append({
                "pdb": pdb.upper(), "chain_req": chain_in, "token_name": code, "token_number": "",
                "alias_used": comp_up, "class": comp_class,
                "status": "used", "found_mode": fm,
                "matched_chain": ch.id, "matched_auth_number": matched_auth_number,
                "matched_names": matched_names
            })
            n_used_nameonly += 1

    # Write outputs
    with open(OUTPUT_JSON, "w") as f:
        json.dump(out_entries, f, indent=2)
    print(f"Saved {OUTPUT_JSON} with {len(out_entries)} chain entries")

    # Audit CSV
    fieldnames = [
        "pdb","chain_req","token_name","token_number","alias_used","class",
        "status","found_mode","matched_chain","matched_auth_number","matched_names"
    ]
    with open(AUDIT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in audit_rows:
            w.writerow(r)
    print(f"Saved audit report {AUDIT_CSV} with {len(audit_rows)} rows")

    # Tiny summary to sanity-check duplication
    print(f"[summary] used classic: {n_used_classic} | used name-only: {n_used_nameonly} "
          f"| name-only skipped (classic present): {n_skipped_nameonly_due_to_classic}")

if __name__ == "__main__":
    main()
