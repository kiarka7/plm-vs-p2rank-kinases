# inspect_in_pymol.pml
#
# Visual inspection of GOLD vs P2Rank vs PLM for a single PDB.
# CONFIG: change these three lines before running:

python
KINASE_FOLDER = "Kinase_Type_III"   # e.g. "Kinase_Type_I", "Kinase_ALLO", ...
PDB_ID        = "1s9j"            # PDB to inspect
PLM_THR       = 0.75              # PLM probability threshold
P2_THR        = 0.50              # P2Rank probability threshold (if labels absent)
python end


# ------------------------------------------------------------
# DO NOT EDIT BELOW (unless you want to tweak behavior)
# ------------------------------------------------------------

python
import os, json
from pymol import cmd

PDB_ID = PDB_ID.upper()

base_dir = "."
gold_path = os.path.join(base_dir, KINASE_FOLDER, "01.chain_gold_labels.json")
p2_path   = os.path.join(base_dir, KINASE_FOLDER, "02.chain_p2rank_labels.json")
plm_path  = os.path.join(base_dir, KINASE_FOLDER, "03.chain_plm_labels.json")
cif_path  = os.path.join(base_dir, KINASE_FOLDER, "structures", PDB_ID.lower() + ".cif")

print(f"[PyMOL] Inspecting {PDB_ID} in {KINASE_FOLDER}")
print(f"        GOLD:   {gold_path}")
print(f"        P2Rank: {p2_path}")
print(f"        PLM:    {plm_path}")
print(f"        CIF:    {cif_path}")

if not os.path.isfile(cif_path):
    print(f"[PyMOL][!] CIF file not found: {cif_path}")
else:
    cmd.load(cif_path, PDB_ID)

# ---- small helpers ----
def up_chain(x):
    s = str(x).strip()
    return (s or "A").upper()

def load_gold_for_pdb(path, pdb_id):
    data = json.load(open(path))
    gold_pos = {}
    gold_lab = {}
    for e in data:
        p = str(e.get("pdb_id","")).strip()[:4].upper()
        if p != pdb_id:
            continue
        c = up_chain(e.get("chain_id",""))
        key = (p,c)
        pos = [int(x) for x in e.get("sequence_residue_numbers",[])]
        lab = [int(x) for x in e.get("labels",[])]
        if len(pos) != len(lab):
            continue
        if key not in gold_pos:
            gold_pos[key] = pos
            gold_lab[key] = lab
        else:
            # OR-merge labels for same positions
            if gold_pos[key] != pos:
                print(f"[PyMOL][!] GOLD positions differ for {p}.{c}, keeping first.")
                continue
            merged = []
            for a,b in zip(gold_lab[key], lab):
                merged.append(1 if (a==1 or b==1) else 0)
            gold_lab[key] = merged
    return gold_pos, gold_lab

def load_p2_for_pdb(path, pdb_id):
    if not os.path.isfile(path):
        print(f"[PyMOL][i] P2Rank JSON not found: {path}")
        return {}, {}, {}
    data = json.load(open(path))
    pos_map, prob_map, lab_map = {}, {}, {}
    for e in data:
        p = str(e.get("pdb_id","")).strip()[:4].upper()
        if p != pdb_id:
            continue
        c = up_chain(e.get("chain_id",""))
        key = (p,c)
        pos_map[key] = [int(x) for x in e.get("sequence_residue_numbers",[])]
        prob_map[key]= [float(x) for x in e.get("pred_probs",[])]
        lab_map[key] = [int(x) for x in e.get("labels",[])] if "labels" in e else []
    return pos_map, prob_map, lab_map

def load_plm_for_pdb(path, pdb_id):
    if not os.path.isfile(path):
        print(f"[PyMOL][i] PLM JSON not found: {path}")
        return {}, {}
    data = json.load(open(path))
    pos_map, prob_map = {}, {}
    for e in data:
        p = str(e.get("pdb_id","")).strip()[:4].upper()
        if p != pdb_id:
            continue
        c = up_chain(e.get("chain_id",""))
        key = (p,c)
        pos_map[key] = [int(x) for x in e.get("sequence_residue_numbers",[])]
        prob_map[key]= [float(x) for x in e.get("pred_probs",[])]
    return pos_map, prob_map

def make_selection_from_residues(object_name, chain_id, residues, sel_name):
    """Create a PyMOL selection like: model PDB_ID and chain A and resi 10+15+23"""
    if not residues:
        return
    resi_str = "+".join(str(r) for r in sorted(set(residues)))
    expr = f"({object_name} and chain {chain_id} and resi {resi_str})"
    cmd.select(sel_name, expr)

# ---- load data ----
gold_pos, gold_lab = load_gold_for_pdb(gold_path, PDB_ID)
p2_pos, p2_probs, p2_labels = load_p2_for_pdb(p2_path, PDB_ID)
plm_pos, plm_probs = load_plm_for_pdb(plm_path, PDB_ID)

if not gold_pos:
    print(f"[PyMOL][!] No GOLD entries found for {PDB_ID} in {gold_path}")
    # still leave the structure loaded if present
else:
    print(f"[PyMOL] GOLD chains for {PDB_ID}: {sorted({ch for (p,ch) in gold_pos.keys()})}")

# Global selections (across chains)
all_gold_res = []
all_p2_res   = []
all_plm_res  = []

# Per-chain selections
for (pdb, chain) in sorted(gold_pos.keys()):
    pos = gold_pos[(pdb, chain)]
    lab = gold_lab[(pdb, chain)]
    L = len(pos)
    if L == 0:
        continue

    gold_res = [r for r,y in zip(pos, lab) if y == 1]
    all_gold_res.extend([(chain,r) for r in gold_res])

    # P2Rank
    p2_res = []
    p2_pos_list = p2_pos.get((pdb, chain), [])
    p2_probs_list = p2_probs.get((pdb, chain), [])
    p2_labels_list = p2_labels.get((pdb, chain), [])

    if p2_pos_list and p2_pos_list == pos:
        if p2_labels_list and len(p2_labels_list) == L:
            p2_res = [r for r,y in zip(pos, p2_labels_list) if y == 1]
        elif p2_probs_list and len(p2_probs_list) == L:
            p2_res = [r for r,pr in zip(pos, p2_probs_list) if pr >= P2_THR]
        else:
            print(f"[PyMOL][!] P2Rank arrays mismatched for {pdb}.{chain}")
    elif p2_pos_list:
        print(f"[PyMOL][!] P2Rank positions not aligned with GOLD for {pdb}.{chain} "
              f"(len GOLD={L}, len P2={len(p2_pos_list)})")

    all_p2_res.extend([(chain,r) for r in p2_res])

    # PLM
    plm_res = []
    plm_pos_list = plm_pos.get((pdb, chain), [])
    plm_probs_list = plm_probs.get((pdb, chain), [])

    if plm_pos_list and plm_pos_list == pos:
        if plm_probs_list and len(plm_probs_list) == L:
            plm_res = [r for r,pr in zip(pos, plm_probs_list) if pr >= PLM_THR]
        else:
            print(f"[PyMOL][!] PLM arrays mismatched for {pdb}.{chain}")
    elif plm_pos_list:
        print(f"[PyMOL][!] PLM positions not aligned with GOLD for {pdb}.{chain} "
              f"(len GOLD={L}, len PLM={len(plm_pos_list)})")

    all_plm_res.extend([(chain,r) for r in plm_res])

    # Per-chain selections
    if gold_res:
        make_selection_from_residues(PDB_ID, chain, gold_res, f"{PDB_ID}_gold_{chain}")
    if p2_res:
        make_selection_from_residues(PDB_ID, chain, p2_res, f"{PDB_ID}_p2_{chain}")
    if plm_res:
        make_selection_from_residues(PDB_ID, chain, plm_res, f"{PDB_ID}_plm_{chain}")

# Global selections (across chains): use a trick with temporary selections
def make_global_selection(tag, items):
    """
    items: list of (chain, resi)
    Creates selection PDBID_tag_all
    """
    if not items:
        return
    parts = []
    for ch, r in items:
        parts.append(f"(model {PDB_ID} and chain {ch} and resi {r})")
    expr = " or ".join(parts)
    cmd.select(f"{PDB_ID}_{tag}_all", expr)

make_global_selection("gold", all_gold_res)
make_global_selection("p2",   all_p2_res)
make_global_selection("plm",  all_plm_res)

# Basic styling
if cmd.count_atoms(PDB_ID) > 0:
    cmd.hide("everything", PDB_ID)
    cmd.show("cartoon", PDB_ID)
    cmd.color("grey70", PDB_ID)

    # color and show sticks for selections
    if cmd.count_atoms(f"{PDB_ID}_gold_all") > 0:
        cmd.show("sticks", f"{PDB_ID}_gold_all")
        cmd.color("red", f"{PDB_ID}_gold_all")
    if cmd.count_atoms(f"{PDB_ID}_p2_all") > 0:
        cmd.show("sticks", f"{PDB_ID}_p2_all")
        cmd.color("cyan", f"{PDB_ID}_p2_all")
    if cmd.count_atoms(f"{PDB_ID}_plm_all") > 0:
        cmd.show("sticks", f"{PDB_ID}_plm_all")
        cmd.color("yellow", f"{PDB_ID}_plm_all")

    # little bit nicer view
    cmd.orient(PDB_ID)
    cmd.zoom(PDB_ID, 5.0)

print("[PyMOL] Done. Use selections like:")
print(f"        {PDB_ID}_gold_all (red)  – GOLD positives")
print(f"        {PDB_ID}_p2_all   (cyan) – P2Rank positives")
print(f"        {PDB_ID}_plm_all  (yellow) – PLM positives (thr={PLM_THR})")
python end
