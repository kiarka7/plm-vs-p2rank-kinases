#!/usr/bin/env python3
"""
99.pymol_batch.py

Generate PyMOL images for GOLD vs P2Rank vs PLM for all chains.
Output: dataset/diag/pymol_img/*.png
"""

import os, json
from pymol import cmd

def load_json(path):
    return json.load(open(path))

def main(folder):
    gold = load_json(os.path.join(folder, "01.chain_gold_labels.json"))
    p2   = load_json(os.path.join(folder, "02.chain_p2rank_labels.json"))
    plm  = load_json(os.path.join(folder, "03.chain_plm_labels.json"))

    idx_p2  = {(e["pdb_id"], e["chain_id"]): e for e in p2}
    idx_plm = {(e["pdb_id"], e["chain_id"]): e for e in plm}

    outdir = os.path.join(folder, "diag", "pymol_img")
    os.makedirs(outdir, exist_ok=True)

    for g in gold:
        pdb, ch = g["pdb_id"], g["chain_id"]
        cif = os.path.join(folder, "structures", f"{pdb.lower()}.cif")
        if not os.path.isfile(cif):
            continue

        cmd.reinitialize()
        cmd.load(cif, pdb)

        # GOLD
        gold_sel = f"{pdb}_gold"
        gold_pos = g["gold_residue_numbers"]
        if gold_pos:
            cmd.select(gold_sel, f"resi { '+'.join(str(x) for x in gold_pos) } and chain {ch} and {pdb}")
            cmd.color("red", gold_sel)
            cmd.show("sticks", gold_sel)

        # P2Rank
        if (pdb,ch) in idx_p2:
            p2_sel = f"{pdb}_p2"
            p2_pos = idx_p2[(pdb,ch)]["pred_labels_fixed"]
            p2_pos = [i+1 for i,v in enumerate(p2_pos) if v==1]
            if p2_pos:
                cmd.select(p2_sel, f"resi { '+'.join(str(x) for x in p2_pos) } and chain {ch} and {pdb}")
                cmd.color("cyan", p2_sel)
                cmd.show("sticks", p2_sel)

        # PLM
        if (pdb,ch) in idx_plm:
            plm_sel = f"{pdb}_plm"
            plm_pos = idx_plm[(pdb,ch)]["pred_labels_fixed"]
            plm_pos = [i+1 for i,v in enumerate(plm_pos) if v==1]
            if plm_pos:
                cmd.select(plm_sel, f"resi { '+'.join(str(x) for x in plm_pos) } and chain {ch} and {pdb}")
                cmd.color("yellow", plm_sel)
                cmd.show("sticks", plm_sel)

        # styling
        cmd.bg_color("white")
        cmd.orient(pdb)
        cmd.zoom(pdb, 3.0)

        outfile = os.path.join(outdir, f"{pdb}_{ch}.png")
        cmd.png(outfile, width=1200, height=1000, dpi=300)
        print(f"[ok] Saved {outfile}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    args = ap.parse_args()

    main(args.folder)
