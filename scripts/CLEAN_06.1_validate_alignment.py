import json

GOLD = "Kinase_Type_ALLO/01.chain_gold_labels_CLEAN.json"
PLM  = "Kinase_Type_ALLO/03.chain_plm_labels.json"

print("=== VALIDATION OF GOLD vs PLM ===")

gold = {(e["pdb_id"], e["chain_id"]): e
        for e in json.load(open(GOLD))}
plm  = {(e["pdb_id"], e["chain_id"]): e
        for e in json.load(open(PLM))}

common = sorted(set(gold) & set(plm))
print(f"Common chains: {len(common)}")

errors = 0

for key in common:
    g = gold[key]
    p = plm[key]

    gpos = g["sequence_residue_numbers"]
    ppos = p["sequence_residue_numbers"]

    if len(gpos) != len(ppos):
        print(f"[LENGTH MISMATCH] {key}: GOLD={len(gpos)} PLM={len(ppos)}")
        errors += 1
        continue

    for i, (a, b) in enumerate(zip(gpos, ppos)):
        if int(a) != int(b):
            print(f"[POSITION MISMATCH] {key} at idx {i}: GOLD={a}, PLM={b}")
            errors += 1
            break

if errors == 0:
    print("Perfect alignment: no mismatches found!")
else:
    print(f"X Errors detected: {errors}")
