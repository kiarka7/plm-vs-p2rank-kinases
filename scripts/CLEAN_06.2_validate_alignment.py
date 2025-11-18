import json

GOLD = "Kinase_Type_I/01.chain_gold_labels_CLEAN.json"
P2   = "Kinase_Type_I/02.chain_p2rank_labels.json"

gold_json = json.load(open(GOLD))
p2_json   = json.load(open(P2))

# 🔥 Normalize both PDB ID and chain ID to uppercase
gold = {(e["pdb_id"].upper(), e["chain_id"].upper())
        for e in gold_json}

p2 = {(e["pdb_id"].upper(), e["chain_id"].upper())
      for e in p2_json}

missing = gold - p2

print(f"Chains in GOLD: {len(gold)}")
print(f"Chains in P2:   {len(p2)}")

if missing:
    print("X P2Rank missing chains:", missing)
else:
    print("P2Rank covers all GOLD chains")
