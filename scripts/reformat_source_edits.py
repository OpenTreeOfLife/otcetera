#!/usr/bin/env python3
import sys
import json
import itertools

def main(edits_fp):
    with open(edits_fp, 'r') as inp:
        eblob = json.load(inp)
    a = {
        "add taxon": [],
        "delete taxon": [],
        "add synonym": [],
        "delete synonym": [],
        "change id": [],
        "change name": [],
        "change rank": {},
        "change flags": {},
    }
    aat = a["add taxon"]
    adt = a["delete taxon"]
    aas = a["add synonym"]
    ads = a["delete synonym"]
    aci = a["change id"]
    acn = a["change name"]
    acr = a["change rank"]
    acf = a["change flags"]
    for el in eblob.get("alpha", []):
        op = el["operation"]
        if op == "add taxon":
            aat.append((el["taxon_id"], el["name"], el["rank"]))
        elif op =="delete taxon":
            adt.append(el["taxon_id"])
        elif op =="add synonym":
            aas.append((el["taxon_id"], el["synonym"], el.get("type", "")))
        elif op =="delete synonym":
            ads.append((el["taxon_id"], el["synonym"], el.get("type", "")))
        elif op =="change id":
            aci.append((el["from"], el["to"]))
        elif op =="change name":
            acn.append((el["taxon_id"], el["to"]))
        elif op =="change rank":
            acr.setdefault(el["from"], {}).setdefault(el["to"], []).append(el["taxon_id"])
        elif op =="change flags":
            acf.setdefault(el["from"], {}).setdefault(el["to"], []).append(el["taxon_id"])
        else:
            raise RuntimeError(f"Unkown operation {op} in alpha edits")
    g = {
        "add taxa": [],
        "deleted grouping": [],
        "new grouping": [],
        "add synonym": [],
        "delete synonym": [],
        "change id": [],
        "change name": [],
        "change rank": {},
        "change flags": {}
    }
    gat = g["add taxa"]
    gdg = g["deleted grouping"]
    gng = g["new grouping"]
    gas = g["add synonym"]
    gds = g["delete synonym"]
    gci = g["change id"]
    gcn = g["change name"]
    gcr = g["change rank"]
    gcf = g["change flags"]
    for el in itertools.chain(eblob.get("alpha_groups", []), eblob.get("higher", [])):
        op = el["operation"]
        if op == "add taxa":
            gat.append((el["taxon_id"], el["added"]))
        elif op =="deleted grouping":
            gdg.append(el["taxon_id"])
        elif op =="new grouping":
            gng.append((el["taxon_id"], el["name"], el.get("rank", ""), el["added"]))
        elif op =="add synonym":
            gas.append((el["taxon_id"], el["synonym"], el.get("type", "")))
        elif op =="delete synonym":
            gds.append((el["taxon_id"], el["synonym"], el.get("type", "")))
        elif op =="change id":
            gci.append((el["from"], el["to"]))
        elif op =="change name":
            gcn.append((el["taxon_id"], el["to"]))
        elif op =="change rank":
            gcr.setdefault(el["from"], {}).setdefault(el["to"], []).append(el["taxon_id"])
        elif op =="change flags":
            gcf.setdefault(el["from"], {}).setdefault(el["to"], []).append(el["taxon_id"])
        else:
            raise RuntimeError(f"Unkown operation {op} in alpha edits")
    ref_blob = {"alpha": a, "groups": g}
    json.dump(ref_blob, sys.stdout, indent=0)
    





# "operation": "add taxon",
# "operation": "delete taxon",

# "operation": "add synonym",
# "operation": "delete synonym",
# "operation": "change id",
# "operation": "change name",
# "operation": "change rank",
# "operation": "change flags",

# "operation": "add taxa",
# "operation": "deleted grouping",
# "operation": "new grouping",

if __name__ == '__main__':
    main(sys.argv[1])


