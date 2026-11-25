#!/usr/bin/env python3
"""Fix Ruslan's gene sets: his mitochondrial list (142, from his workspace) and the KEGG human pathways he used
(OXPHOS hsa00190, Thermogenesis hsa04714) plus his negative control (cell cycle hsa04110) and ribosome hsa03010."""
import hashlib, json, re, urllib.request
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
RES = RUN.parent
raw = RUN / "reference/kegg"
raw.mkdir(parents=True, exist_ok=True)


def fetch(url, name):
    p = raw / name
    if not p.exists():
        p.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    return p.read_text()


sym = {}
for line in fetch("https://rest.kegg.jp/list/hsa", "hsa_list.txt").splitlines():
    parts = line.split("\t")
    if len(parts) >= 4 and parts[1] == "CDS":
        sym[parts[0]] = parts[3].split(";")[0].split(",")[0].strip()
sets = {}
src = {}
for name, pid in [("OXPHOS_hsa00190", "hsa00190"), ("Thermogenesis_hsa04714", "hsa04714"),
                  ("CellCycle_hsa04110_control", "hsa04110"), ("Ribosome_hsa03010", "hsa03010")]:
    txt = fetch(f"https://rest.kegg.jp/link/hsa/path:{pid}", f"{pid}.txt")
    ids = [l.split("\t")[1] for l in txt.splitlines() if "\t" in l]
    sets[name] = sorted({sym[i] for i in ids if i in sym})
    src[name] = {"file": f"reference/kegg/{pid}.txt", "sha256": hashlib.sha256((raw / f"{pid}.txt").read_bytes()).hexdigest(), "n_kegg_ids": len(ids)}
mt = RES / "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt"
sets["Ruslan_MT142"] = sorted({l.strip() for l in mt.read_text().splitlines() if l.strip()})
src["Ruslan_MT142"] = {"file": "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt", "sha256": hashlib.sha256(mt.read_bytes()).hexdigest()}
src["hsa_list"] = {"file": "reference/kegg/hsa_list.txt", "sha256": hashlib.sha256((raw / "hsa_list.txt").read_bytes()).hexdigest()}
(RUN / "config/gene_sets.json").write_text(json.dumps({"sets": sets, "sources": src}, indent=1) + "\n")
print({k: len(v) for k, v in sets.items()})
