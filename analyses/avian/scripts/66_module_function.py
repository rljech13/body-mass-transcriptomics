#!/usr/bin/env python3
"""What the WGCNA modules are, functionally.

Over-representation with Fisher's exact test, and the background is the measured
universe (the 2890 genes of M_fib), not the whole KEGG library — using the library
as background is the standard way to manufacture enrichment.

Pathway definitions: the complete human KEGG pathway-to-gene map, pulled in two
REST calls and cached. Chicken maps exist for metabolism but not for everything
(gga04714 does not exist at all), so hsa is used as the common vocabulary and
gga00190 is checked separately as a sensitivity in the report.

BH is applied within each module across pathways.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, RES, normalize_symbol  # noqa: E402

IN = RES / "wgcna_input"
WG = RES / "wgcna_out"
MV = RES / "modules_vs_mass"
OUT = RES / "module_function"
CACHE = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
MIN_SET = 5
MAX_SET = 500


def kegg_all() -> dict[str, tuple[str, set[str]]]:
    if CACHE.exists():
        d = pd.read_csv(CACHE, sep="\t")
        return {r["pathway"]: (r["name"], set(str(r["genes"]).split(";"))) for _, r in d.iterrows()}
    link = urllib.request.urlopen("https://rest.kegg.jp/link/hsa/pathway", timeout=120).read().decode()
    names = urllib.request.urlopen("https://rest.kegg.jp/list/pathway/hsa", timeout=120).read().decode()
    genes = urllib.request.urlopen("https://rest.kegg.jp/list/hsa", timeout=120).read().decode()
    id2sym = {}
    for line in genes.splitlines():
        p = line.split("\t")
        if len(p) < 2:
            continue
        sym = p[-1].split(";")[0].split(",")[0].strip() if ";" in p[-1] else p[1].split(",")[0].strip()
        s = normalize_symbol(sym)
        if s:
            id2sym[p[0]] = s
    pname = {}
    for line in names.splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            pname[p[0].replace("path:", "")] = p[1].split(" - ")[0].strip()
    out: dict[str, set[str]] = {}
    for line in link.splitlines():
        p = line.split("\t")
        if len(p) < 2:
            continue
        pw = p[0].replace("path:", "")
        g = id2sym.get(p[1])
        if g:
            out.setdefault(pw, set()).add(g)
    rows = [{"pathway": k, "name": pname.get(k, k), "genes": ";".join(sorted(v))} for k, v in out.items()]
    pd.DataFrame(rows).to_csv(CACHE, sep="\t", index=False)
    return {k: (pname.get(k, k), v) for k, v in out.items()}


def ora(members: list[str], universe: set[str], sets: dict[str, tuple[str, set[str]]]) -> pd.DataFrame:
    mem = set(members) & universe
    rows = []
    for pw, (name, g) in sets.items():
        s = g & universe
        if not (MIN_SET <= len(s) <= MAX_SET):
            continue
        a = len(mem & s)
        if a < 3:
            continue
        b = len(mem) - a
        c = len(s) - a
        d = len(universe) - a - b - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(dict(путь=pw, название=name, n_в_пути=len(s), пересечение=a,
                         ожидалось=round(len(mem) * len(s) / len(universe), 1),
                         odds=round(float(odds), 2), p=p))
    t = pd.DataFrame(rows)
    if t.empty:
        return t
    return t.sort_values("p")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = kegg_all()
    print(f"путей KEGG hsa: {len(sets)}")
    gm = pd.read_csv(WG / "gene_modules_fib.tsv", sep="\t")
    universe = set(pd.read_csv(IN / "M_fib_w1000.tsv", sep="\t", index_col=0).index)
    print(f"фон (измеренная вселенная): {len(universe)} генов")
    det = pd.read_csv(MV / "oxphos_modules_detail.tsv", sep="\t")
    trait = pd.read_csv(MV / "module_trait_fib.tsv", sep="\t").set_index("модуль")

    order = trait.sort_values("r_по_видам").index.tolist()
    allrows = []
    for mod in order:
        if mod == "grey":
            continue
        members = gm[gm["module"] == mod]["gene"].tolist()
        t = ora(members, universe, sets)
        if t.empty:
            continue
        t.insert(0, "модуль", mod)
        allrows.append(t)
        top = t.head(4)
        r = trait.loc[mod, "r_по_видам"]
        nox = int(trait.loc[mod, "n_OXPHOS"])
        hit = det[det["модуль"] == mod]
        pn = hit["p_против_нуля"].iloc[0] if len(hit) else None
        print(f"\n--- {mod}: {len(members)} генов, OXPHOS {nox}, r по видам {r:+.3f}"
              + (f", против нуля {pn}" if pn is not None else "") + " ---")
        for _, x in top.iterrows():
            print(f"    {x['название'][:56]:56s} {x['пересечение']:>3}/{x['n_в_пути']:<4} "
                  f"ожид {x['ожидалось']:>5} p={x['p']:.2e}")
    if allrows:
        full = pd.concat(allrows, ignore_index=True)
        full["bh_в_модуле"] = np.nan
        for mod, g in full.groupby("модуль"):
            full.loc[g.index, "bh_в_модуле"] = stats.false_discovery_control(
                g["p"].to_numpy(), method="bh")
        full["bh_глобальный"] = stats.false_discovery_control(full["p"].to_numpy(), method="bh")
        full = full.sort_values(["модуль", "p"])
        full.to_csv(OUT / "module_ora_fib.tsv", sep="\t", index=False)
        sig_mod = full[full["bh_в_модуле"] < 0.05]
        sig_all = full[full["bh_глобальный"] < 0.05]
        print(f"\nBH внутри модуля: {len(sig_mod)} терминов в {sig_mod['модуль'].nunique()} модулях")
        print(f"BH по всем модулям сразу: {len(sig_all)} терминов в {sig_all['модуль'].nunique()} модулях")
        print(f"wrote {OUT / 'module_ora_fib.tsv'}")


if __name__ == "__main__":
    main()
