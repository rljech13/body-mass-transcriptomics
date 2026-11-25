#!/usr/bin/env python3
"""Modules vs mass on the rebuilt panel, plus the centred-within-species network.

Three things, written next to the old run rather than over it:

1. The standard module-trait table of script 65, both library-wise (optimistic)
   and species-wise, plus matched-null p with BH across modules — that last
   correction was missing from the first write-up of green and yellow.

2. The centred-within-species network of script 67, same defaults.

3. A sweep of mergeCutHeight / deepSplit / minModuleSize so the centred network
   is not judged on the four-module degeneration that made yellow's BH of 0.028
   cheaper than green's 0.054. The question is whether a growth-like module
   survives when 15-20 modules are asked of the same data.

The OXPHOS working set is the intersection of the old 55 with the new fibroblast
universe (47 genes). The eight that dropped out are listed, not forced back in.
"""
from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER  # noqa: E402

IN = RES / "wgcna_input_v2"
WG = RES / "wgcna_out_v2"
OUT = RES / "modules_v2"
TAG = "common_only_u5000d1000"
OXFILE = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
OLD_YELLOW = RES / "modules_within_species" / "gene_modules_fib_centered.tsv"
RSCRIPT = "/path/to/user-home/miniconda3/envs/wgcna_r/bin/Rscript"


def bh(p) -> np.ndarray:
    return stats.false_discovery_control(np.nan_to_num(np.asarray(p, float), nan=1.0), method="bh")


R_CODE = r"""
suppressMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads(16)
a <- commandArgs(TRUE)
inf <- a[1]; out <- a[2]
mch <- as.numeric(a[3]); ds <- as.integer(a[4]); minsz <- as.integer(a[5])
m <- as.matrix(read.table(inf, sep="\t", header=TRUE, row.names=1, check.names=FALSE))
d <- t(m)
gsg <- goodSamplesGenes(d, verbose=0); if(!gsg$allOK) d <- d[gsg$goodSamples, gsg$goodGenes, drop=FALSE]
sft <- pickSoftThreshold(d, powerVector=c(1:10, seq(12,30,2)), networkType="signed hybrid", verbose=0)
ok <- sft$fitIndices$SFT.R.sq >= 0.8
p <- if(any(ok)) sft$fitIndices$Power[which(ok)[1]] else 14
cat(sprintf("softPower=%d  mergeCutHeight=%.2f  deepSplit=%d  minModuleSize=%d  (R2>=0.8: %s)\n",
            p, mch, ds, minsz, any(ok)))
net <- blockwiseModules(d, power=p, networkType="signed hybrid", TOMType="signed",
                        minModuleSize=minsz, deepSplit=ds, mergeCutHeight=mch,
                        numericLabels=TRUE, pamRespectsDendro=FALSE, maxBlockSize=20000, verbose=0)
col <- labels2colors(net$colors)
write.table(data.frame(gene=colnames(d), module=col), out, sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("модулей %d (без grey %d), grey %d\n", length(unique(col)),
            length(setdiff(unique(col),"grey")), sum(col=="grey")))
"""


def trait_table(net: str, samples: pd.DataFrame, ox: list[str]) -> pd.DataFrame:
    gm = pd.read_csv(WG / f"gene_modules_{net}.tsv", sep="\t")
    me = pd.read_csv(WG / f"module_eigengenes_{net}.tsv", sep="\t").set_index("sample")
    s = samples.loc[me.index]
    mass_lib = s["log10_mass"].to_numpy(float)
    rows = []
    for mod in me.columns:
        v = me[mod].to_numpy(float)
        r1, p1 = stats.pearsonr(mass_lib, v)
        per = pd.DataFrame({"sp": s["species"].to_numpy(), "v": v}).groupby("sp")["v"].mean()
        per = per.reindex(SPECIES_ORDER)
        r2, p2 = s53.fit(per.to_numpy(float))
        name = mod.replace("ME", "")
        n_in = int((gm["module"] == name).sum())
        n_ox = int(((gm["module"] == name) & gm["gene"].isin(ox)).sum())
        rows.append(dict(сеть=net, модуль=name, n_генов=n_in, n_OXPHOS=n_ox,
                         r_по_библиотекам=round(r1, 3), p_по_библиотекам=p1,
                         r_по_видам=round(r2, 3), p_по_видам=p2))
    t = pd.DataFrame(rows)
    t["bh_по_библиотекам"] = bh(t["p_по_библиотекам"]).round(4)
    t["bh_по_видам"] = bh(t["p_по_видам"]).round(4)
    t["p_по_библиотекам"] = t["p_по_библиотекам"].round(5)
    t["p_по_видам"] = t["p_по_видам"].round(4)
    return t.sort_values("r_по_видам")


def centred_and_test(M: pd.DataFrame, meta: pd.DataFrame, ox: list[str],
                     mch: float, ds: int, minsz: int, tag: str,
                     old_yellow: set[str]) -> pd.DataFrame:
    C = M.copy()
    for sp in SPECIES_ORDER:
        cols = [c for c in M.columns if c.startswith(f"{sp}|")]
        C[cols] = M[cols].sub(M[cols].mean(axis=1), axis=0)
    C = C[C.std(axis=1) > 0]
    fmat = OUT / "M_fib_centered.tsv"
    if not fmat.exists():
        C.to_csv(fmat, sep="\t")
    rf = OUT / "_run.R"
    rf.write_text(R_CODE)
    gmf = OUT / f"gene_modules_centered_{tag}.tsv"
    r = subprocess.run([RSCRIPT, str(rf), str(fmat), str(gmf),
                        str(mch), str(ds), str(minsz)],
                       capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print("R упал:", r.stderr.strip()[-400:])
        return pd.DataFrame()

    gm = pd.read_csv(gmf, sep="\t")
    permean = M.T.groupby(meta["species"].to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
    dist = gm[gm["gene"].isin(ox)]["module"].value_counts()
    rows = []
    for mod in sorted(gm["module"].unique()):
        if mod == "grey":
            continue
        members = [g for g in gm[gm["module"] == mod]["gene"] if g in permean.index]
        if len(members) < 10:
            continue
        r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
        pn, _ = s53.matched_null_p(permean, members, r_set)
        overlap = len(set(members) & old_yellow)
        rows.append(dict(настройка=tag, merge=mch, deepSplit=ds, minSize=minsz,
                         модуль=mod, n_генов=len(members),
                         n_OXPHOS=int(dist.get(mod, 0)),
                         overlap_old_yellow=overlap,
                         r_по_видам=round(r_set, 3), p=round(p_set, 4),
                         p_против_нуля=None if np.isnan(pn) else round(pn, 4)))
    t = pd.DataFrame(rows)
    if t.empty:
        return t
    t["bh_регрессии"] = bh(t["p"]).round(4)
    t["bh_нуля"] = bh(t["p_против_нуля"].fillna(1.0)).round(4)
    t = t.sort_values("r_по_видам")
    n_mod = gm["module"].nunique() - int("grey" in set(gm["module"]))
    print(f"  настроек {tag}: модулей без grey {n_mod}; "
          f"BH нуля < 0.05 у {int((t['bh_нуля'] < 0.05).sum())}")
    return t


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(IN / f"samples_{TAG}.tsv", sep="\t").set_index("sample")
    ox_all = [l.strip() for l in OXFILE.read_text().splitlines() if l.strip()]
    uf = set(pd.read_csv(IN / f"universe_fib_{TAG}.txt", header=None)[0])
    ox = [g for g in ox_all if g in uf]
    lost = [g for g in ox_all if g not in uf]
    print(f"рабочий набор OXPHOS: {len(ox)} из {len(ox_all)} в новой вселенной; нет: {', '.join(lost)}\n")
    (OUT / "working_oxphos_in_universe.txt").write_text("\n".join(ox) + "\n")

    old_yellow: set[str] = set()
    if OLD_YELLOW.exists():
        oy = pd.read_csv(OLD_YELLOW, sep="\t")
        old_yellow = set(oy.loc[oy["module"] == "yellow", "gene"])
        print(f"старый yellow центрированной сети: {len(old_yellow)} генов\n")

    for net in ("fib", "epi", "all_tc"):
        t = trait_table(net, samples, ox)
        t.to_csv(OUT / f"module_trait_{net}.tsv", sep="\t", index=False)
        print("=" * 100)
        print(f"сеть {net}: модулей {len(t) - int('grey' in set(t['модуль']))} без grey")
        print(f"  BH<0.05 по библиотекам: {int(((t['модуль']!='grey') & (t['bh_по_библиотекам']<0.05)).sum())}")
        print(f"  BH<0.05 по видам:       {int(((t['модуль']!='grey') & (t['bh_по_видам']<0.05)).sum())}")
        hit = t[t["n_OXPHOS"] > 0]
        if len(hit):
            print(hit[["модуль", "n_генов", "n_OXPHOS", "r_по_видам", "p_по_видам", "bh_по_видам"]].to_string(index=False))

    gm = pd.read_csv(WG / "gene_modules_fib.tsv", sep="\t")
    M = pd.read_csv(IN / f"M_fib_{TAG}.tsv", sep="\t", index_col=0)
    meta = samples.loc[M.columns]
    permean = M.T.groupby(meta["species"].to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
    inmod = gm[gm["gene"].isin(ox)]
    print("\nкуда попали гены OXPHOS (сеть фибробластов v2):")
    print(inmod["module"].value_counts().to_string())
    det = []
    for mod, k in inmod["module"].value_counts().items():
        members = gm[gm["module"] == mod]["gene"].tolist()
        r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
        pn, _ = s53.matched_null_p(permean, [g for g in members if g in permean.index], r_set)
        det.append(dict(модуль=mod, n_всего=len(members), n_OXPHOS=int(k),
                        r_модуля_по_видам=round(r_set, 3), p=round(p_set, 4),
                        p_против_нуля=None if np.isnan(pn) else round(pn, 4)))
        print(f"  {mod}: n={len(members)} OXPHOS={k} r={r_set:+.3f} нуль "
              f"{'—' if np.isnan(pn) else round(pn, 4)}")
    dett = pd.DataFrame(det)
    dett["bh_нуля"] = bh(dett["p_против_нуля"].fillna(1.0)).round(4)
    dett.to_csv(OUT / "oxphos_modules_detail.tsv", sep="\t", index=False)
    print("BH нуля по модулям с OXPHOS:")
    print(dett.sort_values("p_против_нуля").to_string(index=False))

    print("\n=== центрированная сеть, сетка настроек ===")
    grid = [
        (0.25, 2, 30, "default"),
        (0.15, 2, 30, "merge015"),
        (0.10, 2, 20, "merge010_sz20"),
        (0.10, 3, 20, "merge010_ds3"),
        (0.05, 4, 15, "split"),
    ]
    parts = []
    for mch, ds, minsz, tag in grid:
        print(f"\n--- {tag} ---")
        t = centred_and_test(M, meta, ox, mch, ds, minsz, tag, old_yellow)
        if len(t):
            parts.append(t)
            print(t.to_string(index=False))
    if parts:
        allc = pd.concat(parts, ignore_index=True)
        allc.to_csv(OUT / "centered_sweep.tsv", sep="\t", index=False)
        print(f"\nwrote {OUT}/centered_sweep.tsv")


if __name__ == "__main__":
    main()
