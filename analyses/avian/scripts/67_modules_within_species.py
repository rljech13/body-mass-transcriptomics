#!/usr/bin/env python3
"""Refining the modules: co-expression within species, not between them.

Script 63 showed that inside the fibroblast matrix species identity explains a
median 78% of the variance per gene. So the modules of script 64 are built mostly
on between-species differences — which is the axis body mass lives on, but it also
means the modules partly encode assembly and annotation differences.

The refinement: centre every gene within species first. What is left is
within-species, between-library covariation, i.e. co-regulation in the ordinary
sense. Modules are then defined on that, with no mass information anywhere, and
only afterwards are their UNCENTRED eigengene levels correlated with mass.

Also computes the allometric slope of the working OXPHOS set, so it can be put
next to the published exponents: cellular oxygen consumption in freshly isolated
hepatocytes scales as M^-0.18 (Porter and Brand 1995) and whole-organism Kleiber
predicts M^-0.25 per cell.
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

IN = RES / "wgcna_input"
WG = RES / "wgcna_out"
OUT = RES / "modules_within_species"
W = 1000
OXFILE = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
RSCRIPT = "/path/to/user-home/miniconda3/envs/wgcna_r/bin/Rscript"

R_CODE = r"""
suppressMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads(16)
a <- commandArgs(TRUE); inf <- a[1]; out <- a[2]
m <- as.matrix(read.table(inf, sep="\t", header=TRUE, row.names=1, check.names=FALSE))
d <- t(m)
gsg <- goodSamplesGenes(d, verbose=0); if(!gsg$allOK) d <- d[gsg$goodSamples, gsg$goodGenes, drop=FALSE]
sft <- pickSoftThreshold(d, powerVector=c(1:10, seq(12,30,2)), networkType="signed hybrid", verbose=0)
ok <- sft$fitIndices$SFT.R.sq >= 0.8
p <- if(any(ok)) sft$fitIndices$Power[which(ok)[1]] else 14
cat(sprintf("softPower=%d  (R2>=0.8 достигнут: %s)\n", p, any(ok)))
net <- blockwiseModules(d, power=p, networkType="signed hybrid", TOMType="signed",
                        minModuleSize=30, deepSplit=2, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE, maxBlockSize=20000, verbose=0)
col <- labels2colors(net$colors)
write.table(data.frame(gene=colnames(d), module=col), out, sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("модулей %d (без grey %d), grey %d\n", length(unique(col)),
            length(setdiff(unique(col),"grey")), sum(col=="grey")))
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    M = pd.read_csv(IN / f"M_fib_w{W}.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(IN / f"samples_w{W}.tsv", sep="\t").set_index("sample").loc[M.columns]
    logm = np.log10([SPECIES[s]["mass_g"] for s in SPECIES_ORDER])
    ox = [l.strip() for l in OXFILE.read_text().splitlines() if l.strip()]

    # slope of the working OXPHOS set, for comparison with published exponents
    permean = M.T.groupby(meta["species"].to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
    y = np.log10(permean.reindex([g for g in ox if g in permean.index]).mean(axis=0).to_numpy(float))
    res = stats.linregress(logm, y)
    print("НАКЛОН рабочего набора OXPHOS (55 генов), CAGE фибробластов")
    print(f"  log10(средний уровень) ~ log10(масса):  наклон {res.slope:+.4f}  r {res.rvalue:+.3f}  p {res.pvalue:.4f}")
    print(f"  95% интервал наклона: {res.slope - 1.96*res.stderr:+.4f} … {res.slope + 1.96*res.stderr:+.4f}")
    print("  для сравнения: потребление кислорода свежими гепатоцитами M^-0.18 (Porter, Brand 1995);")
    print("  Клейбер на клетку предсказывает M^-0.25\n")
    pd.DataFrame([dict(набор="kegg_00190 рабочий 55", n=len(y), наклон=round(res.slope, 4),
                       se=round(res.stderr, 4), r=round(res.rvalue, 3), p=round(res.pvalue, 4),
                       ci_low=round(res.slope - 1.96 * res.stderr, 4),
                       ci_high=round(res.slope + 1.96 * res.stderr, 4))]).to_csv(
        OUT / "oxphos_slope.tsv", sep="\t", index=False)

    # centre each gene within species
    C = M.copy()
    for sp in SPECIES_ORDER:
        cols = [c for c in M.columns if c.startswith(f"{sp}|")]
        C[cols] = M[cols].sub(M[cols].mean(axis=1), axis=0)
    keep = C.std(axis=1) > 0
    C = C[keep]
    f = OUT / "M_fib_centered.tsv"
    C.to_csv(f, sep="\t")
    print(f"центрированная матрица: {C.shape[0]} генов x {C.shape[1]} библиотек")

    rf = OUT / "_run.R"
    rf.write_text(R_CODE)
    gmf = OUT / "gene_modules_fib_centered.tsv"
    r = subprocess.run([RSCRIPT, str(rf), str(f), str(gmf)], capture_output=True, text=True)
    print(r.stdout.strip()[-400:])
    if r.returncode != 0:
        print("R упал:", r.stderr.strip()[-400:])
        return

    gm = pd.read_csv(gmf, sep="\t")
    dist = gm[gm["gene"].isin(ox)]["module"].value_counts()
    print("\nкуда попали 55 генов OXPHOS в сети на центрированных данных:")
    print(dist.to_string())

    rows = []
    for mod in sorted(gm["module"].unique()):
        if mod == "grey":
            continue
        members = [g for g in gm[gm["module"] == mod]["gene"] if g in permean.index]
        if len(members) < 10:
            continue
        # uncentred level of the module against mass
        r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
        pn, _ = s53.matched_null_p(permean, members, r_set)
        rows.append(dict(модуль=mod, n_генов=len(members),
                         n_OXPHOS=int(dist.get(mod, 0)),
                         r_по_видам=round(r_set, 3), p=round(p_set, 4),
                         p_против_нуля=None if np.isnan(pn) else round(pn, 4)))
    t = pd.DataFrame(rows).sort_values("r_по_видам")
    t["bh"] = stats.false_discovery_control(t["p"].to_numpy(), method="bh").round(3)
    t.to_csv(OUT / "centered_modules_vs_mass.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\nмодули центрированной сети против массы (уровень нецентрированный):")
    print(t.to_string(index=False))
    print(f"\nпроходят BH 0.05: {int((t['bh'] < 0.05).sum())}")
    print(f"wrote {OUT}/centered_modules_vs_mass.tsv, oxphos_slope.tsv")


if __name__ == "__main__":
    main()
