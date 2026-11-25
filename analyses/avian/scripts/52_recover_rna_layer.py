#!/usr/bin/env python3
"""Carry the gene recovery over to the RNA-seq layer.

The CAGE recovery (scripts 50, 51) works in genome coordinates. The GEO RNA-seq
tables cannot be touched that way: their rows are Augustus transcript ids such as
'g1.t1', not coordinates. The link exists inside the CAGE table, which carries
both the Augustus model id and that model's coordinates (geneChr / geneStart /
geneEnd / geneStrand / geneId). So:

  miniprot locus of gene X  ->  overlapping Augustus model g<N>  ->  RNA rows
  g<N>.t1, g<N>.t2 ...      ->  value for gene X in that species

Only the three Augustus species need this; the four native species already carry
symbols in both layers. The RNA panel is fibroblasts only — there is no epiblast
RNA-seq for these seven species (checked: two norm columns per file, both with
fibroblast C-codes).
"""
from __future__ import annotations

import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    load_geo_rnaseq,
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
    nuclear_oxphos,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
GFFDIR = RES / "gene_recovery"
OUT = RES / "recovered_rna"
AUG = ["Turkey", "Peacock", "Emu"]
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))


def fit(y) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    r = stats.linregress(LOGM[ok], y[ok])
    return float(r.rvalue), float(r.pvalue)


def norm(it) -> list[str]:
    return sorted({normalize_symbol(g) for g in it} - {""})


def miniprot_loci(species: str) -> pd.DataFrame:
    """best locus per gene, from the GFF written by script 50"""
    gff = GFFDIR / f"{species}_miniprot.gff"
    # the query FASTA written by script 50 carries only the protein id in its
    # header; the symbol lives in the source chicken protein file
    src = Path("/path/to/storage/kleiber_birds_gemcat/ref/chicken_prot/query_proteins_all.faa")
    pid2sym = dict(re.findall(r"^>(\S+)\s+(\S+)\s+\[organism=", src.read_text(), re.M))
    rows = []
    for line in gff.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        at = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
        sym = pid2sym.get(at.get("Target", "").split()[0])
        if not sym:
            continue
        seq = f[0].split()[0]
        if seq.startswith("ENA|"):
            seq = seq.split("|")[-1]
        rows.append(dict(gene=sym, seqname=seq, start=int(f[3]), end=int(f[4]),
                         strand=f[6], score=float(f[5]) if f[5] != "." else 0.0))
    d = pd.DataFrame(rows)
    return d.sort_values("score", ascending=False).drop_duplicates("gene")


def augustus_models(species: str) -> pd.DataFrame:
    """Coordinates of the Augustus models, taken from the peak rows.

    geneChr is ChIPseeker's internal chromosome index ('1', '2', ...), not the
    sequence name, so the model's sequence is read from the peak's own seqnames.
    geneStrand is coded 1 / 2.
    """
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    m = d[["geneId", "seqnames", "geneStart", "geneEnd", "geneStrand"]].dropna().copy()
    # Turkey and Peacock carry gene-level ids ('g23687'), Emu carries transcript-level
    # ids ('g2047.t1'). RNA rows are transcript-level everywhere. Strip on both sides.
    m["geneId"] = m["geneId"].astype(str).str.replace(r"\.t\d+$", "", regex=True)
    m = m.drop_duplicates("geneId")
    m["geneStrand"] = m["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna("+")
    return m.rename(columns={"seqnames": "seqname", "geneStart": "start", "geneEnd": "end", "geneStrand": "strand"})


def rna_values(species: str, scale: str = "share") -> pd.Series:
    """Augustus model id (without .tN) -> value, summed over the model's transcripts.

    scale='share' rescales every library to a within-table share x 1e6 before
    averaging, because GEO 'norm' totals are not a common cross-species scale
    (3.96e5 in Zebra against 1.02e6 in Peacock). scale='as_is' keeps the raw
    'norm' values and is kept only for comparison.
    """
    d = pd.read_excel(DATA / SPECIES[species]["geo_rna"])
    cols = [c for c in d.columns if str(c).startswith("norm")]
    key = d.columns[0]
    m = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if scale == "share":
        m = m / m.sum(axis=0) * 1e6
    v = m.mean(axis=1)
    model = d[key].astype(str).str.replace(r"\.t\d+$", "", regex=True)
    return v.groupby(model).sum()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # base matrix on the common scale: within-table share x 1e6, canonical symbols
    cols = {}
    for sp in SPECIES_ORDER:
        g = load_geo_rnaseq(sp)
        m = g["gene_mat"].div(g["lib_total"], axis=1) * 1e6
        cols[sp] = m.mean(axis=1)
        print(f"[{sp}] строк сопоставлено {g['n_mapped_rows']:>6} из {g['n_raw']:>6}, символов {g['n_symbols']:>6}")
    rna = pd.DataFrame(cols)
    before = rna.copy()

    added: dict[str, list[str]] = {}
    for sp in AUG:
        loci = miniprot_loci(sp)
        mods = augustus_models(sp)
        vals = rna_values(sp)
        got = []
        for _, L in loci.iterrows():
            g = L["gene"]
            if g in rna.index and pd.notna(rna.loc[g, sp]):
                continue
            c = mods[(mods["seqname"] == L["seqname"]) & (mods["strand"] == L["strand"])].copy()
            if c.empty:
                continue
            c["ov"] = np.minimum(c["end"], L["end"]) - np.maximum(c["start"], L["start"])
            c = c[c["ov"] > 0]
            if c.empty:
                continue
            mid = c.sort_values("ov", ascending=False).iloc[0]["geneId"]
            if mid not in vals.index:
                continue
            val = float(vals.loc[mid])
            if val <= 0:
                continue
            if g not in rna.index:
                rna.loc[g] = np.nan
            rna.loc[g, sp] = val
            got.append(g)
        added[sp] = sorted(got)
        print(f"[{sp}] РНК: добавлено {len(got)} генов")

    pd.DataFrame([{"species": k, "n_added": len(v), "genes": ";".join(v)} for k, v in added.items()]).to_csv(
        OUT / "rna_genes_added.tsv", sep="\t", index=False
    )

    sets = {
        "KEGG hsa00190": norm(load_geneset("kegg_00190")),
        "его OXPHOS 133": norm(load_geneset("ruslan_oxphos_133")),
        "наш ядерный OXPHOS 77": norm(nuclear_oxphos()),
        "KEGG hsa04714": norm(load_geneset("kegg_04714")),
        "контроль: гликолиз": norm(load_geneset("GLYCOLYSIS")),
    }
    rows = []
    for name, genes in sets.items():
        for tag, w in (("до", before), ("после", rna)):
            hit = w.reindex(genes)
            cnt = hit.notna().sum(axis=0).to_numpy(float)
            mean = np.log10(hit.mean(axis=0).to_numpy(float))
            tot = np.log10(hit.sum(axis=0).to_numpy(float))
            rc, pc = fit(cnt); rm, pm = fit(mean); rs, ps = fit(tot)
            rows.append(dict(gene_set=name, вариант=tag,
                             n_detected=";".join(str(int(v)) for v in cnt),
                             r_count=round(rc, 3), p_count=round(pc, 4),
                             r_mean=round(rm, 3), p_mean=round(pm, 4),
                             r_sum=round(rs, 3), p_sum=round(ps, 4)))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "rna_before_after.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n" + "=" * 110)
    print("RNA-seq фибробластов, GEO, до и после восстановления")
    print("=" * 110)
    print(tab.to_string(index=False))

    # fixed core on the RNA layer alone
    print("\nЛестница по широте на его 133, ядро по RNA:")
    genes = norm(load_geneset("ruslan_oxphos_133"))
    for tag, w in (("до", before), ("после", rna)):
        b = (w.notna() & (w > 0)).sum(axis=1)
        line = []
        for k in (7, 6, 5, 4, 3):
            core = [g for g in genes if b.get(g, 0) >= k]
            r, p = fit(np.log10(w.reindex(core).mean(axis=0).to_numpy(float)))
            line.append(f"k={k} n={len(core):>3} r={r:+.3f} p={p:.3f}")
        print(f"  {tag:>6}: " + " | ".join(line))

    print(f"\nwrote {OUT / 'rna_before_after.tsv'}")


if __name__ == "__main__":
    main()
