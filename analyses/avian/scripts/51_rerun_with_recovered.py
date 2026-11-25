#!/usr/bin/env python3
"""Recount everything after recovering the unannotated genes.

Adds the genes found by script 50 to the three Augustus-annotated species and
asks the two questions the whole project turns on:

  1. does the count term flatten — i.e. does the number of set genes a species
     detects stop falling with body mass;
  2. does the decline of the level survive once the gene list is comparable,
     including on the strictly fixed core.

Per-gene value is built exactly as elsewhere in the project (oxphos_thread_lib
.gene_by_lib_matrix): promoter peaks, the gene's major strand, each peak
thresholded at 1 TPM per library, summed per library, then averaged over the
libraries of the tissue.

One deliberate difference for recovered genes: the 'promoter' class filter is
NOT applied to them. That class was assigned by the Augustus annotation we are
correcting, so 'not promoter' there means 'Augustus had no TSS here', which is
the very thing being fixed. A recovered gene is accepted only if a peak lies
within 1000 bp of the miniprot 5' end on the predicted strand AND carries
non-zero fibroblast signal — the second condition is what keeps the negative
controls empty (Emu NPPA has a peak at 103 bp but 0.0 in fibroblasts).
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geneset,
    normalize_symbol,
    nuclear_oxphos,
    peak_to_symbol_de,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
HITS = RES / "gene_recovery" / "recovery_hits.tsv"
OUT = RES / "recovered_rerun"
AUG = ["Turkey", "Peacock", "Emu"]
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))
WINDOW = 1000
THR = 1.0


def fit(y) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    r = stats.linregress(LOGM[ok], y[ok])
    return float(r.rvalue), float(r.pvalue)


def norm(it) -> list[str]:
    return sorted({normalize_symbol(g) for g in it} - {""})


def tissue_mean(sub: pd.DataFrame, cols: list[str]) -> float:
    m = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    m = m.where(m >= THR, 0.0)
    return float(m.sum(axis=0).mean())


def build(species: str, hits: pd.DataFrame) -> tuple[dict[str, pd.Series], list[str]]:
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    out: dict[str, pd.Series] = {}
    cols = {t: lib_cols(d, t) for t in ("fib", "epi")}
    prom = d[d["class"].astype(str).str.lower().eq("promoter") & d["symbol"].ne("")]
    for t in ("fib", "epi"):
        c = cols[t]
        m = prom[c].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        tot = m.sum(axis=1)
        maj = (
            pd.DataFrame({"symbol": prom["symbol"].values, "strand": prom["strand"].values, "t": tot.values})
            .groupby(["symbol", "strand"], observed=True)["t"]
            .sum()
            .reset_index()
            .sort_values("t", ascending=False)
            .drop_duplicates("symbol")
            .set_index("symbol")["strand"]
        )
        keep = prom[prom["strand"].values == prom["symbol"].map(maj).values]
        mk = keep[c].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        mk = mk.where(mk >= THR, 0.0)
        out[t] = mk.assign(symbol=keep["symbol"].values).groupby("symbol")[c].sum().mean(axis=1)

    added: list[str] = []
    if species in AUG:
        h = hits[(hits.species == species) & (~hits.is_control) & (hits.n_peaks_1000 > 0)]
        for _, r in h.iterrows():
            g = r["gene"]
            if g in out["fib"].index:
                continue
            near = d[
                (d["seqnames"] == r["seqname"])
                & (d["strand"] == r["strand"])
                & (((d["start"] + d["end"]) / 2 - r["five_prime"]).abs() <= WINDOW)
            ]
            if near.empty:
                continue
            fib = tissue_mean(near, cols["fib"])
            if fib <= 0:
                continue  # keeps the controls empty
            out["fib"].loc[g] = fib
            out["epi"].loc[g] = tissue_mean(near, cols["epi"])
            added.append(g)
    return out, sorted(added)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hits = pd.read_csv(HITS, sep="\t")
    base = pickle.load(open(RES / "reverify_confounds" / "gene_matrices.pkl", "rb"))

    mats: dict[str, dict[str, pd.Series]] = {}
    added: dict[str, list[str]] = {}
    for sp in SPECIES_ORDER:
        mats[sp], added[sp] = build(sp, hits)
        if added[sp]:
            print(f"[{sp}] добавлено {len(added[sp])} генов: {', '.join(added[sp])}")
    pd.DataFrame(
        [{"species": sp, "n_added": len(v), "genes": ";".join(v)} for sp, v in added.items()]
    ).to_csv(OUT / "genes_added.tsv", sep="\t", index=False)

    wide_after = {
        t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")
    }
    wide_before = {
        "fib": pd.DataFrame({sp: base["CAGE_fib"][sp][0] for sp in SPECIES_ORDER}),
        "epi": pd.DataFrame({sp: base["CAGE_epi"][sp][0] for sp in SPECIES_ORDER}),
    }

    sets = {
        "KEGG hsa00190": norm(load_geneset("kegg_00190")),
        "KEGG hsa04714": norm(load_geneset("kegg_04714")),
        "его OXPHOS 133": norm(load_geneset("ruslan_oxphos_133")),
        "наш ядерный OXPHOS 77": norm(nuclear_oxphos()),
        "контроль: гликолиз": norm(load_geneset("GLYCOLYSIS")),
    }

    rows = []
    for name, genes in sets.items():
        for tis in ("fib", "epi"):
            for tag, w in (("до", wide_before[tis]), ("после", wide_after[tis])):
                hit = w.reindex(genes)
                cnt = hit.notna().sum(axis=0).to_numpy(float)
                mean = np.log10(hit.mean(axis=0).to_numpy(float))
                tot = np.log10(hit.sum(axis=0).to_numpy(float))
                rc, pc = fit(cnt)
                rm, pm = fit(mean)
                rs, ps = fit(tot)
                rows.append(
                    dict(
                        gene_set=name, tissue=tis, вариант=tag,
                        n_detected=";".join(str(int(v)) for v in cnt),
                        r_count=round(rc, 3), p_count=round(pc, 4),
                        r_mean=round(rm, 3), p_mean=round(pm, 4),
                        r_sum=round(rs, 3), p_sum=round(ps, 4),
                    )
                )
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "before_after.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    for tis in ("fib", "epi"):
        print("\n" + "=" * 120)
        print("CAGE фибробласты" if tis == "fib" else "CAGE эпибласт")
        print("=" * 120)
        print(tab[tab.tissue == tis].drop(columns=["tissue"]).to_string(index=False))

    # fixed core, his list, all layers = fib and epi only (RNA is not recoverable here)
    print("\n" + "=" * 120)
    print("ЛЕСТНИЦА ПО ШИРОТЕ на его 133, ядро по CAGE (fib+epi), до и после")
    print("=" * 120)
    grows = []
    for tag, w in (("до", wide_before), ("после", wide_after)):
        idx = sorted(set(w["fib"].index) | set(w["epi"].index))
        ok = None
        for t in ("fib", "epi"):
            m = w[t].reindex(idx)
            g = m.notna() & (m > 0)
            ok = g if ok is None else (ok & g)
        breadth = ok.sum(axis=1)
        genes = norm(load_geneset("ruslan_oxphos_133"))
        for k in (7, 6, 5, 4, 3):
            core = [g for g in genes if breadth.get(g, 0) >= k]
            for t in ("fib", "epi"):
                r, p = fit(np.log10(w[t].reindex(core).mean(axis=0).to_numpy(float)))
                grows.append(dict(вариант=tag, k=k, n_core=len(core), tissue=t,
                                  r=round(r, 3), p=round(p, 4)))
    g = pd.DataFrame(grows)
    g.to_csv(OUT / "gradient_before_after.tsv", sep="\t", index=False)
    print(g.pivot_table(index=["k"], columns=["вариант", "tissue"], values=["n_core", "r", "p"]).to_string())

    print(f"\nwrote {OUT / 'before_after.tsv'}")
    print(f"wrote {OUT / 'gradient_before_after.tsv'}")


if __name__ == "__main__":
    main()
