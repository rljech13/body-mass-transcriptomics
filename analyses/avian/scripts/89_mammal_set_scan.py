#!/usr/bin/env python3
"""A priori gene-set sums vs log10 mass on the three honest mammal matrices.

Panels are never pooled. Within each panel, BH is within a catalog.
The null is 2000 sets matched on size and on deciles of median log10
abundance across the species of that panel. A gene missing from the
honest universe is dropped, not zeroed.

Frozen r from scripts 42 and 43 are copied beside the honest recomputation.
Those files are not rewritten.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, RES, normalize_symbol  # noqa: E402

MATDIR = RES / "mammal_honest"
OUT = RES / "mammal_size_sets"
CAT = DATA / "genesets" / "mammal_catalogs"
KEGG = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
FROZEN = RES / "file_panels_ruslan"
N_NULL = 2000
MIN_N = 8
GO_MIN, GO_MAX = 8, 50
RNG = np.random.default_rng(20260922)

# Same dual names as script 42. KEGG writes ATP6; matrices carry MT-ATP6.
KEGG_MT = {
    "ATP6": "MT-ATP6", "ATP8": "MT-ATP8",
    "COX1": "MT-CO1", "COX2": "MT-CO2", "COX3": "MT-CO3",
    "CYTB": "MT-CYB",
    "ND1": "MT-ND1", "ND2": "MT-ND2", "ND3": "MT-ND3",
    "ND4": "MT-ND4", "ND4L": "MT-ND4L", "ND5": "MT-ND5", "ND6": "MT-ND6",
}
KEGG_MT_REV = {v: k for k, v in KEGG_MT.items()}

PANELS = ("ma2016", "wagner", "gse105160")


def aliases(symbol: str) -> list[str]:
    raw = str(symbol).upper().strip()
    norm = normalize_symbol(raw)
    out: list[str] = []
    for g in (norm, raw, KEGG_MT.get(raw, ""), KEGG_MT.get(norm, ""), KEGG_MT_REV.get(norm, "")):
        if g and g not in out:
            out.append(g)
    return out


def resolve(genes: list[str], index: set[str]) -> list[str]:
    hit: list[str] = []
    seen: set[str] = set()
    for g in genes:
        for a in aliases(g):
            if a in index and a not in seen:
                hit.append(a)
                seen.add(a)
                break
    return hit


def parse_semi(path: Path) -> list[tuple[str, str, list[str]]]:
    df = pd.read_csv(path, sep="\t")
    rows = []
    for rec in df.itertuples(index=False):
        genes = [x for x in str(rec.genes).split(";") if x.strip()]
        rows.append((str(rec.gs_name), str(rec.gs_name), genes))
    return rows


def load_catalogs() -> dict[str, list[tuple[str, str, list[str]]]]:
    kegg_rows = []
    kegg = pd.read_csv(KEGG, sep="\t")
    for rec in kegg.itertuples(index=False):
        genes = [x for x in str(rec.genes).split(";") if x.strip()]
        kegg_rows.append((str(rec.pathway), str(rec.name), genes))

    go_rows = []
    n_go = 0
    for sid, label, genes in parse_semi(CAT / "gobp.tsv"):
        n_go += 1
        uniq = {normalize_symbol(g) for g in genes}
        uniq.discard("")
        if GO_MIN <= len(uniq) <= GO_MAX:
            go_rows.append((sid, label, genes))

    corum_rows = []
    cor = pd.read_csv(CAT / "corum_humanComplexes.txt", sep="\t", usecols=["complex_id", "complex_name", "subunits_gene_name"])
    for rec in cor.itertuples(index=False):
        genes = [x.strip() for x in str(rec.subunits_gene_name).split(";") if x.strip() and x.strip() not in {"-", "None"}]
        corum_rows.append((f"CORUM:{rec.complex_id}", str(rec.complex_name), genes))

    catalogs = {
        "kegg": kegg_rows,
        "hallmark": parse_semi(CAT / "hallmark.tsv"),
        "reactome": parse_semi(CAT / "reactome.tsv"),
        "corum": corum_rows,
        "gobp": go_rows,
    }
    print(f"GO BP в каталоге {n_go}, в окне {GO_MIN}–{GO_MAX} генов {len(go_rows)}", flush=True)
    for name, rows in catalogs.items():
        print(f"  {name}: {len(rows)} наборов до пересечения со вселенной", flush=True)
    return catalogs


class AbundanceNull:
    def __init__(self, wide: pd.DataFrame, logm: np.ndarray):
        self.genes = np.asarray(wide.index)
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.mat = wide.to_numpy(float)
        self.logm = logm
        yc = logm - logm.mean()
        self.yc = yc
        self.yss = float(yc @ yc)
        med = np.log10(np.median(self.mat, axis=1))
        self.dec = pd.qcut(pd.Series(med).rank(method="first"), 10, labels=False, duplicates="drop").to_numpy()
        self.nbins = int(self.dec.max()) + 1
        self.buckets = [np.flatnonzero(self.dec == k) for k in range(self.nbins)]

    def r_of(self, y: np.ndarray) -> float:
        z = y - y.mean()
        den = np.sqrt(float(z @ z) * self.yss)
        return float((z @ self.yc) / den) if den else np.nan

    def r_batch(self, y: np.ndarray) -> np.ndarray:
        z = y - y.mean(axis=1, keepdims=True)
        den = np.sqrt((z * z).sum(axis=1) * self.yss)
        num = z @ self.yc
        out = np.full(y.shape[0], np.nan)
        ok = den > 0
        out[ok] = num[ok] / den[ok]
        return out

    def score_rows(self, rows: np.ndarray) -> tuple[float, float, np.ndarray]:
        y = np.log10(self.mat[rows].sum(axis=0))
        r = self.r_of(y)
        p = float(stats.pearsonr(self.logm, y)[1]) if np.isfinite(r) else np.nan
        return r, p, y

    def matched(self, rows: np.ndarray, obs_r: float) -> tuple[float, float]:
        if len(rows) < MIN_N or not np.isfinite(obs_r):
            return np.nan, np.nan
        want = np.bincount(self.dec[rows], minlength=self.nbins)
        blocks = []
        for k, need in enumerate(want):
            if need == 0:
                continue
            av = self.buckets[k]
            take = min(int(need), len(av))
            if take == len(av):
                blocks.append(np.broadcast_to(av, (N_NULL, take)).copy())
                continue
            keys = RNG.random((N_NULL, len(av)))
            sel = np.argpartition(keys, take - 1, axis=1)[:, :take]
            blocks.append(av[sel])
        idx = np.concatenate(blocks, axis=1)
        y = np.log10(self.mat[idx].sum(axis=1))
        null = self.r_batch(y)
        null = null[np.isfinite(null)]
        return float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null))


def frozen_anchors() -> pd.DataFrame:
    """Read stored fits. Do not rewrite them. Wagner was not in fits.tsv."""
    fits = pd.read_csv(FROZEN / "fits.tsv", sep="\t")
    rows = []
    for rec in fits.itertuples(index=False):
        panel = {"ma2016": "ma2016", "gse105160_mock": "gse105160"}[rec.panel]
        sid = "hsa00190" if rec.set == "kegg_00190" else "hsa04714"
        rows.append(dict(panel=panel, id=sid, r_frozen=rec.r, p_frozen=rec.p, frozen_from="scripts/42 fits.tsv"))
    wag = pd.read_csv(FROZEN / "wagner_species_scores.tsv", sep="\t")
    for sid, key in (("hsa00190", "kegg_00190"), ("hsa04714", "kegg_04714")):
        g = wag[wag.set == key]
        r = stats.linregress(np.log10(g.mass_g), np.log10(g.sum_tpm))
        rows.append(dict(
            panel="wagner", id=sid, r_frozen=float(r.rvalue), p_frozen=float(r.pvalue),
            frozen_from="scripts/43 wagner_species_scores.tsv",
        ))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    catalogs = load_catalogs()
    frozen = frozen_anchors()
    all_rows = []
    for panel in PANELS:
        wide = pd.read_csv(MATDIR / f"{panel}_M.tsv", sep="\t", index_col=0)
        mass = pd.read_csv(MATDIR / f"{panel}_mass.tsv", sep="\t")
        wide = wide[mass.species.tolist()]
        logm = np.log10(mass.mass_g.to_numpy(float))
        null = AbundanceNull(wide, logm)
        index = set(wide.index)
        print(f"\n{panel}: {wide.shape[0]} генов, {wide.shape[1]} видов", flush=True)
        for cat, sets in catalogs.items():
            rows = []
            n_skip = 0
            for i, (sid, label, genes) in enumerate(sets):
                hit = resolve(genes, index)
                if len(hit) < MIN_N:
                    n_skip += 1
                    continue
                pos = np.array([null.pos[g] for g in hit], dtype=int)
                r, p, _ = null.score_rows(pos)
                pn, nm = null.matched(pos, r)
                rows.append(dict(
                    panel=panel, catalog=cat, id=sid, name=label, n=len(hit),
                    r=r, p_ols=p, p_null=pn, null_median=nm,
                    genes=";".join(hit),
                ))
                if sid in {"hsa00190", "hsa04714"} or (i + 1) % 500 == 0:
                    print(f"  {cat} {i+1}/{len(sets)} {sid} n={len(hit)} r={r:+.3f} нуль {pn:.4f}", flush=True)
            tab = pd.DataFrame(rows)
            if tab.empty:
                print(f"  {cat}: нет наборов с ≥{MIN_N} генами", flush=True)
                continue
            tab["bh"] = stats.false_discovery_control(tab["p_null"].fillna(1.0).to_numpy(), method="bh")
            n_pass = int((tab.bh < 0.05).sum())
            print(f"  {cat}: испытано {len(tab)}, отсеяно по размеру {n_skip}, BH<0.05 {n_pass}", flush=True)
            all_rows.append(tab)
    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv(OUT / "sets_vs_mass.tsv", sep="\t", index=False)

    honest = out[out.id.isin(["hsa00190", "hsa04714"])][
        ["panel", "id", "n", "r", "p_ols", "p_null", "bh"]
    ].rename(columns={"r": "r_honest", "p_ols": "p_honest", "n": "n_honest", "p_null": "p_null_honest", "bh": "bh_honest"})
    anchors = frozen.merge(honest, on=["panel", "id"], how="left")
    anchors.to_csv(OUT / "anchors_frozen_vs_honest.tsv", sep="\t", index=False)
    print("\nякоря, старый счёт рядом с честной вселенной:")
    print(anchors.drop(columns=["frozen_from"], errors="ignore").to_string(index=False))
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
