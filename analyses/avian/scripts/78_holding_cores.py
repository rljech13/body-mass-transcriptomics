#!/usr/bin/env python3
"""Holding cores inside kept KEGG maps and inside each WGCNA module.

The 16 OXPHOS genes were a complete-case set whose sum tracks mass. This
asks the same question inside four maps the user kept, and inside each
fibroblast WGCNA module, without mixing maps or modules.

Rule (frozen before looking): among genes detected in all seven species,
keep those whose own Pearson r with log10 mass has the same sign as the
parent set. Score log10(mean TPM) of that core vs log10 mass against 2000
abundance-matched nulls. That filter uses the mass axis, so it is a
diagnostic, not a new a priori map. Calibration: apply it to hsa00190 and
measure overlap with the known 16.

Literature subsets (official lists or named cassettes, not a peel):
  thermogenesis minus OXPHOS, Hippo cassette, thyroid-specific receptors.

Module hubs: kME_own ≥ 0.7, not selected by mass.

Do not put fraction on the log10 axis. Do not glue maps or modules.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
)

OUT = RES / "holding_cores"
M_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
KEGG = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
MODS = RES / "wgcna_out_v2" / "gene_modules_fib.tsv"
LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
N_NULL = 2000
RNG = np.random.default_rng(20260921)
MIN_N = 8
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]
YC = LOGM - LOGM.mean()
YSS = float(YC @ YC)

HIPPO_CASSETTE = {
    "STK3", "STK4", "SAV1", "LATS1", "LATS2", "MOB1A", "MOB1B",
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "NF2",
    "WWC1", "FRMD6",
}
THYROID_SPECIFIC = {
    "THRA", "THRB", "DIO1", "DIO2", "DIO3",
    "SLC16A2", "SLC16A10", "SLCO1C1",
}
THYROID_COREQ = {
    "NCOR1", "NCOA1", "NCOA2", "NCOA3", "RXRA", "RXRB", "RXRG",
    "HDAC1", "HDAC2", "HDAC3", "CREBBP", "EP300", "SIN3A",
    "KAT2A", "KAT2B", "MED1", "MED4", "MED12", "MED12L", "MED13",
    "MED13L", "MED14", "MED16", "MED17", "MED24", "MED27", "MED30",
}
PI3K_AKT = {
    "AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB", "PIK3CD",
    "PIK3R1", "PIK3R2", "PIK3R3", "MTOR", "RHEB", "TSC2", "PDPK1",
    "FOXO1", "GSK3B", "BAD", "CASP9",
}
MAPK_RAS = {
    "HRAS", "KRAS", "NRAS", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
    "SRC", "STAT1", "MYC", "CCND1", "HIF1A", "TP53", "MDM2",
}


def r_of_y(y: np.ndarray) -> float:
    z = y - y.mean()
    den = np.sqrt(float(z @ z) * YSS)
    return float((z @ YC) / den) if den else np.nan


def fit(y) -> tuple[float, float]:
    r = stats.linregress(LOGM, np.asarray(y, float))
    return float(r.rvalue), float(r.pvalue)


class AbundanceNull:
    def __init__(self, wide: pd.DataFrame):
        self.genes = np.asarray(wide.index)
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.mat = wide.to_numpy(float)
        chick = np.log10(self.mat[:, SPECIES_ORDER.index("Chicken")].clip(1e-12))
        self.dec = pd.qcut(pd.Series(chick).rank(method="first"), 10,
                           labels=False, duplicates="drop").to_numpy()
        self.buckets = [np.flatnonzero(self.dec == k) for k in range(10)]

    def rows_of(self, genes: list[str]) -> np.ndarray:
        return np.array([self.pos[g] for g in genes if g in self.pos], dtype=int)

    def r_mean(self, rows: np.ndarray) -> float:
        return r_of_y(np.log10(self.mat[rows].mean(axis=0)))

    def y_mean(self, genes: list[str]) -> np.ndarray:
        rows = self.rows_of(genes)
        return self.mat[rows].mean(axis=0)

    def matched(self, rows: np.ndarray, obs_r: float) -> tuple[float, float]:
        if len(rows) < MIN_N or not np.isfinite(obs_r):
            return np.nan, np.nan
        want = np.bincount(self.dec[rows], minlength=10)
        null = np.empty(N_NULL)
        for i in range(N_NULL):
            parts = []
            for k, need in enumerate(want):
                if need == 0:
                    continue
                av = self.buckets[k]
                parts.append(RNG.choice(av, size=min(int(need), len(av)), replace=False))
            pick = np.concatenate(parts)
            null[i] = self.r_mean(pick)
        null = null[np.isfinite(null)]
        return float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null))


def kegg_genes(catalog: pd.DataFrame, pid: str) -> list[str]:
    rec = catalog.loc[catalog.pathway == pid].iloc[0]
    genes = [normalize_symbol(x) for x in str(rec.genes).split(";") if x.strip()]
    return [g for g in genes if g]


def thyroid_role(g: str) -> str:
    if g in THYROID_SPECIFIC:
        return "рецептор / дейодиназа / транспорт T3"
    if g in THYROID_COREQ or g.startswith("MED"):
        return "корегулятор транскрипции"
    if g in PI3K_AKT:
        return "PI3K–AKT–mTOR"
    if g in MAPK_RAS:
        return "MAPK / RAS / эффект"
    if g.startswith("ATP1") or g.startswith("ATP2") or g.startswith("SLC") or g.startswith("PLC"):
        return "ионный / метаболический эффектор"
    return "прочее на карте"


def hippo_role(g: str) -> str:
    if g in HIPPO_CASSETTE:
        return "кассета Hippo (MST–LATS–YAP/TAZ–TEAD)"
    if g.startswith("WNT") or g.startswith("FZD") or g.startswith("DVL") or g in {
        "CTNNB1", "AXIN1", "AXIN2", "APC", "APC2", "TCF7", "TCF7L1", "TCF7L2", "LEF1", "GSK3B",
    }:
        return "WNT на карте"
    if g.startswith("TGF") or g.startswith("BMP") or g.startswith("SMAD") or g.startswith("GDF"):
        return "TGF/BMP на карте"
    return "адгезия / прочее"


def gene_r(null: AbundanceNull, g: str) -> tuple[float, float]:
    y = np.log10(null.mat[null.pos[g]])
    return fit(y)


def score_set(null: AbundanceNull, genes: list[str]) -> dict:
    rows = null.rows_of(genes)
    n = int(len(rows))
    if n == 0:
        return dict(n=0, r=np.nan, p=np.nan, p_нуля=np.nan, медиана_нуля=np.nan,
                    y=[np.nan] * 7)
    y = null.mat[rows].mean(axis=0)
    r, p = fit(np.log10(y))
    pn, nm = null.matched(rows, r) if n >= MIN_N else (np.nan, np.nan)
    return dict(n=n, r=r, p=p, p_нуля=pn, медиана_нуля=nm, y=y.tolist())


def same_sign_core(null: AbundanceNull, genes: list[str], parent_r: float) -> list[str]:
    if not np.isfinite(parent_r) or parent_r == 0:
        return []
    want = np.sign(parent_r)
    keep = []
    for g in genes:
        rg, _ = gene_r(null, g)
        if np.isfinite(rg) and np.sign(rg) == want:
            keep.append(g)
    return keep


def loo_table(null: AbundanceNull, genes: list[str], parent_r: float) -> pd.DataFrame:
    rows = []
    parent_rows = null.rows_of(genes)
    parent_y = null.mat[parent_rows].mean(axis=0)
    for g in genes:
        rg, pg = gene_r(null, g)
        rest = [x for x in genes if x != g]
        r_loo, _ = fit(np.log10(null.y_mean(rest)))
        d_abs = abs(parent_r) - abs(r_loo)
        rows.append(dict(
            gene=g, r_гена=rg, p_гена=pg,
            r_без=r_loo, d_abs=d_abs,
            держит=d_abs > 0,
            chick_tpm=float(null.mat[null.pos[g], SPECIES_ORDER.index("Chicken")]),
        ))
    t = pd.DataFrame(rows).sort_values("d_abs", ascending=False)
    t["parent_r"] = parent_r
    t["parent_mean"] = float(parent_y.mean())
    return t


def parse_kegg(catalog: pd.DataFrame) -> dict[str, list[str]]:
    ox = kegg_genes(catalog, "hsa00190")
    th = kegg_genes(catalog, "hsa04714")
    return {
        "hsa00190": ox,
        "hsa04714": th,
        "hsa03008": kegg_genes(catalog, "hsa03008"),
        "hsa04390": kegg_genes(catalog, "hsa04390"),
        "hsa04919": kegg_genes(catalog, "hsa04919"),
        "04714_minus_00190": sorted(set(th) - set(ox)),
        "hippo_cassette": sorted(HIPPO_CASSETTE),
        "thyroid_specific": sorted(THYROID_SPECIFIC),
        "set16": list(SIXTEEN),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(M_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
    det = raw.notna() & (raw > 0)
    wide = raw.loc[det.all(axis=1)]
    print(f"вселенная {len(wide)} генов, CAGE во всех семи", flush=True)
    null = AbundanceNull(wide)
    catalog = pd.read_csv(KEGG, sep="\t")
    lists = parse_kegg(catalog)

    labels = {
        "hsa00190": "OXPHOS (калибровка)",
        "hsa04714": "термогенез",
        "hsa03008": "биогенез рибосом",
        "hsa04390": "Hippo",
        "hsa04919": "щитовидка, сигнал",
        "04714_minus_00190": "термогенез минус OXPHOS",
        "hippo_cassette": "кассета Hippo",
        "thyroid_specific": "рецепторы / дейодиназы / MCT",
        "set16": "известные 16",
    }

    parent_rows = []
    gene_rows = []
    y_rows = []

    map_ids = ["hsa00190", "hsa04714", "hsa03008", "hsa04390", "hsa04919"]
    extra_ids = ["set16", "04714_minus_00190", "hippo_cassette", "thyroid_specific"]

    def present(genes: list[str]) -> list[str]:
        return [g for g in genes if g in null.pos]

    for sid in map_ids + extra_ids:
        hit = present(lists[sid])
        sc = score_set(null, hit)
        parent_rows.append(dict(
            id=sid, набор=labels[sid], kind="parent",
            n=sc["n"], r=sc["r"], p=sc["p"],
            p_нуля=sc["p_нуля"], медиана_нуля=sc["медиана_нуля"],
            overlap16=len(set(hit) & set(SIXTEEN)),
        ))
        for sp, val in zip(SPECIES_ORDER, sc["y"]):
            y_rows.append(dict(id=sid, kind="parent", species=sp, mean_tpm=val))

        if sid in map_ids and sc["n"] >= 3:
            loo = loo_table(null, hit, sc["r"])
            loo["id"] = sid
            loo["набор"] = labels[sid]
            if sid == "hsa04919":
                loo["роль"] = loo["gene"].map(thyroid_role)
            elif sid == "hsa04390":
                loo["роль"] = loo["gene"].map(hippo_role)
            elif sid == "hsa04714":
                loo["роль"] = np.where(loo["gene"].isin(set(lists["hsa00190"])),
                                       "OXPHOS на карте", "не OXPHOS")
            elif sid == "hsa00190":
                loo["роль"] = np.where(loo["gene"].isin(SIXTEEN), "из 16", "остальной OXPHOS")
            else:
                loo["роль"] = "биогенез рибосом"
            gene_rows.append(loo)

            core = same_sign_core(null, hit, sc["r"])
            csc = score_set(null, core)
            parent_rows.append(dict(
                id=sid, набор=labels[sid] + " — ядро того же знака",
                kind="same_sign",
                n=csc["n"], r=csc["r"], p=csc["p"],
                p_нуля=csc["p_нуля"], медиана_нуля=csc["медиана_нуля"],
                overlap16=len(set(core) & set(SIXTEEN)),
            ))
            for sp, val in zip(SPECIES_ORDER, csc["y"]):
                y_rows.append(dict(id=sid, kind="same_sign", species=sp, mean_tpm=val))
            print(
                f"  {sid:12s} parent n={sc['n']:3d} r={sc['r']:+.3f} нуль {sc['p_нуля']:.3f}"
                f"  | core n={csc['n']:3d} r={csc['r']:+.3f} нуль {csc['p_нуля']:.3f}"
                f"  ∩16={len(set(core) & set(SIXTEEN))}",
                flush=True,
            )

    cores = pd.DataFrame(parent_rows)
    genes = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    ys = pd.DataFrame(y_rows)

    ox_core = same_sign_core(null, present(lists["hsa00190"]),
                             float(cores.loc[(cores.id == "hsa00190") & (cores.kind == "parent"), "r"].iloc[0]))
    print(
        f"\nкалибровка 00190 same-sign: {len(ox_core)} генов, "
        f"из 16 попало {len(set(ox_core) & set(SIXTEEN))}/16, "
        f"лишних {len(set(ox_core) - set(SIXTEEN))}",
        flush=True,
    )

    gm = pd.read_csv(MODS, sep="\t")
    mrows = []
    mgenes = []
    for mod, gtab in gm.groupby("module"):
        if mod == "grey":
            continue
        hit = [x for x in gtab["gene"] if x in null.pos]
        if len(hit) < MIN_N:
            continue
        kme = gtab.set_index("gene")["kME_own"]
        sc = score_set(null, hit)
        core = same_sign_core(null, hit, sc["r"])
        csc = score_set(null, core)
        hubs = [g for g in hit if g in kme.index and kme[g] >= 0.7]
        hsc = score_set(null, hubs)
        mrows.append(dict(
            модуль=mod, n_parent=sc["n"], r_parent=sc["r"], p_нуля_parent=sc["p_нуля"],
            n_core=csc["n"], r_core=csc["r"], p_нуля_core=csc["p_нуля"],
            n_hub=hsc["n"], r_hub=hsc["r"], p_нуля_hub=hsc["p_нуля"],
            overlap16_core=len(set(core) & set(SIXTEEN)),
            overlap16_hub=len(set(hubs) & set(SIXTEEN)),
        ))
        for g in hit:
            rg, pg = gene_r(null, g)
            mgenes.append(dict(
                модуль=mod, gene=g, r_гена=rg, p_гена=pg,
                kME=float(kme[g]) if g in kme.index else np.nan,
                same_sign=g in core, hub=g in hubs,
            ))
        if sc["p_нуля"] < 0.05 or csc["p_нуля"] < 0.05 or abs(sc["r"]) >= 0.7 or abs(csc["r"]) >= 0.85:
            print(
                f"  mod {mod:16s} n={sc['n']:3d} r={sc['r']:+.3f} нуль {sc['p_нуля']:.3f}"
                f"  | core n={csc['n']:3d} r={csc['r']:+.3f} нуль {csc['p_нуля']:.3f}"
                f"  | hub n={hsc['n']:3d} r={hsc['r']:+.3f}",
                flush=True,
            )

    mt = pd.DataFrame(mrows).sort_values("r_core")
    mt["bh_core"] = stats.false_discovery_control(
        mt["p_нуля_core"].fillna(1.0).to_numpy(), method="bh")
    mt["bh_hub"] = stats.false_discovery_control(
        mt["p_нуля_hub"].fillna(1.0).to_numpy(), method="bh")
    mt["bh_parent"] = stats.false_discovery_control(
        mt["p_нуля_parent"].fillna(1.0).to_numpy(), method="bh")

    map_cores = cores[cores.kind == "same_sign"].copy()
    map_cores["bh"] = stats.false_discovery_control(
        map_cores["p_нуля"].fillna(1.0).to_numpy(), method="bh")
    cores = cores.merge(map_cores[["id", "kind", "bh"]], on=["id", "kind"], how="left")

    chick = SPECIES_ORDER.index("Chicken")
    extra = []
    extra_y = []
    extra_genes = []
    for sid in map_ids:
        hit = present(lists[sid])
        tpm = {g: float(null.mat[null.pos[g], chick]) for g in hit}
        top = sorted(hit, key=lambda g: -tpm[g])[:16]
        sc = score_set(null, top)
        extra.append(dict(
            id=sid, набор=labels[sid] + " — 16 самых обильных у курицы",
            kind="abund16",
            n=sc["n"], r=sc["r"], p=sc["p"],
            p_нуля=sc["p_нуля"], медиана_нуля=sc["медиана_нуля"],
            overlap16=len(set(top) & set(SIXTEEN)),
        ))
        extra_genes.append(pd.DataFrame(dict(
            id=sid, kind="abund16", gene=top,
            chick_tpm=[tpm[g] for g in top],
        )))
        for sp, val in zip(SPECIES_ORDER, sc["y"]):
            extra_y.append(dict(id=sid, kind="abund16", species=sp, mean_tpm=val))
        if "ACTB" in hit and sid in {"hsa04390", "hsa04919", "hsa04714"}:
            rest = [g for g in hit if g != "ACTB"]
            scb = score_set(null, rest)
            extra.append(dict(
                id=sid, набор=labels[sid] + " без ACTB",
                kind="drop_actb",
                n=scb["n"], r=scb["r"], p=scb["p"],
                p_нуля=scb["p_нуля"], медиана_нуля=scb["медиана_нуля"],
                overlap16=len(set(rest) & set(SIXTEEN)),
            ))
            for sp, val in zip(SPECIES_ORDER, scb["y"]):
                extra_y.append(dict(id=sid, kind="drop_actb", species=sp, mean_tpm=val))

    hub16_rows = []
    for mod, gtab in gm.groupby("module"):
        if mod == "grey":
            continue
        hit = [x for x in gtab["gene"] if x in null.pos]
        if len(hit) < 16:
            continue
        kme = gtab.set_index("gene")["kME_own"]
        top = sorted(hit, key=lambda g: -float(kme.get(g, 0)))[:16]
        sc = score_set(null, top)
        hub16_rows.append(dict(
            модуль=mod, n=16, r=sc["r"], p=sc["p"], p_нуля=sc["p_нуля"],
            overlap16=len(set(top) & set(SIXTEEN)),
            гены=";".join(top),
        ))
    h16 = pd.DataFrame(hub16_rows).sort_values("r")
    h16["bh"] = stats.false_discovery_control(
        h16["p_нуля"].fillna(1.0).to_numpy(), method="bh")

    cores = pd.concat([cores, pd.DataFrame(extra)], ignore_index=True)
    ys = pd.concat([ys, pd.DataFrame(extra_y)], ignore_index=True)
    pd.concat(extra_genes, ignore_index=True).to_csv(
        OUT / "map_abund16_genes.tsv", sep="\t", index=False)
    h16.to_csv(OUT / "module_hub16_vs_mass.tsv", sep="\t", index=False)

    cores.to_csv(OUT / "map_cores_vs_mass.tsv", sep="\t", index=False)
    genes.to_csv(OUT / "map_genes_loo.tsv", sep="\t", index=False)
    ys.to_csv(OUT / "map_species_means.tsv", sep="\t", index=False)
    mt.to_csv(OUT / "module_cores_vs_mass.tsv", sep="\t", index=False)
    pd.DataFrame(mgenes).to_csv(OUT / "module_genes.tsv", sep="\t", index=False)

    thy = genes[genes.id == "hsa04919"].copy()
    print("\nщитовидка, сигнал — гены complete-case:")
    for rec in thy.sort_values("r_гена", ascending=False).itertuples():
        print(f"  {rec.gene:12s} r={rec.r_гена:+.3f}  {rec.роль}  {'держит' if rec.держит else ''}")

    print(f"\nмодули: same-sign BH<0.05: {(mt.bh_core < 0.05).sum()} / {len(mt)}")
    print(f"модули: hub     BH<0.05: {(mt.bh_hub < 0.05).sum()} / {len(mt)}")
    print(f"модули: hub-16  BH<0.05: {(h16.bh < 0.05).sum()} / {len(h16)}")
    print(f"карты:  same-sign BH<0.05: {(map_cores.bh < 0.05).sum()} / {len(map_cores)}")
    print("\nabund16 / без ACTB:")
    for rec in cores[cores.kind.isin(["abund16", "drop_actb"])].itertuples():
        print(f"  {rec.kind:10s} {rec.id:12s} n={rec.n:3d} r={rec.r:+.3f} нуль {rec.p_нуля:.3f} ∩16={rec.overlap16}")
    print("\nhub-16 модулей |r|≥0.65:")
    for rec in h16[h16.r.abs() >= 0.65].itertuples():
        print(f"  {rec.модуль:16s} r={rec.r:+.3f} нуль {rec.p_нуля:.3f} BH={rec.bh:.3f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
