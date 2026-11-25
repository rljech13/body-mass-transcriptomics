#!/usr/bin/env python3
"""Build 7-species gene correspondence and recompute CAGE fib+epi.

Genes come from what we can actually download/read:
  - RefSeq GTF on the assemblies we already have (turkey, emu)
  - Ensembl GTF peacock AIIM (mostly IDs, few names)
  - NCBI gene_info + gene_orthologs (chicken hub)
  - Ensembl Compara REST for OXPHOS/module/glycolysis
  - Ruslan DE + Dropbox CAGE chicken_genes / geneId (native tables)

CAGE phenotype uses the same recipe as GEO RNA: transcriptome fraction,
plus Y1, soft-core (≥5/7) and hard-core (7/7). Frozen set/module files
are not rewritten; aliases are collapsed only at match time.
"""
from __future__ import annotations

import gzip
import json
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    EXP,
    FIG,
    PREFERRED_ALIAS,
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    is_augustus_id,
    lib_cols,
    load_geneset,
    normalize_symbol,
    nuclear_oxphos,
    ols_logmass,
    peak_to_symbol_de,
    spearman_logmass,
)

NCBI = DATA / "ncbi"
ANNOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/ref/annot")
OUT = RES / "gene_correspondence"
CAGE_OUT = RES / "cage_geo_style_recompute"

TAX = {
    "Human": 9606,
    "Chicken": 9031,
    "Turkey": 9103,
    "Emu": 8790,
    "Peacock": 9049,
    "Duck": 8839,
    "Quail": 93934,
    "Zebra": 59729,
}
SP_TO_TAX = {k: TAX[k] for k in SPECIES_ORDER}
ENSEMBL_SP = {
    "Zebra": "taeniopygia_guttata",
    "Quail": "coturnix_japonica",
    "Chicken": "gallus_gallus",
    "Duck": "anas_platyrhynchos",
    "Peacock": "pavo_cristatus",
    "Turkey": "meleagris_gallopavo",
    "Emu": "dromaius_novaehollandiae",
}
GTF_PATHS = {
    "Turkey": ANNOT / "turkey_GCF_905368555.1/ncbi_dataset/data/GCF_905368555.1/genomic.gtf",
    "Emu": ANNOT / "emu_GCF_036370855.1/ncbi_dataset/data/GCF_036370855.1/genomic.gtf",
    "Peacock": ANNOT / "Pavo_cristatus.AIIM_Pcri_1.0.gtf.gz",
}


def canonical(sym: str) -> str:
    return normalize_symbol(sym)


def collapse_set(genes: set[str]) -> set[str]:
    return {canonical(g) for g in genes if canonical(g)}


def parse_gtf_genes(path: Path) -> pd.DataFrame:
    openf = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with openf(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attr = p[8]
            d = dict(re.findall(r'(\w+)\s+"([^"]+)"', attr))
            d.update({k: v.strip() for k, v in re.findall(r"(\w+)=([^;]+)", attr)})
            gid = d.get("gene_id") or d.get("ID") or ""
            raw = d.get("gene_name") or d.get("Name") or d.get("gene") or ""
            rows.append(
                {
                    "native_id": gid,
                    "symbol_raw": raw,
                    "symbol": canonical(raw) if raw else "",
                    "biotype": d.get("gene_biotype") or d.get("biotype") or "",
                    "chrom": p[0],
                    "start": int(p[3]),
                    "end": int(p[4]),
                    "strand": p[6],
                }
            )
    return pd.DataFrame(rows)


def load_gene_info_birds() -> pd.DataFrame:
    path = NCBI / "gene_info_birds8.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=["#tax_id", "GeneID", "Symbol", "Synonyms", "type_of_gene"])
    df = df.rename(columns={"#tax_id": "tax_id"})
    df["tax_id"] = df["tax_id"].astype(int)
    df["symbol"] = df["Symbol"].map(canonical)
    return df


def load_ncbi_orthologs(info: pd.DataFrame) -> pd.DataFrame:
    path = NCBI / "gene_orthologs.gz"
    if not path.exists() or path.stat().st_size < 1000:
        return pd.DataFrame()
    bird = {TAX[s] for s in SPECIES_ORDER}
    hub = {TAX["Chicken"], TAX["Human"]}
    id2sym = {}
    if not info.empty:
        id2sym = dict(zip(info["GeneID"].astype(str), info["symbol"]))
    rows = []
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                t1, g1, rel, t2, g2 = int(p[0]), p[1], p[2], int(p[3]), p[4]
            except ValueError:
                continue
            if "ortholog" not in rel.lower():
                continue
            # chicken/human hub ↔ a panel bird (skip human–human noise)
            if not ((t1 in hub and t2 in bird) or (t2 in hub and t1 in bird)):
                continue
            rows.append((t1, g1, rel, t2, g2))
    if not rows:
        return pd.DataFrame()
    odf = pd.DataFrame(rows, columns=["tax1", "gene1", "relationship", "tax2", "gene2"])
    tax2sp = {v: k for k, v in TAX.items()}
    odf["sp1"] = odf["tax1"].map(tax2sp)
    odf["sp2"] = odf["tax2"].map(tax2sp)
    odf["sym1"] = odf["gene1"].map(id2sym).fillna("")
    odf["sym2"] = odf["gene2"].map(id2sym).fillna("")
    return odf


def fetch_compara(symbols: list[str]) -> pd.DataFrame:
    """Chicken-symbol → per-species Compara orthologs (one REST call per gene)."""
    rows = []
    for i, sym in enumerate(symbols):
        url = (
            "https://rest.ensembl.org/homology/symbol/gallus_gallus/"
            f"{sym}?content-type=application/json;type=orthologues"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "kleiber-birds/22"})
        data = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=40) as resp:
                    data = json.loads(resp.read().decode())
                break
            except urllib.error.HTTPError as e:
                if e.code == 400:
                    data = {"data": []}
                    break
                time.sleep(1.5 * (attempt + 1))
            except Exception:
                time.sleep(1.5 * (attempt + 1))
        if data is None:
            rows.append({"chicken_symbol": sym, "compara_ok": 0})
            continue
        hits = []
        for block in data.get("data", []):
            for h in block.get("homologies", []):
                tgt = h.get("target") or {}
                src = h.get("source") or {}
                sp = (tgt.get("species") or "").lower()
                hits.append(
                    {
                        "chicken_symbol": src.get("symbol") or sym,
                        "chicken_id": src.get("id"),
                        "ensembl_species": sp,
                        "target_symbol": (tgt.get("symbol") or "").upper(),
                        "target_id": tgt.get("id"),
                        "homology_type": h.get("type"),
                    }
                )
        if not hits:
            rows.append({"chicken_symbol": sym, "compara_ok": 1, "n_homologs": 0})
        else:
            for h in hits:
                h["compara_ok"] = 1
                h["n_homologs"] = len(hits)
                rows.append(h)
        if (i + 1) % 20 == 0:
            print(f"  Compara {i+1}/{len(symbols)}", flush=True)
        time.sleep(0.08)
    cdf = pd.DataFrame(rows)
    inv = {v: k for k, v in ENSEMBL_SP.items()}
    if "ensembl_species" in cdf.columns:
        cdf["species"] = cdf["ensembl_species"].map(inv)
    return cdf


def de_and_cage_symbols() -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    """Observed symbols per species + native_id → symbol maps."""
    present: dict[str, set[str]] = {}
    native_maps: dict[str, dict[str, str]] = {}
    for sp, meta in SPECIES.items():
        de = pd.read_excel(RUSLAN / meta["de"], usecols=lambda c: c in {"id", "genes", "geneId", "transcriptId", "chicken_genes"})
        cage = pd.read_csv(EXP / meta["file"], sep="\t", usecols=lambda c: c in {"id", "geneId", "transcriptId", "chicken_genes"})
        syms: set[str] = set()
        nmap: dict[str, str] = {}

        def add(native: str, raw: str) -> None:
            s = canonical(raw)
            if not s:
                return
            syms.add(s)
            key = str(native).strip()
            if key and key.lower() not in {"nan", "none", ""}:
                nmap[key.upper()] = s

        for _, r in de.iterrows():
            raw = r["genes"] if "genes" in r and not is_augustus_id(str(r["genes"])) else ""
            if not raw and "chicken_genes" in r:
                raw = r["chicken_genes"]
            if not raw and "geneId" in r and not is_augustus_id(str(r["geneId"])):
                raw = r["geneId"]
            for keycol in ("transcriptId", "geneId", "id"):
                if keycol in r:
                    add(r[keycol], raw)
            if raw:
                syms.add(canonical(str(raw)))
        if "chicken_genes" in cage.columns:
            for native, raw in zip(cage.get("geneId", []), cage["chicken_genes"]):
                add(native, raw)
        if "geneId" in cage.columns:
            for g in cage["geneId"].astype(str):
                s = canonical(g)
                if s:
                    syms.add(s)
                    nmap[g.upper()] = s
        present[sp] = {s for s in syms if s}
        native_maps[sp] = nmap
        print(f"  observed {sp}: {len(present[sp])} symbols, {len(nmap)} native IDs", flush=True)
    return present, native_maps


def fit_pair(y: pd.Series, masses: pd.Series) -> dict:
    ols = ols_logmass(y, masses)
    spr = spearman_logmass(y, masses)
    return {**{f"ols_{k}": v for k, v in ols.items()}, "spearman_rho": spr["rho"], "spearman_p": spr["p"]}


def promoter_gene_matrix(df: pd.DataFrame, tissue: str) -> tuple[pd.DataFrame, pd.Series]:
    cols = lib_cols(df, tissue)
    if not cols:
        return pd.DataFrame(), pd.Series(dtype=float)
    sub = df[df["class"].astype(str).str.lower().eq("promoter")].copy()
    mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib_total = mat.sum(axis=0)
    keep = sub["symbol"].astype(str).str.len() > 0
    tmp = mat.loc[keep].copy()
    tmp["_sym"] = sub.loc[keep, "symbol"].map(canonical).values
    tmp = tmp[tmp["_sym"] != ""]
    gene_mat = tmp.groupby("_sym", observed=True)[cols].sum()
    return gene_mat, lib_total


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CAGE_OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    frozen_set = nuclear_oxphos()
    frozen_mod = set(pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"].str.upper())
    defs_frozen = {
        "set": frozen_set,
        "module": {g for g in frozen_mod if g},
        "GLYCOLYSIS": load_geneset("GLYCOLYSIS"),
    }
    defs = {k: collapse_set(v) for k, v in defs_frozen.items()}
    interest = sorted(set.union(*defs.values()))
    print(f"canonical set {len(defs['set'])} (frozen nuclear {len(frozen_set)}); "
          f"module {len(defs['module'])}; glycolysis {len(defs['GLYCOLYSIS'])}", flush=True)

    # --- official gene lists from GTF ---
    print("=== GTF gene lists ===", flush=True)
    gtf_syms: dict[str, set[str]] = {}
    for sp, path in GTF_PATHS.items():
        gdf = parse_gtf_genes(path)
        gdf.to_csv(OUT / f"genes_{sp.lower()}_from_gtf.tsv", sep="\t", index=False)
        named = set(gdf.loc[gdf["symbol"] != "", "symbol"])
        gtf_syms[sp] = named
        print(f"  {sp} GTF genes={len(gdf)} named={len(named)}", flush=True)

    # --- observed symbols in the tables we actually quantify ---
    print("=== DE + CAGE observed symbols ===", flush=True)
    observed, native_maps = de_and_cage_symbols()
    for sp, nmap in native_maps.items():
        pd.Series(nmap, name="symbol").rename_axis("native_id").reset_index().to_csv(
            OUT / f"native_to_symbol_{sp.lower()}.tsv", sep="\t", index=False
        )

    # --- NCBI gene_info / orthologs if present ---
    print("=== NCBI ===", flush=True)
    info = load_gene_info_birds()
    ncbi_syms: dict[str, set[str]] = {}
    if not info.empty:
        tax2sp = {v: k for k, v in TAX.items()}
        info["species"] = info["tax_id"].map(tax2sp)
        for sp, sub in info.groupby("species"):
            ncbi_syms[sp] = set(sub.loc[sub["symbol"] != "", "symbol"])
            print(f"  gene_info {sp}: {len(ncbi_syms[sp])} symbols", flush=True)
    odf = load_ncbi_orthologs(info)
    if odf.empty:
        print("  gene_orthologs not ready — correspondence uses symbols + Compara", flush=True)
    else:
        odf.to_csv(OUT / "ncbi_orthologs_birds.tsv", sep="\t", index=False)
        print(f"  ortholog rows kept: {len(odf)}", flush=True)

    # --- Ensembl Compara for the genes we test ---
    print("=== Ensembl Compara (interest genes) ===", flush=True)
    compara = fetch_compara(interest)
    compara.to_csv(OUT / "ensembl_compara_interest.tsv", sep="\t", index=False)
    compara_by_sp: dict[str, set[str]] = {sp: set() for sp in SPECIES_ORDER}
    if "species" in compara.columns:
        one2one = compara[compara["homology_type"].astype(str).str.contains("ortholog_one2one", na=False)]
        for sp, sub in one2one.groupby("species"):
            if sp in compara_by_sp:
                compara_by_sp[sp] = set(sub["chicken_symbol"].map(canonical)) | set(
                    sub["target_symbol"].map(canonical)
                )
                # keep chicken symbols that have a 1:1 in that species
                compara_by_sp[sp] = set(sub["chicken_symbol"].map(canonical))
        print("  Compara 1:1 interest genes:", {sp: len(compara_by_sp[sp]) for sp in SPECIES_ORDER}, flush=True)

    # --- chicken-hub correspondence: one row per canonical symbol ---
    print("=== correspondence table ===", flush=True)
    universe = set()
    for d in (gtf_syms, observed, ncbi_syms):
        for s in d.values():
            universe |= s
    universe |= set(interest)
    # chicken gene_info as hub backbone if present
    if "Chicken" in ncbi_syms:
        universe |= ncbi_syms["Chicken"]
    universe = sorted(g for g in universe if g)

    recs = []
    for g in universe:
        rec = {"symbol": g, "in_oxphos_set": int(g in defs["set"]), "in_module": int(g in defs["module"]),
               "in_glycolysis": int(g in defs["GLYCOLYSIS"])}
        n = 0
        sources = []
        for sp in SPECIES_ORDER:
            flags = []
            if g in observed.get(sp, set()):
                flags.append("table")
            if g in gtf_syms.get(sp, set()):
                flags.append("gtf")
            if g in ncbi_syms.get(sp, set()):
                flags.append("ncbi")
            if g in compara_by_sp.get(sp, set()):
                flags.append("compara1to1")
            rec[f"{sp}_present"] = int(bool(flags))
            rec[f"{sp}_source"] = ",".join(flags)
            if flags:
                n += 1
                sources.extend(flags)
        rec["n_species"] = n
        recs.append(rec)
    corr = pd.DataFrame(recs).sort_values(["in_oxphos_set", "in_module", "n_species", "symbol"], ascending=[False, False, False, True])
    corr.to_csv(OUT / "gene_correspondence_7spp.tsv", sep="\t", index=False)

    ox = corr[corr.in_oxphos_set == 1]
    ox.to_csv(OUT / "oxphos_correspondence_7spp.tsv", sep="\t", index=False)
    print(f"  full table {len(corr)} symbols; OXPHOS canonical {len(ox)}", flush=True)
    print(ox[["symbol", "n_species"] + [f"{sp}_present" for sp in SPECIES_ORDER]].to_string(index=False), flush=True)

    cov_rows = []
    for dname, genes in defs.items():
        for sp in SPECIES_ORDER:
            have = genes & observed.get(sp, set())
            cov_rows.append(
                dict(
                    definition=dname,
                    species=sp,
                    n_canonical=len(genes),
                    n_in_table=len(have),
                    n_in_gtf=len(genes & gtf_syms.get(sp, set())),
                    n_in_ncbi=len(genes & ncbi_syms.get(sp, set())),
                    n_compara1to1=len(genes & compara_by_sp.get(sp, set())),
                    n_any=len(genes & (observed.get(sp, set()) | gtf_syms.get(sp, set()) | ncbi_syms.get(sp, set()))),
                )
            )
    pd.DataFrame(cov_rows).to_csv(OUT / "interest_coverage_by_source.tsv", sep="\t", index=False)

    # --- CAGE recompute ---
    print("=== CAGE fib + epi ===", flush=True)
    masses = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})
    tables = {}
    present_cage = {tiss: {} for tiss in ("fib", "epi")}
    gene_mats = {tiss: {} for tiss in ("fib", "epi")}
    lib_totals = {tiss: {} for tiss in ("fib", "epi")}
    for sp, meta in SPECIES.items():
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        id2 = peak_to_symbol_de(RUSLAN / meta["de"])
        df = annotate_table(df, id2)
        df["symbol"] = df["symbol"].map(canonical)
        tables[sp] = df
        for tiss in ("fib", "epi"):
            gmat, tot = promoter_gene_matrix(df, tiss)
            gene_mats[tiss][sp] = gmat
            lib_totals[tiss][sp] = tot
            present_cage[tiss][sp] = set(gmat.index)
            print(f"  {sp} {tiss}: {gmat.shape[0]} genes, {len(tot)} libs", flush=True)

    cores = {}
    softs = {}
    for tiss in ("fib", "epi"):
        for dname, genes in defs.items():
            counts = defaultdict(int)
            for sp in SPECIES_ORDER:
                for g in present_cage[tiss][sp] & genes:
                    counts[g] += 1
            cores[(tiss, dname)] = {g for g, n in counts.items() if n == 7}
            softs[(tiss, dname)] = {g for g, n in counts.items() if n >= 5}

    phen_rows, lib_rows, cov2 = [], [], []
    for tiss in ("fib", "epi"):
        for dname, genes in defs.items():
            core = cores[(tiss, dname)]
            soft = softs[(tiss, dname)]
            for sp in SPECIES_ORDER:
                gmat = gene_mats[tiss][sp]
                tot = lib_totals[tiss][sp]
                detected = sorted(present_cage[tiss][sp] & genes)
                soft_g = sorted(present_cage[tiss][sp] & soft)
                core_g = sorted(core)
                gmean = gmat.mean(axis=1) if len(gmat) else pd.Series(dtype=float)
                den = float(tot.mean()) if len(tot) else np.nan

                def frac(names: list[str]) -> float:
                    if not names or not den:
                        return float("nan")
                    return float(gmean.reindex(names).fillna(0).sum()) / den

                def y1(names: list[str]) -> float:
                    if not names:
                        return float("nan")
                    return float(gmean.reindex(names).dropna().mean())

                phen_rows.append(
                    dict(
                        tissue=tiss,
                        definition=dname,
                        species=sp,
                        mass_g=masses[sp],
                        n_detected=len(detected),
                        n_soft=len(soft_g),
                        n_core=len(core_g),
                        n_soft_def=len(soft),
                        n_core_def=len(core),
                        frac_detected=frac(detected),
                        frac_soft=frac(soft_g),
                        frac_core=frac(core_g),
                        Y1_detected=y1(detected),
                        Y1_soft=y1(soft_g),
                        Y1_core=y1(core_g),
                        lib_total=den,
                    )
                )
                cov2.append(
                    dict(
                        tissue=tiss,
                        definition=dname,
                        species=sp,
                        n_canonical=len(genes),
                        n_detected=len(detected),
                        n_soft_def=len(soft),
                        n_core_def=len(core),
                    )
                )
                for col in gmat.columns:
                    ser = gmat[col]
                    lib_den = float(tot[col])
                    lib_rows.append(
                        dict(
                            tissue=tiss,
                            definition=dname,
                            species=sp,
                            library=col,
                            mass_g=masses[sp],
                            frac_detected=float(ser.reindex(detected).fillna(0).sum()) / lib_den if lib_den else np.nan,
                            frac_soft=float(ser.reindex(soft_g).fillna(0).sum()) / lib_den if lib_den else np.nan,
                            frac_core=float(ser.reindex(core_g).fillna(0).sum()) / lib_den if lib_den and core_g else np.nan,
                            Y1_detected=float(ser.reindex(detected).dropna().mean()) if detected else np.nan,
                        )
                    )

    phen = pd.DataFrame(phen_rows)
    phen.to_csv(CAGE_OUT / "phenotype_by_species.tsv", sep="\t", index=False)
    lib = pd.DataFrame(lib_rows)
    lib.to_csv(CAGE_OUT / "phenotype_by_library.tsv", sep="\t", index=False)
    pd.DataFrame(cov2).to_csv(CAGE_OUT / "gene_coverage.tsv", sep="\t", index=False)

    fit_rows = []
    for tiss in ("fib", "epi"):
        for dname in defs:
            sub = phen[(phen.tissue == tiss) & (phen.definition == dname)].set_index("species").reindex(SPECIES_ORDER)
            for metric in ("frac_detected", "frac_soft", "frac_core", "Y1_detected", "Y1_soft", "Y1_core"):
                fit = fit_pair(sub[metric], masses)
                fit_rows.append(dict(tissue=tiss, definition=dname, unit="species_mean_n7", metric=metric, **fit))
            subl = lib[(lib.tissue == tiss) & (lib.definition == dname)]
            xmass = pd.Series(subl["mass_g"].values)
            for metric in ("frac_detected", "frac_soft", "frac_core", "Y1_detected"):
                y = pd.Series(subl[metric].values)
                y.index = xmass.index = range(len(y))
                fit = fit_pair(y, xmass)
                fit_rows.append(dict(tissue=tiss, definition=dname, unit=f"library_n{len(y)}", metric=metric, **fit))
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(CAGE_OUT / "ols_spearman_vs_log10mass.tsv", sep="\t", index=False)

    # figures
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.0), constrained_layout=True)
    x = np.log10(masses.reindex(SPECIES_ORDER).values)
    for i, tiss in enumerate(("fib", "epi")):
        for j, dname in enumerate(("set", "module", "GLYCOLYSIS")):
            ax = axes[i, j]
            sub = phen[(phen.tissue == tiss) & (phen.definition == dname)].set_index("species").reindex(SPECIES_ORDER)
            y = sub["frac_detected"]
            ax.scatter(x, y.values, c="#1d3557" if tiss == "fib" else "#2a9d8f", s=50)
            fit = ols_logmass(y, masses)
            if np.isfinite(fit["slope"]):
                xx = np.linspace(x.min(), x.max(), 40)
                ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#1d3557" if tiss == "fib" else "#2a9d8f", lw=1.4)
            for sp, xi, yi in zip(SPECIES_ORDER, x, y.values):
                ax.text(xi, yi, sp[:3], fontsize=7)
            ax.set_title(f"{tiss} {dname} frac\nOLS p={fit['p']:.3g} ρ={spearman_logmass(y, masses)['rho']:.2f}")
            ax.set_xlabel("log10 mass (g)")
            if j == 0:
                ax.set_ylabel("gene-set / all promoters")
    fig.suptitle("CAGE native tables: OXPHOS promoter fraction vs body mass", fontsize=11)
    fig.savefig(FIG / "fig_cage_fraction_fib_epi_vs_mass.png", dpi=160)
    fig.savefig(CAGE_OUT / "fig_cage_fraction_fib_epi_vs_mass.png", dpi=160)
    plt.close()

    show = fits[
        (fits.unit == "species_mean_n7")
        & (fits.metric.isin(["frac_detected", "frac_soft", "frac_core", "Y1_detected"]))
    ]
    print("\n=== fits species n=7 ===", flush=True)
    print(show.to_string(index=False), flush=True)
    print("\n=== coverage fib set ===", flush=True)
    print(phen[(phen.tissue == "fib") & (phen.definition == "set")][
        ["species", "n_detected", "n_soft", "n_core", "frac_detected", "Y1_detected"]
    ].to_string(index=False), flush=True)

    summary = {
        "date": "2026-08-14",
        "canonical_set_n": len(defs["set"]),
        "frozen_set_n": len(frozen_set),
        "alias_note": "ATP5A1/ATP5F1A etc collapsed at match time; frozen tsv not rewritten",
        "correspondence_n": int(len(corr)),
        "oxphos_n_species_ge5": int((ox.n_species >= 5).sum()),
        "oxphos_n_species_eq7": int((ox.n_species == 7).sum()),
        "cores": {f"{t}_{d}": sorted(cores[(t, d)]) for t, d in cores},
        "softs_n": {f"{t}_{d}": len(softs[(t, d)]) for t, d in softs},
        "primary_fits": show.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (CAGE_OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote {OUT} and {CAGE_OUT}", flush=True)


if __name__ == "__main__":
    main()
