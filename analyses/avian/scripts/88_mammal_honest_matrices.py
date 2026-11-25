#!/usr/bin/env python3
"""Honest gene x species matrices for three mammal fibroblast panels.

Same inputs as scripts 42 and 43. A gene stays only if its species score is
> 0 in every species of that panel. A missing gene is dropped, not zeroed.
Species score is the mean across libraries. Libraries are not the regression
unit. Panels are written separately and never joined.

Ma and GSE105160: per library, 1e6 * gene counts / library size (library size
is the sum of raw counts, mapped or not). Wagner quant.sf is already TPM:
mean TPM across libraries, not divided by the TPM total again.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util

from oxphos_thread_lib import DATA, RES, normalize_symbol  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = Path(__file__).resolve().parent
s42 = _load("s42", HERE / "42_file_panels_ruslan.py")
s37 = _load("s37", HERE / "37_wagner_sf_oxphos.py")

OUT = RES / "mammal_honest"


def strip_ver(x: str) -> str:
    return str(x).split(".", 1)[0]


def species_mean(lib_frames: list[tuple[str, pd.Series]]) -> pd.DataFrame:
    """lib_frames: (species, gene -> score). Mean across libraries, genes as rows."""
    wide = {}
    for i, (sp, s) in enumerate(lib_frames):
        wide[f"{sp}__{i}"] = s
    df = pd.DataFrame(wide)
    sp = pd.Series([c.split("__", 1)[0] for c in df.columns], index=df.columns)
    return df.T.groupby(sp.to_numpy()).mean().T


def honest(mean: pd.DataFrame, species_order: list[str]) -> pd.DataFrame:
    m = mean.reindex(columns=species_order)
    keep = m.notna().all(axis=1) & (m > 0).all(axis=1)
    out = m.loc[keep].sort_index()
    out.index.name = "gene"
    return out


def write_panel(name: str, mat: pd.DataFrame, mass: pd.DataFrame, n_before: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mat.to_csv(OUT / f"{name}_M.tsv", sep="\t")
    mass.to_csv(OUT / f"{name}_mass.tsv", sep="\t", index=False)
    print(
        f"{name}: библиотек сведено в {mass.shape[0]} видов, "
        f"генов с символом {n_before}, честных (во всех видах > 0) {mat.shape[0]}",
        flush=True,
    )


def matrix_ma() -> None:
    id2h = s42.mouse_ens_to_human()
    id2h = {strip_ver(k): v for k, v in id2h.items()}
    aw = s42.ma_paper_aw()
    raw = pd.read_excel(s42.MA_XLSX, sheet_name="A. raw counts")
    raw = raw.rename(columns={raw.columns[0]: "gene"})
    raw["gene"] = raw["gene"].map(strip_ver)
    frames = []
    mass_rows = []
    seen_sp = set()
    for col in raw.columns[1:]:
        m = __import__("re").match(r"^(.+)\.\d", str(col))
        if not m:
            continue
        key = m.group(1)
        paper = s42.MA_PAPER_NAME[key]
        sp = s42.MA_LABEL[key]
        counts = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
        counts.index = raw["gene"].astype(str)
        counts = counts.groupby(level=0).sum()
        lib = float(counts.sum())
        mapped = counts.groupby(counts.index.map(lambda i: id2h.get(i, ""))).sum()
        mapped = mapped[mapped.index != ""]
        score = mapped * (1e6 / lib if lib else np.nan)
        frames.append((sp, score))
        if sp not in seen_sp:
            seen_sp.add(sp)
            mass_rows.append(dict(
                species=sp, mass_g=aw[paper], mass_source="Ma eLife fig1 AW",
                organism=paper,
            ))
    mean = species_mean(frames)
    mass = pd.DataFrame(mass_rows).sort_values("mass_g")
    mat = honest(mean, mass.species.tolist())
    write_panel("ma2016", mat, mass, int(mean.shape[0]))


def matrix_gse() -> None:
    mock = s42.parse_soft_mock()
    ens_h = {strip_ver(k): v for k, v in s42.human_ens_to_symbol().items()}
    ens_m = {strip_ver(k): v for k, v in s42.mouse_ens_to_human().items()}
    files = {}
    for p in s42.COUNTS.glob("*"):
        hit = __import__("re").search(r"GSM\d+", p.name)
        if hit:
            files[hit.group()] = p
    frames = []
    mass_rows = []
    seen = set()
    n_missing = 0
    for rec in mock.itertuples(index=False):
        path = files.get(rec.gsm)
        if path is None:
            n_missing += 1
            print("MISSING", rec.gsm, rec.organism, flush=True)
            continue
        counts = s42.read_count_file(path)
        counts.index = counts.index.map(strip_ver)
        counts = counts.groupby(level=0).sum()
        id2sym = ens_m if rec.organism == "Mus musculus" else ens_h
        lib = float(counts.sum())
        mapped = counts.groupby(counts.index.map(lambda i: id2sym.get(i, ""))).sum()
        mapped = mapped[mapped.index != ""]
        score = mapped * (1e6 / lib if lib else np.nan)
        sp = s42.GSE_LABEL[rec.organism]
        frames.append((sp, score))
        if sp not in seen:
            seen.add(sp)
            mass, src = s42.GSE_MASS[rec.organism]
            mass_rows.append(dict(
                species=sp, mass_g=mass, mass_source=src, organism=rec.organism,
            ))
    if n_missing:
        raise SystemExit(f"GSE105160: нет файлов у {n_missing} mock-библиотек")
    mean = species_mean(frames)
    mass = pd.DataFrame(mass_rows).sort_values("mass_g")
    mat = honest(mean, mass.species.tolist())
    write_panel("gse105160", mat, mass, int(mean.shape[0]))


def matrix_wagner() -> None:
    frames = []
    mass_rows = []
    seen = set()
    for sp, sample, _srr in s37.SAMPLES:
        qf = s37.QUANT / sample / "quant.sf"
        if not qf.exists():
            # Recorded gap: horse has two libraries, Ec_SF_3 was never quantified.
            if sample == "Ec_SF_3":
                print("skip Ec_SF_3: no quant.sf (horse stays at two libraries)", flush=True)
                continue
            raise SystemExit(f"MISSING {qf}")
        t2g = pd.read_csv(s37.build_tx2gene(sp), sep="\t")
        q = pd.read_csv(qf, sep="\t", usecols=["Name", "TPM"])
        q = q.merge(t2g, left_on="Name", right_on="tx", how="left")
        sym = q["symbol"].fillna("").map(lambda s: normalize_symbol(s) if s else "")
        if "raw_symbol" in q.columns:
            miss = sym.eq("")
            sym.loc[miss] = q.loc[miss, "raw_symbol"].fillna("").map(normalize_symbol)
        q = q.assign(symbol=sym)
        q = q[q["symbol"] != ""]
        score = q.groupby("symbol", observed=True)["TPM"].sum()
        frames.append((sp, score))
        if sp not in seen:
            seen.add(sp)
            mass_rows.append(dict(
                species=sp, mass_g=s37.MASS_G[sp],
                mass_source="script 37 MASS_G (Wagner panel)",
                organism=s37.LATIN[sp],
            ))
    mean = species_mean(frames)
    mass = pd.DataFrame(mass_rows).sort_values("mass_g")
    mat = honest(mean, mass.species.tolist())
    write_panel("wagner", mat, mass, int(mean.shape[0]))


def main() -> None:
    matrix_ma()
    matrix_gse()
    matrix_wagner()
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
