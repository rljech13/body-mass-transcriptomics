#!/usr/bin/env python3
"""Wagner / Ba skin-fibroblast panel (PRJNA564062 + PRJNA674672 SF only).

Same frozen 29-gene module and within-table fraction as GSE167580.
Does not pool with GSE167580 or GSE213253.
"""
from __future__ import annotations

import gzip
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    PREFERRED_ALIAS,
    RES,
    normalize_symbol,
    ols_logmass,
    spearman_logmass,
)

ROOT = Path("/path/to/storage/kleiber_birds_gemcat/wagner_sf")
CDNA = ROOT / "cdna"
IDX = ROOT / "salmon_index"
FQ = ROOT / "fastq"
QUANT = ROOT / "quant"
OUT = RES / "wagner_sf_module"
SALMON = "salmon"

MASS_G = {
    "opossum": 104.0,
    "rat": 300.0,
    "guinea_pig": 728.0,
    "rabbit": 1800.0,
    "cat": 4000.0,
    "rhesus": 8000.0,
    "dog": 40000.0,
    "human": 62035.0,
    "sheep": 80000.0,
    "pig": 130000.0,
    "horse": 400000.0,
    "cow": 750000.0,
}
LATIN = {
    "opossum": "M. domestica",
    "rat": "R. norvegicus",
    "guinea_pig": "C. porcellus",
    "rabbit": "O. cuniculus",
    "cat": "F. catus",
    "rhesus": "M. mulatta",
    "dog": "C. familiaris",
    "human": "H. sapiens",
    "sheep": "O. aries",
    "pig": "S. scrofa",
    "horse": "E. caballus",
    "cow": "B. taurus",
}

# Skin fibroblasts only (not ESF). One lab, DMEM high + 10% FBS (Wagner).
SAMPLES = [
    # PRJNA564062
    ("human", "Hs_SF_1", "SRR10075991"),
    ("human", "Hs_SF_2", "SRR10075990"),
    ("human", "Hs_SF_3", "SRR10075989"),
    ("rat", "Rn_SF_1", "SRR10072255"),
    ("rat", "Rn_SF_2", "SRR10072254"),
    ("rat", "Rn_SF_3", "SRR10072253"),
    ("rabbit", "Oc_SF_1", "SRR10072272"),
    ("rabbit", "Oc_SF_2", "SRR10072271"),
    ("rabbit", "Oc_SF_3", "SRR10072260"),
    ("cat", "Fc_SF_1", "SRR10072262"),
    ("cat", "Fc_SF_2", "SRR10072261"),
    ("cat", "Fc_SF_3", "SRR10072259"),
    ("sheep", "Oa_SF_1", "SRR10072265"),
    ("sheep", "Oa_SF_2", "SRR10072264"),
    ("sheep", "Oa_SF_3", "SRR10072263"),
    ("cow", "Bt_SF_1", "SRR10072258"),
    ("cow", "Bt_SF_2", "SRR10072257"),
    ("cow", "Bt_SF_3", "SRR10072256"),
    ("horse", "Ec_SF_1", "SRR10072268"),
    ("horse", "Ec_SF_2", "SRR10072267"),
    ("horse", "Ec_SF_3", "SRR10072266"),
    ("guinea_pig", "Cp_SF_1", "SRR10072252"),
    ("guinea_pig", "Cp_SF_2", "SRR10072270"),
    ("guinea_pig", "Cp_SF_3", "SRR10072269"),
    # PRJNA674672
    ("dog", "Cl_SF_1", "SRR12984020"),
    ("dog", "Cl_SF_2", "SRR12984019"),
    ("dog", "Cl_SF_3", "SRR12984018"),
    ("pig", "Ss_SF_1", "SRR12984014"),
    ("pig", "Ss_SF_2", "SRR12984013"),
    ("pig", "Ss_SF_3", "SRR12984012"),
    ("rhesus", "Mmul_SF_1", "SRR12984010"),
    ("rhesus", "Mmul_SF_2", "SRR12984009"),
    ("rhesus", "Mmul_SF_3", "SRR12984008"),
    ("opossum", "Md_SF_1", "SRR12984017"),
    ("opossum", "Md_SF_2", "SRR12984016"),
    ("opossum", "Md_SF_3", "SRR12984015"),
]


def module_genes() -> set[str]:
    raw = pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"]
    return {normalize_symbol(g) for g in raw if normalize_symbol(g)}


def cdna_fa(species: str) -> Path:
    hits = sorted(CDNA.glob(f"{species}.*cdna.all.fa.gz"))
    if not hits:
        raise FileNotFoundError(species)
    return hits[0]


def parse_header(h: str) -> tuple[str, str, str]:
    tid = h.split()[0].lstrip(">")
    gene = ""
    sym = ""
    m = re.search(r"gene:(\S+)", h)
    if m:
        gene = m.group(1)
    m = re.search(r"gene_symbol:(\S+)", h)
    if m:
        sym = m.group(1)
    return tid, gene, sym


def build_tx2gene(species: str) -> Path:
    dest = ROOT / "tx2gene" / f"{species}.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    rows = []
    with gzip.open(cdna_fa(species), "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            tid, gene, sym = parse_header(line)
            rows.append((tid, gene, normalize_symbol(sym) or "", sym))
    pd.DataFrame(rows, columns=["tx", "gene", "symbol", "raw_symbol"]).to_csv(
        dest, sep="\t", index=False
    )
    return dest


def ensure_index(species: str) -> Path:
    dest = IDX / species
    ok = dest / "complete_ref_lens.bin"
    if dest.exists() and ok.exists():
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    fa = ROOT / "fa" / f"{species}.fa"
    fa.parent.mkdir(parents=True, exist_ok=True)
    if not fa.exists():
        subprocess.check_call(["gzip", "-cd", str(cdna_fa(species))], stdout=open(fa, "wb"))
    log = dest / "index.log"
    with open(log, "w") as lf:
        subprocess.check_call(
            [SALMON, "index", "-t", str(fa), "-i", str(dest), "-k", "23", "-p", "8"],
            stdout=lf,
            stderr=subprocess.STDOUT,
        )
    return dest


def ena_fastq_url(srr: str) -> str:
    """ENA vol1 layout depends on accession length (see ena-docs file-download)."""
    prefix = srr[:6]
    n = len(srr)
    if n == 9:
        return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{srr}/{srr}.fastq.gz"
    if n == 10:
        mid = f"00{srr[-1]}"
    elif n == 11:
        mid = f"0{srr[-2:]}"
    else:
        mid = srr[-3:]
    return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{mid}/{srr}/{srr}.fastq.gz"


def ensure_fastq(srr: str) -> Path:
    FQ.mkdir(parents=True, exist_ok=True)
    dest = FQ / f"{srr}.fastq.gz"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest
    url = ena_fastq_url(srr)
    tmp = dest.with_suffix(".tmp")
    subprocess.check_call(["wget", "-q", "-O", str(tmp), url])
    tmp.rename(dest)
    return dest


def quant_one(species: str, sample: str, srr: str) -> Path:
    out = QUANT / sample
    if (out / "quant.sf").exists() and (out / "quant.sf").stat().st_size > 1000:
        return out
    out.mkdir(parents=True, exist_ok=True)
    idx = ensure_index(species)
    fq = ensure_fastq(srr)
    log = out / "salmon.log"
    with open(log, "w") as lf:
        subprocess.check_call(
            [
                SALMON,
                "quant",
                "-i",
                str(idx),
                "-l",
                "A",
                "-r",
                str(fq),
                "-o",
                str(out),
                "-p",
                "8",
                "--validateMappings",
            ],
            stdout=lf,
            stderr=subprocess.STDOUT,
        )
    try:
        fq.unlink()
    except OSError:
        pass
    return out


def score_quant(species: str, sample: str, genes: set[str]) -> dict:
    t2g = pd.read_csv(build_tx2gene(species), sep="\t")
    q = pd.read_csv(QUANT / sample / "quant.sf", sep="\t")
    q = q.merge(t2g, left_on="Name", right_on="tx", how="left")
    q["symbol"] = q["symbol"].fillna("").map(lambda s: normalize_symbol(s) if s else "")
    # fallback: raw_symbol
    if "raw_symbol" in q.columns:
        miss = q["symbol"].eq("")
        q.loc[miss, "symbol"] = q.loc[miss, "raw_symbol"].fillna("").map(normalize_symbol)
    lib = float(q["TPM"].sum())
    sub = q[q["symbol"].isin(genes)]
    collapsed = sub.groupby("symbol", observed=True)["TPM"].sum()
    present = [g for g in genes if g in collapsed.index]
    num = float(collapsed.reindex(present).fillna(0).sum()) if present else 0.0
    return {
        "module_n_mapped": len(present),
        "module_n_expr": int((collapsed.reindex(present).fillna(0) > 0).sum()) if present else 0,
        "module_frac": num / lib if lib else np.nan,
        "lib_total": lib,
        "n_tx": int(len(q)),
        "n_tx_with_symbol": int((q["symbol"] != "").sum()),
    }


def fit_pair(y: pd.Series, masses: pd.Series) -> dict:
    ols = ols_logmass(y, masses)
    spr = spearman_logmass(y, masses)
    return {**{f"ols_{k}": v for k, v in ols.items()}, "spearman_rho": spr["rho"], "spearman_p": spr["p"]}


def write_results() -> None:
    genes = module_genes()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp, sample, srr in SAMPLES:
        qf = QUANT / sample / "quant.sf"
        if not qf.exists():
            continue
        sc = score_quant(sp, sample, genes)
        rows.append(
            dict(
                study="Wagner_SF",
                bioproject="PRJNA564062" if srr.startswith("SRR100") else "PRJNA674672",
                species=sp,
                latin=LATIN[sp],
                sample=sample,
                srr=srr,
                mass_g=MASS_G[sp],
                **sc,
            )
        )
    lib = pd.DataFrame(rows)
    lib.to_csv(OUT / "phenotype_by_library.tsv", sep="\t", index=False)
    if lib.empty:
        print("no quants yet")
        return
    phen = []
    for sp, g in lib.groupby("species"):
        phen.append(
            dict(
                study="Wagner_SF",
                species=sp,
                latin=g["latin"].iloc[0],
                mass_g=float(g["mass_g"].iloc[0]),
                n_samples=int(len(g)),
                samples=";".join(g["sample"]),
                module_n_mapped=int(g["module_n_mapped"].median()),
                module_frac=float(g["module_frac"].mean()),
            )
        )
    phen = pd.DataFrame(phen).sort_values("mass_g")
    phen.to_csv(OUT / "phenotype_by_species.tsv", sep="\t", index=False)
    y = pd.Series(phen["module_frac"].values)
    m = pd.Series(phen["mass_g"].values)
    fit = fit_pair(y, m)
    (OUT / "ols_spearman_vs_log10mass.tsv")
    pd.DataFrame([{"unit": "wagner_sf_species", "metric": "module_frac", "n_species": int(len(phen)), **fit}]).to_csv(
        OUT / "ols_spearman_vs_log10mass.tsv", sep="\t", index=False
    )
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "date": time.strftime("%Y-%m-%d"),
                "n_quantified": int(len(lib)),
                "n_species": int(phen.species.nunique()),
                "fit": fit,
                "species": phen.to_dict(orient="records"),
            },
            indent=2,
            default=str,
        )
    )
    print(phen.to_string(index=False))
    print(fit)


def cmd_prep() -> None:
    for sp in MASS_G:
        print("tx2gene", sp)
        t = build_tx2gene(sp)
        n = sum(1 for _ in open(t)) - 1
        print(" ", t.name, n)
        print("index", sp)
        ensure_index(sp)


def cmd_quant() -> None:
    for sp, sample, srr in SAMPLES:
        print(time.strftime("%H:%M:%S"), "quant", sample, srr, flush=True)
        try:
            quant_one(sp, sample, srr)
        except Exception as e:
            print("FAIL", sample, e, flush=True)
            continue
        write_results()


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    if arg == "prep":
        cmd_prep()
    elif arg == "quant":
        cmd_quant()
    elif arg == "score":
        write_results()
    else:
        cmd_prep()
        cmd_quant()
        write_results()


if __name__ == "__main__":
    main()
