#!/usr/bin/env python3
"""A5: GEMCAT on ESC/fib — iES1300 primary, Recon3D baseline (recomputed, not Ruslan's outputs)."""
from __future__ import annotations

import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
RES = ROOT / "results"
FIG = ROOT / "figures"
IES = DATA / "ies1300/S1_File/S1_File in sbml format.xml"
RECON = DATA / "Recon3D.xml"
ENV = {**os.environ, "LD_PRELOAD": "/usr/lib/x86_64-linux-gnu/libsqlite3.so.0"}


def load_esc_fib_logfc() -> pd.Series:
    df = pd.read_csv(RES / "chicken_esc_vs_fib_gene_tpm.tsv", sep="\t", index_col=0)
    s = df["logFC_ESC_over_fib"].astype(float)
    s.index = s.index.astype(str).str.upper()
    return s.replace([np.inf, -np.inf], np.nan).dropna()


def run_gemcat(model: Path | str, de_csv: Path, out_csv: Path) -> None:
    cmd = ["gemcat", str(model), str(de_csv), "-e", "logFC", "-o", str(out_csv)]
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True, env=ENV)


def sort_gemcat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 2:
        df.columns = ["metabolite"]
        df["score"] = np.nan
    else:
        df = df.iloc[:, :2]
        df.columns = ["metabolite", "score"]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.dropna(subset=["score"]).sort_values("score", ascending=False)


def expand_to_recon_at(entrez_fc: pd.DataFrame) -> pd.DataFrame:
    import cobra

    model = cobra.io.read_sbml_model(str(RECON))
    entrez_to_at = defaultdict(list)
    for g in model.genes:
        m = re.match(r"^(\d+)_AT\d+$", g.id)
        if m:
            entrez_to_at[m.group(1)].append(g.id)
    rows = []
    for e, v in zip(entrez_fc["gene"].astype(str), entrez_fc["logFC"]):
        for gid in entrez_to_at.get(e, []):
            rows.append((gid, v))
    return pd.DataFrame(rows, columns=["gene", "logFC"])


def main():
    logfc = load_esc_fib_logfc()
    sym2chk = pd.read_csv(DATA / "symbol_to_chicken_entrez.tsv", sep="\t")
    sym2chk["symbol"] = sym2chk["symbol"].astype(str).str.upper()
    sym2hum = pd.read_csv(DATA / "symbol_to_human_entrez_ies_overlap.tsv", sep="\t")
    sym2hum["symbol"] = sym2hum["symbol"].astype(str).str.upper()

    # --- iES1300 primary: only model genes ---
    ies_rows = []
    for _, r in sym2chk.iterrows():
        sym = r["symbol"]
        if sym in logfc.index:
            ies_rows.append((str(int(float(r["chicken_entrez"]))), float(logfc.loc[sym])))
    ies_de = pd.DataFrame(ies_rows, columns=["gene", "logFC"]).drop_duplicates("gene")
    ies_de_path = RES / "esc_fib_logFC_ies1300_chicken_entrez.csv"
    ies_de.to_csv(ies_de_path, index=False)
    print(f"iES1300 DE genes mapped: {len(ies_de)}")

    ies_out = RES / "gemcat_esc_fib_ies1300.csv"
    run_gemcat(IES, ies_de_path, ies_out)
    ies = sort_gemcat(ies_out)
    ies.to_csv(RES / "gemcat_esc_fib_ies1300_sorted.tsv", sep="\t", index=False)
    ies.head(50).to_csv(RES / "gemcat_esc_fib_ies1300_top50.tsv", sep="\t", index=False)
    ies.tail(50).to_csv(RES / "gemcat_esc_fib_ies1300_bottom50.tsv", sep="\t", index=False)

    # --- Recon3D baseline (recomputed) ---
    hum_rows = []
    for _, r in sym2hum.iterrows():
        sym = r["symbol"]
        if sym in logfc.index:
            hum_rows.append((str(int(float(r["human_entrez"]))), float(logfc.loc[sym])))
    hum_de = pd.DataFrame(hum_rows, columns=["gene", "logFC"]).drop_duplicates("gene")
    fixed = expand_to_recon_at(hum_de)
    fixed_path = RES / "esc_fib_logFC_recon3d_fixedID.csv"
    fixed.to_csv(fixed_path, index=False)
    print(f"Recon3D DE Entrez: {len(hum_de)}; AT rows: {len(fixed)}")

    recon_out = RES / "gemcat_esc_fib_recon3d.csv"
    run_gemcat(RECON, fixed_path, recon_out)
    recon = sort_gemcat(recon_out)
    recon.to_csv(RES / "gemcat_esc_fib_recon3d_sorted.tsv", sep="\t", index=False)
    recon.head(50).to_csv(RES / "gemcat_esc_fib_recon3d_top50.tsv", sep="\t", index=False)

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, df, title in [
        (axes[0], ies, "iES1300 PRIMARY"),
        (axes[1], recon, "Recon3D baseline"),
    ]:
        top = pd.concat([df.head(12), df.tail(12)])
        ax.barh(top.metabolite.astype(str), top.score)
        ax.set_title(title)
        ax.set_xlabel("GEMCAT score (ESC/fib logFC)")
    fig.suptitle("ESC vs fibroblast GEMCAT — chicken iES1300 vs human Recon3D")
    fig.savefig(FIG / "fig_A5_gemcat_ies_vs_recon.png", dpi=160)
    plt.close(fig)

    # overlap of top metabolite names (string match; imperfect across models)
    top_ies = set(ies.head(50).metabolite.astype(str))
    top_rec = set(recon.head(50).metabolite.astype(str))
    pd.Series(sorted(top_ies & top_rec)).to_csv(RES / "gemcat_top50_name_overlap_ies_recon.txt", index=False, header=False)

    with open(RES / "A5_model_policy_runtime.json", "w") as f:
        import json

        json.dump(
            {
                "primary": "iES1300",
                "baseline": "Recon3D",
                "ruslan_outputs_reused": False,
                "ies1300_input_genes": len(ies_de),
                "recon3d_input_entrez": len(hum_de),
                "recon3d_input_AT": len(fixed),
                "top50_metabolite_name_overlap": len(top_ies & top_rec),
            },
            f,
            indent=2,
        )
    print("Done A5 dual GEMCAT.")


if __name__ == "__main__":
    main()
