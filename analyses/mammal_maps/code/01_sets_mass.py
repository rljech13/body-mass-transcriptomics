#!/usr/bin/env python3
"""Step 1 of Ruslan's approach in mammalian fibroblasts: his gene sets versus adult mass.
Fixes relative to the bird version: identical gene composition across species within a panel,
phylogenetic GLS, and an expression-matched random-set null. Plan frozen in config/plan.json."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg, stats

RUN = Path(__file__).resolve().parents[1]
RES = RUN.parent
GI = RES / "growth_regulators_2026-10-21_v1/inputs"
plan = json.loads((RUN / "config/plan.json").read_text())
sets = json.loads((RUN / "config/gene_sets.json").read_text())["sets"]
adm = pd.read_csv(RES / "nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv", sep="\t")
rng = np.random.default_rng(plan["seed"])


def gls(x, y, C):
    X = np.column_stack([np.ones(len(x)), x])
    L = np.linalg.cholesky(C)
    W = linalg.solve_triangular(L, X, lower=True)
    v = linalg.solve_triangular(L, y, lower=True)
    b = np.linalg.lstsq(W, v, rcond=None)[0]
    df = len(x) - 2
    e = v - W @ b
    se = np.sqrt((e @ e / df) * np.linalg.inv(W.T @ W)[1, 1])
    t = b[1] / se
    return float(b[1]), float(se), float(t), float(2 * stats.t.sf(abs(t), df))


def gls_t_many(x, Y, C):
    """t statistics of the mass slope for many response columns at once (same design)."""
    X = np.column_stack([np.ones(len(x)), x])
    L = np.linalg.cholesky(C)
    W = linalg.solve_triangular(L, X, lower=True)
    V = linalg.solve_triangular(L, Y, lower=True)
    B = np.linalg.lstsq(W, V, rcond=None)[0]
    E = V - W @ B
    df = len(x) - 2
    se = np.sqrt((E * E).sum(0) / df * np.linalg.inv(W.T @ W)[1, 1])
    return B[1] / se


def load(panel):
    s = pd.read_csv(GI / f"{panel}_full.samples.tsv", sep="\t")
    s = s.rename(columns={s.columns[0]: "sample"})
    f = pd.read_csv(GI / f"{panel}_full.features.tsv", sep="\t")
    f = f.rename(columns={f.columns[0]: "feature"}).set_index("feature")
    lin = pd.read_csv(GI / f"{panel}_full.linear.tsv.gz", sep="\t", index_col=0)
    lg = pd.read_csv(GI / f"{panel}_full.log.tsv", sep="\t", index_col=0)
    cov = np.load(GI / f"{panel}_full.cov.npz", allow_pickle=True)
    sp = [str(x) for x in cov["species"]]
    a = adm[(adm.panel == panel) & (adm.status == "admitted")].set_index("label")
    keep = [x for x in sp if x in a.index]
    ix = [sp.index(x) for x in keep]
    C = cov["C"][np.ix_(ix, ix)]
    C = C / np.mean(np.diag(C))
    smap = s.set_index("sample").species
    def by_species(M):
        M = M[[c for c in M.columns if c in smap.index]]
        return M.T.groupby(smap).mean().T[keep]
    fs = f[f.symbols.notna() & ~f.symbols.astype(str).str.contains(";")]
    fs = fs[~fs.symbols.duplicated(keep=False)]
    lin_s = by_species(lin).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
    log_s = by_species(lg).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
    # identical composition: measured with positive linear expression in every admitted species
    ok = (lin_s > 0).all(axis=1) & np.isfinite(log_s).all(axis=1)
    lin_s, log_s = lin_s[ok], log_s.loc[lin_s.index[ok]]
    mass = np.log10(a.loc[keep, "accepted_mass_g"].to_numpy(float))
    return dict(species=keep, C=C, mass=mass, lin=lin_s, log=log_s)


rows, species_rows, gene_rows = [], [], []
for panel in plan["panels"]:
    d = load(panel)
    Z = d["log"].sub(d["log"].mean(1), axis=0).div(d["log"].std(1), axis=0)
    Z = Z[np.isfinite(Z).all(1)]
    genes_all = Z.index.to_numpy()
    mean_expr = d["log"].loc[genes_all].mean(1).to_numpy()
    dec = pd.qcut(mean_expr, plan["null_expression_bins"], labels=False, duplicates="drop")
    bins = {b: np.flatnonzero(dec == b) for b in np.unique(dec)}
    for set_name, members in sets.items():
        g = [x for x in members if x in Z.index]
        if len(g) < plan["min_genes"]:
            rows.append(dict(panel=panel, set=set_name, n_species=len(d["species"]), n_genes_list=len(members), n_genes_used=len(g), status="too_few_genes"))
            continue
        z_score = Z.loc[g].mean(0).to_numpy()
        log_sum = np.log10(d["lin"].loc[g].sum(0).to_numpy())
        r_z, p_r_z = stats.pearsonr(d["mass"], z_score)
        r_s, p_r_s = stats.pearsonr(d["mass"], log_sum)
        b, se, t, p = gls(d["mass"], z_score, d["C"])
        bs, ses, ts, ps = gls(d["mass"], log_sum, d["C"])
        # expression-matched null: same number of genes drawn from the same expression bins
        gi = np.flatnonzero(np.isin(genes_all, g))
        need = pd.Series(dec[gi]).value_counts()
        draws = np.empty((plan["n_null"], len(gi)), dtype=int)
        col = 0
        for bnum, k in need.items():
            pool = bins[bnum]
            draws[:, col:col + k] = np.stack([rng.choice(pool, k, replace=False) for _ in range(plan["n_null"])])
            col += k
        Zv = Z.to_numpy()
        null_scores = Zv[draws].mean(1).T  # species x draws
        t_null = gls_t_many(d["mass"], null_scores, d["C"])
        p_emp = (1 + np.sum(np.abs(t_null) >= abs(t))) / (1 + len(t_null))
        rows.append(dict(panel=panel, set=set_name, n_species=len(d["species"]), n_genes_list=len(members), n_genes_used=len(g),
                         status="tested", pearson_r_z=r_z, pearson_p_z=p_r_z, PGLS_beta_z=b, PGLS_se_z=se, PGLS_t_z=t, PGLS_p_z=p,
                         null_emp_p_z=p_emp, null_median_t=float(np.median(t_null)),
                         pearson_r_logsum=r_s, pearson_p_logsum=p_r_s, PGLS_beta_logsum=bs, PGLS_p_logsum=ps))
        for s_, m_, zz, ls in zip(d["species"], d["mass"], z_score, log_sum):
            species_rows.append(dict(panel=panel, set=set_name, species=s_, log10_mass=m_, mean_z=zz, log10_sum_fixed=ls))
        tg = gls_t_many(d["mass"], Z.loc[g].to_numpy().T, d["C"])
        dfree = len(d["species"]) - 2
        for gene, tt in zip(g, tg):
            gene_rows.append(dict(panel=panel, set=set_name, gene=gene, PGLS_t=float(tt), PGLS_p=float(2 * stats.t.sf(abs(tt), dfree))))
    print(panel, len(d["species"]), "species", Z.shape[0], "genes in fixed background")

out = RUN / "results"
out.mkdir(exist_ok=True)
R = pd.DataFrame(rows)
tested = R.status.eq("tested")
from statsmodels.stats.multitest import multipletests
R.loc[tested, "BH_q_emp_all"] = multipletests(R.loc[tested, "null_emp_p_z"], method="fdr_bh")[1]
R.to_csv(out / "01_sets_mass.tsv", sep="\t", index=False)
pd.DataFrame(species_rows).to_csv(out / "01_species_scores.tsv", sep="\t", index=False)
pd.DataFrame(gene_rows).to_csv(out / "01_gene_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(R[["panel", "set", "n_species", "n_genes_used", "pearson_r_z", "PGLS_beta_z", "PGLS_p_z", "null_emp_p_z", "BH_q_emp_all", "pearson_r_logsum", "PGLS_p_logsum"]].round(4).to_string())
print("DONE01")
