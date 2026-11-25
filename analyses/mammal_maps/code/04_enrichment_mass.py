#!/usr/bin/env python3
"""Step 2b: per species x set x motif Mantel-Haenszel log odds ratio (GC quintiles, Haldane 0.5),
then PGLS of log OR on log10 adult mass per panel, BH over everything, cross-panel replication (plan_step2.json)."""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN / "code"))
import common_rm  # noqa: E402

plan = json.loads((RUN / "config/plan_step2.json").read_text())
sets = json.loads((RUN / "config/gene_sets.json").read_text())["sets"]
info = json.loads((RUN / "results/02_panel_backgrounds.json").read_text())
OUT = RUN / "results"


def mh_logor(pos, fg, strata):
    """pos: genes x motifs bool; fg: genes bool; strata: genes int. Returns log MH OR per motif."""
    num = np.zeros(pos.shape[1]); den = np.zeros(pos.shape[1])
    for s in np.unique(strata):
        m = strata == s
        f, b = fg & m, ~fg & m
        if f.sum() == 0 or b.sum() == 0:
            continue
        a = pos[f].sum(0) + 0.5; c = f.sum() - pos[f].sum(0) + 0.5
        bb = pos[b].sum(0) + 0.5; d = b.sum() - pos[b].sum(0) + 0.5
        n = a + bb + c + d
        num += a * d / n; den += bb * c / n
    return np.log(num / den)


def gls_many(x, Y, C):
    from scipy import linalg
    X = np.column_stack([np.ones(len(x)), x])
    L = np.linalg.cholesky(C)
    W = linalg.solve_triangular(L, X, lower=True)
    V = linalg.solve_triangular(L, Y, lower=True)
    B = np.linalg.lstsq(W, V, rcond=None)[0]
    E = V - W @ B
    df = len(x) - 2
    se = np.sqrt((E * E).sum(0) / df * np.linalg.inv(W.T @ W)[1, 1])
    t = B[1] / se
    return B[1], t, 2 * stats.t.sf(np.abs(t), df)


lor_rows, test_rows = [], []
motifs = None
null_inputs = {}
for panel in ["ma2016", "gse105160", "wagner"]:
    d = common_rm.load(panel)
    labs = d["species"]
    sci = info["species"][panel]
    have = [l for l in labs if (OUT / "motifs" / f"{sci[l]}.npz").exists()]
    ix = [labs.index(l) for l in have]
    C = d["C"][np.ix_(ix, ix)]
    mass = d["mass"][ix]
    bg = set(info["background"][panel])
    L = {s: [] for s in sets}
    for lab in have:
        z = np.load(OUT / "motifs" / f"{sci[lab]}.npz", allow_pickle=True)
        motifs = [str(m) for m in z["motifs"]]
        g = np.array([str(x) for x in z["genes"]])
        keep = np.isin(g, list(bg))
        P, g, gc = z["presence"][keep], g[keep], z["gc"][keep]
        strata = pd.qcut(gc, 5, labels=False, duplicates="drop")
        for s, members in sets.items():
            fg = np.isin(g, members)
            lo = mh_logor(P, fg, strata)
            L[s].append(lo)
            lor_rows.append(pd.DataFrame({"panel": panel, "species": lab, "set": s, "n_fg": int(fg.sum()), "n_bg": int((~fg).sum()), "motif": motifs, "log_MH_OR": lo}))
    null_inputs[panel] = (mass, C, {s: np.array(L[s]) for s in sets})
    for s in sets:
        Y = np.array(L[s])  # species x motifs
        beta, t, p = gls_many(mass, Y, C)
        r = np.array([stats.pearsonr(mass, Y[:, j])[0] if np.std(Y[:, j]) > 0 else np.nan for j in range(Y.shape[1])])
        test_rows.append(pd.DataFrame({"panel": panel, "set": s, "n_species": len(have), "motif": motifs, "PGLS_beta": beta, "PGLS_t": t, "PGLS_p": p, "pearson_r": r}))
    print(panel, "species with genomes", len(have), "of", len(labs))

T = pd.concat(test_rows, ignore_index=True)
ok = np.isfinite(T.PGLS_p)
T.loc[ok, "BH_q_all"] = multipletests(T.loc[ok, "PGLS_p"], method="fdr_bh")[1]
T["tf"] = T.motif.str.split("|").str[1].str.upper()
pd.concat(lor_rows, ignore_index=True).to_csv(OUT / "04_species_logOR.tsv.gz", sep="\t", index=False)
T.to_csv(OUT / "04_motif_mass_tests.tsv.gz", sep="\t", index=False)

rep = []
for (s, m), g in T.groupby(["set", "motif"]):
    sig = g[g.PGLS_p < 0.05]
    for sign in (1, -1):
        k = sig[np.sign(sig.PGLS_beta) == sign]
        if len(k) >= 2:
            rep.append(dict(set=s, motif=m, tf=g.tf.iloc[0], direction="up_with_mass" if sign > 0 else "down_with_mass",
                            panels=",".join(k.panel), min_q=float(g.BH_q_all.min()),
                            **{f"p_{r.panel}": r.PGLS_p for r in g.itertuples()}, **{f"beta_{r.panel}": r.PGLS_beta for r in g.itertuples()}))
R = pd.DataFrame(rep)
R.to_csv(OUT / "04_replicated_motifs.tsv", sep="\t", index=False)

rng = np.random.default_rng(20261022)
set_names = list(sets)
null_counts = []
for _ in range(1000):
    up = down = None
    for panel, (mass, C, Ys) in null_inputs.items():
        sim = rng.multivariate_normal(np.zeros(len(mass)), C)
        sim = (sim - sim.mean()) / sim.std() * mass.std() + mass.mean()
        b = np.concatenate([gls_many(sim, Ys[s], C)[0] for s in set_names])
        p = np.concatenate([gls_many(sim, Ys[s], C)[2] for s in set_names])
        u = ((p < 0.05) & (b > 0)).astype(int); dn = ((p < 0.05) & (b < 0)).astype(int)
        up = u if up is None else up + u
        down = dn if down is None else down + dn
    null_counts.append(int(((up >= 2) | (down >= 2)).sum()))
null_counts = np.array(null_counts)
obs = len(R.drop_duplicates(["set", "motif"])) if len(R) else 0
rep_null = dict(observed_replicated_pairs=obs, null_median=float(np.median(null_counts)), null_q95=float(np.quantile(null_counts, 0.95)),
                empirical_p=float((1 + np.sum(null_counts >= obs)) / (1 + len(null_counts))))
(OUT / "04_replication_null.json").write_text(json.dumps(rep_null, indent=1))
print("replication null:", rep_null)
pre = T[T.tf.isin(plan["prespecified_from_ruslan"])]
pre.to_csv(OUT / "04_ruslan_prespecified.tsv", sep="\t", index=False)
print("tests", len(T), "q<0.05:", int((T.BH_q_all < 0.05).sum()), "min q", float(T.BH_q_all.min()))
print("replicated (same sign, p<0.05 in >=2 panels):", len(R))
expected = sum(1 for _ in T.groupby(["set", "motif"]))
print("set x motif pairs:", expected)
if len(R):
    print(R.sort_values("min_q").head(40).to_string())
print(pre.sort_values("PGLS_p").head(30)[["panel", "set", "tf", "n_species", "PGLS_beta", "PGLS_p", "BH_q_all", "pearson_r"]].to_string())
print("DONE04")
