#!/usr/bin/env python3
"""Phylogenetic calibration of the OXPHOS working-set decline, n = 7.

Reuses the numpy PGLS of script 10 (Newick, Grafen lengths, Brownian VCV, GLS).
The dated tree is the same TimeTree-informed ultrametric Newick already used
there. No ape/phytools.

What is calibrated: the fibroblast CAGE level of the working OXPHOS set against
log10 mass. Two measurements, reported side by side:

  старая панель   55 genes, matrices of script 62 — this is the number the rest
                  of the project quotes (script 67: slope −0.0233, r = −0.806)
  новая панель    47 genes that survive the symmetric re-annotation of script 68

Five numbers per measurement:

1. OLS (the number we already have)
2. PGLS on the dated tree
3. PGLS on Grafen lengths of the same topology (sensitivity to dates)
4. Pagel's λ, profiled 0..1 — description of residual phylogenetic signal
5. Brownian simulation: 10 000 traits on the dated tree, empirical P(|r| >= |r_obs|)
6. Felsenstein independent contrasts as a cross-check (regression through origin)

OU and EVE are not fitted: seven tips do not identify α.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s10 = importlib.import_module("10_pgls_simplicity_S")
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER  # noqa: E402

OUT = RES / "phylo_oxphos"
OXFILE = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
TIPS = ["zebra", "quail", "chicken", "duck", "peacock", "turkey", "emu"]
SP_OF_TIP = {t: t.capitalize() if t != "zebra" else "Zebra" for t in TIPS}
N_SIM = 10_000
RNG = np.random.default_rng(20260921)
TIMETREE = (
    "(emu:110,(zebra:98,(duck:80,(peacock:45,"
    "(turkey:37,(chicken:33,quail:33):4):8):35):18):12);"
)


def species_score(M: pd.DataFrame, genes: list[str]) -> np.ndarray:
    sp = pd.Series([c.split("|")[0] for c in M.columns], index=M.columns)
    per = M.T.groupby(sp.to_numpy()).mean().T
    hit = [g for g in genes if g in per.index]
    y = np.log10(per.reindex(hit)[[SP_OF_TIP[t] for t in TIPS]].mean(axis=0).to_numpy(float))
    return y, hit


def pagel_lambda(y: np.ndarray, x: np.ndarray, V: np.ndarray) -> tuple[float, pd.DataFrame]:
    """Profile ML of Pagel's λ for y ~ 1 + x under V_λ."""
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    rows = []
    grid = np.linspace(0.0, 1.0, 101)
    for lam in grid:
        Vlam = lam * V.copy()
        np.fill_diagonal(Vlam, np.diag(V))
        Vlam = Vlam + np.eye(n) * 1e-9 * np.trace(V) / n
        sign, logdet = np.linalg.slogdet(Vlam)
        if sign <= 0:
            rows.append(dict(lambda_=lam, ll=np.nan))
            continue
        try:
            Vi = np.linalg.inv(Vlam)
        except np.linalg.LinAlgError:
            Vi = np.linalg.pinv(Vlam)
        beta = np.linalg.pinv(X.T @ Vi @ X) @ (X.T @ Vi @ y)
        resid = y - X @ beta
        sigma2 = float(resid @ Vi @ resid) / n
        if sigma2 <= 0:
            rows.append(dict(lambda_=lam, ll=np.nan))
            continue
        ll = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
        rows.append(dict(lambda_=lam, ll=ll, slope=float(beta[1]), sigma2=sigma2))
    tab = pd.DataFrame(rows)
    i = tab["ll"].idxmax()
    return float(tab.loc[i, "lambda_"]), tab


def pic_slope(y: np.ndarray, x: np.ndarray, root: s10.Node, tips: list[str]) -> dict:
    """Felsenstein independent contrasts, regression through the origin."""
    yv = {t: float(v) for t, v in zip(tips, y)}
    xv = {t: float(v) for t, v in zip(tips, x)}
    cy, cx = [], []

    def walk(n: s10.Node) -> tuple[float, float, float]:
        """Return (y, x, extra_length) at this node after reducing children."""
        if not n.children:
            return yv[n.name], xv[n.name], 0.0
        reduced = []
        for c in n.children:
            yc, xc, extra = walk(c)
            reduced.append((yc, xc, c.length + extra))
        while len(reduced) > 1:
            (y1, x1, v1), (y2, x2, v2) = reduced[0], reduced[1]
            cy.append((y1 - y2) / np.sqrt(v1 + v2))
            cx.append((x1 - x2) / np.sqrt(v1 + v2))
            w = 1 / v1 + 1 / v2
            yp = (y1 / v1 + y2 / v2) / w
            xp = (x1 / v1 + x2 / v2) / w
            reduced = [(yp, xp, 1 / w)] + reduced[2:]
        return reduced[0]

    walk(root)
    cy = np.asarray(cy, float)
    cx = np.asarray(cx, float)
    slope = float(np.sum(cx * cy) / np.sum(cx * cx))
    resid = cy - slope * cx
    df = max(len(cx) - 1, 1)
    se = float(np.sqrt(np.sum(resid ** 2) / df / np.sum(cx * cx)))
    t = slope / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan
    r = float(np.corrcoef(cx, cy)[0, 1]) if len(cx) > 1 else np.nan
    return dict(n_contrasts=len(cx), slope=slope, se=se, t=t, p=p, r=r)


def bm_empirical_p(x: np.ndarray, V: np.ndarray, r_obs: float) -> float:
    L = np.linalg.cholesky(V + np.eye(len(V)) * 1e-9 * np.trace(V) / len(V))
    null = np.empty(N_SIM)
    for i in range(N_SIM):
        z = RNG.normal(size=len(x))
        y = L @ z
        null[i] = stats.pearsonr(x, y)[0]
    return float(np.mean(np.abs(null) >= abs(r_obs)))


def fit_one(label: str, y: np.ndarray, x: np.ndarray, V_tt: np.ndarray, V_gr: np.ndarray,
            root_tt: s10.Node) -> dict:
    ols = s10.ols(y, x)
    X = np.column_stack([np.ones(len(y)), x])
    gls_tt = s10.pgls(y, X, V_tt)
    gls_gr = s10.pgls(y, X, V_gr)
    lam, _ = pagel_lambda(y, x, V_tt)
    p_bm = bm_empirical_p(x, V_tt, ols["r"])
    pic = pic_slope(y, x, root_tt, TIPS)
    return dict(
        измерение=label, n=len(y),
        ols_slope=round(ols["slope"], 4), ols_r=round(ols["r"], 3), ols_p=round(ols["p"], 4),
        pgls_tt_slope=round(float(gls_tt["beta"][1]), 4),
        pgls_tt_p=round(float(gls_tt["p"][1]), 4),
        pgls_gr_slope=round(float(gls_gr["beta"][1]), 4),
        pgls_gr_p=round(float(gls_gr["p"][1]), 4),
        lambda_pegel=round(lam, 2),
        p_броуновский=round(p_bm, 4),
        pic_slope=round(pic["slope"], 4), pic_r=round(pic["r"], 3), pic_p=round(pic["p"], 4),
        pic_n=pic["n_contrasts"],
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ox = [l.strip() for l in OXFILE.read_text().splitlines() if l.strip()]
    x = np.log10(np.array([SPECIES[SP_OF_TIP[t]]["mass_g"] for t in TIPS], float))

    root_tt = s10.parse_newick(TIMETREE)
    V_tt = s10.vcv_brownian(root_tt, TIPS)
    pd.DataFrame(V_tt, index=TIPS, columns=TIPS).to_csv(OUT / "vcv_timetree.tsv", sep="\t")
    root_gr = s10.parse_newick(TIMETREE)
    s10.apply_grafen(root_gr)
    V_gr = s10.vcv_brownian(root_gr, TIPS)

    old = pd.read_csv(RES / "wgcna_input" / "M_fib_w1000.tsv", sep="\t", index_col=0)
    y_old, hit_old = species_score(old, ox)
    new = pd.read_csv(RES / "wgcna_input_v2" / "M_fib_common_only_u5000d1000.tsv", sep="\t", index_col=0)
    uf = [l.strip() for l in (RES / "wgcna_input_v2" / "universe_fib_common_only_u5000d1000.txt").read_text().splitlines() if l.strip()]
    y_new, hit_new = species_score(new, [g for g in ox if g in uf])

    print("виды:", ", ".join(TIPS))
    print("x = log10 массы:", np.round(x, 3))
    print(f"старая панель: {len(hit_old)} генов, y =", np.round(y_old, 3))
    print(f"новая панель:  {len(hit_new)} генов, y =", np.round(y_new, 3))

    rows = [
        fit_one(f"старая панель, {len(hit_old)} генов", y_old, x, V_tt, V_gr, root_tt),
        fit_one(f"новая панель, {len(hit_new)} генов", y_new, x, V_tt, V_gr, root_tt),
    ]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "phylo_oxphos.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n", t.to_string(index=False))
    print(f"\nwrote {OUT}/phylo_oxphos.tsv")


if __name__ == "__main__":
    main()
