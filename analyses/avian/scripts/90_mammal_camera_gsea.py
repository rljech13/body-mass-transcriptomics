#!/usr/bin/env python3
"""CAMERA and rank GSEA on the honest mammal matrices.

CAMERA: species means, covariate log10 mass, inter-gene correlation estimated
inside the set and floored at 0 (Wu and Smyth 2012; limma's strict mode when
inter.gene.cor is left empty). Libraries are not samples. Moderated t follows
Smyth 2004 (no trend). limma itself is not installed for this R, so the
statistic is computed here to that specification.

GSEA: gene rank is the correlation of log expression with log10 mass.
The null is 10 000 permutations of the mass vector, re-ranking each time.
Set-membership permutation is not the reported null.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy import stats
from scipy.optimize import brentq
from scipy.special import digamma, polygamma

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES  # noqa: E402

MATDIR = RES / "mammal_honest"
OUT = RES / "mammal_size_sets"
N_PERM = 10_000
RNG = np.random.default_rng(20260923)
PANELS = ("ma2016", "wagner", "gse105160")


def fit_f_dist(s2: np.ndarray, df: float) -> tuple[float, float]:
    """Prior scale and prior df for residual variances. Smyth 2004, no covariate."""
    x = np.asarray(s2, float)
    x = x[np.isfinite(x) & (x > 0)]
    z = np.log(x)
    e = float(z.mean())
    v = float(np.var(z, ddof=1))
    noise = float(polygamma(1, df / 2.0))
    resid = v - noise
    if resid <= 1e-8:
        d0 = 1e6
    else:
        def f(d0: float) -> float:
            return float(polygamma(1, d0 / 2.0) - resid)

        lo, hi = 1e-3, 1e8
        if f(lo) <= 0:
            d0 = lo
        elif f(hi) >= 0:
            d0 = hi
        else:
            d0 = float(brentq(f, lo, hi))
    log_s20 = e - digamma(df / 2.0) + np.log(df / 2.0) + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    return float(np.exp(log_s20)), float(d0)


def moderated_t(y: np.ndarray, logm: np.ndarray) -> tuple[np.ndarray, float]:
    """y: genes x species, already log2. Design is intercept + log10 mass."""
    n = y.shape[1]
    x = np.column_stack([np.ones(n), logm])
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y.T  # 2 x genes
    fitted = (x @ beta).T
    resid = y - fitted
    df = float(n - 2)
    s2 = (resid * resid).sum(axis=1) / df
    s20, d0 = fit_f_dist(s2, df)
    s2_post = (d0 * s20 + df * s2) / (d0 + df)
    unscaled = np.sqrt(xtx_inv[1, 1])
    tstat = beta[1] / (np.sqrt(s2_post) * unscaled)
    return tstat, df + d0


def set_rho(resid: np.ndarray, rows: np.ndarray) -> float:
    """Mean pairwise Pearson correlation of residuals inside the set, floored at 0."""
    r = resid[rows]
    if r.shape[0] < 2:
        return 0.0
    c = np.corrcoef(r)
    m = c.shape[0]
    rho = float((c.sum() - m) / (m * (m - 1)))
    if not np.isfinite(rho):
        return 0.0
    return max(0.0, rho)


def camera_one(tstat: np.ndarray, rows: np.ndarray, rho: float, df: float) -> tuple[float, float]:
    g = len(tstat)
    m = len(rows)
    m2 = g - m
    if m < 1 or m2 < 1:
        return np.nan, np.nan
    mean_in = float(tstat[rows].mean())
    mean_out = float(np.delete(tstat, rows).mean())
    var_t = float(tstat.var(ddof=1))
    vif = 1.0 + (m - 1) * rho
    se = np.sqrt(var_t * (vif / m + 1.0 / m2))
    if se <= 0 or not np.isfinite(se):
        return np.nan, np.nan
    t = (mean_in - mean_out) / se
    p = float(2 * stats.t.sf(abs(t), df))
    return float(t), p


@njit
def es_sorted(pos: np.ndarray, w: np.ndarray, n: int) -> float:
    """Weighted ES from hit positions already sorted along the ranking."""
    nh = pos.shape[0]
    if nh == 0 or nh == n:
        return 0.0
    wsum = 0.0
    for i in range(nh):
        wsum += w[i]
    if wsum <= 0:
        return 0.0
    miss = 1.0 / (n - nh)
    acc = 0.0
    peak = 0.0
    prev = -1
    for i in range(nh):
        gap = pos[i] - prev - 1
        if gap > 0:
            acc -= gap * miss
            av = abs(acc)
            if av > abs(peak):
                peak = acc
        acc += w[i] / wsum
        av = abs(acc)
        if av > abs(peak):
            peak = acc
        prev = pos[i]
    tail = n - 1 - prev
    if tail > 0:
        acc -= tail * miss
        if abs(acc) > abs(peak):
            peak = acc
    return peak


@njit
def gsea_batch(inv_rank: np.ndarray, weights: np.ndarray, hit: np.ndarray) -> np.ndarray:
    """inv_rank, weights: n_perm x n_genes. hit: gene indices of the set."""
    n_perm, n = inv_rank.shape
    nh = hit.shape[0]
    out = np.empty(n_perm)
    pos = np.empty(nh, dtype=np.int64)
    w = np.empty(nh)
    for p in range(n_perm):
        for i in range(nh):
            g = hit[i]
            pos[i] = inv_rank[p, g]
            w[i] = weights[p, g]
        for i in range(1, nh):
            keyp = pos[i]
            keyw = w[i]
            j = i - 1
            while j >= 0 and pos[j] > keyp:
                pos[j + 1] = pos[j]
                w[j + 1] = w[j]
                j -= 1
            pos[j + 1] = keyp
            w[j + 1] = keyw
        out[p] = es_sorted(pos, w, n)
    return out


def gene_corr(logy: np.ndarray, logm: np.ndarray) -> np.ndarray:
    yc = logy - logy.mean(axis=1, keepdims=True)
    xc = logm - logm.mean()
    num = yc @ xc
    den = np.sqrt((yc * yc).sum(axis=1) * float(xc @ xc))
    r = np.zeros(logy.shape[0])
    ok = den > 0
    r[ok] = num[ok] / den[ok]
    return r


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = pd.read_csv(OUT / "sets_vs_mass.tsv", sep="\t")
    cam_rows = []
    gsea_rows = []
    for panel in PANELS:
        wide = pd.read_csv(MATDIR / f"{panel}_M.tsv", sep="\t", index_col=0)
        mass = pd.read_csv(MATDIR / f"{panel}_mass.tsv", sep="\t")
        wide = wide[mass.species.tolist()]
        logm = np.log10(mass.mass_g.to_numpy(float))
        logy = np.log2(wide.to_numpy(float))
        pos = {g: i for i, g in enumerate(wide.index)}
        tstat, df = moderated_t(logy, logm)
        x = np.column_stack([np.ones(len(logm)), logm])
        beta = np.linalg.lstsq(x, logy.T, rcond=None)[0]
        resid = logy - (x @ beta).T
        print(f"{panel}: moderated df {df:.1f}", flush=True)

        sub = sets[sets.panel == panel]
        # observed correlations and a bank of mass permutations
        obs_r = gene_corr(logy, logm)
        # permute masses
        perms = np.empty((N_PERM, len(logm)))
        for i in range(N_PERM):
            perms[i] = RNG.permutation(logm)
        # correlations for all perms, in batches, then ES per set
        # store orders and |r| weights: too big to hold 10000 x genes x 2 of int+float
        # genes ~ 12k, 10000 * 12000 * 8 = 960MB for r, plus argsort. Do batches of 500.
        print(f"  GSEA permutations {N_PERM}", flush=True)
        batch = 250
        # per-set hit indices
        hit_list = []
        meta = []
        for rec in sub.itertuples(index=False):
            genes = [g for g in str(rec.genes).split(";") if g]
            rows = np.array([pos[g] for g in genes if g in pos], dtype=np.int64)
            hit_list.append(rows)
            meta.append(rec)
        n_sets = len(hit_list)
        n_genes = logy.shape[0]
        es_obs = np.empty(n_sets)
        ge = np.zeros(n_sets, dtype=np.int64)
        order = np.argsort(-obs_r).astype(np.int64)
        inv = np.empty(n_genes, dtype=np.int64)
        inv[order] = np.arange(n_genes, dtype=np.int64)
        abs_r = np.abs(obs_r)
        for s, rows in enumerate(hit_list):
            pos = np.sort(inv[rows])
            w = abs_r[rows][np.argsort(inv[rows])]
            es_obs[s] = es_sorted(pos.astype(np.int64), w.astype(float), n_genes)
        yc = logy - logy.mean(axis=1, keepdims=True)
        yc_ss = (yc * yc).sum(axis=1, keepdims=True)
        for start in range(0, N_PERM, batch):
            b = perms[start:start + batch]
            xc = b - b.mean(axis=1, keepdims=True)
            num = yc @ xc.T
            den = np.sqrt(yc_ss * (xc * xc).sum(axis=1)[None, :])
            r = np.zeros_like(num)
            np.divide(num, den, out=r, where=den > 0)
            nb = r.shape[1]
            inv_rank = np.empty((nb, n_genes), dtype=np.int64)
            weights = np.abs(r).T.copy()
            for j in range(nb):
                order = np.argsort(-r[:, j])
                inv_rank[j, order] = np.arange(n_genes, dtype=np.int64)
            for s, rows in enumerate(hit_list):
                null_es = gsea_batch(inv_rank, weights, rows)
                ge[s] += int(np.sum(np.abs(null_es) >= abs(es_obs[s])))
            print(f"  perm {start + nb}/{N_PERM}", flush=True)
        for s, rec in enumerate(meta):
            rows = hit_list[s]
            rho = set_rho(resid, rows)
            t, p = camera_one(tstat, rows, rho, df)
            cam_rows.append(dict(
                panel=panel, catalog=rec.catalog, id=rec.id, name=rec.name, n=rec.n,
                rho=rho, t_camera=t, p_camera=p,
            ))
            gsea_rows.append(dict(
                panel=panel, catalog=rec.catalog, id=rec.id, name=rec.name, n=rec.n,
                es=float(es_obs[s]), p_gsea=float(ge[s] / N_PERM),
            ))
        print(f"  {panel} done", flush=True)

    cam = pd.DataFrame(cam_rows)
    gsea = pd.DataFrame(gsea_rows)
    cam["bh_camera"] = cam.groupby(["panel", "catalog"], group_keys=False)["p_camera"].transform(
        lambda s: stats.false_discovery_control(s.fillna(1.0).to_numpy(), method="bh")
    )
    gsea["bh_gsea"] = gsea.groupby(["panel", "catalog"], group_keys=False)["p_gsea"].transform(
        lambda s: stats.false_discovery_control(s.fillna(1.0).to_numpy(), method="bh")
    )
    cam.to_csv(OUT / "camera.tsv", sep="\t", index=False)
    gsea.to_csv(OUT / "gsea.tsv", sep="\t", index=False)
    print(f"CAMERA BH<0.05 {(cam.bh_camera < 0.05).sum()} / {len(cam)}")
    print(f"GSEA BH<0.05 {(gsea.bh_gsea < 0.05).sum()} / {len(gsea)}")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
