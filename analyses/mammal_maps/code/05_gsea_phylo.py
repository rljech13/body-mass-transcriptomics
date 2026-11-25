#!/usr/bin/env python3
"""Step 3: genome-wide phylogenetic GSEA in three mammalian fibroblast panels with cross-panel concordance (plan_step3.json)."""
import json, os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
from pathlib import Path
import numpy as np
import pandas as pd
from numba import njit, prange
from scipy import linalg, stats
from statsmodels.stats.multitest import multipletests

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN / "code"))
import common_rm  # noqa: E402

plan = json.loads((RUN / "config/plan_step3.json").read_text())
GS = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/genesets")
N = 1000
rng = np.random.default_rng(plan["seed"])

catalog = {}
for name, f in [("GO_BP", "mammal_catalogs/gobp.tsv"), ("Reactome", "mammal_catalogs/reactome.tsv"), ("Hallmark", "mammal_catalogs/hallmark.tsv")]:
    for r in pd.read_csv(GS / f, sep="\t").itertuples():
        catalog[f"{name}:{r.gs_name}"] = r.genes.split(";")
for r in pd.read_csv(GS / "kegg_hsa_all_pathways.tsv", sep="\t").itertuples():
    catalog[f"KEGG:{r.pathway} {r.name}"] = str(r.genes).split(";")
catalog["Ruslan:OXPHOS_133"] = (GS / "ruslan_oxphos_133.txt").read_text().split()
catalog["Ruslan:Thermogenesis_229"] = (GS / "ruslan_thermo_229.txt").read_text().split()
catalog["Ruslan:MT142"] = (RUN.parent / "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt").read_text().split()


def gene_t(Zt, x, Linv):
    """Zt: species x genes. PGLS t of each gene on x."""
    X = Linv @ np.column_stack([np.ones(len(x)), x])
    V = Linv @ Zt
    B, *_ = np.linalg.lstsq(X, V, rcond=None)
    E = V - X @ B
    df = len(x) - 2
    se = np.sqrt((E * E).sum(0) / df * np.linalg.inv(X.T @ X)[1, 1])
    return B[1] / se


@njit(parallel=True)
def all_es(order_rank, w, set_ptr, set_idx, n):
    """order_rank[g] = rank of gene g (0 = most positive t); w[g] = |t|. Returns ES per set."""
    S = len(set_ptr) - 1
    out = np.empty(S)
    for s in prange(S):
        idx = set_idx[set_ptr[s]:set_ptr[s + 1]]
        k = len(idx)
        pos = np.empty(k, np.int64); ww = np.empty(k)
        for j in range(k):
            pos[j] = order_rank[idx[j]]; ww[j] = w[idx[j]]
        o = np.argsort(pos); pos = pos[o]; ww = ww[o]
        total = ww.sum(); miss = 1.0 / (n - k)
        acc = 0.0; peak = 0.0; prev = -1
        for j in range(k):
            acc -= (pos[j] - prev - 1) * miss
            if abs(acc) > abs(peak):
                peak = acc
            acc += ww[j] / total
            if abs(acc) > abs(peak):
                peak = acc
            prev = pos[j]
        out[s] = peak
    return out


def es_for(t, set_ptr, set_idx):
    order = np.argsort(-t)
    rank = np.empty(len(t), np.int64); rank[order] = np.arange(len(t))
    return all_es(rank, np.abs(t), set_ptr, set_idx, len(t))


def nes_p(obs, null):
    """obs: sets; null: draws x sets. Same-sign NES and p."""
    pos_mean = np.where(null > 0, null, np.nan); neg_mean = np.where(null < 0, null, np.nan)
    mp = np.nanmean(pos_mean, 0); mn = np.nanmean(neg_mean, 0)
    nes = np.where(obs >= 0, obs / mp, -obs / mn)
    cnt_pos = np.sum(null >= obs[None, :], 0); n_pos = np.sum(null >= 0, 0)
    cnt_neg = np.sum(null <= obs[None, :], 0); n_neg = np.sum(null < 0, 0)
    p = np.where(obs >= 0, (cnt_pos + 1) / (n_pos + 1), (cnt_neg + 1) / (n_neg + 1))
    return nes, np.minimum(p, 1.0)


panel_res, null_sig = {}, {}
for panel in plan["panels"]:
    d = common_rm.load(panel)
    Z = common_rm.background(d)
    genes = list(Z.index); gi = {g: i for i, g in enumerate(genes)}
    names, ptr, idx = [], [0], []
    for nm, members in catalog.items():
        h = sorted({gi[g] for g in members if g in gi})
        if 15 <= len(h) <= 500:
            names.append(nm); idx += h; ptr.append(len(idx))
    ptr = np.array(ptr, np.int64); idx = np.array(idx, np.int64)
    C = d["C"]; Linv = linalg.inv(np.linalg.cholesky(C)); x = d["mass"]; Zt = Z.to_numpy().T
    t_obs = gene_t(Zt, x, Linv)
    es_obs = es_for(t_obs, ptr, idx)
    null = np.empty((N, len(names)))
    for b in range(N):
        sim = rng.multivariate_normal(np.zeros(len(x)), C)
        sim = (sim - sim.mean()) / sim.std() * x.std() + x.mean()
        null[b] = es_for(gene_t(Zt, sim, Linv), ptr, idx)
    nes, p = nes_p(es_obs, null)
    q = multipletests(p, method="fdr_bh")[1]
    # each null draw scored against the null (leave-in), for the concordance null
    rk = np.empty_like(null)
    for s in range(null.shape[1]):
        col = null[:, s]
        pos = col >= 0
        sp = np.sort(col[pos]); sn = np.sort(col[~pos])
        r = np.empty(N)
        r[pos] = (len(sp) - np.searchsorted(sp, col[pos], "left") + 0) / max(1, len(sp))
        r[~pos] = (np.searchsorted(sn, col[~pos], "right")) / max(1, len(sn))
        rk[:, s] = r
    null_sig[panel] = pd.DataFrame(np.sign(null) * (rk < 0.05), columns=names)
    panel_res[panel] = pd.DataFrame({"set": names, "size": np.diff(ptr), "ES": es_obs, "NES": nes, "p": p, "q_BH_panel": q})
    pd.DataFrame({"gene": genes, "PGLS_t": t_obs}).to_csv(RUN / f"results/05_gene_rank_{panel}.tsv", sep="\t", index=False)
    print(panel, len(genes), "genes", len(names), "sets", "q<0.05:", int((q < 0.05).sum()), "min q", float(q.min()))

common = sorted(set.intersection(*[set(v.set) for v in panel_res.values()]))
W = pd.DataFrame({"set": common})
for panel, r in panel_res.items():
    r = r.set_index("set").loc[common]
    for c in ["size", "NES", "p", "q_BH_panel"]:
        W[f"{c}_{panel}"] = r[c].to_numpy()
sig = np.stack([np.sign(W[f"NES_{p}"]) * (W[f"p_{p}"] < 0.05) for p in plan["panels"]])
W["concordant_3"] = (np.abs(sig.sum(0)) == 3)
W["concordant_2plus"] = (np.maximum((sig > 0).sum(0), (sig < 0).sum(0)) >= 2)
z = np.stack([np.sign(W[f"NES_{p}"]) * stats.norm.isf(W[f"p_{p}"] / 2) for p in plan["panels"]])
W["stouffer_z"] = z.sum(0) / np.sqrt(3)
W["stouffer_p"] = 2 * stats.norm.sf(np.abs(W.stouffer_z))
W["stouffer_q"] = multipletests(W.stouffer_p, method="fdr_bh")[1]
nsig = np.stack([null_sig[p][common].to_numpy() for p in plan["panels"]])  # panels x draws x sets
null3 = (np.abs(nsig.sum(0)) == 3).sum(1)
null2 = (np.maximum((nsig > 0).sum(0), (nsig < 0).sum(0)) >= 2).sum(1)
summ = dict(n_common_sets=len(common), observed_concordant_3=int(W.concordant_3.sum()), null3_median=float(np.median(null3)), null3_q95=float(np.quantile(null3, .95)),
            p_concordant_3=float((1 + (null3 >= W.concordant_3.sum()).sum()) / (1 + N)),
            observed_concordant_2plus=int(W.concordant_2plus.sum()), null2_median=float(np.median(null2)), null2_q95=float(np.quantile(null2, .95)),
            p_concordant_2plus=float((1 + (null2 >= W.concordant_2plus.sum()).sum()) / (1 + N)),
            stouffer_q_lt_0_05=int((W.stouffer_q < 0.05).sum()),
            per_panel_q_lt_0_05={p: int((r.q_BH_panel < 0.05).sum()) for p, r in panel_res.items()})
for p, r in panel_res.items():
    r.to_csv(RUN / f"results/05_gsea_{p}.tsv", sep="\t", index=False)
W.sort_values("stouffer_p").to_csv(RUN / "results/05_gsea_cross_panel.tsv", sep="\t", index=False)
(RUN / "results/05_summary.json").write_text(json.dumps(summ, indent=1))
print(json.dumps(summ, indent=1))
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
cols = ["set"] + [f"NES_{p}" for p in plan["panels"]] + [f"p_{p}" for p in plan["panels"]] + ["stouffer_q"]
print("concordant in all three:\n", W[W.concordant_3][cols].round(4).to_string())
print("top Stouffer:\n", W.sort_values("stouffer_p").head(25)[cols].round(4).to_string())
print("Ruslan sets:\n", W[W.set.str.startswith("Ruslan")][cols].round(4).to_string())
for p, r in panel_res.items():
    print(p, "top by q:\n", r.sort_values("p").head(10).round(4).to_string())
print("DONE05")
