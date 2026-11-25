#!/usr/bin/env python3
"""Gene groups vs body mass, beyond KEGG.

Same score as scripts 75 and 76: promoter CAGE within 2 kb of the CDS,
present in all seven species, missing genes dropped rather than zeroed,
log10(mean TPM) against log10 mass. Fibroblast and epiblast stay separate.

Catalogs are pre-specified (MSigDB 2026.1 Hallmark, Reactome, GO BP;
CORUM 5.3 human complexes). The sum is compared with 2000 abundance-matched
draws, the same null as script 76. CAMERA estimates the within-set
correlation (Wu and Smyth 2012, strict mode). Ranked GSEA permutes the
seven masses (7! = 5040), not the genes.

De novo groups come only from average-linkage clustering of the
within-species centred fibroblast matrix. Descent is Meinshausen 2008:
a node is tested only if its parent was rejected, and the p-value is
multiplied by (genes in the tree) / (genes in the node). The tree is not
built from the seven species means.

The nine locked sets are not rewritten.
"""
from __future__ import annotations

import sys
import time
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES,
    SPECIES,
    SPECIES_ORDER,
    normalize_symbol,
)

CAT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/genesets/catalog_2026-09-22")
NINE_DISJOINT = Path(
    "/path/to/user-home/chicken/kleiber_birds_gemcat/corpus_kleiber_2026-08-06"
    "/02_база_знаний/working_sets_disjoint_2026-09-22.tsv"
)
NINE_OVERLAP = Path(
    "/path/to/user-home/chicken/kleiber_birds_gemcat/corpus_kleiber_2026-08-06"
    "/02_база_знаний/working_sets_overlap_ok_2026-09-22.tsv"
)
OUT = RES / "size_gene_groups"
M_FIB = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
M_EPI = RES / "story_epi_traits" / "M_epi_nearest2kb.tsv"
M_CENTERED = RES / "modules_v2" / "M_fib_centered.tsv"

LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
YC = LOGM - LOGM.mean()
YSS = float(YC @ YC)
N_NULL = 2000
ALPHA = 0.05
RNG = np.random.default_rng(20260922)
PERMS = np.array(list(permutations(range(7))), dtype=np.int16)
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]


def bh(p: np.ndarray) -> np.ndarray:
    return stats.false_discovery_control(np.nan_to_num(np.asarray(p, float), nan=1.0), method="bh")


def load_universe(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", index_col=0)
    raw.index = pd.Index([normalize_symbol(g) for g in raw.index])
    raw = raw[raw.index != ""]
    raw = raw[~raw.index.duplicated(keep="first")]
    raw = raw.reindex(columns=SPECIES_ORDER)
    det = raw.notna() & (raw > 0)
    wide = raw.loc[det.all(axis=1)]
    keep = wide.std(axis=1) > 0
    return wide.loc[keep]


def r_of(y: np.ndarray) -> float:
    z = y - y.mean()
    den = np.sqrt(float(z @ z) * YSS)
    return float(z @ YC / den) if den else np.nan


def fit(y: np.ndarray) -> tuple[float, float]:
    r = stats.linregress(LOGM, np.asarray(y, float))
    return float(r.rvalue), float(r.pvalue)


def perm_p(y: np.ndarray, obs_r: float) -> float:
    z = y - y.mean()
    den = np.sqrt(float(z @ z) * YSS)
    if not den or not np.isfinite(obs_r):
        return np.nan
    rperm = (YC[PERMS] @ z) / den
    return float(np.mean(np.abs(rperm) >= abs(obs_r) - 1e-12))


class AbundanceNull:
    def __init__(self, wide: pd.DataFrame):
        self.genes = np.asarray(wide.index)
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.mat = wide.to_numpy(float)
        chick = np.log10(self.mat[:, SPECIES_ORDER.index("Chicken")].clip(1e-12))
        self.dec = pd.qcut(
            pd.Series(chick).rank(method="first"), 10, labels=False, duplicates="drop"
        ).to_numpy()
        n_dec = int(self.dec.max()) + 1
        self.buckets = [np.flatnonzero(self.dec == k) for k in range(n_dec)]
        self._cache: dict[tuple[int, ...], np.ndarray] = {}

    def rows_of(self, genes: list[str]) -> np.ndarray:
        return np.array([self.pos[g] for g in genes if g in self.pos], dtype=int)

    def r_mean(self, rows: np.ndarray) -> float:
        return r_of(np.log10(self.mat[rows].mean(axis=0)))

    def _null_rs(self, want: np.ndarray) -> np.ndarray:
        key = tuple(int(x) for x in want)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        parts = []
        for k, need in enumerate(want):
            need = int(need)
            if need <= 0:
                continue
            av = self.buckets[k]
            if need >= len(av):
                parts.append(np.broadcast_to(av, (N_NULL, len(av))).copy())
                continue
            keys = RNG.random((N_NULL, len(av)))
            sel = np.argpartition(keys, need - 1, axis=1)[:, :need]
            parts.append(av[sel])
        pick = np.concatenate(parts, axis=1)
        y = np.log10(self.mat[pick].mean(axis=1))
        z = y - y.mean(axis=1, keepdims=True)
        den = np.sqrt((z * z).sum(axis=1) * YSS)
        rs = (z @ YC) / den
        self._cache[key] = rs
        return rs

    def matched(self, rows: np.ndarray, obs_r: float) -> tuple[float, float]:
        if len(rows) < 3 or not np.isfinite(obs_r):
            return np.nan, np.nan
        want = np.bincount(self.dec[rows], minlength=len(self.buckets))
        null = self._null_rs(want)
        null = null[np.isfinite(null)]
        # k/N as in script 76. A complete miss is still only resolved to 1/N;
        # storing 0 would make BH identically 0.
        p = float(np.mean(np.abs(null) >= abs(obs_r)))
        if p == 0.0:
            p = 1.0 / N_NULL
        return p, float(np.median(null))


def load_msig(path: Path) -> list[tuple[str, str, list[str]]]:
    df = pd.read_csv(path, sep="\t", usecols=["gs_id", "gs_name", "db_gene_symbol"])
    out = []
    for (sid, name), g in df.groupby(["gs_id", "gs_name"], sort=False):
        genes = sorted({normalize_symbol(x) for x in g["db_gene_symbol"] if normalize_symbol(x)})
        out.append((str(sid), str(name), genes))
    return out


def load_corum(path: Path) -> list[tuple[str, str, list[str]]]:
    df = pd.read_csv(path, sep="\t", usecols=["complex_id", "complex_name", "subunits_gene_name"])
    out = []
    for rec in df.itertuples(index=False):
        genes = sorted({
            normalize_symbol(x)
            for x in str(rec.subunits_gene_name).split(";")
            if normalize_symbol(x)
        })
        out.append((f"CORUM:{rec.complex_id}", str(rec.complex_name), genes))
    return out


def prepare(catalog: list[tuple[str, str, list[str]]], universe: set[str],
            lo: int, hi: int | None) -> list[tuple[str, str, list[str]]]:
    kept = []
    seen: set[str] = set()
    for sid, name, genes in catalog:
        if sid in seen:
            continue
        seen.add(sid)
        hit = [g for g in genes if g in universe]
        if len(hit) < lo:
            continue
        if hi is not None and len(hit) > hi:
            continue
        kept.append((sid, name, hit))
    return kept


def score_sums(null: AbundanceNull, sets: list[tuple[str, str, list[str]]]) -> pd.DataFrame:
    rows = []
    t0 = time.time()
    for i, (sid, name, genes) in enumerate(sets, start=1):
        idx = null.rows_of(genes)
        y = np.log10(null.mat[idx].mean(axis=0))
        r, p = fit(y)
        pn, med = null.matched(idx, r)
        rows.append(dict(id=sid, name=name, n=int(len(idx)), r=r, p=p, p_null=pn, null_median=med))
        if i % 400 == 0:
            print(f"    суммы {i}/{len(sets)}  {time.time() - t0:.0f} с", flush=True)
    t = pd.DataFrame(rows)
    t["bh_null"] = bh(t["p_null"].to_numpy())
    return t


def camera(wide: pd.DataFrame, sets: list[tuple[str, str, list[str]]]) -> pd.DataFrame:
    """Competitive test, correlation estimated per set. limma camera, inter.gene.cor=NA."""
    x = np.log10(wide.to_numpy(float))
    genes = list(wide.index)
    pos = {g: i for i, g in enumerate(genes)}
    n = x.shape[1]
    design = np.column_stack([np.ones(n), LOGM])
    q, rfac = np.linalg.qr(design, mode="complete")
    p_rank = design.shape[1]
    df = n - p_rank
    q2 = q[:, p_rank:]
    u = x @ q2
    sigma = np.sqrt(np.sum(u * u, axis=1) / df)
    sigma = np.where(sigma > 0, sigma, np.nan)
    us = u / sigma[:, None]
    qp = q[:, :p_rank]
    coef = np.linalg.solve(rfac[:p_rank], qp.T @ x.T)
    xtx_inv = np.linalg.inv(rfac[:p_rank].T @ rfac[:p_rank])
    v = float(xtx_inv[1, 1])
    tstat = coef[1] / (sigma * np.sqrt(v))
    g = len(tstat)
    mean_all = float(np.nanmean(tstat))
    var_stat = float(np.nanvar(tstat, ddof=1))
    rows = []
    for sid, name, geneset in sets:
        idx = np.array([pos[g] for g in geneset], dtype=int)
        m = len(idx)
        m2 = g - m
        uset = us[idx]
        cm = np.nanmean(uset, axis=0)
        vif = float(m * np.nanmean(cm * cm))
        vif = float(np.clip(vif, 1.0, m))
        rho = (vif - 1.0) / (m - 1.0) if m > 1 else np.nan
        mean_set = float(np.nanmean(tstat[idx]))
        delta = g / m2 * (mean_set - mean_all)
        var_pooled = ((g - 1) * var_stat - delta ** 2 * m * m2 / g) / (g - 2)
        if var_pooled <= 0 or not np.isfinite(var_pooled):
            p = np.nan
            tt = np.nan
        else:
            tt = delta / np.sqrt(var_pooled * (vif / m + 1.0 / m2))
            p = float(2 * stats.t.sf(abs(tt), df))
        rows.append(dict(
            id=sid, camera_t=tt, camera_vif=vif, camera_rho=rho, camera_p=p,
        ))
    t = pd.DataFrame(rows)
    t["bh_camera"] = bh(t["camera_p"].to_numpy())
    return t


def es_from_positions(pos: np.ndarray, weight: np.ndarray, n_genes: int) -> np.ndarray:
    """Weighted running-sum ES. pos and weight are (m, n_perm), pos sorted."""
    m, n_perm = pos.shape
    n_miss = n_genes - m
    nr = weight.sum(axis=0)
    hit_cum = np.cumsum(weight, axis=0) / nr
    misses = (pos - np.arange(m)[:, None]) / n_miss
    after = hit_cum - misses
    before = np.empty_like(after)
    before[0] = -pos[0] / n_miss
    if m > 1:
        before[1:] = hit_cum[:-1] - (pos[1:] - np.arange(1, m)[:, None]) / n_miss
    both = np.concatenate([before, after], axis=0)
    pick = np.argmax(np.abs(both), axis=0)
    return both[pick, np.arange(n_perm)]


def es_slow(correl: np.ndarray, in_set: np.ndarray) -> float:
    order = np.argsort(-correl)
    c = correl[order]
    hit = in_set[order]
    weight = np.abs(c)
    nr = float(weight[hit].sum())
    n_miss = int((~hit).sum())
    step = np.where(hit, weight / nr, -1.0 / n_miss)
    walk = np.cumsum(step)
    return float(walk[np.argmax(np.abs(walk))])


def check_es() -> None:
    rng = np.random.default_rng(1)
    correl = rng.normal(size=40)
    in_set = np.zeros(40, dtype=bool)
    in_set[:7] = True
    slow = es_slow(correl, in_set)
    order = np.argsort(-correl)
    rank = np.empty(40, dtype=int)
    rank[order] = np.arange(40)
    idx = np.flatnonzero(in_set)
    pos = rank[idx]
    ord_in = np.argsort(pos)
    fast = es_from_positions(pos[ord_in][:, None], np.abs(correl[idx][ord_in])[:, None], 40)[0]
    if abs(slow - fast) > 1e-8:
        raise SystemExit(f"GSEA расходится: slow {slow} fast {fast}")


def gsea(wide: pd.DataFrame, sets: list[tuple[str, str, list[str]]]) -> pd.DataFrame:
    x = np.log10(wide.to_numpy(float))
    genes = list(wide.index)
    pos = {g: i for i, g in enumerate(genes)}
    xc = x - x.mean(axis=1, keepdims=True)
    xn = np.sqrt((xc * xc).sum(axis=1))
    corr = (xc @ YC[PERMS].T) / (xn[:, None] * np.sqrt(YSS))
    g, n_perm = corr.shape
    order = np.argsort(-corr, axis=0)
    rank = np.empty_like(order)
    rank[order, np.arange(n_perm)] = np.arange(g)[:, None]
    ident = int(np.where((PERMS == np.arange(7)).all(axis=1))[0][0])
    out = []
    t0 = time.time()
    for i, (sid, _name, geneset) in enumerate(sets, start=1):
        idx = np.array([pos[gene] for gene in geneset], dtype=int)
        positions = rank[idx]
        ord_in = np.argsort(positions, axis=0)
        pos_s = np.take_along_axis(positions, ord_in, axis=0)
        weight = np.take_along_axis(np.abs(corr[idx]), ord_in, axis=0)
        es = es_from_positions(pos_s, weight, g)
        obs = float(es[ident])
        p = float(np.mean(np.abs(es) >= abs(obs) - 1e-12))
        out.append(dict(id=sid, gsea_es=obs, gsea_p=p))
        if i % 400 == 0:
            print(f"    GSEA {i}/{len(sets)}  {time.time() - t0:.0f} с", flush=True)
    t = pd.DataFrame(out)
    t["bh_gsea"] = bh(t["gsea_p"].to_numpy())
    return t


def check_controls(fib: pd.DataFrame, epi: pd.DataFrame) -> None:
    for label, wide, expect in (
        ("фибробласт", fib, -0.860),
        ("эпибласт", epi, -0.220),
    ):
        rows = [g for g in SIXTEEN if g in wide.index]
        if len(rows) != 16:
            raise SystemExit(f"OXPHOS 16 в {label}: нашлось {len(rows)}")
        y = np.log10(wide.loc[rows].to_numpy(float).mean(axis=0))
        r, p = fit(y)
        print(f"  контроль 16, {label}: r={r:+.3f} p={p:.4f}", flush=True)
        if abs(r - expect) > 0.015:
            raise SystemExit(f"контроль 16 не сошёлся в {label}: {r:+.3f}, ждали {expect}")


def nine_sets() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for path in (NINE_DISJOINT, NINE_OVERLAP):
        df = pd.read_csv(path, sep="\t")
        for sid, g in df.groupby("set_id"):
            out[str(sid)] = {normalize_symbol(x) for x in g["gene"]}
    return out


def overlap_columns(sets: list[tuple[str, str, list[str]]], nine: dict[str, set[str]]) -> pd.DataFrame:
    union = set().union(*nine.values())
    rows = []
    for sid, _name, genes in sets:
        gs = set(genes)
        inter_all = gs & union
        best_name, best_j = "", 0.0
        for nid, ng in nine.items():
            inter = len(gs & ng)
            j = inter / len(gs | ng) if gs or ng else 0.0
            if j > best_j:
                best_name, best_j = nid, j
        rows.append(dict(
            id=sid,
            n_in_nine=len(inter_all),
            frac_in_nine=len(inter_all) / len(gs),
            n_new=len(gs - union),
            best_nine=best_name,
            jaccard_best=best_j,
        ))
    return pd.DataFrame(rows)


def jaccard_pairs(passed: list[tuple[str, list[str]]], path: Path) -> None:
    rows = []
    for i, (a, ga) in enumerate(passed):
        sa = set(ga)
        for b, gb in passed[i + 1:]:
            sb = set(gb)
            inter = len(sa & sb)
            if inter == 0:
                continue
            j = inter / len(sa | sb)
            if j >= 0.5:
                rows.append(dict(a=a, b=b, n_inter=inter, jaccard=j))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def run_catalog(name: str, sets: list[tuple[str, str, list[str]]],
                fib: AbundanceNull, epi: AbundanceNull,
                wide_fib: pd.DataFrame, wide_epi: pd.DataFrame,
                nine: dict[str, set[str]]) -> pd.DataFrame:
    print(f"\n{name}: {len(sets)} наборов", flush=True)
    print("  суммы, фибробласт", flush=True)
    sf = score_sums(fib, sets).rename(columns=lambda c: c if c in {"id", "name", "n"} else f"{c}_fib")
    print("  суммы, эпибласт", flush=True)
    se = score_sums(epi, sets)
    se = se.drop(columns=["name", "n"]).rename(columns=lambda c: c if c == "id" else f"{c}_epi")
    print("  CAMERA", flush=True)
    cf = camera(wide_fib, sets).rename(columns=lambda c: c if c == "id" else f"{c}_fib")
    ce = camera(wide_epi, sets).rename(columns=lambda c: c if c == "id" else f"{c}_epi")
    print("  GSEA", flush=True)
    gf = gsea(wide_fib, sets).rename(columns=lambda c: c if c == "id" else f"{c}_fib")
    ge = gsea(wide_epi, sets).rename(columns=lambda c: c if c == "id" else f"{c}_epi")
    ov = overlap_columns(sets, nine)
    t = sf.merge(se, on="id").merge(cf, on="id").merge(ce, on="id").merge(gf, on="id").merge(ge, on="id").merge(ov, on="id")
    t.insert(0, "catalog", name)
    t = t.sort_values("p_null_fib")
    n_null = int((t["bh_null_fib"] < ALPHA).sum() + (t["bh_null_epi"] < ALPHA).sum())
    n_cam = int((t["bh_camera_fib"] < ALPHA).sum() + (t["bh_camera_epi"] < ALPHA).sum())
    n_g = int((t["bh_gsea_fib"] < ALPHA).sum() + (t["bh_gsea_epi"] < ALPHA).sum())
    print(f"  BH<0.05 строк (ткань считается отдельно): нуль {n_null}, CAMERA {n_cam}, GSEA {n_g}", flush=True)
    return t


def descend_tree(leaf_rows: np.ndarray, mat: pd.DataFrame, epi: pd.DataFrame) -> pd.DataFrame:
    """Average linkage on the centred matrix. Meinshausen on the species sum."""
    print("\nдерево: связь профилей внутри вида", flush=True)
    # mat and epi are aligned to the same gene order as leaf_rows indexes into them
    x = mat.to_numpy(float)
    xe = epi.to_numpy(float)
    root = to_tree(linkage(pdist(leaf_rows, metric="correlation"), method="average"))
    nodes: dict[int, dict] = {}
    sys.setrecursionlimit(max(20000, len(leaf_rows) + 100))

    def rec(node) -> int:
        if node.is_leaf():
            i = int(node.id)
            nodes[i] = {"genes": [i], "left": None, "right": None}
            return i
        left = rec(node.get_left())
        right = rec(node.get_right())
        i = 10_000_000 + int(node.id)
        nodes[i] = {"genes": nodes[left]["genes"] + nodes[right]["genes"], "left": left, "right": right}
        return i

    root_id = rec(root)
    m_total = len(nodes[root_id]["genes"])
    print(f"  листьев {m_total}", flush=True)
    tested = []
    stack = [(root_id, True, 0.0, -1)]
    while stack:
        nid, parent_ok, parent_h, parent_id = stack.pop()
        if not parent_ok:
            continue
        node = nodes[nid]
        idx = np.asarray(node["genes"], dtype=int)
        size = int(len(idx))
        y = np.log10(x[idx].mean(axis=0))
        r, p_an = fit(y)
        p_perm = perm_p(y, r)
        p_adj = p_perm * m_total / size
        p_h = max(p_adj, parent_h)
        rejected = bool(np.isfinite(p_h) and p_h <= ALPHA)
        ye = np.log10(xe[idx].mean(axis=0))
        r_e, p_e = fit(ye)
        tested.append(dict(
            node=nid, parent=parent_id, n=size, r_fib=r, p_analytic_fib=p_an,
            p_perm_fib=p_perm, p_adj=p_adj, p_hier=p_h, rejected=int(rejected),
            r_epi=r_e, p_epi=p_e, is_root=int(nid == root_id),
            genes=";".join(mat.index[idx]),
        ))
        if rejected and node["left"] is not None:
            stack.append((node["left"], True, p_h, nid))
            stack.append((node["right"], True, p_h, nid))
    t = pd.DataFrame(tested)
    parents = set(t.loc[t.rejected == 1, "parent"].astype(int))
    t["minimal"] = ((t.rejected == 1) & ~t.node.isin(parents)).astype(int)
    print(f"  проверено узлов {len(t)}, отвергнуто {(t.rejected == 1).sum()}, "
          f"из них крайних {(t.minimal == 1).sum()}", flush=True)
    return t


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    check_es()
    print("вселенные", flush=True)
    fib = load_universe(M_FIB)
    epi = load_universe(M_EPI)
    shared = sorted(set(fib.index) & set(epi.index))
    fib = fib.loc[shared]
    epi = epi.loc[shared]
    print(f"  общая вселенная {len(shared)}", flush=True)
    check_controls(fib, epi)
    universe = set(shared)
    nine = nine_sets()

    catalogs = {
        "hallmark": prepare(load_msig(CAT / "hallmark_genes.tsv"), universe, 8, None),
        "reactome": prepare(load_msig(CAT / "reactome_genes.tsv"), universe, 8, None),
        "go_bp": prepare(load_msig(CAT / "gobp_genes.tsv"), universe, 8, 50),
        "corum": prepare(load_corum(CAT / "corum_human.txt"), universe, 3, 50),
    }
    for name, sets in catalogs.items():
        print(f"  {name}: к тесту {len(sets)}", flush=True)

    null_fib = AbundanceNull(fib)
    null_epi = AbundanceNull(epi)
    frames = []
    gene_of = {}
    for name, sets in catalogs.items():
        frames.append(run_catalog(name, sets, null_fib, null_epi, fib, epi, nine))
        for sid, _n, genes in sets:
            gene_of[sid] = genes
    table = pd.concat(frames, ignore_index=True)
    table.to_csv(OUT / "catalogs_vs_mass.tsv", sep="\t", index=False)

    pass_mask = (
        (table["bh_null_fib"] < ALPHA) | (table["bh_null_epi"] < ALPHA)
        | (table["bh_camera_fib"] < ALPHA) | (table["bh_camera_epi"] < ALPHA)
        | (table["bh_gsea_fib"] < ALPHA) | (table["bh_gsea_epi"] < ALPHA)
    )
    passed = table.loc[pass_mask].copy()
    passed.to_csv(OUT / "catalogs_bh05.tsv", sep="\t", index=False)
    jaccard_pairs([(r.id, gene_of[r.id]) for r in passed.itertuples()], OUT / "pass_jaccard.tsv")
    gene_rows = []
    for rec in passed.itertuples():
        for g in gene_of[rec.id]:
            gene_rows.append(dict(catalog=rec.catalog, id=rec.id, name=rec.name, gene=g))
    pd.DataFrame(gene_rows).to_csv(OUT / "pass_genes.tsv", sep="\t", index=False)
    print(f"\nпрошли хотя бы одну поправку: {len(passed)}", flush=True)

    # Tree. Profiles are the centred library matrix; the test uses species means.
    centered = pd.read_csv(M_CENTERED, sep="\t", index_col=0)
    centered.index = pd.Index([normalize_symbol(g) for g in centered.index])
    centered = centered[centered.index != ""]
    centered = centered[~centered.index.duplicated(keep="first")]
    centered = centered.loc[centered.index.isin(universe)]
    centered = centered.loc[centered.std(axis=1) > 0]
    centered = centered.dropna(axis=0)
    order = list(centered.index)
    profiles = centered.to_numpy(float)
    fib_al = fib.loc[order]
    epi_al = epi.loc[order]
    # descend_tree clusters `profiles`, tests rows of fib_al in the same order
    tree = descend_tree(profiles, fib_al, epi_al)
    tree.to_csv(OUT / "tree_tested.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
