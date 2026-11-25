#!/usr/bin/env python3
"""PGLS for mammal gene sets that passed the abundance-null BH.

OLS rows are kept. Dated-tree PGLS is a second reading: if it removes the
slope, both numbers stay and the row is marked as a clade effect.
Trees are approximate ultrametric topologies with published-order node ages,
not a fresh TimeTree download. Grafen lengths are the date-free check.
Bird OXPHOS phylogeny is not reopened.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = Path(__file__).resolve().parent
s10 = _load("s10_pgls", HERE / "10_pgls_simplicity_S.py")
s73 = _load("s73_phylo", HERE / "73_phylo_oxphos.py")

OUT = RES / "mammal_size_sets"
MATDIR = RES / "mammal_honest"
N_SIM = 10_000
RNG = np.random.default_rng(20260924)

# Node age is MYA. A pair is (age, left, right). A tip is a string.
# Ages are approximate and chosen so every tip is the same distance from the root.
TREES = {
    "ma2016": (
        90,
        (81, "shrew", (20, "little_brown_bat", "big_brown_bat")),
        (
            73,
            (45, "porcupine", (33, "guinea_pig", "chinchilla")),
            (
                66,
                (40, "chipmunk", (9, "red_squirrel", "fox_squirrel")),
                (
                    62,
                    "beaver",
                    (
                        27,
                        "grass_rat",
                        (20, "cotton_rat", (12, "meadow_vole", (3, "deer_mouse", "white_footed_mouse"))),
                    ),
                ),
            ),
        ),
    ),
    "wagner": (
        160,
        "opossum",
        (
            96,
            (90, (29, "human", "rhesus"), (80, "rabbit", (70, "rat", "guinea_pig"))),
            (
                85,
                (55, "cat", "dog"),
                (78, "horse", (60, "pig", (25, "cow", "sheep"))),
            ),
        ),
    ),
    "gse105160": (
        90,
        "mouse",
        (
            43,
            "squirrel_monkey",
            (
                29,
                (12, "olive_baboon", (5, "pigtail_macaque", "rhesus")),
                (15, "orangutan", (9, "gorilla", (6, "human", (2, "bonobo", "chimp")))),
            ),
        ),
    ),
}


def tip_name(species: str) -> str:
    return species.replace("-", "_").replace(" ", "_")


def to_newick(node, parent_age: float) -> str:
    if isinstance(node, str):
        return f"{node}:{parent_age}"
    age, left, right = node
    inner = f"({to_newick(left, age)},{to_newick(right, age)})"
    if parent_age == age:
        return inner
    return f"{inner}:{parent_age - age}"


def tip_depths(root: s10.Node) -> dict[str, float]:
    out = {}

    def walk(n: s10.Node, dist: float):
        dist = dist + (0.0 if n is root else n.length)
        if not n.children:
            out[n.name] = dist
            return
        for c in n.children:
            walk(c, dist)

    walk(root, 0.0)
    return out


def bm_p(x: np.ndarray, V: np.ndarray, r_obs: float) -> float:
    n = len(x)
    L = np.linalg.cholesky(V + np.eye(n) * 1e-9 * np.trace(V) / n)
    null = np.empty(N_SIM)
    for i in range(N_SIM):
        y = L @ RNG.normal(size=n)
        null[i] = np.corrcoef(x, y)[0, 1]
    return float(np.mean(np.abs(null) >= abs(r_obs)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = pd.read_csv(OUT / "sets_vs_mass.tsv", sep="\t")
    passed = sets[sets.bh < 0.05].copy()
    tree_rows = []
    for panel, node in TREES.items():
        nwk = to_newick(node, node[0]) + ";"
        (OUT / f"tree_{panel}.nwk").write_text(nwk + "\n")
        root = s10.parse_newick(nwk)
        depths = tip_depths(root)
        span = max(depths.values()) - min(depths.values())
        if span > 1e-6:
            raise SystemExit(f"{panel} tree is not ultrametric: {depths}")
        tree_rows.append(dict(panel=panel, root_mya=node[0], n_tips=len(depths), depth_span=span, note="approximate"))
    pd.DataFrame(tree_rows).to_csv(OUT / "trees_note.tsv", sep="\t", index=False)

    if passed.empty:
        pd.DataFrame(columns=[
            "panel", "catalog", "id", "name", "n", "ols_r", "ols_p",
            "pgls_tt_slope", "pgls_tt_p", "pgls_gr_slope", "pgls_gr_p",
            "lambda_pegel", "p_brownian", "pic_r", "pic_p", "note",
        ]).to_csv(OUT / "pgls.tsv", sep="\t", index=False)
        print("BH суммы не прошёл ни один набор. PGLS не считался.", flush=True)
        return

    rows = []
    for panel, g in passed.groupby("panel"):
        mass = pd.read_csv(MATDIR / f"{panel}_mass.tsv", sep="\t")
        wide = pd.read_csv(MATDIR / f"{panel}_M.tsv", sep="\t", index_col=0)[mass.species.tolist()]
        tips = [tip_name(s) for s in mass.species]
        x = np.log10(mass.mass_g.to_numpy(float))
        nwk = (OUT / f"tree_{panel}.nwk").read_text().strip()
        root_tt = s10.parse_newick(nwk)
        # parse_newick tip order follows the newick, not our species order.
        v_order = list(s10.path_from_root(root_tt))
        # reorder V to mass-file order
        V_tt = s10.vcv_brownian(root_tt, tips)
        root_gr = s10.parse_newick(nwk)
        s10.apply_grafen(root_gr)
        V_gr = s10.vcv_brownian(root_gr, tips)
        X = np.column_stack([np.ones(len(x)), x])
        print(f"{panel}: PGLS для {len(g)} наборов", flush=True)
        for rec in g.itertuples(index=False):
            genes = [z for z in str(rec.genes).split(";") if z in wide.index]
            y = np.log10(wide.loc[genes].sum(axis=0).to_numpy(float))
            ols = s10.ols(y, x)
            gls_tt = s10.pgls(y, X, V_tt)
            gls_gr = s10.pgls(y, X, V_gr)
            lam, _ = s73.pagel_lambda(y, x, V_tt)
            p_bm = bm_p(x, V_tt, ols["r"])
            pic = s73.pic_slope(y, x, root_tt, tips)
            note = ""
            if float(gls_tt["p"][1]) >= 0.05:
                note = "клад, не размер"
            rows.append(dict(
                panel=panel, catalog=rec.catalog, id=rec.id, name=rec.name, n=len(genes),
                ols_r=ols["r"], ols_p=ols["p"],
                pgls_tt_slope=float(gls_tt["beta"][1]), pgls_tt_p=float(gls_tt["p"][1]),
                pgls_gr_slope=float(gls_gr["beta"][1]), pgls_gr_p=float(gls_gr["p"][1]),
                lambda_pegel=lam, p_brownian=p_bm,
                pic_r=pic["r"], pic_p=pic["p"], note=note,
            ))
            print(f"  {rec.id} ols r={ols['r']:+.3f} pgls p={float(gls_tt['p'][1]):.3f} {note}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "pgls.tsv", sep="\t", index=False)
    print(f"wrote {OUT / 'pgls.tsv'}", flush=True)


if __name__ == "__main__":
    main()
