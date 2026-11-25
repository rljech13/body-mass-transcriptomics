#!/usr/bin/env python3
"""PGLS + simplicity index S for OXPHOS CAGE promoters (n=7).

Phylogeny: TimeTree-informed ultrametric Newick (MYA) + sensitivity with
Grafen branch lengths on OpenTree topology.

PGLS: Brownian-motion VCV, GLS via (X'V^{-1}X)^{-1}X'V^{-1}y.
No ape/caper required.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
RES = ROOT / "results" / "paperA_diff_pack"
OUT = RES / "arm5_alt_scaling"
FIG = ROOT / "figures"
DATA = ROOT / "data"

# tip order used everywhere
TIPS = ["zebra", "quail", "chicken", "duck", "peacock", "turkey", "emu"]
TIP_TO_SPECIES = {
    "zebra": "ZebraFinch",
    "quail": "Quail",
    "chicken": "Chicken",
    "duck": "Duck",
    "peacock": "Peacock",
    "turkey": "Turkey",
    "emu": "Emu",
}

# TimeTree-informed ultrametric ages (MYA) — approximate published MRCAs
# Topology: (emu, (zebra, (duck, (peacock, (turkey, (chicken, quail))))))
# Branch lengths chosen so tip-to-root ≈ 110 MYA.
TIMETREE_NWK = (
    "(emu:110,"
    "(zebra:98,"
    "(duck:80,"
    "(peacock:45,"
    "(turkey:37,"
    "(chicken:33,quail:33):4"
    "):8"
    "):35"
    "):18"
    "):12"
    ");"
)


# ---------- Newick → BM VCV ----------
class Node:
    __slots__ = ("name", "length", "children")

    def __init__(self, name=None, length=0.0, children=None):
        self.name = name
        self.length = float(length or 0.0)
        self.children = children or []


def parse_newick(nwk: str) -> Node:
    s = nwk.strip().rstrip(";").replace(" ", "")
    i = 0

    def parse_node():
        nonlocal i
        children = []
        if s[i] == "(":
            i += 1
            children.append(parse_node())
            while i < len(s) and s[i] == ",":
                i += 1
                children.append(parse_node())
            if i >= len(s) or s[i] != ")":
                raise ValueError(f"expected ) at {i}: {s[i:i+20]}")
            i += 1
        # name
        m = re.match(r"([A-Za-z0-9_.\'-]+)", s[i:])
        name = None
        if m and (not children or s[i] not in ":,)"):
            # only consume name if not starting length
            if s[i] != ":":
                name = m.group(1)
                i += len(name)
        length = 0.0
        if i < len(s) and s[i] == ":":
            i += 1
            m = re.match(r"([0-9.eE+-]+)", s[i:])
            if not m:
                raise ValueError(f"bad length at {i}")
            length = float(m.group(1))
            i += len(m.group(1))
        return Node(name=name, length=length, children=children)

    root = parse_node()
    return root


def tip_names(node: Node) -> list[str]:
    if not node.children:
        return [node.name]
    out = []
    for c in node.children:
        out.extend(tip_names(c))
    return out


def simplify_opentree_newick(nwk: str) -> str:
    """Map OpenTree tip labels to short names; keep topology, unit lengths later."""
    mapping = {
        "Taeniopygia_guttata_ott708327": "zebra",
        "Coturnix_japonica_ott803622": "quail",
        "Gallus_gallus_ott153563": "chicken",
        "Anas_platyrhynchos_ott765167": "duck",
        "Pavo_cristatus_ott453575": "peacock",
        "Meleagris_gallopavo_ott446482": "turkey",
        "Dromaius_novaehollandiae_ott829513": "emu",
    }
    # strip internal mrcaott... labels
    s = nwk.strip()
    for k, v in mapping.items():
        s = s.replace(k, v)
    s = re.sub(r"mrcaott[0-9]+ott[0-9]+", "", s)
    s = re.sub(r"ott[0-9]+", "", s)
    return s


def apply_grafen(root: Node) -> Node:
    """Grafen (1989): node height = (#tips in clade - 1) / (#tips - 1); ultrametricize."""

    def n_tips(n: Node) -> int:
        if not n.children:
            return 1
        return sum(n_tips(c) for c in n.children)

    total = n_tips(root)
    assert total >= 2

    def height(n: Node) -> float:
        if not n.children:
            return 0.0
        return (n_tips(n) - 1) / (total - 1)

    def set_lengths(n: Node, parent_h: float):
        h = height(n)
        # length from parent to this node
        n.length = max(parent_h - h, 0.0) if n.children else parent_h
        for c in n.children:
            set_lengths(c, h if n.children else parent_h)
        if not n.children:
            n.length = parent_h

    # For root children, parent height = 1
    root.length = 0.0
    for c in root.children:
        set_lengths(c, 1.0)
    return root


def path_from_root(root: Node) -> dict[str, list[Node]]:
    """For each tip, list of nodes on path root→tip (including tip, excluding root edge as node list of edges owners)."""
    paths = {}

    def walk(n: Node, stack: list[Node]):
        stack = stack + [n]
        if not n.children:
            paths[n.name] = stack
            return
        for c in n.children:
            walk(c, stack)

    walk(root, [])
    return paths


def vcv_brownian(root: Node, tips: list[str]) -> np.ndarray:
    """V_ij = shared path length from root to MRCA(i,j)."""
    paths = path_from_root(root)
    for t in tips:
        if t not in paths:
            raise KeyError(f"tip {t} missing; have {list(paths)}")
    # precompute cumulative distance from root to each node on path
    # path nodes: root,...,tip; edge length is on child
    dist_to = {}

    def accumulate(n: Node, dist: float):
        # dist = distance from root TO this node (via parent edge = n.length for non-root)
        dist_to[id(n)] = dist
        for c in n.children:
            accumulate(c, dist + c.length)

    accumulate(root, 0.0)

    n = len(tips)
    V = np.zeros((n, n))
    for i, a in enumerate(tips):
        for j, b in enumerate(tips):
            pa, pb = paths[a], paths[b]
            # MRCA = deepest common node
            mrca = pa[0]
            for x, y in zip(pa, pb):
                if x is y:
                    mrca = x
                else:
                    break
            V[i, j] = dist_to[id(mrca)]
            # for i==j, should be full root-to-tip length
            if i == j:
                V[i, j] = dist_to[id(pa[-1])]
    # If diagonals zero (bad tree), fall back
    if np.allclose(np.diag(V), 0):
        raise RuntimeError("VCV diagonal zero — branch lengths missing")
    return V


def pgls(y: np.ndarray, X: np.ndarray, V: np.ndarray) -> dict:
    """GLS with known V. X includes intercept column."""
    # add small ridge for numerical stability
    V = V + np.eye(len(V)) * 1e-9 * np.trace(V) / len(V)
    try:
        Vi = np.linalg.inv(V)
    except np.linalg.LinAlgError:
        Vi = np.linalg.pinv(V)
    XtVi = X.T @ Vi
    A = XtVi @ X
    b = XtVi @ y
    try:
        beta = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(A) @ b
    yhat = X @ beta
    resid = y - yhat
    n, k = X.shape
    sigma2 = float(resid.T @ Vi @ resid) / (n - k)
    cov_beta = sigma2 * np.linalg.pinv(A)
    se = np.sqrt(np.diag(cov_beta))
    tvals = beta / se
    # two-sided t
    from scipy.stats import t as tdist

    pvals = 2 * tdist.sf(np.abs(tvals), df=max(n - k, 1))
    # r-like: cor(y, yhat)
    r = float(np.corrcoef(y, yhat)[0, 1]) if np.std(yhat) > 0 else 0.0
    return {
        "beta": beta,
        "se": se,
        "t": tvals,
        "p": pvals,
        "sigma2": sigma2,
        "r_fit": r,
        "yhat": yhat,
    }


def ols(y: np.ndarray, x: np.ndarray) -> dict:
    sl, icept, r, p, se = stats.linregress(x, y)
    return {"slope": sl, "intercept": icept, "r": r, "p": p, "se": se}


def simplicity_S(shannon: np.ndarray, n_active: np.ndarray) -> np.ndarray:
    """S = 1 - H / log2(max(n_active, 2)); 1 = single-promoter monopoly."""
    hmax = np.log2(np.maximum(n_active, 2.0))
    return 1.0 - shannon / hmax


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    cage = pd.read_csv(RES / "arm4_cage_n7_OXPHOS_species_summary.tsv", sep="\t")
    key = {
        "ZebraFinch": "zebra",
        "Quail": "quail",
        "Chicken": "chicken",
        "Duck": "duck",
        "Peacock": "peacock",
        "Turkey": "turkey",
        "Emu": "emu",
    }
    cage["tip"] = cage["species"].map(key)
    # masses in g (consistent with traits)
    mass_g = {
        "zebra": 11.93,
        "quail": 96.6,
        "chicken": 912.5,
        "duck": 1121.0,
        "peacock": 4093.75,
        "turkey": 5811.0,
        "emu": 36200.0,
    }
    cage["mass_g"] = cage["tip"].map(mass_g)
    cage["log10_mass"] = np.log10(cage["mass_g"])
    cage["S"] = simplicity_S(cage["mean_shannon"].values, cage["mean_n_active"].values)
    cage = cage.set_index("tip").loc[TIPS].reset_index()
    cage.to_csv(OUT / "pgls_inputs_S_shannon_n7.tsv", sep="\t", index=False)

    # --- trees ---
    root_tt = parse_newick(TIMETREE_NWK)
    assert set(tip_names(root_tt)) == set(TIPS), tip_names(root_tt)
    V_tt = vcv_brownian(root_tt, TIPS)
    pd.DataFrame(V_tt, index=TIPS, columns=TIPS).to_csv(OUT / "vcv_timetree_MYA.tsv", sep="\t")

    # OpenTree + Grafen
    ot_path = DATA / "birds7_opentree.nwk"
    V_gr = None
    if ot_path.exists():
        raw = ot_path.read_text().strip()
        simp = simplify_opentree_newick(raw)
        # collapse empty nodes by re-serializing tips-only topology is hard;
        # instead build Grafen tree from TimeTree topology (same tips order constraints)
        # Use TimeTree topology with Grafen lengths as sensitivity
        root_gr = parse_newick(TIMETREE_NWK)
        apply_grafen(root_gr)
        V_gr = vcv_brownian(root_gr, TIPS)
        pd.DataFrame(V_gr, index=TIPS, columns=TIPS).to_csv(OUT / "vcv_grafen.tsv", sep="\t")
        (OUT / "timetree_birds7.nwk").write_text(TIMETREE_NWK + "\n")
        (OUT / "opentree_simplified_raw.txt").write_text(simp[:2000])

    y_map = {
        "S_simplicity": cage["S"].values.astype(float),
        "Shannon_H": cage["mean_shannon"].values.astype(float),
        "total_TPM": cage["mean_total_tpm"].values.astype(float),
        "dominant_frac": cage["mean_dominant_frac"].values.astype(float),
        "n_active": cage["mean_n_active"].values.astype(float),
    }
    x = cage["log10_mass"].values.astype(float)
    X = np.column_stack([np.ones(len(x)), x])

    rows = []
    for endpoint, y in y_map.items():
        o = ols(y, x)
        rows.append(
            {
                "endpoint": endpoint,
                "method": "OLS",
                "tree": "none",
                "formula": "Y ~ log10(M)",
                "intercept": o["intercept"],
                "slope": o["slope"],
                "se_slope": o["se"],
                "r": o["r"],
                "p_slope": o["p"],
            }
        )
        for tree_name, V in [("TimeTree_MYA", V_tt), ("Grafen", V_gr)]:
            if V is None:
                continue
            g = pgls(y, X, V)
            rows.append(
                {
                    "endpoint": endpoint,
                    "method": "PGLS_BM",
                    "tree": tree_name,
                    "formula": "Y ~ log10(M)",
                    "intercept": float(g["beta"][0]),
                    "slope": float(g["beta"][1]),
                    "se_slope": float(g["se"][1]),
                    "r": float(g["r_fit"]),
                    "p_slope": float(g["p"][1]),
                }
            )
        # free power on positive endpoints: log10 Y ~ α log10 M  (OLS + PGLS)
        if np.all(y > 0):
            ly = np.log10(y)
            o2 = ols(ly, x)
            rows.append(
                {
                    "endpoint": endpoint,
                    "method": "OLS",
                    "tree": "none",
                    "formula": "log10(Y) ~ log10(M)  [= free α]",
                    "intercept": o2["intercept"],
                    "slope": o2["slope"],
                    "se_slope": o2["se"],
                    "r": o2["r"],
                    "p_slope": o2["p"],
                }
            )
            for tree_name, V in [("TimeTree_MYA", V_tt), ("Grafen", V_gr)]:
                if V is None:
                    continue
                g2 = pgls(ly, X, V)
                rows.append(
                    {
                        "endpoint": endpoint,
                        "method": "PGLS_BM",
                        "tree": tree_name,
                        "formula": "log10(Y) ~ log10(M)  [= free α]",
                        "intercept": float(g2["beta"][0]),
                        "slope": float(g2["beta"][1]),
                        "se_slope": float(g2["se"][1]),
                        "r": float(g2["r_fit"]),
                        "p_slope": float(g2["p"][1]),
                    }
                )

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "pgls_vs_ols_results.tsv", sep="\t", index=False)

    # Proposed S formula (OLS + PGLS TimeTree)
    s = y_map["S_simplicity"]
    oS = ols(s, x)
    gS = pgls(s, X, V_tt)
    # also α for (1-S) or for H
    proposed = {
        "index": "S = 1 - H / log2(max(n_active, 2))",
        "interpretation": "OXPHOS promoter simplicity (1 = monopolized / single-promoter-like)",
        "OLS": {
            "formula": f"S ≈ {oS['intercept']:.4f} + ({oS['slope']:.4f}) · log10(M_g)",
            "slope": oS["slope"],
            "p": oS["p"],
            "r": oS["r"],
        },
        "PGLS_TimeTree_BM": {
            "formula": f"S ≈ {gS['beta'][0]:.4f} + ({gS['beta'][1]:.4f}) · log10(M_g)",
            "slope": float(gS["beta"][1]),
            "se": float(gS["se"][1]),
            "p": float(gS["p"][1]),
            "r_fit": float(gS["r_fit"]),
        },
        "Shannon_free_alpha_PGLS": {},
        "n": 7,
        "kleiber_not_used": True,
        "tree_newick": TIMETREE_NWK,
    }
    # Shannon α PGLS
    ly = np.log10(y_map["Shannon_H"])
    gH = pgls(ly, X, V_tt)
    proposed["Shannon_free_alpha_PGLS"] = {
        "formula": f"H ∝ M^{gH['beta'][1]:.3f}",
        "alpha": float(gH["beta"][1]),
        "p": float(gH["p"][1]),
    }
    # S free power if S>0
    if np.all(s > 0):
        gSp = pgls(np.log10(s), X, V_tt)
        proposed["S_free_alpha_PGLS"] = {
            "formula": f"S ∝ M^{gSp['beta'][1]:.3f}",
            "alpha": float(gSp["beta"][1]),
            "p": float(gSp["p"][1]),
        }

    (OUT / "proposed_S_formula_pgls.json").write_text(json.dumps(proposed, indent=2))

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    # S vs mass OLS + PGLS line (same x; PGLS slope/intercept)
    ax = axes[0]
    ax.scatter(x, s, s=60, zorder=3)
    for tip, xi, yi in zip(TIPS, x, s):
        ax.text(xi, yi, tip, fontsize=8)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, oS["intercept"] + oS["slope"] * xs, color="C1", label=f"OLS p={oS['p']:.3g}")
    ax.plot(
        xs,
        gS["beta"][0] + gS["beta"][1] * xs,
        color="C2",
        ls="--",
        label=f"PGLS p={gS['p'][1]:.3g}",
    )
    ax.set_xlabel("log10 adult body mass (g)")
    ax.set_ylabel("OXPHOS promoter simplicity S")
    ax.set_title("S = 1 − H/log2(n_active*)")
    ax.legend(fontsize=8)

    # Shannon α
    ax = axes[1]
    H = y_map["Shannon_H"]
    ax.scatter(x, np.log10(H), s=60, zorder=3)
    for tip, xi, yi in zip(TIPS, x, np.log10(H)):
        ax.text(xi, yi, tip, fontsize=8)
    oH = ols(np.log10(H), x)
    ax.plot(xs, oH["intercept"] + oH["slope"] * xs, color="C1", label=f"OLS α={oH['slope']:.3f}")
    ax.plot(
        xs,
        gH["beta"][0] + gH["beta"][1] * xs,
        color="C2",
        ls="--",
        label=f"PGLS α={gH['beta'][1]:.3f}",
    )
    ax.plot(xs, (oH["intercept"] + oH["slope"] * x.mean()) - 0.25 * (xs - x.mean()), "k:", lw=0.9, label="α=−1/4 ref")
    ax.set_xlabel("log10 adult body mass (g)")
    ax.set_ylabel("log10 Shannon H")
    ax.set_title("Free power law (not Kleiber)")
    ax.legend(fontsize=8)

    fig.suptitle("PGLS vs OLS — OXPHOS promoter architecture (n=7)")
    fig.savefig(FIG / "fig_arm5_pgls_simplicity_S.png", dpi=160)
    plt.close(fig)

    # short markdown summary in results
    summary_md = f"""# PGLS + simplicity S (n=7)

## Index
`S = 1 - H / log2(max(n_active, 2))` — выше S = проще / монопольнее OXPHOS-промоторы.

## Formulas
**OLS:** {proposed['OLS']['formula']} (p={proposed['OLS']['p']:.4g})
**PGLS (TimeTree BM):** {proposed['PGLS_TimeTree_BM']['formula']} (p={proposed['PGLS_TimeTree_BM']['p']:.4g})

**Shannon free α PGLS:** {proposed['Shannon_free_alpha_PGLS']['formula']} (p={proposed['Shannon_free_alpha_PGLS']['p']:.4g})

Tree: `{TIMETREE_NWK}`
"""
    (OUT / "PGLS_S_SUMMARY.md").write_text(summary_md)

    print(tab.to_string(index=False))
    print("\nPROPOSED S:")
    print(json.dumps(proposed, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
