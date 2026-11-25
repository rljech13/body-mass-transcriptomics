"""Existing gene-group, pathway and promoter-regulation analyses.

Calculation bodies are retained from the source commit; configuration and execution scope are explicit.
"""

from __future__ import annotations
from types import SimpleNamespace


def _ruslan_style_scatters(context):
    """Ruslan-style log–log scatters: theme_light, lm, species labels, 'cor. ; p ='.

    Y is log10 of the OXPHOS *sum in TPM*, not log10 of the library fraction.
    birds - ESC.R: fibroblast = mean(colSums(norm.C of oxphos_genes)); y = log10(fibroblast).
    RNA fraction f is the same number as sum/1e6, so log10(f) = log10(sum) - 6.
    Plotting log10(f) made every RNA panel sit on a negative axis; that was scale,
    not negative expression, and not his plotted values.

    Source: analyses/avian/scripts/40_ruslan_style_scatters.py at d230e56ab942."""
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from adjustText import adjust_text
    from scipy import stats

    ROOT = Path(str(context.resolve_path("/path/to/user-home/chicken/kleiber_birds_gemcat")))
    FIG = ROOT / "figures"
    OUT = ROOT / "results" / "oxphos_species_thread" / "ruslan_metric_fib_panels"
    BIRDS_CAGE = [
        ("zebra", 11.93, 29098.4),
        ("quail", 96.6, 18002.8),
        ("chicken", 912.5, 17047.4),
        ("duck", 1121.0, 17494.5),
        ("peacock", 4093.75, 7897.5),
        ("turkey", 5811.0, 13194.1),
        ("emu", 36200.0, 8055.4),
    ]
    BIRDS_RNA = [
        ("zebra", 11.93, 16649.0),
        ("quail", 96.6, 12127.0),
        ("chicken", 912.5, 10624.0),
        ("duck", 1121.0, 15297.0),
        ("peacock", 4093.75, 8194.0),
        ("turkey", 5811.0, 9684.0),
        ("emu", 36200.0, 6636.0),
    ]
    WAGNER = [
        ("opossum", 104.0, 7373.0),
        ("rat", 300.0, 12512.0),
        ("guinea pig", 728.0, 7501.0),
        ("rabbit", 1800.0, 7043.0),
        ("cat", 4000.0, 8070.0),
        ("rhesus", 8000.0, 13097.0),
        ("dog", 40000.0, 10164.0),
        ("human", 62035.0, 13781.0),
        ("sheep", 80000.0, 12811.0),
        ("pig", 130000.0, 5853.0),
        ("horse", 400000.0, 8464.0),
        ("cow", 750000.0, 13170.0),
    ]
    GSE167580 = [
        ("mouse", 20.5, 12949.0),
        ("chicken", 912.5, 9894.0),
        ("rabbit", 1800.0, 13059.0),
        ("cat", 4000.0, 8383.0),
        ("rhesus", 8000.0, 16221.0),
        ("dog", 40000.0, 8065.0),
        ("human", 62035.0, 12182.0),
        ("sheep", 80000.0, 18170.0),
        ("pig", 130000.0, 8582.0),
        ("cow", 750000.0, 10559.0),
    ]
    GSE167580_MAMMALS = [r for r in GSE167580 if r[0] != "chicken"]
    LINE = "#3366FF"
    POINT = "#1a1a1a"

    def theme_light(ax) -> None:
        ax.set_facecolor("white")
        ax.grid(True, which="major", color="#D0D0D0", lw=0.6)
        ax.grid(True, which="minor", color="#ECECEC", lw=0.4)
        for sp in ax.spines.values():
            sp.set_color("#B0B0B0")
            sp.set_linewidth(0.8)
        ax.tick_params(colors="#4A4A4A", labelsize=8)
        ax.xaxis.label.set_color("#4A4A4A")
        ax.yaxis.label.set_color("#4A4A4A")
        ax.set_axisbelow(True)

    def panel(ax, rows, ylabel: str) -> dict:
        names = [r[0] for r in rows]
        x = np.log10(np.array([r[1] for r in rows], float))
        y = np.log10(np.array([r[2] for r in rows], float))
        fit = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 80)
        ax.plot(xs, fit.intercept + fit.slope * xs, color=LINE, lw=1.6, zorder=2)
        ax.scatter(x, y, s=22, c=POINT, zorder=3, linewidths=0)
        texts = [
            ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
            for name, xi, yi in zip(names, x, y)
        ]
        adjust_text(
            texts,
            ax=ax,
            x=x,
            y=y,
            expand=(1.15, 1.25),
            force_text=(0.25, 0.4),
            force_static=(0.15, 0.2),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )
        theme_light(ax)
        ax.set_xlabel("log10_body_mass")
        ax.set_ylabel(ylabel or "log10_fib")
        ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
        ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
        return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(len(rows))}

    def save(fig, stem: str) -> list[Path]:
        FIG.mkdir(parents=True, exist_ok=True)
        OUT.mkdir(parents=True, exist_ok=True)
        paths = []
        for folder in (FIG, OUT):
            png = folder / f"{stem}.png"
            pdf = folder / f"{stem}.pdf"
            fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
            fig.savefig(pdf, bbox_inches="tight", facecolor="white")
            paths.extend([png, pdf])
        return paths

    def single(rows, ylabel: str, stem: str) -> dict:
        fig, ax = plt.subplots(figsize=(4.6, 3.8))
        fit = panel(ax, rows, ylabel)
        fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
        save(fig, stem)
        plt.close(fig)
        return fit

    def species_means(blocks: list[list[tuple]]) -> list[tuple]:
        acc: dict[str, list] = {}
        mass: dict[str, float] = {}
        for rows in blocks:
            for name, m, frac in rows:
                acc.setdefault(name, []).append(frac)
                mass[name] = m
        return [(n, mass[n], float(np.mean(acc[n]))) for n in mass]

    def glued() -> dict:
        all_pts = [("wagner", *r) for r in WAGNER] + [("gse167580", *r) for r in GSE167580_MAMMALS]
        x = np.log10(np.array([r[2] for r in all_pts], float))
        y = np.log10(np.array([r[3] for r in all_pts], float))
        fit_all = stats.linregress(x, y)
        means = species_means([WAGNER, GSE167580_MAMMALS])
        xm = np.log10(np.array([r[1] for r in means], float))
        ym = np.log10(np.array([r[2] for r in means], float))
        fit_mean = stats.linregress(xm, ym)
        fig, ax = plt.subplots(figsize=(5.2, 4.2))
        xs = np.linspace(x.min(), x.max(), 80)
        ax.plot(xs, fit_all.intercept + fit_all.slope * xs, color=LINE, lw=1.6, zorder=2)
        styles = {
            "wagner": dict(c=POINT, marker="o", s=28),
            "gse167580": dict(c="white", marker="s", s=26, edgecolors=POINT, linewidths=1.1),
        }
        texts = []
        for study, name, mass, frac in all_pts:
            xi, yi = (np.log10(mass), np.log10(frac))
            ax.scatter([xi], [yi], zorder=3, **styles[study])
            texts.append(ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4))
        ax.scatter([], [], c=POINT, marker="o", s=28, label="Wagner SF")
        ax.scatter(
            [], [], c="white", marker="s", s=26, edgecolors=POINT, linewidths=1.1, label="GSE167580"
        )
        adjust_text(
            texts,
            ax=ax,
            x=x,
            y=y,
            expand=(1.12, 1.2),
            force_text=(0.22, 0.35),
            force_static=(0.12, 0.18),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )
        theme_light(ax)
        ax.set_xlabel("log10_body_mass")
        ax.set_ylabel("log10_fib")
        ax.legend(frameon=False, fontsize=8, loc="lower left")
        ptxt = f"{fit_all.pvalue:.3f}" if fit_all.pvalue >= 0.001 else f"{fit_all.pvalue:.3g}"
        ax.set_title(f"cor. {fit_all.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
        fig.subplots_adjust(left=0.15, right=0.97, top=0.88, bottom=0.15)
        save(fig, "fig_ruslan_style_wagner_gse")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(4.8, 4.0))
        panel(ax, means, "log10_fib")
        fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
        save(fig, "fig_ruslan_style_wagner_gse_mean")
        plt.close(fig)
        return {
            "all_n": int(len(all_pts)),
            "all_r": float(fit_all.rvalue),
            "all_p": float(fit_all.pvalue),
            "mean_n": int(len(means)),
            "mean_r": float(fit_mean.rvalue),
            "mean_p": float(fit_mean.pvalue),
            "means": means,
        }

    def main() -> None:
        w = single(WAGNER, "log10_fib", "fig_ruslan_style_wagner")
        g = single(GSE167580_MAMMALS, "log10_fib", "fig_ruslan_style_gse167580")
        print("Wagner standalone", w)
        print("GSE167580 mammals-only standalone", g)
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
        wt = panel(axes[0], WAGNER, "log10_fib")
        gt = panel(axes[1], GSE167580_MAMMALS, "log10_fib")
        axes[0].set_title(f"Wagner SF   cor. {wt['r']:.3f}; p = {wt['p']:.3f}", fontsize=10, pad=8)
        axes[1].set_title(
            f"GSE167580 mammals   cor. {gt['r']:.3f}; p = {gt['p']:.3f}", fontsize=10, pad=8
        )
        fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.16)
        save(fig, "fig_ruslan_style_mammals")
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
        c = panel(axes[0], BIRDS_CAGE, "log10_fib")
        r = panel(axes[1], BIRDS_RNA, "log10_fib")
        axes[0].set_title(
            f"GSE213253 CAGE   cor. {c['r']:.3f}; p = {c['p']:.3f}", fontsize=10, pad=8
        )
        axes[1].set_title(
            f"GSE213253 RNA   cor. {r['r']:.3f}; p = {r['p']:.3f}", fontsize=10, pad=8
        )
        fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.16)
        save(fig, "fig_ruslan_style_birds")
        plt.close(fig)
        fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.6))
        specs = [
            (axes[0, 0], BIRDS_CAGE, "log10_fib", "birds CAGE"),
            (axes[0, 1], BIRDS_RNA, "log10_fib", "birds RNA"),
            (axes[1, 0], WAGNER, "log10_fib", "Wagner SF"),
            (axes[1, 1], GSE167580_MAMMALS, "log10_fib", "GSE167580 mammals"),
        ]
        for ax, rows, ylab, name in specs:
            fit = panel(ax, rows, ylab)
            ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            ax.set_title(f"{name}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
        fig.subplots_adjust(wspace=0.32, hspace=0.38, left=0.08, right=0.99, top=0.95, bottom=0.08)
        save(fig, "fig_ruslan_style_four")
        plt.close(fig)
        print("Wagner two-up", wt)
        print("GSE167580 mammals two-up", gt)
        print("birds CAGE log10(sum TPM)~log10(mass)", c)
        print("birds RNA log10(sum TPM)~log10(mass)", r)

    return SimpleNamespace(**locals())


def _holding_cores(context):
    """Holding cores inside kept KEGG maps and inside each WGCNA module.

    The 16 OXPHOS genes were a complete-case set whose sum tracks mass. This
    asks the same question inside four maps the user kept, and inside each
    fibroblast WGCNA module, without mixing maps or modules.

    Rule (frozen before looking): among genes detected in all seven species,
    keep those whose own Pearson r with log10 mass has the same sign as the
    parent set. Score log10(mean TPM) of that core vs log10 mass against 2000
    abundance-matched nulls. That filter uses the mass axis, so it is a
    diagnostic, not a new a priori map. Calibration: apply it to hsa00190 and
    measure overlap with the known 16.

    Literature subsets (official lists or named cassettes, not a peel):
      thermogenesis minus OXPHOS, Hippo cassette, thyroid-specific receptors.

    Module hubs: kME_own ≥ 0.7, not selected by mass.

    Do not put fraction on the log10 axis. Do not glue maps or modules.

    Source: analyses/avian/scripts/78_holding_cores.py at d230e56ab942."""
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    DATA = context.unit("avian.oxphos_thread_lib").DATA
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    load_geneset = context.unit("avian.oxphos_thread_lib").load_geneset
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    OUT = RES / "holding_cores"
    M_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
    KEGG = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
    MODS = RES / "wgcna_out_v2" / "gene_modules_fib.tsv"
    LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
    N_NULL = 2000
    RNG = np.random.default_rng(20260921)
    MIN_N = 8
    SIXTEEN = [
        "ATP5F1E",
        "ATP6V0E1",
        "ATP6V0E2",
        "COX6A1",
        "COX6C",
        "NDUFA5",
        "NDUFA7",
        "NDUFAB1",
        "NDUFB2",
        "NDUFB5",
        "NDUFB6",
        "NDUFB8",
        "NDUFC2",
        "SDHD",
        "UQCR10",
        "UQCRC1",
    ]
    YC = LOGM - LOGM.mean()
    YSS = float(YC @ YC)
    HIPPO_CASSETTE = {
        "STK3",
        "STK4",
        "SAV1",
        "LATS1",
        "LATS2",
        "MOB1A",
        "MOB1B",
        "YAP1",
        "WWTR1",
        "TEAD1",
        "TEAD2",
        "TEAD3",
        "TEAD4",
        "NF2",
        "WWC1",
        "FRMD6",
    }
    THYROID_SPECIFIC = {"THRA", "THRB", "DIO1", "DIO2", "DIO3", "SLC16A2", "SLC16A10", "SLCO1C1"}
    THYROID_COREQ = {
        "NCOR1",
        "NCOA1",
        "NCOA2",
        "NCOA3",
        "RXRA",
        "RXRB",
        "RXRG",
        "HDAC1",
        "HDAC2",
        "HDAC3",
        "CREBBP",
        "EP300",
        "SIN3A",
        "KAT2A",
        "KAT2B",
        "MED1",
        "MED4",
        "MED12",
        "MED12L",
        "MED13",
        "MED13L",
        "MED14",
        "MED16",
        "MED17",
        "MED24",
        "MED27",
        "MED30",
    }
    PI3K_AKT = {
        "AKT1",
        "AKT2",
        "AKT3",
        "PIK3CA",
        "PIK3CB",
        "PIK3CD",
        "PIK3R1",
        "PIK3R2",
        "PIK3R3",
        "MTOR",
        "RHEB",
        "TSC2",
        "PDPK1",
        "FOXO1",
        "GSK3B",
        "BAD",
        "CASP9",
    }
    MAPK_RAS = {
        "HRAS",
        "KRAS",
        "NRAS",
        "RAF1",
        "MAP2K1",
        "MAP2K2",
        "MAPK1",
        "MAPK3",
        "SRC",
        "STAT1",
        "MYC",
        "CCND1",
        "HIF1A",
        "TP53",
        "MDM2",
    }

    def r_of_y(y: np.ndarray) -> float:
        z = y - y.mean()
        den = np.sqrt(float(z @ z) * YSS)
        return float(z @ YC / den) if den else np.nan

    def fit(y) -> tuple[float, float]:
        r = stats.linregress(LOGM, np.asarray(y, float))
        return (float(r.rvalue), float(r.pvalue))

    class AbundanceNull:

        def __init__(self, wide: pd.DataFrame):
            self.genes = np.asarray(wide.index)
            self.pos = {g: i for i, g in enumerate(self.genes)}
            self.mat = wide.to_numpy(float)
            chick = np.log10(self.mat[:, SPECIES_ORDER.index("Chicken")].clip(1e-12))
            self.dec = pd.qcut(
                pd.Series(chick).rank(method="first"), 10, labels=False, duplicates="drop"
            ).to_numpy()
            self.buckets = [np.flatnonzero(self.dec == k) for k in range(10)]

        def rows_of(self, genes: list[str]) -> np.ndarray:
            return np.array([self.pos[g] for g in genes if g in self.pos], dtype=int)

        def r_mean(self, rows: np.ndarray) -> float:
            return r_of_y(np.log10(self.mat[rows].mean(axis=0)))

        def y_mean(self, genes: list[str]) -> np.ndarray:
            rows = self.rows_of(genes)
            return self.mat[rows].mean(axis=0)

        def matched(self, rows: np.ndarray, obs_r: float) -> tuple[float, float]:
            if len(rows) < MIN_N or not np.isfinite(obs_r):
                return (np.nan, np.nan)
            want = np.bincount(self.dec[rows], minlength=10)
            null = np.empty(N_NULL)
            for i in range(N_NULL):
                parts = []
                for k, need in enumerate(want):
                    if need == 0:
                        continue
                    av = self.buckets[k]
                    parts.append(RNG.choice(av, size=min(int(need), len(av)), replace=False))
                pick = np.concatenate(parts)
                null[i] = self.r_mean(pick)
            null = null[np.isfinite(null)]
            return (float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null)))

    def kegg_genes(catalog: pd.DataFrame, pid: str) -> list[str]:
        rec = catalog.loc[catalog.pathway == pid].iloc[0]
        genes = [normalize_symbol(x) for x in str(rec.genes).split(";") if x.strip()]
        return [g for g in genes if g]

    def thyroid_role(g: str) -> str:
        if g in THYROID_SPECIFIC:
            return "рецептор / дейодиназа / транспорт T3"
        if g in THYROID_COREQ or g.startswith("MED"):
            return "корегулятор транскрипции"
        if g in PI3K_AKT:
            return "PI3K–AKT–mTOR"
        if g in MAPK_RAS:
            return "MAPK / RAS / эффект"
        if (
            g.startswith("ATP1")
            or g.startswith("ATP2")
            or g.startswith("SLC")
            or g.startswith("PLC")
        ):
            return "ионный / метаболический эффектор"
        return "прочее на карте"

    def hippo_role(g: str) -> str:
        if g in HIPPO_CASSETTE:
            return "кассета Hippo (MST–LATS–YAP/TAZ–TEAD)"
        if (
            g.startswith("WNT")
            or g.startswith("FZD")
            or g.startswith("DVL")
            or (
                g
                in {
                    "CTNNB1",
                    "AXIN1",
                    "AXIN2",
                    "APC",
                    "APC2",
                    "TCF7",
                    "TCF7L1",
                    "TCF7L2",
                    "LEF1",
                    "GSK3B",
                }
            )
        ):
            return "WNT на карте"
        if (
            g.startswith("TGF")
            or g.startswith("BMP")
            or g.startswith("SMAD")
            or g.startswith("GDF")
        ):
            return "TGF/BMP на карте"
        return "адгезия / прочее"

    def gene_r(null: AbundanceNull, g: str) -> tuple[float, float]:
        y = np.log10(null.mat[null.pos[g]])
        return fit(y)

    def score_set(null: AbundanceNull, genes: list[str]) -> dict:
        rows = null.rows_of(genes)
        n = int(len(rows))
        if n == 0:
            return dict(n=0, r=np.nan, p=np.nan, p_нуля=np.nan, медиана_нуля=np.nan, y=[np.nan] * 7)
        y = null.mat[rows].mean(axis=0)
        r, p = fit(np.log10(y))
        pn, nm = null.matched(rows, r) if n >= MIN_N else (np.nan, np.nan)
        return dict(n=n, r=r, p=p, p_нуля=pn, медиана_нуля=nm, y=y.tolist())

    def same_sign_core(null: AbundanceNull, genes: list[str], parent_r: float) -> list[str]:
        if not np.isfinite(parent_r) or parent_r == 0:
            return []
        want = np.sign(parent_r)
        keep = []
        for g in genes:
            rg, _ = gene_r(null, g)
            if np.isfinite(rg) and np.sign(rg) == want:
                keep.append(g)
        return keep

    def loo_table(null: AbundanceNull, genes: list[str], parent_r: float) -> pd.DataFrame:
        rows = []
        parent_rows = null.rows_of(genes)
        parent_y = null.mat[parent_rows].mean(axis=0)
        for g in genes:
            rg, pg = gene_r(null, g)
            rest = [x for x in genes if x != g]
            r_loo, _ = fit(np.log10(null.y_mean(rest)))
            d_abs = abs(parent_r) - abs(r_loo)
            rows.append(
                dict(
                    gene=g,
                    r_гена=rg,
                    p_гена=pg,
                    r_без=r_loo,
                    d_abs=d_abs,
                    держит=d_abs > 0,
                    chick_tpm=float(null.mat[null.pos[g], SPECIES_ORDER.index("Chicken")]),
                )
            )
        t = pd.DataFrame(rows).sort_values("d_abs", ascending=False)
        t["parent_r"] = parent_r
        t["parent_mean"] = float(parent_y.mean())
        return t

    def parse_kegg(catalog: pd.DataFrame) -> dict[str, list[str]]:
        ox = kegg_genes(catalog, "hsa00190")
        th = kegg_genes(catalog, "hsa04714")
        return {
            "hsa00190": ox,
            "hsa04714": th,
            "hsa03008": kegg_genes(catalog, "hsa03008"),
            "hsa04390": kegg_genes(catalog, "hsa04390"),
            "hsa04919": kegg_genes(catalog, "hsa04919"),
            "04714_minus_00190": sorted(set(th) - set(ox)),
            "hippo_cassette": sorted(HIPPO_CASSETTE),
            "thyroid_specific": sorted(THYROID_SPECIFIC),
            "set16": list(SIXTEEN),
        }

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        raw = pd.read_csv(M_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
        det = raw.notna() & (raw > 0)
        wide = raw.loc[det.all(axis=1)]
        print(f"вселенная {len(wide)} генов, CAGE во всех семи", flush=True)
        null = AbundanceNull(wide)
        catalog = pd.read_csv(KEGG, sep="\t")
        lists = parse_kegg(catalog)
        labels = {
            "hsa00190": "OXPHOS (калибровка)",
            "hsa04714": "термогенез",
            "hsa03008": "биогенез рибосом",
            "hsa04390": "Hippo",
            "hsa04919": "щитовидка, сигнал",
            "04714_minus_00190": "термогенез минус OXPHOS",
            "hippo_cassette": "кассета Hippo",
            "thyroid_specific": "рецепторы / дейодиназы / MCT",
            "set16": "известные 16",
        }
        parent_rows = []
        gene_rows = []
        y_rows = []
        map_ids = ["hsa00190", "hsa04714", "hsa03008", "hsa04390", "hsa04919"]
        extra_ids = ["set16", "04714_minus_00190", "hippo_cassette", "thyroid_specific"]

        def present(genes: list[str]) -> list[str]:
            return [g for g in genes if g in null.pos]

        for sid in map_ids + extra_ids:
            hit = present(lists[sid])
            sc = score_set(null, hit)
            parent_rows.append(
                dict(
                    id=sid,
                    набор=labels[sid],
                    kind="parent",
                    n=sc["n"],
                    r=sc["r"],
                    p=sc["p"],
                    p_нуля=sc["p_нуля"],
                    медиана_нуля=sc["медиана_нуля"],
                    overlap16=len(set(hit) & set(SIXTEEN)),
                )
            )
            for sp, val in zip(SPECIES_ORDER, sc["y"]):
                y_rows.append(dict(id=sid, kind="parent", species=sp, mean_tpm=val))
            if sid in map_ids and sc["n"] >= 3:
                loo = loo_table(null, hit, sc["r"])
                loo["id"] = sid
                loo["набор"] = labels[sid]
                if sid == "hsa04919":
                    loo["роль"] = loo["gene"].map(thyroid_role)
                elif sid == "hsa04390":
                    loo["роль"] = loo["gene"].map(hippo_role)
                elif sid == "hsa04714":
                    loo["роль"] = np.where(
                        loo["gene"].isin(set(lists["hsa00190"])), "OXPHOS на карте", "не OXPHOS"
                    )
                elif sid == "hsa00190":
                    loo["роль"] = np.where(loo["gene"].isin(SIXTEEN), "из 16", "остальной OXPHOS")
                else:
                    loo["роль"] = "биогенез рибосом"
                gene_rows.append(loo)
                core = same_sign_core(null, hit, sc["r"])
                csc = score_set(null, core)
                parent_rows.append(
                    dict(
                        id=sid,
                        набор=labels[sid] + " — ядро того же знака",
                        kind="same_sign",
                        n=csc["n"],
                        r=csc["r"],
                        p=csc["p"],
                        p_нуля=csc["p_нуля"],
                        медиана_нуля=csc["медиана_нуля"],
                        overlap16=len(set(core) & set(SIXTEEN)),
                    )
                )
                for sp, val in zip(SPECIES_ORDER, csc["y"]):
                    y_rows.append(dict(id=sid, kind="same_sign", species=sp, mean_tpm=val))
                print(
                    f"  {sid:12s} parent n={sc['n']:3d} r={sc['r']:+.3f} нуль {sc['p_нуля']:.3f}  | core n={csc['n']:3d} r={csc['r']:+.3f} нуль {csc['p_нуля']:.3f}  ∩16={len(set(core) & set(SIXTEEN))}",
                    flush=True,
                )
        cores = pd.DataFrame(parent_rows)
        genes = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
        ys = pd.DataFrame(y_rows)
        ox_core = same_sign_core(
            null,
            present(lists["hsa00190"]),
            float(cores.loc[(cores.id == "hsa00190") & (cores.kind == "parent"), "r"].iloc[0]),
        )
        print(
            f"\nкалибровка 00190 same-sign: {len(ox_core)} генов, из 16 попало {len(set(ox_core) & set(SIXTEEN))}/16, лишних {len(set(ox_core) - set(SIXTEEN))}",
            flush=True,
        )
        gm = pd.read_csv(MODS, sep="\t")
        mrows = []
        mgenes = []
        for mod, gtab in gm.groupby("module"):
            if mod == "grey":
                continue
            hit = [x for x in gtab["gene"] if x in null.pos]
            if len(hit) < MIN_N:
                continue
            kme = gtab.set_index("gene")["kME_own"]
            sc = score_set(null, hit)
            core = same_sign_core(null, hit, sc["r"])
            csc = score_set(null, core)
            hubs = [g for g in hit if g in kme.index and kme[g] >= 0.7]
            hsc = score_set(null, hubs)
            mrows.append(
                dict(
                    модуль=mod,
                    n_parent=sc["n"],
                    r_parent=sc["r"],
                    p_нуля_parent=sc["p_нуля"],
                    n_core=csc["n"],
                    r_core=csc["r"],
                    p_нуля_core=csc["p_нуля"],
                    n_hub=hsc["n"],
                    r_hub=hsc["r"],
                    p_нуля_hub=hsc["p_нуля"],
                    overlap16_core=len(set(core) & set(SIXTEEN)),
                    overlap16_hub=len(set(hubs) & set(SIXTEEN)),
                )
            )
            for g in hit:
                rg, pg = gene_r(null, g)
                mgenes.append(
                    dict(
                        модуль=mod,
                        gene=g,
                        r_гена=rg,
                        p_гена=pg,
                        kME=float(kme[g]) if g in kme.index else np.nan,
                        same_sign=g in core,
                        hub=g in hubs,
                    )
                )
            if (
                sc["p_нуля"] < 0.05
                or csc["p_нуля"] < 0.05
                or abs(sc["r"]) >= 0.7
                or (abs(csc["r"]) >= 0.85)
            ):
                print(
                    f"  mod {mod:16s} n={sc['n']:3d} r={sc['r']:+.3f} нуль {sc['p_нуля']:.3f}  | core n={csc['n']:3d} r={csc['r']:+.3f} нуль {csc['p_нуля']:.3f}  | hub n={hsc['n']:3d} r={hsc['r']:+.3f}",
                    flush=True,
                )
        mt = pd.DataFrame(mrows).sort_values("r_core")
        mt["bh_core"] = stats.false_discovery_control(
            mt["p_нуля_core"].fillna(1.0).to_numpy(), method="bh"
        )
        mt["bh_hub"] = stats.false_discovery_control(
            mt["p_нуля_hub"].fillna(1.0).to_numpy(), method="bh"
        )
        mt["bh_parent"] = stats.false_discovery_control(
            mt["p_нуля_parent"].fillna(1.0).to_numpy(), method="bh"
        )
        map_cores = cores[cores.kind == "same_sign"].copy()
        map_cores["bh"] = stats.false_discovery_control(
            map_cores["p_нуля"].fillna(1.0).to_numpy(), method="bh"
        )
        cores = cores.merge(map_cores[["id", "kind", "bh"]], on=["id", "kind"], how="left")
        chick = SPECIES_ORDER.index("Chicken")
        extra = []
        extra_y = []
        extra_genes = []
        for sid in map_ids:
            hit = present(lists[sid])
            tpm = {g: float(null.mat[null.pos[g], chick]) for g in hit}
            top = sorted(hit, key=lambda g: -tpm[g])[:16]
            sc = score_set(null, top)
            extra.append(
                dict(
                    id=sid,
                    набор=labels[sid] + " — 16 самых обильных у курицы",
                    kind="abund16",
                    n=sc["n"],
                    r=sc["r"],
                    p=sc["p"],
                    p_нуля=sc["p_нуля"],
                    медиана_нуля=sc["медиана_нуля"],
                    overlap16=len(set(top) & set(SIXTEEN)),
                )
            )
            extra_genes.append(
                pd.DataFrame(
                    dict(id=sid, kind="abund16", gene=top, chick_tpm=[tpm[g] for g in top])
                )
            )
            for sp, val in zip(SPECIES_ORDER, sc["y"]):
                extra_y.append(dict(id=sid, kind="abund16", species=sp, mean_tpm=val))
            if "ACTB" in hit and sid in {"hsa04390", "hsa04919", "hsa04714"}:
                rest = [g for g in hit if g != "ACTB"]
                scb = score_set(null, rest)
                extra.append(
                    dict(
                        id=sid,
                        набор=labels[sid] + " без ACTB",
                        kind="drop_actb",
                        n=scb["n"],
                        r=scb["r"],
                        p=scb["p"],
                        p_нуля=scb["p_нуля"],
                        медиана_нуля=scb["медиана_нуля"],
                        overlap16=len(set(rest) & set(SIXTEEN)),
                    )
                )
                for sp, val in zip(SPECIES_ORDER, scb["y"]):
                    extra_y.append(dict(id=sid, kind="drop_actb", species=sp, mean_tpm=val))
        hub16_rows = []
        for mod, gtab in gm.groupby("module"):
            if mod == "grey":
                continue
            hit = [x for x in gtab["gene"] if x in null.pos]
            if len(hit) < 16:
                continue
            kme = gtab.set_index("gene")["kME_own"]
            top = sorted(hit, key=lambda g: -float(kme.get(g, 0)))[:16]
            sc = score_set(null, top)
            hub16_rows.append(
                dict(
                    модуль=mod,
                    n=16,
                    r=sc["r"],
                    p=sc["p"],
                    p_нуля=sc["p_нуля"],
                    overlap16=len(set(top) & set(SIXTEEN)),
                    гены=";".join(top),
                )
            )
        h16 = pd.DataFrame(hub16_rows).sort_values("r")
        h16["bh"] = stats.false_discovery_control(h16["p_нуля"].fillna(1.0).to_numpy(), method="bh")
        cores = pd.concat([cores, pd.DataFrame(extra)], ignore_index=True)
        ys = pd.concat([ys, pd.DataFrame(extra_y)], ignore_index=True)
        pd.concat(extra_genes, ignore_index=True).to_csv(
            OUT / "map_abund16_genes.tsv", sep="\t", index=False
        )
        h16.to_csv(OUT / "module_hub16_vs_mass.tsv", sep="\t", index=False)
        cores.to_csv(OUT / "map_cores_vs_mass.tsv", sep="\t", index=False)
        genes.to_csv(OUT / "map_genes_loo.tsv", sep="\t", index=False)
        ys.to_csv(OUT / "map_species_means.tsv", sep="\t", index=False)
        mt.to_csv(OUT / "module_cores_vs_mass.tsv", sep="\t", index=False)
        pd.DataFrame(mgenes).to_csv(OUT / "module_genes.tsv", sep="\t", index=False)
        thy = genes[genes.id == "hsa04919"].copy()
        print("\nщитовидка, сигнал — гены complete-case:")
        for rec in thy.sort_values("r_гена", ascending=False).itertuples():
            print(
                f"  {rec.gene:12s} r={rec.r_гена:+.3f}  {rec.роль}  {('держит' if rec.держит else '')}"
            )
        print(f"\nмодули: same-sign BH<0.05: {(mt.bh_core < 0.05).sum()} / {len(mt)}")
        print(f"модули: hub     BH<0.05: {(mt.bh_hub < 0.05).sum()} / {len(mt)}")
        print(f"модули: hub-16  BH<0.05: {(h16.bh < 0.05).sum()} / {len(h16)}")
        print(f"карты:  same-sign BH<0.05: {(map_cores.bh < 0.05).sum()} / {len(map_cores)}")
        print("\nabund16 / без ACTB:")
        for rec in cores[cores.kind.isin(["abund16", "drop_actb"])].itertuples():
            print(
                f"  {rec.kind:10s} {rec.id:12s} n={rec.n:3d} r={rec.r:+.3f} нуль {rec.p_нуля:.3f} ∩16={rec.overlap16}"
            )
        print("\nhub-16 модулей |r|≥0.65:")
        for rec in h16[h16.r.abs() >= 0.65].itertuples():
            print(f"  {rec.модуль:16s} r={rec.r:+.3f} нуль {rec.p_нуля:.3f} BH={rec.bh:.3f}")
        print(f"\nwrote {OUT}")

    return SimpleNamespace(**locals())


def _story_epi_and_traits(context):
    """Story sets in epiblast, and the same sums vs Ruslan's other three traits.

    Same homologous-CAGE assignment as scripts 75/79 (nearest chicken-protein CDS
    ≤ 2 kb, peak kept if fibroblast mean > 0). Epiblast is the mean of those same
    peaks' norm.E libraries, not a new promoter call. Missing ≠ 0 at assignment;
    a gene already in the fibroblast matrix may have epi = 0 at that promoter.

    y = log10(sum TPM of the set). Predictors from birds-ESC.R species_traits:
    log10 mass, incubation_d, max_longevity_y, log10 genome_size_gb. Tissues
    and traits are separate panels. Do not glue sets or studies.

    Source: analyses/avian/scripts/80_story_epi_and_traits.py at d230e56ab942."""
    import importlib.util
    import json
    import sys
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import stats

    s74 = context.unit("avian.74_nineteen_via_cage")
    s75 = context.unit("avian.75_honest_gene_mass_scan")
    DATA = context.unit("avian.oxphos_thread_lib").DATA
    FIG = context.unit("avian.oxphos_thread_lib").FIG
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    OUT = RES / "story_epi_traits"
    M_FIB_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
    THR = 1.0
    SETS = [
        ("oxphos_16", "OXPHOS 16", "fig_story_oxphos_16"),
        ("folate_3", "folate 3", "fig_story_folate_3"),
        ("ribo_5", "ribo 5", "fig_story_ribo_5"),
        ("green_hubs_7", "green hubs 7", "fig_story_green_hubs_7"),
        ("green_sig_51", "green 51", "fig_story_green_51"),
    ]
    TRAITS = {
        "Zebra": dict(mass_g=11.93, incubation_d=14.15, max_longevity_y=14.5, genome_size_gb=1.223),
        "Quail": dict(mass_g=96.6, incubation_d=17.0, max_longevity_y=14.58, genome_size_gb=1.379),
        "Chicken": dict(
            mass_g=912.5, incubation_d=20.5, max_longevity_y=30.0, genome_size_gb=1.223
        ),
        "Duck": dict(mass_g=1121.0, incubation_d=26.5, max_longevity_y=29.1, genome_size_gb=1.428),
        "Peacock": dict(
            mass_g=4093.75, incubation_d=28.0, max_longevity_y=25.0, genome_size_gb=1.13
        ),
        "Turkey": dict(
            mass_g=5811.0, incubation_d=28.0, max_longevity_y=13.0, genome_size_gb=1.281
        ),
        "Emu": dict(mass_g=36200.0, incubation_d=56.0, max_longevity_y=16.6, genome_size_gb=1.594),
    }
    PREDICTORS = [
        ("mass", "log10_body_mass", lambda t: np.log10(t["mass_g"])),
        ("incubation", "incubation_d", lambda t: float(t["incubation_d"])),
        ("longevity", "max_longevity_y", lambda t: float(t["max_longevity_y"])),
        ("genome", "log10_genome_size", lambda t: np.log10(t["genome_size_gb"])),
    ]

    def load40():
        mod = context.unit("avian.40_ruslan_style_scatters")
        return mod

    def genes_of(stem: str) -> list[str]:
        return [
            ln.strip()
            for ln in (DATA / "genesets" / f"story_{stem}.txt").read_text().splitlines()
            if ln.strip()
        ]

    def load_peaks_both(species: str) -> pd.DataFrame:
        """Same fibroblast filter as script 74, plus mean of E columns on those peaks."""
        d = pd.read_csv(s74.s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        fcols = lib_cols(d, "fib")
        ecols = lib_cols(d, "epi")
        fnum = d[fcols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        fnum = fnum.where(fnum >= THR, 0.0)
        if ecols:
            enum = d[ecols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
            enum = enum.where(enum >= THR, 0.0)
            epi = enum.mean(axis=1)
        else:
            epi = pd.Series(0.0, index=d.index)
        d = d.assign(
            mid=(d["start"] + d["end"]) / 2.0,
            fib=fnum.mean(axis=1),
            epi=epi,
            seq=d["seqnames"].map(s74.contig),
            n_fib=len(fcols),
            n_epi=len(ecols),
        )
        return d[d["fib"] > 0].copy()

    def assign_both(peaks: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
        """Identical catchment to script 75; also sums epiblast on the same peaks."""
        if loci.empty or peaks.empty:
            return pd.DataFrame(columns=["gene", "fib", "epi", "n_peaks", "dist"])
        index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
        for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
            g = g.sort_values("five")
            index[seq, st] = (g["five"].to_numpy(float), g["symbol"].tolist())
        acc: dict[str, list] = {}
        for r in peaks.itertuples():
            packed = index.get((r.seq, r.strand))
            if packed is None:
                continue
            fives, syms = packed
            mid = float(r.mid)
            i = int(np.searchsorted(fives, mid))
            cands: list[int] = []
            if i > 0:
                cands.append(i - 1)
            if i < len(fives):
                cands.append(i)
            if not cands:
                continue
            ds = [abs(mid - fives[j]) for j in cands]
            jbest = cands[int(np.argmin(ds))]
            d = abs(mid - fives[jbest])
            if d > s75.MAX_BP:
                continue
            if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-06:
                continue
            gene = syms[jbest]
            acc.setdefault(gene, []).append((float(r.fib), float(r.epi), d))
        rows = []
        for gene, hits in acc.items():
            fib = float(sum((h[0] for h in hits)))
            if fib <= 0:
                continue
            rows.append(
                dict(
                    gene=gene,
                    fib=fib,
                    epi=float(sum((h[1] for h in hits))),
                    n_peaks=len(hits),
                    dist=min((h[2] for h in hits)),
                )
            )
        return pd.DataFrame(rows)

    def scatter(ax, names, x, y, xlabel: str, ylabel: str, s40) -> dict:
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        fit = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 80)
        ax.plot(xs, fit.intercept + fit.slope * xs, color=s40.LINE, lw=1.6, zorder=2)
        ax.scatter(x, y, s=22, c=s40.POINT, zorder=3, linewidths=0)
        texts = [
            ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
            for name, xi, yi in zip(names, x, y)
        ]
        s40.adjust_text(
            texts,
            ax=ax,
            x=x,
            y=y,
            expand=(1.15, 1.25),
            force_text=(0.25, 0.4),
            force_static=(0.15, 0.2),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )
        s40.theme_light(ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
        ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
        return {
            "r": float(fit.rvalue),
            "p": float(fit.pvalue),
            "slope": float(fit.slope),
            "intercept": float(fit.intercept),
            "n": int(len(names)),
        }

    def save(fig, stem: str, s40) -> None:
        FIG.mkdir(parents=True, exist_ok=True)
        OUT.mkdir(parents=True, exist_ok=True)
        for folder in (FIG, OUT):
            fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
            fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")

    def set_sums(wide: pd.DataFrame, genes: list[str]) -> np.ndarray:
        hit = [g for g in genes if g in wide.index]
        if len(hit) != len(genes):
            raise SystemExit(f"missing from matrix: {set(genes) - set(hit)}")
        return wide.loc[hit, SPECIES_ORDER].sum(axis=0).to_numpy(float)

    def build_matrices() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        print("локусы куриного белка + те же пики, сумма E", flush=True)
        fib_s, epi_s, lib_rows = ({}, {}, [])
        for sp in SPECIES_ORDER:
            print(f"  {sp}", flush=True)
            loci = s75.all_chicken_loci(s75.ASM / sp / "reanno_named.gff")
            peaks = load_peaks_both(sp)
            n_fib = int(peaks["n_fib"].iloc[0]) if len(peaks) else 0
            n_epi = int(peaks["n_epi"].iloc[0]) if len(peaks) else 0
            assigned = assign_both(peaks, loci)
            fib_s[sp] = (
                assigned.set_index("gene")["fib"] if len(assigned) else pd.Series(dtype=float)
            )
            epi_s[sp] = (
                assigned.set_index("gene")["epi"] if len(assigned) else pd.Series(dtype=float)
            )
            lib_rows.append(
                dict(species=sp, n_fib_libs=n_fib, n_epi_libs=n_epi, n_genes=len(assigned))
            )
            print(f"    библиотек C={n_fib} E={n_epi}  генов {len(assigned)}", flush=True)
        Mf = pd.DataFrame(fib_s).reindex(columns=SPECIES_ORDER)
        Me = pd.DataFrame(epi_s).reindex(index=Mf.index, columns=SPECIES_ORDER)
        return (Mf, Me, pd.DataFrame(lib_rows))

    def five_grid(panels, ylabel: str, xlabel: str, stem: str, s40) -> None:
        fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
        axes = axes.ravel()
        for i, (title, names, x, y, fit) in enumerate(panels):
            scatter(axes[i], names, x, y, xlabel, ylabel, s40)
            ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            axes[i].set_title(f"{title}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
        axes[5].set_visible(False)
        fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.06, right=0.99, top=0.94, bottom=0.08)
        save(fig, stem, s40)
        plt.close(fig)

    def main() -> None:
        s40 = load40()
        OUT.mkdir(parents=True, exist_ok=True)
        names = [sp.lower() for sp in SPECIES_ORDER]
        Mf, Me, libs = build_matrices()
        old = pd.read_csv(M_FIB_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
        story_genes = []
        for stem, _, _ in SETS:
            story_genes.extend(genes_of(stem))
        story_genes = sorted(set(story_genes))
        rels = []
        for g in story_genes:
            a = Mf.loc[g].to_numpy(float)
            b = old.loc[g].to_numpy(float)
            rels.append(np.nanmax(np.abs(a - b) / np.clip(b, 1e-06, None)))
        max_rel = float(np.max(rels))
        print(
            f"сверка M_fib с скриптом 75 по {len(story_genes)} генам истории: max отн. {max_rel:.4g}"
        )
        if max_rel > 1e-06:
            raise SystemExit("fibroblast assignment drifted; stop")
        Mf.to_csv(OUT / "M_fib_nearest2kb_check.tsv", sep="\t")
        Me.to_csv(OUT / "M_epi_nearest2kb.tsv", sep="\t")
        libs.to_csv(OUT / "library_counts.tsv", sep="\t", index=False)
        X = {
            key: np.array([fn(TRAITS[sp]) for sp in SPECIES_ORDER], float)
            for key, _, fn in PREDICTORS
        }
        xlab = {key: lab for key, lab, _ in PREDICTORS}
        pred_r = []
        keys = [k for k, _, _ in PREDICTORS]
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                rr = stats.linregress(X[a], X[b])
                pred_r.append(dict(a=a, b=b, r=float(rr.rvalue), p=float(rr.pvalue)))
                print(f"коллинеарность {a} vs {b}: r={rr.rvalue:+.3f} p={rr.pvalue:.3f}")
        pd.DataFrame(pred_r).to_csv(OUT / "predictor_collinearity.tsv", sep="\t", index=False)
        recs, canvas = ([], [])
        tissues = [("fib", Mf, "log10_fib"), ("epi", Me, "log10_epi")]
        for stem, title, figstem in SETS:
            genes = genes_of(stem)
            for tiss, wide, ylab in tissues:
                sums = set_sums(wide, genes)
                n_zero = int((wide.loc[genes] == 0).to_numpy().sum())
                n_na = int(wide.loc[genes].isna().to_numpy().sum())
                for key, xlab_i, _ in PREDICTORS:
                    y = np.log10(sums)
                    if not np.isfinite(y).all():
                        raise SystemExit(f"non-finite log10 sum: {stem} {tiss} {key} {sums}")
                    fit = stats.linregress(X[key], y)
                    recs.append(
                        dict(
                            set=stem,
                            title=title,
                            n=len(genes),
                            tissue=tiss,
                            predictor=key,
                            xlabel=xlab_i,
                            ylabel=ylab,
                            r=float(fit.rvalue),
                            p=float(fit.pvalue),
                            slope=float(fit.slope),
                            intercept=float(fit.intercept),
                            n_zero_cells=n_zero,
                            n_na_cells=n_na,
                            **{
                                f"sum_{sp.lower()}": float(sums[i])
                                for i, sp in enumerate(SPECIES_ORDER)
                            },
                            **{
                                f"y_{sp.lower()}": float(y[i]) for i, sp in enumerate(SPECIES_ORDER)
                            },
                        )
                    )
                    canvas.append(
                        dict(
                            set=stem,
                            title=title,
                            tissue=tiss,
                            predictor=key,
                            xlabel=xlab_i,
                            ylabel=ylab,
                            r=float(fit.rvalue),
                            p=float(fit.pvalue),
                            slope=float(fit.slope),
                            intercept=float(fit.intercept),
                            x=X[key].tolist(),
                            y=y.tolist(),
                            names=names,
                        )
                    )
                print(f"{stem:16s} {tiss} zeros={n_zero} NA={n_na}", flush=True)
        fits = pd.DataFrame(recs)
        fits.to_csv(OUT / "story_set_trait_fits.tsv", sep="\t", index=False)
        (OUT / "canvas_panels.json").write_text(json.dumps(canvas, indent=2))
        epi_mass = []
        for stem, title, figstem in SETS:
            genes = genes_of(stem)
            sums = set_sums(Me, genes)
            y = np.log10(sums)
            x = X["mass"]
            fig, ax = plt.subplots(figsize=(4.6, 3.8))
            fit = scatter(ax, names, x, y, xlab["mass"], "log10_epi", s40)
            fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
            save(fig, f"{figstem}_epi", s40)
            plt.close(fig)
            epi_mass.append((title, names, x, y, fit))
            print(
                f"{stem:16s} epi vs mass  r={fit['r']:+.3f} p={fit['p']:.4f}  sums={[round(s, 1) for s in sums]}"
            )
        five_grid(epi_mass, "log10_epi", xlab["mass"], "fig_story_five_epi", s40)
        fig, axes = plt.subplots(5, 2, figsize=(8.4, 16.2))
        for i, (stem, title, figstem) in enumerate(SETS):
            genes = genes_of(stem)
            for j, (tiss, wide, ylab) in enumerate(tissues):
                y = np.log10(set_sums(wide, genes))
                fit = scatter(axes[i, j], names, X["mass"], y, xlab["mass"], ylab, s40)
                ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
                axes[i, j].set_title(
                    f"{title} {tiss}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=9, pad=6
                )
        fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.08, right=0.99, top=0.98, bottom=0.03)
        save(fig, "fig_story_five_fib_epi", s40)
        plt.close(fig)
        for key in ("incubation", "longevity", "genome"):
            for tiss, wide, ylab in tissues:
                panels = []
                for stem, title, _ in SETS:
                    genes = genes_of(stem)
                    y = np.log10(set_sums(wide, genes))
                    fit = dict(
                        r=float(stats.linregress(X[key], y).rvalue),
                        p=float(stats.linregress(X[key], y).pvalue),
                    )
                    panels.append((title, names, X[key], y, fit))
                five_grid(panels, ylab, xlab[key], f"fig_story_five_{tiss}_{key}", s40)
        fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.6))
        pred_labs = [lab for _, lab, _ in PREDICTORS]
        set_labs = [t for _, t, _ in SETS]
        for ax, tiss, ttl in zip(axes, ("fib", "epi"), ("фибробласты", "эпибласт")):
            mat = np.array(
                [
                    [
                        float(
                            fits[
                                (fits.set == stem) & (fits.tissue == tiss) & (fits.predictor == key)
                            ].r.iloc[0]
                        )
                        for key, _, _ in PREDICTORS
                    ]
                    for stem, _, _ in SETS
                ]
            )
            im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
            ax.set_xticks(range(4), pred_labs, rotation=25, ha="right", fontsize=8)
            ax.set_yticks(range(5), set_labs, fontsize=9)
            ax.set_title(ttl, fontsize=11)
            for i in range(5):
                for j in range(4):
                    ax.text(
                        j,
                        i,
                        f"{mat[i, j]:+.2f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="#111" if abs(mat[i, j]) < 0.55 else "white",
                    )
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.subplots_adjust(wspace=0.35, left=0.12, right=0.98, top=0.88, bottom=0.22)
        save(fig, "fig_story_r_heatmap_traits", s40)
        plt.close(fig)
        print(f"wrote {OUT}")
        show = fits[fits.predictor == "mass"][["set", "tissue", "r", "p"]]
        print(show.to_string(index=False))

    return SimpleNamespace(**locals())


def _story_cre_motifs(context):
    """Story-set CRE: proximal vs distal homologous CAGE, and motif enrichment.

    Promoter = CAGE whose midpoint is nearest this chicken-protein CDS 5' and
    ≤ 2 kb (the same peaks as scripts 75/79). Distal CAGE = same nearest-gene
    rule, 2–100 kb — not the chicken CRE-atlas enhancer class and not Ruslan's
    correlation links.

    Motifs: FIMO JASPAR 2024 CORE vertebrates on proximal peak sequence
    (≥400 bp window). Fisher: peaks of one story set vs other proximal peaks
    in that species. BH inside set × species. Do not take top-10 then correlate
    with mass (that is the step that made his TFBS-to-size plots).

    Sets stay separate. Tissues stay separate.

    Source: analyses/avian/scripts/81_story_cre_motifs.py at d230e56ab942."""
    import importlib.util
    import json
    import subprocess
    import sys
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import stats

    s74 = context.unit("avian.74_nineteen_via_cage")
    s75 = context.unit("avian.75_honest_gene_mass_scan")
    DATA = context.unit("avian.oxphos_thread_lib").DATA
    FIG = context.unit("avian.oxphos_thread_lib").FIG
    FIMO = context.unit("avian.oxphos_thread_lib").FIMO
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    OUT = RES / "story_cre_motifs"
    MEME = RES / "module_consistency" / "fimo_chicken" / "JASPAR2024_vertebrates.meme"
    PROX_BP = 2000
    DISTAL_BP = 100000
    MIN_WIN = 400
    FIMO_THRESH = "1e-4"
    SETS = [
        ("oxphos_16", "OXPHOS 16"),
        ("folate_3", "folate 3"),
        ("ribo_5", "ribo 5"),
        ("green_hubs_7", "green hubs 7"),
        ("green_sig_51", "green 51"),
    ]
    WATCH = ["SNAI3", "PHOX2B", "ZEB1", "SOX14", "MEF2B", "HEY2", "NRF1", "GABPA", "ESRRA", "YY1"]

    def load_py(stem: str, fname: str):
        return context.unit("avian." + Path(fname).stem)

    def genes_of(stem: str) -> list[str]:
        return [
            ln.strip()
            for ln in (DATA / "genesets" / f"story_{stem}.txt").read_text().splitlines()
            if ln.strip()
        ]

    def assign_peaks(peaks: pd.DataFrame, loci: pd.DataFrame, max_bp: int) -> pd.DataFrame:
        if loci.empty or peaks.empty:
            return pd.DataFrame()
        index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
        for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
            g = g.sort_values("five")
            index[seq, st] = (g["five"].to_numpy(float), g["symbol"].tolist())
        rows = []
        for r in peaks.itertuples():
            packed = index.get((r.seq, r.strand))
            if packed is None:
                continue
            fives, syms = packed
            mid = float(r.mid)
            i = int(np.searchsorted(fives, mid))
            cands: list[int] = []
            if i > 0:
                cands.append(i - 1)
            if i < len(fives):
                cands.append(i)
            if not cands:
                continue
            ds = [abs(mid - fives[j]) for j in cands]
            jbest = cands[int(np.argmin(ds))]
            d = abs(mid - fives[jbest])
            if d > max_bp:
                continue
            if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-06:
                continue
            rows.append(
                dict(
                    gene=syms[jbest],
                    peak_id=str(r.id),
                    seq=str(r.seq),
                    start=int(r.start),
                    end=int(r.end),
                    strand=str(r.strand),
                    fib=float(r.fib),
                    epi=float(r.epi),
                    dist=float(d),
                )
            )
        return pd.DataFrame(rows)

    def klass(dist: float) -> str:
        return "proximal" if dist <= PROX_BP else "distal"

    def rc(seq: str) -> str:
        tbl = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(tbl)[::-1]

    def fasta_alias(fa) -> dict[str, str]:
        """Map CAGE seqnames onto FASTA keys.

        Turkey and emu headers are ENA|HG999680|HG999680.1 while exp_tab
        stores HG999680.1. Exact key wins; otherwise the last pipe field.
        """
        alias: dict[str, str] = {}
        for key in fa.keys():
            alias.setdefault(key, key)
            alias.setdefault(key.split("|")[-1], key)
        return alias

    def write_fasta(peaks: pd.DataFrame, fa_path: Path, out_fa: Path) -> int:
        from pyfaidx import Fasta

        fa = Fasta(str(fa_path))
        alias = fasta_alias(fa)
        n = 0
        missed = 0
        with out_fa.open("w") as out:
            for r in peaks.itertuples():
                lo, hi = (int(r.start), int(r.end))
                if hi < lo:
                    lo, hi = (hi, lo)
                mid = (lo + hi) // 2
                half = max(MIN_WIN, hi - lo) // 2
                a = max(0, mid - half)
                b = mid + half
                name = alias.get(str(r.seq))
                if name is None:
                    missed += 1
                    continue
                b = min(b, len(fa[name]))
                if b - a < 20:
                    continue
                seq = str(fa[name][a:b])
                if str(r.strand) == "-":
                    seq = rc(seq)
                sid = f"{r.gene}|{r.peak_id}|{r.species}|{r.klass}"
                out.write(f">{sid}\n{seq}\n")
                n += 1
        if missed:
            print(f"  без последовательности {missed}", flush=True)
        return n

    def run_fimo(fasta: Path, out_tsv: Path) -> pd.DataFrame:
        if out_tsv.exists() and out_tsv.stat().st_size > 100:
            return pd.read_csv(out_tsv, sep="\t", comment="#")
        cmd = [
            str(FIMO),
            "--text",
            "--no-qvalue",
            "--verbosity",
            "1",
            "--thresh",
            FIMO_THRESH,
            str(MEME),
            str(fasta),
        ]
        with out_tsv.open("w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.DEVNULL)
        if out_tsv.stat().st_size < 50:
            return pd.DataFrame()
        return pd.read_csv(out_tsv, sep="\t", comment="#")

    def fisher_block(hit_seq: set[str], pos: set[str], neg: set[str]) -> dict:
        pp = len(pos & hit_seq)
        np_ = len(pos) - pp
        pn = len(neg & hit_seq)
        nn = len(neg) - pn
        table = np.array([[pp, pn], [np_, nn]], int)
        if table.min() < 0:
            return dict(odds=np.nan, p=1.0, pp=pp, pn=pn, np=np_, nn=nn)
        odds, p = stats.fisher_exact(table, alternative="two-sided")
        return dict(odds=float(odds), p=float(p), pp=int(pp), pn=int(pn), np=int(np_), nn=int(nn))

    def scatter(ax, names, x, y, xlabel, ylabel, s40) -> dict:
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        ok = np.isfinite(x) & np.isfinite(y)
        x, y, names = (x[ok], y[ok], [n for n, k in zip(names, ok) if k])
        if len(x) < 5:
            ax.set_visible(False)
            return {"r": np.nan, "p": np.nan, "n": int(len(x))}
        fit = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 80)
        ax.plot(xs, fit.intercept + fit.slope * xs, color=s40.LINE, lw=1.6, zorder=2)
        ax.scatter(x, y, s=22, c=s40.POINT, zorder=3, linewidths=0)
        texts = [
            ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
            for name, xi, yi in zip(names, x, y)
        ]
        s40.adjust_text(
            texts,
            ax=ax,
            x=x,
            y=y,
            expand=(1.15, 1.25),
            force_text=(0.25, 0.4),
            force_static=(0.15, 0.2),
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
        )
        s40.theme_light(ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
        ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
        return {
            "r": float(fit.rvalue),
            "p": float(fit.pvalue),
            "n": int(len(names)),
            "slope": float(fit.slope),
            "intercept": float(fit.intercept),
        }

    def save(fig, stem: str) -> None:
        FIG.mkdir(parents=True, exist_ok=True)
        OUT.mkdir(parents=True, exist_ok=True)
        for folder in (FIG, OUT):
            fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
            fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")

    def main() -> None:
        s40 = load_py("s40", "40_ruslan_style_scatters.py")
        s80 = load_py("s80", "80_story_epi_and_traits.py")
        OUT.mkdir(parents=True, exist_ok=True)
        story = {stem: genes_of(stem) for stem, _ in SETS}
        all_story = sorted({g for gs in story.values() for g in gs})
        names = [sp.lower() for sp in SPECIES_ORDER]
        logm = np.log10([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER])
        peak_path = OUT / "assigned_peaks.tsv"
        if peak_path.exists():
            peaks = pd.read_csv(peak_path, sep="\t")
            print(f"пики из кэша {peak_path}  n={len(peaks)}", flush=True)
        else:
            chunks = []
            for sp in SPECIES_ORDER:
                print(f"присвоение {sp}", flush=True)
                loci = s75.all_chicken_loci(s75.ASM / sp / "reanno_named.gff")
                raw = s80.load_peaks_both(sp)
                a = assign_peaks(raw, loci, DISTAL_BP)
                if a.empty:
                    continue
                a["species"] = sp
                a["klass"] = a["dist"].map(klass)
                chunks.append(a)
                n_pr = int((a.klass == "proximal").sum())
                n_ds = int((a.klass == "distal").sum())
                print(
                    f"  пиков proximal {n_pr}  distal {n_ds}  генов {a.gene.nunique()}", flush=True
                )
            peaks = pd.concat(chunks, ignore_index=True)
            peaks.to_csv(peak_path, sep="\t", index=False)
        old = pd.read_csv(RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv", sep="\t", index_col=0)
        prox = peaks[peaks.klass == "proximal"]
        rels = []
        for g in all_story:
            for sp in SPECIES_ORDER:
                got = float(prox.loc[(prox.gene == g) & (prox.species == sp), "fib"].sum())
                exp = float(old.loc[g, sp]) if g in old.index else np.nan
                if exp and np.isfinite(exp) and (exp > 0):
                    rels.append(abs(got - exp) / exp)
        print(
            f"сверка proximal vs M_fib: max отн. {(max(rels) if rels else np.nan):.4g}  n={len(rels)}"
        )
        recs = []
        for stem, title in SETS:
            genes = story[stem]
            for klass_i in ("proximal", "distal"):
                sub = peaks[(peaks.klass == klass_i) & peaks.gene.isin(genes)]
                for tissue, col in (("fib", "fib"), ("epi", "epi")):
                    sums, npeaks = ([], [])
                    for sp in SPECIES_ORDER:
                        sl = sub[sub.species == sp]
                        sums.append(float(sl[col].sum()))
                        npeaks.append(int(sl.peak_id.nunique()))
                    y = np.array([np.log10(s) if s > 0 else np.nan for s in sums], float)
                    ok = np.isfinite(y)
                    if ok.sum() >= 5:
                        fit = stats.linregress(logm[ok], y[ok])
                        rr, pp, sl, ice = (fit.rvalue, fit.pvalue, fit.slope, fit.intercept)
                    else:
                        rr = pp = sl = ice = np.nan
                    recs.append(
                        dict(
                            set=stem,
                            title=title,
                            klass=klass_i,
                            tissue=tissue,
                            n_genes=len(genes),
                            n_species=int(ok.sum()),
                            r=rr,
                            p=pp,
                            slope=sl,
                            intercept=ice,
                            **{f"sum_{sp.lower()}": sums[i] for i, sp in enumerate(SPECIES_ORDER)},
                            **{
                                f"npeaks_{sp.lower()}": npeaks[i]
                                for i, sp in enumerate(SPECIES_ORDER)
                            },
                        )
                    )
                    print(
                        f"{stem:16s} {klass_i:9s} {tissue}  r={rr:+.3f} p={pp:.3f}  peaks={npeaks} sums={[round(s, 1) for s in sums]}",
                        flush=True,
                    )
        cage = pd.DataFrame(recs)
        cage.to_csv(OUT / "cage_prox_distal_fits.tsv", sep="\t", index=False)
        distal_panels = []
        fig, axes = plt.subplots(5, 2, figsize=(8.4, 16.2))
        for i, (stem, title) in enumerate(SETS):
            for j, klass_i in enumerate(("proximal", "distal")):
                row = cage[
                    (cage.set == stem) & (cage.klass == klass_i) & (cage.tissue == "fib")
                ].iloc[0]
                sums = [float(row[f"sum_{sp.lower()}"]) for sp in SPECIES_ORDER]
                y = np.array([np.log10(s) if s > 0 else np.nan for s in sums])
                ylab = "log10_fib_proximal" if klass_i == "proximal" else "log10_fib_distal"
                fit = scatter(axes[i, j], names, logm, y, "log10_body_mass", ylab, s40)
                ptxt = (
                    "NA"
                    if not np.isfinite(fit["p"])
                    else f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
                )
                rtxt = "NA" if not np.isfinite(fit["r"]) else f"{fit['r']:.3f}"
                axes[i, j].set_title(
                    f"{title} {klass_i}   cor. {rtxt}; p = {ptxt}", fontsize=9, pad=6
                )
                if klass_i == "distal":
                    distal_panels.append((title, names, logm, y, fit))
        fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.08, right=0.99, top=0.98, bottom=0.03)
        save(fig, "fig_story_prox_distal_fib")
        plt.close(fig)
        fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
        axes = axes.ravel()
        for i, (title, nm, x, y, fit) in enumerate(distal_panels):
            scatter(axes[i], nm, x, y, "log10_body_mass", "log10_fib_distal", s40)
            ptxt = (
                "NA"
                if not np.isfinite(fit["p"])
                else f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            )
            rtxt = "NA" if not np.isfinite(fit["r"]) else f"{fit['r']:.3f}"
            axes[i].set_title(f"{title}   cor. {rtxt}; p = {ptxt}", fontsize=10, pad=7)
        axes[5].set_visible(False)
        fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.06, right=0.99, top=0.94, bottom=0.08)
        save(fig, "fig_story_five_distal_fib")
        plt.close(fig)
        if not MEME.exists() or not FIMO.exists():
            raise SystemExit(f"need FIMO={FIMO.exists()} MEME={MEME.exists()}")
        hit_rows = []
        fish_rows = []
        watch_rows = []
        pending = [
            sp
            for sp in SPECIES_ORDER
            if not (OUT / f"{sp}_fimo_proximal.tsv").exists()
            or (OUT / f"{sp}_fimo_proximal.tsv").stat().st_size < 100
        ]
        species_loop = pending + [sp for sp in SPECIES_ORDER if sp not in pending]
        for sp in species_loop:
            print(f"FIMO {sp}", flush=True)
            sp_peaks = prox[prox.species == sp].copy()
            fa = OUT / f"{sp}_proximal.fa"
            nseq = write_fasta(sp_peaks, s74.GENOME[sp], fa)
            print(f"  fasta {nseq} seq", flush=True)
            hits = run_fimo(fa, OUT / f"{sp}_fimo_proximal.tsv")
            if hits.empty:
                print("  нет хитов")
                continue
            name_col = "motif_alt_id" if "motif_alt_id" in hits.columns else "motif_id"
            seq_col = "sequence_name"
            hits["gene"] = hits[seq_col].astype(str).str.split("|").str[0]
            hits["tf"] = hits[name_col].astype(str)
            sp_peaks = sp_peaks.copy()
            sp_peaks["sid"] = (
                sp_peaks.gene.astype(str) + "|" + sp_peaks.peak_id.astype(str) + f"|{sp}|proximal"
            )
            all_sid = set(sp_peaks["sid"])
            by_tf: dict[str, set[str]] = {}
            for tf, seq in zip(hits["tf"], hits[seq_col].astype(str)):
                by_tf.setdefault(tf.upper(), set()).add(seq)
            for stem, title in SETS:
                genes = set(story[stem])
                pos = set(sp_peaks.loc[sp_peaks.gene.isin(genes), "sid"])
                neg = all_sid - pos
                if len(pos) < 2 or len(neg) < 20:
                    continue
                ps = []
                for tf, hit_seq in by_tf.items():
                    rec = fisher_block(hit_seq, pos, neg)
                    rec.update(
                        set=stem, title=title, species=sp, tf=tf, n_pos=len(pos), n_neg=len(neg)
                    )
                    fish_rows.append(rec)
                    ps.append(rec["p"])
                print(f"  {stem} pos={len(pos)} neg={len(neg)} motifs={len(by_tf)}", flush=True)
            for stem, title in SETS:
                genes = set(story[stem])
                pos = set(sp_peaks.loc[sp_peaks.gene.isin(genes), "sid"])
                neg = all_sid - pos
                for w in WATCH:
                    hit_seq = set()
                    for tf, seqs in by_tf.items():
                        parts = {p.strip() for p in tf.replace("::", " ").split()}
                        if tf == w or w in parts:
                            hit_seq |= seqs
                    rec = fisher_block(hit_seq, pos, neg)
                    rec.update(
                        set=stem, title=title, species=sp, tf=w, n_pos=len(pos), n_neg=len(neg)
                    )
                    watch_rows.append(rec)
        fish = pd.DataFrame(fish_rows)
        if not fish.empty:
            fish["bh"] = np.nan
            for (stem, sp), idx in fish.groupby(["set", "species"]).groups.items():
                p = fish.loc[idx, "p"].to_numpy(float)
                fish.loc[idx, "bh"] = stats.false_discovery_control(
                    np.nan_to_num(p, nan=1.0), method="bh"
                )
            fish.to_csv(OUT / "motif_fisher_proximal.tsv", sep="\t", index=False)
            sig = fish[fish.bh < 0.05]
            print("\nBH<0.05 proximal:")
            print(sig.groupby(["set", "species"]).size().to_string() if len(sig) else "  нет")
            keep = []
            for stem, g in sig.groupby("set"):
                vc = g.groupby("tf")["species"].nunique()
                for tf, k in vc.items():
                    if k >= 3:
                        keep.append(dict(set=stem, tf=tf, n_species=int(k)))
            pd.DataFrame(keep).to_csv(OUT / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
            print("мотивы BH в ≥3 видах:", keep)
        watch = pd.DataFrame(watch_rows)
        if not watch.empty:
            watch.to_csv(OUT / "motif_watch_fisher.tsv", sep="\t", index=False)
            wsum = []
            for (stem, tf), g in watch.groupby(["set", "tf"]):
                g = g.set_index("species").reindex(SPECIES_ORDER)
                y = -np.log10(np.clip(g["p"].to_numpy(float), 1e-12, 1))
                if np.isfinite(y).sum() >= 5:
                    fit = stats.linregress(logm[np.isfinite(y)], y[np.isfinite(y)])
                    wsum.append(
                        dict(
                            set=stem,
                            tf=tf,
                            r=float(fit.rvalue),
                            p=float(fit.pvalue),
                            **{f"p_{sp.lower()}": float(g.loc[sp, "p"]) for sp in SPECIES_ORDER},
                            **{
                                f"odds_{sp.lower()}": float(g.loc[sp, "odds"])
                                for sp in SPECIES_ORDER
                            },
                        )
                    )
            pd.DataFrame(wsum).to_csv(
                OUT / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
            )
            print("\nwatch TF, −log10 p vs mass:")
            ws = pd.DataFrame(wsum).sort_values("p")
            print(ws[["set", "tf", "r", "p"]].head(20).to_string(index=False))
        print(f"wrote {OUT}")

    return SimpleNamespace(**locals())


def _motif_method_compare(context):
    """Three lightweight motif tests next to the frozen JASPAR FIMO.

    Does not rewrite story_cre_motifs/*_fimo_proximal.tsv or the frozen tables.

      hocomoco  FIMO 5.3.0, threshold 1e-4, two-sided Fisher, BH within
                set x species. Catalog is HOCOMOCO v12 H12CORE, one model
                per factor (index 0). Names are HGNC symbols (ERR1 -> ESRRA).
      sea       SEA from MEME 5.5 against the other proximal promoters.
                Score threshold is chosen on a hold-out of at least 6
                sequences when the set is large enough for that.
      monalisa  monaLisa::calcBinnedMotifEnrR, bins set vs rest, weights
                by GC and 3-mers. Calls the R script.
      compare   Gene-symbol lists. No averaged score.

    Source: analyses/avian/scripts/82_motif_method_compare.py at d230e56ab942."""
    import importlib.util
    import math
    import re
    import subprocess
    import sys
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    ROOT = Path(str(context.resolve_path("/path/to/user-home/chicken/kleiber_birds_gemcat")))
    SCRIPTS = ROOT / "scripts"
    STORY = ROOT / "results" / "oxphos_species_thread" / "story_cre_motifs"
    CMP = STORY / "method_compare"
    FROZEN = CMP / "fimo_jaspar_frozen"
    HOCO = CMP / "hocomoco"
    SEA_DIR = CMP / "sea"
    MONA = CMP / "monalisa"
    REPORT = CMP / "compare"
    FIMO = Path(str(context.resolve_path("/path/to/user-home/miniconda3/envs/meme/bin/fimo")))
    SEA_BIN = CMP / "meme55" / "bin" / "sea"
    R_SCRIPT = context.root("repository") / "R/regulation.R"
    FIMO_THRESH = "1e-4"
    MOTIF_RE = re.compile("^(.+)\\.H12CORE\\.(\\d+)\\.([^.]+)\\.([A-D])$")
    JASPAR_MOTIF_RE = re.compile("^MOTIF\\s+(\\S+)(?:\\s+(\\S+))?")

    def load_s81():
        mod = context.unit("avian.81_story_cre_motifs")
        return mod

    def split_meme(text: str) -> tuple[list[str], list[tuple[str, list[str]]]]:
        header: list[str] = []
        blocks: list[tuple[str, list[str]]] = []
        cur_name = None
        cur: list[str] = []
        for line in text.splitlines(keepends=True):
            if line.startswith("MOTIF "):
                if cur_name is not None:
                    blocks.append((cur_name, cur))
                cur_name = line.split()[1]
                cur = [line]
            elif cur_name is None:
                header.append(line)
            else:
                cur.append(line)
        if cur_name is not None:
            blocks.append((cur_name, cur))
        return (header, blocks)

    def write_blocks(path: Path, header: list[str], blocks: list[list[str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as fh:
            fh.writelines(header)
            for block in blocks:
                fh.writelines(block)
                if not block[-1].endswith("\n"):
                    fh.write("\n")

    def hocomoco_symbol_table() -> pd.DataFrame:
        master = pd.read_csv(HOCO / "hocomoco_v12_masterlist.tsv", sep="\t")
        mapping = pd.read_csv(HOCO / "hocomoco_benchmarking_ID_mapping.tsv", sep="\t")
        human = master[master["curated:uniprot_id"].astype(str).str.endswith("_HUMAN")]
        u2g = dict(zip(human["curated:uniprot_id"], human["auto:gene_symbol"]))
        mapping = mapping.copy()
        mapping["prefix"] = mapping["HOCOMOCO_ID"].str.split(".").str[0]
        rows = []
        seen = set()
        for rec in mapping.itertuples(index=False):
            if rec.prefix in seen:
                continue
            if ".H12CORE.0." not in rec.HOCOMOCO_ID:
                continue
            seen.add(rec.prefix)
            gene = u2g.get(rec.Uniprot_HUMAN)
            if gene is None or (isinstance(gene, float) and math.isnan(gene)):
                raise SystemExit(f"no HGNC for {rec.prefix} / {rec.Uniprot_HUMAN}")
            rows.append(dict(prefix=rec.prefix, uniprot=rec.Uniprot_HUMAN, tf=str(gene).upper()))
        tab = pd.DataFrame(rows)
        if tab["tf"].duplicated().any():
            bad = tab[tab["tf"].duplicated(keep=False)]
            raise SystemExit(f"two HOCOMOCO factors map to one symbol\n{bad}")
        return tab

    def build_hocomoco_meme() -> Path:
        src = HOCO / "H12CORE_meme_format.meme"
        dst = HOCO / "H12CORE_one_per_tf.meme"
        symbols = hocomoco_symbol_table().set_index("prefix")
        header, blocks = split_meme(src.read_text())
        kept = []
        meta = []
        for name, block in blocks:
            m = MOTIF_RE.match(name)
            if m is None or int(m.group(2)) != 0:
                continue
            prefix = m.group(1)
            if prefix not in symbols.index:
                raise SystemExit(f"motif {name} has no HGNC map")
            kept.append(block)
            meta.append(
                dict(
                    motif_id=name,
                    prefix=prefix,
                    origin=m.group(3),
                    grade=m.group(4),
                    tf=symbols.loc[prefix, "tf"],
                    uniprot=symbols.loc[prefix, "uniprot"],
                )
            )
        write_blocks(dst, header, kept)
        pd.DataFrame(meta).to_csv(HOCO / "models_kept.tsv", sep="\t", index=False)
        print(f"HOCOMOCO одна модель на фактор: {len(kept)}  -> {dst}", flush=True)
        return dst

    def jaspar_nsites(block: list[str]) -> int:
        for line in block:
            m = re.search("nsites=\\s*([0-9]+)", line)
            if m:
                return int(m.group(1))
        return 0

    def jaspar_version(motif_id: str) -> int:
        m = re.search("\\.(\\d+)$", motif_id)
        return int(m.group(1)) if m else 0

    def build_jaspar_one_per_symbol(meme: Path) -> Path:
        """One JASPAR matrix per uppercase alt id.

        The frozen FIMO ORs every version that shares an alt id. Versions
        that share a name are 53 alt ids; the OXPHOS six and the watch list
        each have a single matrix. Here the kept matrix is the one with the
        most sites, then the higher JASPAR version.
        """
        dst = CMP / "JASPAR2024_one_per_symbol.meme"
        header, blocks = split_meme(meme.read_text())
        best: dict[str, tuple[int, int, str, list[str]]] = {}
        for name, block in blocks:
            m = JASPAR_MOTIF_RE.match(block[0].rstrip("\n"))
            if m is None:
                continue
            motif_id, alt = (m.group(1), m.group(2) or m.group(1))
            symbol = alt.upper()
            key = (jaspar_nsites(block), jaspar_version(motif_id))
            prev = best.get(symbol)
            if prev is None or key > (prev[0], prev[1]):
                best[symbol] = (key[0], key[1], motif_id, block)
        chosen = list(best.values())
        write_blocks(dst, header, [c[3] for c in chosen])
        pd.DataFrame(
            [dict(tf=sym, motif_id=rec[2], nsites=rec[0]) for sym, rec in best.items()]
        ).to_csv(CMP / "jaspar_one_per_symbol.tsv", sep="\t", index=False)
        print(f"JASPAR одна матрица на имя: {len(chosen)} из {len(blocks)}", flush=True)
        return dst

    def fasta_records(path: Path) -> list[tuple[str, str, str]]:
        """Return (header_without_gt, gene, raw_record_including_header)."""
        out = []
        header = None
        buf: list[str] = []
        with path.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        out.append((header, header.split("|", 1)[0], "".join(buf)))
                    header = line[1:].strip()
                    buf = [line]
                else:
                    buf.append(line)
        if header is not None:
            out.append((header, header.split("|", 1)[0], "".join(buf)))
        return out

    def bh_within(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["bh"] = np.nan
        for (_, _), idx in df.groupby(["set", "species"]).groups.items():
            p = df.loc[idx, "p"].to_numpy(float)
            ok = np.isfinite(p)
            if ok.sum() == 0:
                continue
            adj = np.full(len(p), np.nan)
            adj[ok] = stats.false_discovery_control(p[ok], method="bh")
            df.loc[idx, "bh"] = adj
        return df

    def ge3_table(df: pd.DataFrame) -> pd.DataFrame:
        sig = df[df["bh"] < 0.05]
        rows = []
        if sig.empty:
            return pd.DataFrame(columns=["set", "tf", "n_species"])
        for (stem, tf), g in sig.groupby(["set", "tf"]):
            k = int(g["species"].nunique())
            if k >= 3:
                rows.append(dict(set=stem, tf=tf, n_species=k))
        out = pd.DataFrame(rows)
        if out.empty:
            return pd.DataFrame(columns=["set", "tf", "n_species"])
        return out.sort_values(["set", "n_species", "tf"], ascending=[True, False, True])

    def watch_vs_mass(
        watch: pd.DataFrame, species_order: list[str], logm: np.ndarray
    ) -> pd.DataFrame:
        rows = []
        for (stem, tf), g in watch.groupby(["set", "tf"]):
            g = g.set_index("species").reindex(species_order)
            p = g["p"].to_numpy(float)
            y = -np.log10(np.clip(p, 1e-300, 1))
            ok = np.isfinite(y)
            rec = dict(set=stem, tf=tf, n_points=int(ok.sum()), r=np.nan, p_mass=np.nan)
            if ok.sum() >= 5:
                fit = stats.linregress(logm[ok], y[ok])
                rec["r"] = float(fit.rvalue)
                rec["p_mass"] = float(fit.pvalue)
            bh = g["bh"] if "bh" in g.columns else pd.Series(np.nan, index=g.index)
            rec["n_species_bh"] = int((bh < 0.05).sum())
            for sp in species_order:
                rec[f"p_{sp.lower()}"] = (
                    float(g.loc[sp, "p"])
                    if sp in g.index and np.isfinite(g.loc[sp, "p"])
                    else np.nan
                )
            rows.append(rec)
        return pd.DataFrame(rows)

    def attach_watch_bh(watch: pd.DataFrame, fisher: pd.DataFrame) -> pd.DataFrame:
        key = fisher.set_index(["set", "species", "tf"])["bh"]
        watch = watch.copy()
        watch["bh"] = [
            key.get((r.set, r.species, r.tf), np.nan) for r in watch.itertuples(index=False)
        ]
        return watch

    def run_hocomoco(s81) -> None:
        meme = build_hocomoco_meme()
        symbol = (
            pd.read_csv(HOCO / "models_kept.tsv", sep="\t").set_index("motif_id")["tf"].to_dict()
        )
        HOCO.mkdir(parents=True, exist_ok=True)
        jobs = []
        for sp in s81.SPECIES_ORDER:
            fasta = STORY / f"{sp}_proximal.fa"
            if not fasta.exists():
                raise SystemExit(f"нет fasta {fasta}")
            out = HOCO / f"{sp}_fimo.tsv"
            log = HOCO / f"{sp}_fimo.log"
            jobs.append((sp, fasta, meme, out, log))

        def one(job):
            sp, fasta, meme_path, out, log = job
            if (
                out.exists()
                and out.stat().st_size > 1000
                and log.exists()
                and ("fimo failed" not in log.read_text(errors="replace"))
            ):
                if (HOCO / f"{sp}_fimo.ok").exists():
                    return (sp, "cached", out.stat().st_size)
            cmd = [
                str(FIMO),
                "--text",
                "--no-qvalue",
                "--verbosity",
                "1",
                "--thresh",
                FIMO_THRESH,
                str(meme_path),
                str(fasta),
            ]
            with out.open("w") as fh, log.open("w") as err:
                proc = subprocess.run(cmd, stdout=fh, stderr=err)
            if proc.returncode != 0:
                log.write_text(log.read_text() + f"\nfimo failed {proc.returncode}\n")
                raise SystemExit(f"FIMO {sp} rc={proc.returncode}")
            (HOCO / f"{sp}_fimo.ok").write_text("ok\n")
            return (sp, "ran", out.stat().st_size)

        print("FIMO HOCOMOCO, семь видов", flush=True)
        with ThreadPoolExecutor(max_workers=7) as pool:
            futs = [pool.submit(one, job) for job in jobs]
            for fut in as_completed(futs):
                sp, how, size = fut.result()
                print(f"  {sp} {how} {size}", flush=True)
        story = {stem: set(s81.genes_of(stem)) for stem, _ in s81.SETS}
        fish_rows = []
        watch_rows = []
        for sp in s81.SPECIES_ORDER:
            recs = fasta_records(STORY / f"{sp}_proximal.fa")
            by_gene: dict[str, set[str]] = {}
            for sid, gene, _ in recs:
                by_gene.setdefault(gene, set()).add(sid)
            all_sid = {sid for sid, _, _ in recs}
            hits = pd.read_csv(
                HOCO / f"{sp}_fimo.tsv",
                sep="\t",
                comment="#",
                usecols=lambda c: c in {"motif_id", "motif_alt_id", "sequence_name"},
            )
            if hits.empty:
                print(f"  {sp} нет хитов", flush=True)
                continue
            id_col = "motif_id"
            hits["tf"] = hits[id_col].map(symbol)
            missing = hits["tf"].isna().sum()
            if missing:
                raise SystemExit(f"{sp}: {missing} хитов без символа")
            by_tf: dict[str, set[str]] = {}
            for tf, seq in zip(hits["tf"], hits["sequence_name"].astype(str)):
                by_tf.setdefault(str(tf), set()).add(seq)
            print(f"  {sp} seq {len(all_sid)} motifs {len(by_tf)}", flush=True)
            for stem, title in s81.SETS:
                pos = (
                    set().union(*(by_gene.get(g, set()) for g in story[stem]))
                    if story[stem]
                    else set()
                )
                neg = all_sid - pos
                if len(pos) < 2 or len(neg) < 20:
                    print(f"    {stem} пропуск pos={len(pos)} neg={len(neg)}", flush=True)
                    continue
                for tf, hit_seq in by_tf.items():
                    rec = s81.fisher_block(hit_seq, pos, neg)
                    rec.update(
                        set=stem, title=title, species=sp, tf=tf, n_pos=len(pos), n_neg=len(neg)
                    )
                    fish_rows.append(rec)
                for w in s81.WATCH:
                    rec = s81.fisher_block(by_tf.get(w, set()), pos, neg)
                    rec.update(
                        set=stem, title=title, species=sp, tf=w, n_pos=len(pos), n_neg=len(neg)
                    )
                    watch_rows.append(rec)
        fish = bh_within(pd.DataFrame(fish_rows))
        fish.to_csv(HOCO / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = ge3_table(fish)
        ge3.to_csv(HOCO / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        watch = attach_watch_bh(pd.DataFrame(watch_rows), fish)
        watch.to_csv(HOCO / "motif_watch_fisher.tsv", sep="\t", index=False)
        logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
        watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
            HOCO / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
        )
        print(
            "HOCOMOCO BH<0.05",
            int((fish.bh < 0.05).sum()),
            "ge3",
            ge3.to_dict("records"),
            flush=True,
        )

    def sea_hofract(n_pos: int) -> tuple[float, str]:
        """SEA drops the hold-out when it would contain fewer than 5 sequences.

        Default 0.1 is kept when it already yields at least 6. Otherwise the
        fraction is raised so the hold-out has 6, and only if at least 5
        sequences remain to be scored. Smaller sets cannot host that hold-out.
        """
        if n_pos < 11:
            return (0.1, "too_few_for_holdout")
        if 0.1 * n_pos >= 6:
            return (0.1, "default_0.1")
        return (6 / n_pos, "raised_holdout_6")

    def split_set_fasta(recs, genes: set[str], pos_path: Path, neg_path: Path) -> tuple[int, int]:
        pos_path.parent.mkdir(parents=True, exist_ok=True)
        n_pos = n_neg = 0
        with pos_path.open("w") as fp, neg_path.open("w") as fn:
            for _sid, gene, raw in recs:
                if gene in genes:
                    fp.write(raw)
                    n_pos += 1
                else:
                    fn.write(raw)
                    n_neg += 1
        return (n_pos, n_neg)

    def run_sea(s81) -> None:
        meme = build_jaspar_one_per_symbol(s81.MEME)
        SEA_DIR.mkdir(parents=True, exist_ok=True)
        story = {stem: set(s81.genes_of(stem)) for stem, _ in s81.SETS}
        jobs = []
        for sp in s81.SPECIES_ORDER:
            recs = fasta_records(STORY / f"{sp}_proximal.fa")
            for stem, title in s81.SETS:
                pos = SEA_DIR / "fasta" / f"{sp}_{stem}_pos.fa"
                neg = SEA_DIR / "fasta" / f"{sp}_{stem}_neg.fa"
                n_pos, n_neg = split_set_fasta(recs, story[stem], pos, neg)
                h, note = sea_hofract(n_pos)
                outdir = SEA_DIR / "runs" / sp / stem
                jobs.append(
                    dict(
                        sp=sp,
                        stem=stem,
                        title=title,
                        pos=pos,
                        neg=neg,
                        n_pos=n_pos,
                        n_neg=n_neg,
                        hofract=h,
                        note=note,
                        outdir=outdir,
                    )
                )

        def one(job):
            tsv = job["outdir"] / "sea.tsv"
            err_path = job["outdir"] / "sea.err"
            if tsv.exists() and tsv.stat().st_size > 50 and (job["outdir"] / "sea.ok").exists():
                return (
                    job["sp"],
                    job["stem"],
                    "cached",
                    err_path.read_text(errors="replace") if err_path.exists() else "",
                )
            if job["n_pos"] < 2 or job["n_neg"] < 20:
                return (job["sp"], job["stem"], "skip", "")
            job["outdir"].mkdir(parents=True, exist_ok=True)
            cmd = [
                str(SEA_BIN),
                "--p",
                str(job["pos"]),
                "--n",
                str(job["neg"]),
                "--m",
                str(meme),
                "--text",
                "--noseqs",
                "--notrim",
                "--hofract",
                f"{job['hofract']:.6f}",
                "--seed",
                "1",
                "--pvalue",
                "--thresh",
                "1",
                "--verbosity",
                "1",
            ]
            with tsv.open("w") as out, err_path.open("w") as err:
                proc = subprocess.run(cmd, stdout=out, stderr=err)
            if proc.returncode != 0 or tsv.stat().st_size < 200:
                raise SystemExit(f"SEA {job['sp']} {job['stem']} rc={proc.returncode}")
            (job["outdir"] / "sea.ok").write_text(job["note"] + "\n")
            return (job["sp"], job["stem"], "ran", err_path.read_text(errors="replace"))

        print(f"SEA, заданий {len(jobs)}", flush=True)
        notes = []
        with ThreadPoolExecutor(max_workers=12) as pool:
            futs = [pool.submit(one, job) for job in jobs]
            done = {(j["sp"], j["stem"]): j for j in jobs}
            for fut in as_completed(futs):
                sp, stem, how, err = fut.result()
                job = done[sp, stem]
                holdout_used = "No hold-out set was created" not in err
                notes.append(
                    dict(
                        species=sp,
                        set=stem,
                        n_pos=job["n_pos"],
                        n_neg=job["n_neg"],
                        hofract=job["hofract"],
                        plan=job["note"],
                        status=how,
                        holdout_used=holdout_used and how != "skip",
                        length_warning="different length distributions" in err,
                    )
                )
                print(f"  {sp} {stem} {how} holdout={holdout_used}", flush=True)
        pd.DataFrame(notes).to_csv(SEA_DIR / "run_notes.tsv", sep="\t", index=False)
        fish_rows = []
        for job in jobs:
            tsv = job["outdir"] / "sea.tsv"
            if not tsv.exists():
                continue
            tab = pd.read_csv(tsv, sep="\t")
            if "ALT_ID" not in tab.columns or "PVALUE" not in tab.columns:
                raise SystemExit(f"непонятный sea.tsv {tsv}: {tab.columns.tolist()}")
            tab = tab.dropna(subset=["ALT_ID", "ID", "PVALUE"])
            tab["tf"] = tab["ALT_ID"].astype(str).str.upper()
            if tab["tf"].duplicated().any():
                raise SystemExit(f"SEA {job['sp']} {job['stem']}: повтор символа")
            for rec in tab.itertuples(index=False):
                fish_rows.append(
                    dict(
                        set=job["stem"],
                        title=job["title"],
                        species=job["sp"],
                        tf=rec.tf,
                        motif_id=getattr(rec, "ID"),
                        p=float(rec.PVALUE),
                        enr_ratio=float(rec.ENR_RATIO) if np.isfinite(rec.ENR_RATIO) else np.nan,
                        score_thr=(
                            float(rec.SCORE_THR)
                            if np.isfinite(getattr(rec, "SCORE_THR", np.nan))
                            else np.nan
                        ),
                        tp=float(rec.TP),
                        fp=float(rec.FP),
                        n_pos=job["n_pos"],
                        n_neg=job["n_neg"],
                        sea_q=float(rec.QVALUE) if np.isfinite(rec.QVALUE) else np.nan,
                    )
                )
        fish = bh_within(pd.DataFrame(fish_rows))
        fish.to_csv(SEA_DIR / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = ge3_table(fish)
        ge3.to_csv(SEA_DIR / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        watch_rows = []
        have = fish.set_index(["set", "species", "tf"])
        for job in jobs:
            if not (job["outdir"] / "sea.tsv").exists():
                continue
            for w in s81.WATCH:
                key = (job["stem"], job["sp"], w)
                if key in have.index:
                    row = have.loc[key]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    p, bh = (float(row["p"]), float(row["bh"]))
                else:
                    p, bh = (1.0, np.nan)
                watch_rows.append(dict(set=job["stem"], species=job["sp"], tf=w, p=p, bh=bh))
        watch = pd.DataFrame(watch_rows)
        watch.to_csv(SEA_DIR / "motif_watch_fisher.tsv", sep="\t", index=False)
        logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
        watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
            SEA_DIR / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
        )
        print(
            "SEA BH<0.05",
            int((fish.bh < 0.05).sum()) if len(fish) else 0,
            "ge3",
            ge3.to_dict("records"),
            flush=True,
        )

    def run_monalisa(s81) -> None:
        meme = build_jaspar_one_per_symbol(s81.MEME)
        MONA.mkdir(parents=True, exist_ok=True)
        geneset_dir = s81.DATA / "genesets"

        def one(sp):
            out = MONA / f"{sp}_binned.tsv"
            log = MONA / f"{sp}.log"
            if (MONA / f"{sp}.ok").exists() and out.exists() and (out.stat().st_size > 50):
                return (sp, "cached")
            cmd = [
                "Rscript",
                "--vanilla",
                str(R_SCRIPT),
                "--meme",
                str(meme),
                "--fasta",
                str(STORY / f"{sp}_proximal.fa"),
                "--species",
                sp,
                "--geneset-dir",
                str(geneset_dir),
                "--out",
                str(out),
                "--cores",
                "12",
            ]
            with log.open("w") as fh:
                proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                raise SystemExit(f"monaLisa {sp} rc={proc.returncode} см. {log}")
            (MONA / f"{sp}.ok").write_text("ok\n")
            return (sp, "ran")

        print("monaLisa, семь видов", flush=True)
        with ThreadPoolExecutor(max_workers=7) as pool:
            futs = [pool.submit(one, sp) for sp in s81.SPECIES_ORDER]
            for fut in as_completed(futs):
                sp, how = fut.result()
                print(f"  {sp} {how}", flush=True)
        parts = [pd.read_csv(MONA / f"{sp}_binned.tsv", sep="\t") for sp in s81.SPECIES_ORDER]
        fish = bh_within(pd.concat(parts, ignore_index=True))
        fish.to_csv(MONA / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = ge3_table(fish)
        ge3.to_csv(MONA / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        have = fish.set_index(["set", "species", "tf"])
        watch_rows = []
        for sp in s81.SPECIES_ORDER:
            for stem, _title in s81.SETS:
                for w in s81.WATCH:
                    key = (stem, sp, w)
                    if key in have.index:
                        row = have.loc[key]
                        if isinstance(row, pd.DataFrame):
                            row = row.iloc[0]
                        p, bh = (float(row["p"]), float(row["bh"]))
                    else:
                        p, bh = (1.0, np.nan)
                    watch_rows.append(dict(set=stem, species=sp, tf=w, p=p, bh=bh))
        watch = pd.DataFrame(watch_rows)
        watch.to_csv(MONA / "motif_watch_fisher.tsv", sep="\t", index=False)
        logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
        watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
            MONA / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
        )
        print(
            "monaLisa BH<0.05",
            int((fish.bh < 0.05).sum()),
            "ge3",
            ge3.to_dict("records"),
            flush=True,
        )

    def _ge3_sets(path: Path) -> dict[str, set[str]]:
        if not path.exists():
            return {}
        df = pd.read_csv(path, sep="\t")
        if df.empty or "tf" not in df.columns:
            return {}
        out: dict[str, set[str]] = {}
        for stem, g in df.groupby("set"):
            out[str(stem)] = set(g["tf"].astype(str))
        return out

    def uniform_watch(method: str, path: Path) -> pd.DataFrame:
        w = pd.read_csv(path, sep="\t")
        out = w[["set", "tf", "n_points", "n_species_bh", "r", "p_mass"]].copy()
        out.insert(0, "method", method)
        return out

    def fimo_watch_uniform() -> pd.DataFrame:
        """Frozen watch file stores the mass-correlation p in column `p` and has no BH count."""
        w = pd.read_csv(FROZEN / "watch_tf_enrichment_vs_mass.tsv", sep="\t")
        fish = pd.read_csv(FROZEN / "motif_fisher_proximal.tsv", sep="\t")
        counts = fish.loc[fish["bh"] < 0.05].groupby(["set", "tf"])["species"].nunique()
        rows = []
        for rec in w.itertuples(index=False):
            rows.append(
                dict(
                    method="fimo_jaspar",
                    set=rec.set,
                    tf=rec.tf,
                    n_points=7,
                    n_species_bh=int(counts.get((rec.set, rec.tf), 0)),
                    r=float(rec.r),
                    p_mass=float(rec.p),
                )
            )
        return pd.DataFrame(rows)

    def jaccard(a: set[str], b: set[str]) -> float:
        if not a and (not b):
            return float("nan")
        return len(a & b) / len(a | b)

    def run_compare(s81) -> None:
        methods = {
            "fimo_jaspar": FROZEN / "motifs_bh_in_ge3_species.tsv",
            "hocomoco": HOCO / "motifs_bh_in_ge3_species.tsv",
            "monalisa": MONA / "motifs_bh_in_ge3_species.tsv",
            "sea": SEA_DIR / "motifs_bh_in_ge3_species.tsv",
        }
        for name, path in methods.items():
            if not path.exists():
                raise SystemExit(f"нет таблицы {name}: {path}")
        live = (STORY / "motifs_bh_in_ge3_species.tsv").read_bytes()
        frozen = (FROZEN / "motifs_bh_in_ge3_species.tsv").read_bytes()
        if live != frozen:
            raise SystemExit("живая таблица FIMO разъехалась с заморозкой")
        lists = {name: _ge3_sets(path) for name, path in methods.items()}
        stems = [stem for stem, _ in s81.SETS]
        jac_rows = []
        ret_rows = []
        long_rows = []
        for stem in stems:
            base = lists["fimo_jaspar"].get(stem, set())
            for name, grouped in lists.items():
                got = grouped.get(stem, set())
                for tf in sorted(got):
                    src = pd.read_csv(methods[name], sep="\t")
                    n = int(src[(src["set"] == stem) & (src["tf"] == tf)]["n_species"].iloc[0])
                    long_rows.append(dict(set=stem, method=name, tf=tf, n_species=n))
                if name == "fimo_jaspar":
                    continue
                jac_rows.append(
                    dict(
                        set=stem,
                        method=name,
                        n_fimo=len(base),
                        n_method=len(got),
                        n_both=len(base & got),
                        jaccard=jaccard(base, got),
                        only_fimo=";".join(sorted(base - got)),
                        only_method=";".join(sorted(got - base)),
                    )
                )
            for tf in sorted(base):
                rec = dict(set=stem, tf=tf)
                for name in ("hocomoco", "monalisa", "sea"):
                    src = pd.read_csv(methods[name], sep="\t")
                    hit = src[(src["set"] == stem) & (src["tf"] == tf)]
                    rec[f"n_{name}"] = int(hit["n_species"].iloc[0]) if len(hit) else 0
                    rec[f"kept_{name}"] = bool(len(hit))
                ret_rows.append(rec)
        REPORT.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(long_rows).to_csv(REPORT / "ge3_by_method.tsv", sep="\t", index=False)
        pd.DataFrame(jac_rows).to_csv(REPORT / "jaccard_vs_fimo.tsv", sep="\t", index=False)
        pd.DataFrame(ret_rows).to_csv(REPORT / "fimo_ge3_retention.tsv", sep="\t", index=False)
        watch = pd.concat(
            [
                fimo_watch_uniform(),
                uniform_watch("hocomoco", HOCO / "watch_tf_enrichment_vs_mass.tsv"),
                uniform_watch("monalisa", MONA / "watch_tf_enrichment_vs_mass.tsv"),
                uniform_watch("sea", SEA_DIR / "watch_tf_enrichment_vs_mass.tsv"),
            ],
            ignore_index=True,
        )
        watch.to_csv(REPORT / "watch_vs_mass.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                dict(
                    method="fimo_jaspar",
                    catalog="JASPAR2024 CORE vertebrates",
                    scan="FIMO 5.3.0, thresh 1e-4",
                    test="Fisher two-sided",
                    background="other proximal peaks",
                    composition="none",
                    matrices="879, versions of one name OR-ed",
                ),
                dict(
                    method="hocomoco",
                    catalog="HOCOMOCO v12 H12CORE, index 0, HGNC",
                    scan="FIMO 5.3.0, thresh 1e-4",
                    test="Fisher two-sided",
                    background="other proximal peaks",
                    composition="none",
                    matrices="949, one model per symbol",
                ),
                dict(
                    method="monalisa",
                    catalog="JASPAR2024, one matrix per symbol",
                    scan="monaLisa 1.18.0, min.score 10",
                    test="Fisher one-sided enrichment",
                    background="other bin",
                    composition="GC bins and k-mers to 3",
                    matrices="816",
                ),
                dict(
                    method="sea",
                    catalog="JASPAR2024, one matrix per symbol",
                    scan="SEA 5.5.9, --pvalue --thresh 1, --notrim",
                    test="SEA score, BH recomputed inside set x species",
                    background="other proximal peaks",
                    composition="none; hold-out raises to 6 seq when n>=11, else refused",
                    matrices="816",
                ),
            ]
        ).to_csv(REPORT / "method_notes.tsv", sep="\t", index=False)
        print(pd.DataFrame(jac_rows).to_string(index=False), flush=True)
        print(pd.DataFrame(ret_rows).to_string(index=False), flush=True)

    def main() -> None:
        stage = (["bodymass"] + context.argv)[1] if len(["bodymass"] + context.argv) > 1 else "all"
        s81 = load_s81()
        if stage in ("hocomoco", "all"):
            run_hocomoco(s81)
        if stage in ("sea", "all"):
            run_sea(s81)
        if stage in ("monalisa", "all"):
            run_monalisa(s81)
        if stage in ("compare", "all"):
            run_compare(s81)

    return SimpleNamespace(**locals())


def _nine_set_motifs(context):
    """Motif enrichment for the nine locked bird sets.

    Does not rewrite story_cre_motifs/*_fimo_proximal.tsv, the HOCOMOCO
    hit tables, or method_compare/. Those scans already cover every
    proximal peak. This script only changes which genes are the foreground.

    JASPAR FIMO and HOCOMOCO FIMO: two-sided Fisher, BH within set x species.
    SEA and monaLisa: same settings as script 82, new output directory.
    A motif is reported when BH < 0.05 in at least three species.

    Source: analyses/avian/scripts/94_nine_set_motifs.py at d230e56ab942."""
    import importlib.util
    import subprocess
    import sys
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    ROOT = Path(str(context.resolve_path("/path/to/user-home/chicken/kleiber_birds_gemcat")))
    SCRIPTS = ROOT / "scripts"
    STORY = ROOT / "results" / "oxphos_species_thread" / "story_cre_motifs"
    OUT = STORY / "nine_sets"
    KB = ROOT / "corpus_kleiber_2026-08-06" / "02_база_знаний"
    HOCO = STORY / "method_compare" / "hocomoco"
    JASPAR_MEME = STORY / "method_compare" / "JASPAR2024_one_per_symbol.meme"
    SEA_BIN = STORY / "method_compare" / "meme55" / "bin" / "sea"
    R_SCRIPT = context.root("repository") / "R/regulation.R"
    SETS = [
        ("oxphos_16", "OXPHOS 16"),
        ("folate_12", "фолат 12"),
        ("ribo_5", "ribo 5"),
        ("green_hubs_7", "green hubs 7"),
        ("green_51_rest", "green, остальные 41"),
        ("carbon_7", "углерод 7"),
        ("amino_5", "аминокислоты 5"),
        ("nafld_4", "NAFLD 4"),
        ("epi_lipid_3", "липиды 3"),
    ]

    def load_s82():
        mod = context.unit("avian.82_motif_method_compare")
        return mod

    def load_sets() -> dict[str, set[str]]:
        dis = pd.read_csv(KB / "working_sets_disjoint_2026-09-22.tsv", sep="\t")
        ov = pd.read_csv(KB / "working_sets_overlap_ok_2026-09-22.tsv", sep="\t")
        genes = {stem: set(dis.loc[dis.set_id == stem, "gene"].astype(str)) for stem, _ in SETS}
        for stem, _ in SETS:
            extra = set(ov.loc[ov.set_id == stem, "gene"].astype(str))
            if extra:
                genes[stem] = extra
        gene_dir = OUT / "genes"
        gene_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for stem, title in SETS:
            if not genes[stem]:
                raise SystemExit(f"пустой набор {stem}")
            path = gene_dir / f"{stem}.txt"
            path.write_text("".join((g + "\n" for g in sorted(genes[stem]))))
            rows.append(dict(set=stem, title=title, n=len(genes[stem])))
        pd.DataFrame(rows).to_csv(OUT / "sets.tsv", sep="\t", index=False)
        return genes

    def species_order(s82) -> list[str]:
        s81 = s82.load_s81()
        return list(s81.SPECIES_ORDER)

    def sequences_by_gene(s82, sp: str) -> tuple[dict[str, set[str]], set[str]]:
        recs = s82.fasta_records(STORY / f"{sp}_proximal.fa")
        by_gene: dict[str, set[str]] = {}
        for sid, gene, _raw in recs:
            by_gene.setdefault(gene, set()).add(sid)
        return (by_gene, {sid for sid, _, _ in recs})

    def fisher_catalog(s82, genes: dict[str, set[str]], catalog: str, hit_path, tf_of) -> None:
        """tf_of(dataframe) -> Series of TF symbols aligned to hits."""
        fisher_block = s82.load_s81().fisher_block
        sps = species_order(s82)
        rows = []
        for sp in sps:
            by_gene, all_sid = sequences_by_gene(s82, sp)
            hits = pd.read_csv(hit_path(sp), sep="\t", comment="#")
            hits = hits[hits["sequence_name"].notna()]
            tf = tf_of(hits)
            by_tf: dict[str, set[str]] = {}
            for name, seq in zip(tf, hits["sequence_name"].astype(str)):
                if not isinstance(name, str) or not name:
                    continue
                by_tf.setdefault(name, set()).add(seq)
            print(f"  {catalog} {sp} seq {len(all_sid)} motifs {len(by_tf)}", flush=True)
            for stem, title in SETS:
                pos = set().union(*(by_gene.get(g, set()) for g in genes[stem]))
                neg = all_sid - pos
                if len(pos) < 2 or len(neg) < 20:
                    print(f"    {stem} пропуск pos={len(pos)}", flush=True)
                    continue
                for name, hit_seq in by_tf.items():
                    rec = fisher_block(hit_seq, pos, neg)
                    rec.update(
                        set=stem, title=title, species=sp, tf=name, n_pos=len(pos), n_neg=len(neg)
                    )
                    rows.append(rec)
        fish = s82.bh_within(pd.DataFrame(rows))
        dest = OUT / catalog
        dest.mkdir(parents=True, exist_ok=True)
        fish.to_csv(dest / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = s82.ge3_table(fish)
        ge3.to_csv(dest / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        n_bh = int((fish.bh < 0.05).sum())
        print(f"{catalog} BH<0.05 {n_bh}  ge3 {len(ge3)}", flush=True)
        if len(ge3):
            print(ge3.to_string(index=False), flush=True)

    def jaspar_tf(hits: pd.DataFrame) -> pd.Series:
        return hits["motif_alt_id"].astype(str).str.upper()

    def hoco_tf(hits: pd.DataFrame) -> pd.Series:
        symbol = pd.read_csv(HOCO / "models_kept.tsv", sep="\t").set_index("motif_id")["tf"]
        return hits["motif_id"].map(symbol)

    def run_fisher(s82, genes) -> None:
        print("JASPAR FIMO, девять наборов", flush=True)
        fisher_catalog(
            s82, genes, "fimo_jaspar", lambda sp: STORY / f"{sp}_fimo_proximal.tsv", jaspar_tf
        )
        print("HOCOMOCO FIMO, девять наборов", flush=True)
        fisher_catalog(s82, genes, "hocomoco", lambda sp: HOCO / f"{sp}_fimo.tsv", hoco_tf)

    def run_sea(s82, genes) -> None:
        s81 = s82.load_s81()
        sea_dir = OUT / "sea"
        jobs = []
        for sp in s81.SPECIES_ORDER:
            recs = s82.fasta_records(STORY / f"{sp}_proximal.fa")
            for stem, title in SETS:
                pos = sea_dir / "fasta" / f"{sp}_{stem}_pos.fa"
                neg = sea_dir / "fasta" / f"{sp}_{stem}_neg.fa"
                n_pos, n_neg = s82.split_set_fasta(recs, genes[stem], pos, neg)
                h, note = s82.sea_hofract(n_pos)
                jobs.append(
                    dict(
                        sp=sp,
                        stem=stem,
                        title=title,
                        pos=pos,
                        neg=neg,
                        n_pos=n_pos,
                        n_neg=n_neg,
                        hofract=h,
                        note=note,
                        outdir=sea_dir / "runs" / sp / stem,
                    )
                )

        def one(job):
            tsv = job["outdir"] / "sea.tsv"
            err_path = job["outdir"] / "sea.err"
            if tsv.exists() and tsv.stat().st_size > 50 and (job["outdir"] / "sea.ok").exists():
                return (job["sp"], job["stem"], "cached")
            if job["n_pos"] < 2 or job["n_neg"] < 20:
                return (job["sp"], job["stem"], "skip")
            job["outdir"].mkdir(parents=True, exist_ok=True)
            cmd = [
                str(SEA_BIN),
                "--p",
                str(job["pos"]),
                "--n",
                str(job["neg"]),
                "--m",
                str(JASPAR_MEME),
                "--text",
                "--noseqs",
                "--notrim",
                "--hofract",
                f"{job['hofract']:.6f}",
                "--seed",
                "1",
                "--pvalue",
                "--thresh",
                "1",
                "--verbosity",
                "1",
            ]
            with tsv.open("w") as out, err_path.open("w") as err:
                proc = subprocess.run(cmd, stdout=out, stderr=err)
            if proc.returncode != 0 or tsv.stat().st_size < 200:
                raise SystemExit(f"SEA {job['sp']} {job['stem']} rc={proc.returncode}")
            (job["outdir"] / "sea.ok").write_text(job["note"] + "\n")
            return (job["sp"], job["stem"], "ran")

        print(f"SEA, заданий {len(jobs)}", flush=True)
        with ThreadPoolExecutor(max_workers=12) as pool:
            futs = [pool.submit(one, job) for job in jobs]
            for fut in as_completed(futs):
                sp, stem, how = fut.result()
                print(f"  {sp} {stem} {how}", flush=True)
        fish_rows = []
        for job in jobs:
            tsv = job["outdir"] / "sea.tsv"
            if not tsv.exists() or tsv.stat().st_size < 50:
                continue
            tab = pd.read_csv(tsv, sep="\t")
            if "ALT_ID" not in tab.columns or "PVALUE" not in tab.columns:
                raise SystemExit(f"непонятный sea.tsv {tsv}")
            tab = tab.dropna(subset=["ALT_ID", "ID", "PVALUE"])
            tab["tf"] = tab["ALT_ID"].astype(str).str.upper()
            if tab["tf"].duplicated().any():
                raise SystemExit(f"SEA {job['sp']} {job['stem']}: повтор символа")
            for rec in tab.itertuples(index=False):
                fish_rows.append(
                    dict(
                        set=job["stem"],
                        title=job["title"],
                        species=job["sp"],
                        tf=rec.tf,
                        p=float(rec.PVALUE),
                        n_pos=job["n_pos"],
                        n_neg=job["n_neg"],
                    )
                )
        if not fish_rows:
            print("SEA: таблица p не собрана, сырые sea.tsv лежат в runs/", flush=True)
            return
        fish = s82.bh_within(pd.DataFrame(fish_rows))
        fish.to_csv(sea_dir / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = s82.ge3_table(fish)
        ge3.to_csv(sea_dir / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        print(f"SEA BH<0.05 {int((fish.bh < 0.05).sum())}  ge3 {len(ge3)}", flush=True)

    def run_monalisa(s82) -> None:
        s81 = s82.load_s81()
        mona = OUT / "monalisa"
        mona.mkdir(parents=True, exist_ok=True)

        def one(sp):
            out = mona / f"{sp}_binned.tsv"
            log = mona / f"{sp}.log"
            if (mona / f"{sp}.ok").exists() and out.exists() and (out.stat().st_size > 50):
                return (sp, "cached")
            cmd = [
                "Rscript",
                "--vanilla",
                str(R_SCRIPT),
                "--meme",
                str(JASPAR_MEME),
                "--fasta",
                str(STORY / f"{sp}_proximal.fa"),
                "--species",
                sp,
                "--geneset-dir",
                str(OUT / "genes"),
                "--sets",
                str(OUT / "sets.tsv"),
                "--out",
                str(out),
                "--cores",
                "16",
            ]
            with log.open("w") as fh:
                proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                raise SystemExit(f"monaLisa {sp} rc={proc.returncode} см. {log}")
            (mona / f"{sp}.ok").write_text("ok\n")
            return (sp, "ran")

        print("monaLisa, семь видов", flush=True)
        with ThreadPoolExecutor(max_workers=4) as pool:
            futs = [pool.submit(one, sp) for sp in s81.SPECIES_ORDER]
            for fut in as_completed(futs):
                sp, how = fut.result()
                print(f"  {sp} {how}", flush=True)
        parts = [pd.read_csv(mona / f"{sp}_binned.tsv", sep="\t") for sp in s81.SPECIES_ORDER]
        fish = s82.bh_within(pd.concat(parts, ignore_index=True))
        fish.to_csv(mona / "motif_fisher_proximal.tsv", sep="\t", index=False)
        ge3 = s82.ge3_table(fish)
        ge3.to_csv(mona / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
        print(f"monaLisa BH<0.05 {int((fish.bh < 0.05).sum())}  ge3 {len(ge3)}", flush=True)

    def main() -> None:
        stage = (
            (["bodymass"] + context.argv)[1] if len(["bodymass"] + context.argv) > 1 else "fisher"
        )
        s82 = load_s82()
        genes = load_sets()
        print("наборы", {k: len(v) for k, v in genes.items()}, flush=True)
        if stage in ("fisher", "all"):
            run_fisher(s82, genes)
        if stage in ("sea", "all"):
            run_sea(s82, genes)
        if stage in ("monalisa", "all"):
            run_monalisa(s82)

    return SimpleNamespace(**locals())


def _nine_set_slopes(context):
    """Nine locked sets, fibroblast and epiblast, Ruslan-style log-log.

    Membership is the frozen 22.09 files. Score is the one in
    Рабочие_наборы_без_пересечений: promoter CAGE already in
    M_fib_nearest2kb / M_epi_nearest2kb, missing gene dropped, y = log10(sum TPM),
    x = log10(adult mass). One grid per tissue. Sets are not pooled.

    Source: analyses/avian/scripts/95_nine_set_slopes.py at d230e56ab942."""
    import importlib.util
    import sys
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import stats

    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    KB = Path(
        str(
            context.resolve_path(
                "/path/to/user-home/chicken/kleiber_birds_gemcat/corpus_kleiber_2026-08-06/02_база_знаний"
            )
        )
    )
    M_FIB = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
    M_EPI = RES / "story_epi_traits" / "M_epi_nearest2kb.tsv"
    OUT = RES / "nine_set_slopes"
    PRES = RES / "presentation"
    EXPECT = {
        "oxphos_16": ((-0.86, 0.013), (-0.22, 0.64)),
        "folate_12": ((0.98, 0.0001), (0.903, 0.005)),
        "ribo_5": ((0.899, 0.006), (0.089, 0.85)),
        "green_hubs_7": ((0.845, 0.017), (0.804, 0.029)),
        "green_51_rest": ((0.926, 0.003), (0.727, 0.064)),
        "carbon_7": ((0.921, 0.003), (0.962, 0.0005)),
        "amino_5": ((0.921, 0.003), (0.964, 0.0005)),
        "nafld_4": ((-0.966, 0.0004), (-0.565, 0.19)),
        "epi_lipid_3": ((0.62, 0.14), (0.969, 0.0003)),
    }
    ORDER = list(EXPECT)

    def load40():
        mod = context.unit("avian.40_ruslan_style_scatters")
        return mod

    def load_sets() -> dict[str, tuple[str, list[str]]]:
        frames = [
            pd.read_csv(KB / "working_sets_disjoint_2026-09-22.tsv", sep="\t"),
            pd.read_csv(KB / "working_sets_overlap_ok_2026-09-22.tsv", sep="\t"),
        ]
        out = {}
        for df in frames:
            for sid, g in df.groupby("set_id", sort=False):
                genes = [normalize_symbol(x) for x in g["gene"]]
                genes = [x for x in genes if x]
                label = str(g["set_label"].iloc[0])
                out[str(sid)] = (label, genes)
        missing = [s for s in ORDER if s not in out]
        if missing:
            raise SystemExit(f"sets missing from the lock files: {missing}")
        return out

    def rows_of(wide: pd.DataFrame, genes: list[str]) -> tuple[list[tuple], list[str]]:
        hit = [g for g in genes if g in wide.index]
        dropped = [g for g in genes if g not in wide.index]
        sums = wide.loc[hit].sum(axis=0)
        rows = []
        for sp in SPECIES_ORDER:
            rows.append((sp.lower(), SPECIES[sp]["mass_g"], float(sums[sp])))
        return (rows, dropped)

    def main() -> None:
        s40 = load40()
        sets = load_sets()
        fib = pd.read_csv(M_FIB, sep="\t", index_col=0)[SPECIES_ORDER]
        epi = pd.read_csv(M_EPI, sep="\t", index_col=0)[SPECIES_ORDER]
        OUT.mkdir(parents=True, exist_ok=True)
        recs = []
        panels = {"fib": [], "epi": []}
        for sid in ORDER:
            label, genes = sets[sid]
            for tissue, wide, slot in (("fib", fib, 0), ("epi", epi, 1)):
                rows, dropped = rows_of(wide, genes)
                if any((r[2] <= 0 for r in rows)):
                    raise SystemExit(f"{sid} {tissue}: non-positive sum, log10 refused")
                x = [r[1] for r in rows]
                y = [r[2] for r in rows]
                fit = stats.pearsonr(np.log10(x), np.log10(y))
                exp_r, exp_p = EXPECT[sid][slot]
                if abs(fit.statistic - exp_r) > 0.015:
                    raise SystemExit(
                        f"{sid} {tissue}: r={fit.statistic:+.3f}, lock {exp_r:+.3f}; dropped {dropped}"
                    )
                recs.append(
                    dict(
                        set=sid,
                        label=label,
                        tissue=tissue,
                        n_listed=len(genes),
                        n_used=len(genes) - len(dropped),
                        dropped=",".join(dropped),
                        r=fit.statistic,
                        p=fit.pvalue,
                        lock_r=exp_r,
                        lock_p=exp_p,
                    )
                )
                panels[tissue].append((label, len(genes) - len(dropped), rows, float(fit.pvalue)))
                print(
                    f"{sid:16s} {tissue:3s} n={len(genes) - len(dropped):3d} r={fit.statistic:+.3f} p={fit.pvalue:.4g}"
                )
        for tissue in panels:
            ps = np.array([p for _, _, _, p in panels[tissue]], float)
            qs = stats.false_discovery_control(ps, method="bh")
            panels[tissue] = [
                (label, n, rows, q) for (label, n, rows, _), q in zip(panels[tissue], qs)
            ]
        for tissue, ylabel, stem in (
            ("fib", "log10_fib", "fig_nine_fib"),
            ("epi", "log10_epi", "fig_nine_epi"),
        ):
            fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.2))
            for ax, (label, n, rows, q) in zip(axes.ravel(), panels[tissue]):
                fit = s40.panel(ax, rows, ylabel)
                qtxt = f"{q:.3f}" if q >= 0.001 else f"{q:.3g}"
                ax.set_title(f"{label}  n={n}\ncor. {fit['r']:.3f}; BH = {qtxt}", fontsize=9, pad=6)
            fig.subplots_adjust(
                wspace=0.38, hspace=0.55, left=0.06, right=0.99, top=0.96, bottom=0.05
            )
            for folder in (OUT, PRES):
                folder.mkdir(parents=True, exist_ok=True)
                fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
                fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
            plt.close(fig)
        table = pd.DataFrame(recs)
        bh_map = {
            (tissue, label): q for tissue, items in panels.items() for label, _, _, q in items
        }
        table["bh"] = [bh_map[r.tissue, r.label] for r in table.itertuples()]
        table.to_csv(OUT / "nine_set_slopes.tsv", sep="\t", index=False)
        print(f"wrote {OUT}")

    return SimpleNamespace(**locals())


def _gene_sets(context):
    """Fix Ruslan's gene sets: his mitochondrial list (142, from his workspace) and the KEGG human pathways he used
    (OXPHOS hsa00190, Thermogenesis hsa04714) plus his negative control (cell cycle hsa04110) and ribosome hsa03010.

    Source: analyses/mammal_maps/code/00_gene_sets.py at d230e56ab942."""

    def main():
        import hashlib, json, re, urllib.request
        from pathlib import Path

        RUN = context.run_dir("maps")
        RES = RUN.parent
        raw = RUN / "reference/kegg"
        raw.mkdir(parents=True, exist_ok=True)

        def fetch(url, name):
            p = raw / name
            if not p.exists():
                p.write_bytes(urllib.request.urlopen(url, timeout=60).read())
            return p.read_text()

        sym = {}
        for line in fetch("https://rest.kegg.jp/list/hsa", "hsa_list.txt").splitlines():
            parts = line.split("\t")
            if len(parts) >= 4 and parts[1] == "CDS":
                sym[parts[0]] = parts[3].split(";")[0].split(",")[0].strip()
        sets = {}
        src = {}
        for name, pid in [
            ("OXPHOS_hsa00190", "hsa00190"),
            ("Thermogenesis_hsa04714", "hsa04714"),
            ("CellCycle_hsa04110_control", "hsa04110"),
            ("Ribosome_hsa03010", "hsa03010"),
        ]:
            txt = fetch(f"https://rest.kegg.jp/link/hsa/path:{pid}", f"{pid}.txt")
            ids = [l.split("\t")[1] for l in txt.splitlines() if "\t" in l]
            sets[name] = sorted({sym[i] for i in ids if i in sym})
            src[name] = {
                "file": f"reference/kegg/{pid}.txt",
                "sha256": hashlib.sha256((raw / f"{pid}.txt").read_bytes()).hexdigest(),
                "n_kegg_ids": len(ids),
            }
        mt = RES / "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt"
        sets["Ruslan_MT142"] = sorted({l.strip() for l in mt.read_text().splitlines() if l.strip()})
        src["Ruslan_MT142"] = {
            "file": "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt",
            "sha256": hashlib.sha256(mt.read_bytes()).hexdigest(),
        }
        src["hsa_list"] = {
            "file": "reference/kegg/hsa_list.txt",
            "sha256": hashlib.sha256((raw / "hsa_list.txt").read_bytes()).hexdigest(),
        }
        (RUN / "config/gene_sets.json").write_text(
            json.dumps({"sets": sets, "sources": src}, indent=1) + "\n"
        )
        print({k: len(v) for k, v in sets.items()})

    return SimpleNamespace(**locals())


def _sets_mass(context):
    """Step 1 of Ruslan's approach in mammalian fibroblasts: his gene sets versus adult mass.
    Fixes relative to the bird version: identical gene composition across species within a panel,
    phylogenetic GLS, and an expression-matched random-set null. Plan frozen in config/plan.json.

    Source: analyses/mammal_maps/code/01_sets_mass.py at d230e56ab942."""

    def main():
        import json
        from pathlib import Path
        import numpy as np
        import pandas as pd
        from scipy import linalg, stats

        RUN = context.run_dir("maps")
        RES = RUN.parent
        GI = RES / "growth_regulators_2026-10-21_v1/inputs"
        plan = json.loads((RUN / "config/plan.json").read_text())
        sets = json.loads((RUN / "config/gene_sets.json").read_text())["sets"]
        adm = pd.read_csv(
            RES / "nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv",
            sep="\t",
        )
        rng = np.random.default_rng(plan["seed"])

        def gls(x, y, C):
            X = np.column_stack([np.ones(len(x)), x])
            L = np.linalg.cholesky(C)
            W = linalg.solve_triangular(L, X, lower=True)
            v = linalg.solve_triangular(L, y, lower=True)
            b = np.linalg.lstsq(W, v, rcond=None)[0]
            df = len(x) - 2
            e = v - W @ b
            se = np.sqrt(e @ e / df * np.linalg.inv(W.T @ W)[1, 1])
            t = b[1] / se
            return (float(b[1]), float(se), float(t), float(2 * stats.t.sf(abs(t), df)))

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
            lin_s = (
                by_species(lin)
                .loc[lambda d: d.index.isin(fs.index)]
                .rename(index=fs.symbols.to_dict())
            )
            log_s = (
                by_species(lg)
                .loc[lambda d: d.index.isin(fs.index)]
                .rename(index=fs.symbols.to_dict())
            )
            ok = (lin_s > 0).all(axis=1) & np.isfinite(log_s).all(axis=1)
            lin_s, log_s = (lin_s[ok], log_s.loc[lin_s.index[ok]])
            mass = np.log10(a.loc[keep, "accepted_mass_g"].to_numpy(float))
            return dict(species=keep, C=C, mass=mass, lin=lin_s, log=log_s)

        rows, species_rows, gene_rows = ([], [], [])
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
                    rows.append(
                        dict(
                            panel=panel,
                            set=set_name,
                            n_species=len(d["species"]),
                            n_genes_list=len(members),
                            n_genes_used=len(g),
                            status="too_few_genes",
                        )
                    )
                    continue
                z_score = Z.loc[g].mean(0).to_numpy()
                log_sum = np.log10(d["lin"].loc[g].sum(0).to_numpy())
                r_z, p_r_z = stats.pearsonr(d["mass"], z_score)
                r_s, p_r_s = stats.pearsonr(d["mass"], log_sum)
                b, se, t, p = gls(d["mass"], z_score, d["C"])
                bs, ses, ts, ps = gls(d["mass"], log_sum, d["C"])
                gi = np.flatnonzero(np.isin(genes_all, g))
                need = pd.Series(dec[gi]).value_counts()
                draws = np.empty((plan["n_null"], len(gi)), dtype=int)
                col = 0
                for bnum, k in need.items():
                    pool = bins[bnum]
                    draws[:, col : col + k] = np.stack(
                        [rng.choice(pool, k, replace=False) for _ in range(plan["n_null"])]
                    )
                    col += k
                Zv = Z.to_numpy()
                null_scores = Zv[draws].mean(1).T
                t_null = gls_t_many(d["mass"], null_scores, d["C"])
                p_emp = (1 + np.sum(np.abs(t_null) >= abs(t))) / (1 + len(t_null))
                rows.append(
                    dict(
                        panel=panel,
                        set=set_name,
                        n_species=len(d["species"]),
                        n_genes_list=len(members),
                        n_genes_used=len(g),
                        status="tested",
                        pearson_r_z=r_z,
                        pearson_p_z=p_r_z,
                        PGLS_beta_z=b,
                        PGLS_se_z=se,
                        PGLS_t_z=t,
                        PGLS_p_z=p,
                        null_emp_p_z=p_emp,
                        null_median_t=float(np.median(t_null)),
                        pearson_r_logsum=r_s,
                        pearson_p_logsum=p_r_s,
                        PGLS_beta_logsum=bs,
                        PGLS_p_logsum=ps,
                    )
                )
                for s_, m_, zz, ls in zip(d["species"], d["mass"], z_score, log_sum):
                    species_rows.append(
                        dict(
                            panel=panel,
                            set=set_name,
                            species=s_,
                            log10_mass=m_,
                            mean_z=zz,
                            log10_sum_fixed=ls,
                        )
                    )
                tg = gls_t_many(d["mass"], Z.loc[g].to_numpy().T, d["C"])
                dfree = len(d["species"]) - 2
                for gene, tt in zip(g, tg):
                    gene_rows.append(
                        dict(
                            panel=panel,
                            set=set_name,
                            gene=gene,
                            PGLS_t=float(tt),
                            PGLS_p=float(2 * stats.t.sf(abs(tt), dfree)),
                        )
                    )
            print(panel, len(d["species"]), "species", Z.shape[0], "genes in fixed background")
        out = RUN / "results"
        out.mkdir(exist_ok=True)
        R = pd.DataFrame(rows)
        tested = R.status.eq("tested")
        from statsmodels.stats.multitest import multipletests

        R.loc[tested, "BH_q_emp_all"] = multipletests(
            R.loc[tested, "null_emp_p_z"], method="fdr_bh"
        )[1]
        R.to_csv(out / "01_sets_mass.tsv", sep="\t", index=False)
        pd.DataFrame(species_rows).to_csv(out / "01_species_scores.tsv", sep="\t", index=False)
        pd.DataFrame(gene_rows).to_csv(out / "01_gene_tests.tsv", sep="\t", index=False)
        pd.set_option("display.width", 250)
        print(
            R[
                [
                    "panel",
                    "set",
                    "n_species",
                    "n_genes_used",
                    "pearson_r_z",
                    "PGLS_beta_z",
                    "PGLS_p_z",
                    "null_emp_p_z",
                    "BH_q_emp_all",
                    "pearson_r_logsum",
                    "PGLS_p_logsum",
                ]
            ]
            .round(4)
            .to_string()
        )
        print("DONE01")

    return SimpleNamespace(**locals())


def _promoters_fimo(context):
    """Step 2a: native promoter windows and FIMO motif presence per species (plan_step2.json).
    Usage: 03_promoters_fimo.py prep            -> results/02_panel_backgrounds.json
           03_promoters_fimo.py <Species_name>  -> results/promoters/<species>.{fa,tsv}, results/motifs/<species>.npz

    Source: analyses/mammal_maps/code/03_promoters_fimo.py at d230e56ab942."""

    def main():
        import gzip, json, re, shutil, subprocess, sys, zipfile
        from collections import defaultdict
        from pathlib import Path
        import numpy as np
        import pandas as pd

        RUN = context.run_dir("maps")
        plan = json.loads((RUN / "config/plan_step2.json").read_text())
        MEME = (
            RUN.parent
            / "size_program_regulons_2026-10-21_v1/reference/JASPAR2024_one_per_symbol.meme"
        )
        FIMO = str(context.resolve_path("/path/to/user-home/miniconda3/envs/meme/bin/fimo"))
        UP, DOWN = (1000, 200)
        OUT = RUN / "results"
        if (["bodymass"] + context.argv)[1] == "prep":
            common_rm = context.unit("mammal_maps.common_rm")
            bg, members = ({}, {})
            for panel in ["ma2016", "gse105160", "wagner"]:
                d = common_rm.load(panel)
                Z = common_rm.background(d)
                a = common_rm.adm[
                    (common_rm.adm.panel == panel) & (common_rm.adm.status == "admitted")
                ].set_index("label")
                bg[panel] = sorted(Z.index)
                members[panel] = {lab: a.loc[lab, "scientific_name"] for lab in d["species"]}
            (OUT / "02_panel_backgrounds.json").write_text(
                json.dumps({"background": bg, "species": members}, indent=1)
            )
            print({p: len(v) for p, v in bg.items()})
            sys.exit()
        species = (["bodymass"] + context.argv)[1]
        info = json.loads((OUT / "02_panel_backgrounds.json").read_text())
        need = sorted(
            {
                g
                for p, labs in info["species"].items()
                if species in labs.values()
                for g in info["background"][p]
            }
        )
        work = RUN / "tmp" / species
        work.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(RUN / "reference/genomes" / f"{species}.zip") as z:
            names = z.namelist()
            fna = [n for n in names if n.endswith(".fna")]
            gff = [n for n in names if n.endswith(".gff")]
            assert len(fna) == 1 and len(gff) == 1, (species, fna, gff)
            for n in fna + gff:
                target = work / Path(n).name
                if not target.exists():
                    with z.open(n) as src, open(target, "wb") as dst:
                        shutil.copyfileobj(src, dst, 1 << 24)
        fna_p, gff_p = (work / Path(fna[0]).name, work / Path(gff[0]).name)
        attr_re = re.compile("([^=;]+)=([^;]*)")
        genes, mrna, exlen = ({}, {}, defaultdict(int))
        with open(gff_p) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9:
                    continue
                a = dict(attr_re.findall(c[8]))
                if c[2] == "gene" and a.get("gene_biotype", "protein_coding") == "protein_coding":
                    nm = a.get("Name") or a.get("gene")
                    if nm:
                        genes[a["ID"]] = nm.upper()
                elif c[2] == "mRNA" and "Parent" in a:
                    mrna[a["ID"]] = (a["Parent"], c[0], int(c[3]), int(c[4]), c[6])
                elif c[2] == "exon" and "Parent" in a:
                    exlen[a["Parent"]] += int(c[4]) - int(c[3]) + 1
        count = defaultdict(int)
        for nm in genes.values():
            count[nm] += 1
        best = {}
        for tid, (gid, chrom, s, e, strand) in mrna.items():
            nm = genes.get(gid)
            if nm is None or count[nm] != 1:
                continue
            L = exlen.get(tid, e - s + 1)
            if nm not in best or L > best[nm][0]:
                best[nm] = (L, chrom, s, e, strand)
        from Bio import SeqIO
        from Bio.Seq import Seq

        by_chrom = defaultdict(list)
        for nm in need:
            if nm in best:
                by_chrom[best[nm][1]].append(nm)
        rows = []
        pdir, mdir = (OUT / "promoters", OUT / "motifs")
        pdir.mkdir(exist_ok=True)
        mdir.mkdir(exist_ok=True)
        fa = pdir / f"{species}.fa"
        with open(fa, "w") as out:
            for rec in SeqIO.parse(str(fna_p), "fasta"):
                if rec.id not in by_chrom:
                    continue
                chrom_seq = str(rec.seq)
                for nm in by_chrom[rec.id]:
                    _, chrom, s, e, strand = best[nm]
                    if strand == "+":
                        lo, hi = (s - 1 - UP, s - 1 + DOWN)
                    else:
                        lo, hi = (e - DOWN, e + UP)
                    if lo < 0 or hi > len(chrom_seq):
                        continue
                    seq = chrom_seq[lo:hi]
                    if strand == "-":
                        seq = str(Seq(seq).reverse_complement())
                    if len(seq) != UP + DOWN or seq.upper().count("N") > 0.2 * len(seq):
                        continue
                    up = seq.upper()
                    gc = (up.count("G") + up.count("C")) / max(1, len(up) - up.count("N"))
                    rows.append(
                        dict(
                            gene=nm,
                            chrom=chrom,
                            strand=strand,
                            window_start0=lo,
                            window_end=hi,
                            gc=gc,
                        )
                    )
                    out.write(f">{nm}\n{seq}\n")
        pd.DataFrame(rows).to_csv(pdir / f"{species}.tsv", sep="\t", index=False)
        motif_ids = [
            l.split()[1] + "|" + l.split()[2]
            for l in MEME.read_text().splitlines()
            if l.startswith("MOTIF")
        ]
        mpos = {m.split("|")[0]: i for i, m in enumerate(motif_ids)}
        gpos = {r["gene"]: i for i, r in enumerate(rows)}
        P = np.zeros((len(rows), len(motif_ids)), dtype=bool)
        lines = fa.read_text().splitlines()
        chunks = []
        for k in range(4):
            cf = work / f"chunk{k}.fa"
            cf.write_text(
                "\n".join((lines[2 * i + j] for i in range(k, len(lines) // 2, 4) for j in (0, 1)))
                + "\n"
            )
            chunks.append((cf, open(work / f"chunk{k}.fimo", "w")))
        procs = [
            subprocess.Popen(
                [FIMO, "--text", "--thresh", "1e-4", "--verbosity", "1", str(MEME), str(cf)],
                stdout=fh,
            )
            for cf, fh in chunks
        ]
        assert all((p.wait() == 0 for p in procs))
        for cf, fh in chunks:
            fh.close()
            with open(work / (cf.stem + ".fimo")) as f:
                for line in f:
                    if line.startswith("motif_id") or line.startswith("#"):
                        continue
                    c = line.split("\t")
                    if len(c) > 3 and c[2] in gpos:
                        P[gpos[c[2]], mpos[c[0]]] = True
        np.savez_compressed(
            mdir / f"{species}.npz",
            presence=P,
            genes=np.array([r["gene"] for r in rows]),
            motifs=np.array(motif_ids),
            gc=np.array([r["gc"] for r in rows]),
        )
        shutil.rmtree(work)
        print(
            species,
            "needed",
            len(need),
            "unique_pc_genes",
            len(best),
            "windows",
            len(rows),
            "positive_pairs",
            int(P.sum()),
        )

    return SimpleNamespace(**locals())


def _enrichment_mass(context):
    """Step 2b: per species x set x motif Mantel-Haenszel log odds ratio (GC quintiles, Haldane 0.5),
    then PGLS of log OR on log10 adult mass per panel, BH over everything, cross-panel replication (plan_step2.json).

    Source: analyses/mammal_maps/code/04_enrichment_mass.py at d230e56ab942."""

    def main():
        import json, sys
        from pathlib import Path
        import numpy as np
        import pandas as pd
        from scipy import stats
        from statsmodels.stats.multitest import multipletests

        RUN = context.run_dir("maps")
        common_rm = context.unit("mammal_maps.common_rm")
        plan = json.loads((RUN / "config/plan_step2.json").read_text())
        sets = json.loads((RUN / "config/gene_sets.json").read_text())["sets"]
        info = json.loads((RUN / "results/02_panel_backgrounds.json").read_text())
        OUT = RUN / "results"

        def mh_logor(pos, fg, strata):
            """pos: genes x motifs bool; fg: genes bool; strata: genes int. Returns log MH OR per motif."""
            num = np.zeros(pos.shape[1])
            den = np.zeros(pos.shape[1])
            for s in np.unique(strata):
                m = strata == s
                f, b = (fg & m, ~fg & m)
                if f.sum() == 0 or b.sum() == 0:
                    continue
                a = pos[f].sum(0) + 0.5
                c = f.sum() - pos[f].sum(0) + 0.5
                bb = pos[b].sum(0) + 0.5
                d = b.sum() - pos[b].sum(0) + 0.5
                n = a + bb + c + d
                num += a * d / n
                den += bb * c / n
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
            return (B[1], t, 2 * stats.t.sf(np.abs(t), df))

        lor_rows, test_rows = ([], [])
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
                P, g, gc = (z["presence"][keep], g[keep], z["gc"][keep])
                strata = pd.qcut(gc, 5, labels=False, duplicates="drop")
                for s, members in sets.items():
                    fg = np.isin(g, members)
                    lo = mh_logor(P, fg, strata)
                    L[s].append(lo)
                    lor_rows.append(
                        pd.DataFrame(
                            {
                                "panel": panel,
                                "species": lab,
                                "set": s,
                                "n_fg": int(fg.sum()),
                                "n_bg": int((~fg).sum()),
                                "motif": motifs,
                                "log_MH_OR": lo,
                            }
                        )
                    )
            null_inputs[panel] = (mass, C, {s: np.array(L[s]) for s in sets})
            for s in sets:
                Y = np.array(L[s])
                beta, t, p = gls_many(mass, Y, C)
                r = np.array(
                    [
                        stats.pearsonr(mass, Y[:, j])[0] if np.std(Y[:, j]) > 0 else np.nan
                        for j in range(Y.shape[1])
                    ]
                )
                test_rows.append(
                    pd.DataFrame(
                        {
                            "panel": panel,
                            "set": s,
                            "n_species": len(have),
                            "motif": motifs,
                            "PGLS_beta": beta,
                            "PGLS_t": t,
                            "PGLS_p": p,
                            "pearson_r": r,
                        }
                    )
                )
            print(panel, "species with genomes", len(have), "of", len(labs))
        T = pd.concat(test_rows, ignore_index=True)
        ok = np.isfinite(T.PGLS_p)
        T.loc[ok, "BH_q_all"] = multipletests(T.loc[ok, "PGLS_p"], method="fdr_bh")[1]
        T["tf"] = T.motif.str.split("|").str[1].str.upper()
        pd.concat(lor_rows, ignore_index=True).to_csv(
            OUT / "04_species_logOR.tsv.gz", sep="\t", index=False
        )
        T.to_csv(OUT / "04_motif_mass_tests.tsv.gz", sep="\t", index=False)
        rep = []
        for (s, m), g in T.groupby(["set", "motif"]):
            sig = g[g.PGLS_p < 0.05]
            for sign in (1, -1):
                k = sig[np.sign(sig.PGLS_beta) == sign]
                if len(k) >= 2:
                    rep.append(
                        dict(
                            set=s,
                            motif=m,
                            tf=g.tf.iloc[0],
                            direction="up_with_mass" if sign > 0 else "down_with_mass",
                            panels=",".join(k.panel),
                            min_q=float(g.BH_q_all.min()),
                            **{f"p_{r.panel}": r.PGLS_p for r in g.itertuples()},
                            **{f"beta_{r.panel}": r.PGLS_beta for r in g.itertuples()},
                        )
                    )
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
                u = ((p < 0.05) & (b > 0)).astype(int)
                dn = ((p < 0.05) & (b < 0)).astype(int)
                up = u if up is None else up + u
                down = dn if down is None else down + dn
            null_counts.append(int(((up >= 2) | (down >= 2)).sum()))
        null_counts = np.array(null_counts)
        obs = len(R.drop_duplicates(["set", "motif"])) if len(R) else 0
        rep_null = dict(
            observed_replicated_pairs=obs,
            null_median=float(np.median(null_counts)),
            null_q95=float(np.quantile(null_counts, 0.95)),
            empirical_p=float((1 + np.sum(null_counts >= obs)) / (1 + len(null_counts))),
        )
        (OUT / "04_replication_null.json").write_text(json.dumps(rep_null, indent=1))
        print("replication null:", rep_null)
        pre = T[T.tf.isin(plan["prespecified_from_ruslan"])]
        pre.to_csv(OUT / "04_ruslan_prespecified.tsv", sep="\t", index=False)
        print(
            "tests",
            len(T),
            "q<0.05:",
            int((T.BH_q_all < 0.05).sum()),
            "min q",
            float(T.BH_q_all.min()),
        )
        print("replicated (same sign, p<0.05 in >=2 panels):", len(R))
        expected = sum((1 for _ in T.groupby(["set", "motif"])))
        print("set x motif pairs:", expected)
        if len(R):
            print(R.sort_values("min_q").head(40).to_string())
        print(
            pre.sort_values("PGLS_p")
            .head(30)[
                ["panel", "set", "tf", "n_species", "PGLS_beta", "PGLS_p", "BH_q_all", "pearson_r"]
            ]
            .to_string()
        )
        print("DONE04")

    return SimpleNamespace(**locals())


def _gsea_phylo(context):
    """Step 3: genome-wide phylogenetic GSEA in three mammalian fibroblast panels with cross-panel concordance (plan_step3.json).

    Source: analyses/mammal_maps/code/05_gsea_phylo.py at d230e56ab942."""

    def main():
        import json, os, sys

        os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
        from pathlib import Path
        import numpy as np
        import pandas as pd
        from numba import njit, prange
        from scipy import linalg, stats
        from statsmodels.stats.multitest import multipletests

        RUN = context.run_dir("maps")
        common_rm = context.unit("mammal_maps.common_rm")
        plan = json.loads((RUN / "config/plan_step3.json").read_text())
        GS = Path(
            str(
                context.resolve_path(
                    "/path/to/user-home/chicken/kleiber_birds_gemcat/data/genesets"
                )
            )
        )
        N = 1000
        rng = np.random.default_rng(plan["seed"])
        catalog = {}
        for name, f in [
            ("GO_BP", "mammal_catalogs/gobp.tsv"),
            ("Reactome", "mammal_catalogs/reactome.tsv"),
            ("Hallmark", "mammal_catalogs/hallmark.tsv"),
        ]:
            for r in pd.read_csv(GS / f, sep="\t").itertuples():
                catalog[f"{name}:{r.gs_name}"] = r.genes.split(";")
        for r in pd.read_csv(GS / "kegg_hsa_all_pathways.tsv", sep="\t").itertuples():
            catalog[f"KEGG:{r.pathway} {r.name}"] = str(r.genes).split(";")
        catalog["Ruslan:OXPHOS_133"] = (GS / "ruslan_oxphos_133.txt").read_text().split()
        catalog["Ruslan:Thermogenesis_229"] = (GS / "ruslan_thermo_229.txt").read_text().split()
        catalog["Ruslan:MT142"] = (
            (RUN.parent / "seven_ruslan_sums_2026-10-20_v1/inputs/data_pc__mt_genes.txt")
            .read_text()
            .split()
        )

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
                idx = set_idx[set_ptr[s] : set_ptr[s + 1]]
                k = len(idx)
                pos = np.empty(k, np.int64)
                ww = np.empty(k)
                for j in range(k):
                    pos[j] = order_rank[idx[j]]
                    ww[j] = w[idx[j]]
                o = np.argsort(pos)
                pos = pos[o]
                ww = ww[o]
                total = ww.sum()
                miss = 1.0 / (n - k)
                acc = 0.0
                peak = 0.0
                prev = -1
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
            rank = np.empty(len(t), np.int64)
            rank[order] = np.arange(len(t))
            return all_es(rank, np.abs(t), set_ptr, set_idx, len(t))

        def nes_p(obs, null):
            """obs: sets; null: draws x sets. Same-sign NES and p."""
            pos_mean = np.where(null > 0, null, np.nan)
            neg_mean = np.where(null < 0, null, np.nan)
            mp = np.nanmean(pos_mean, 0)
            mn = np.nanmean(neg_mean, 0)
            nes = np.where(obs >= 0, obs / mp, -obs / mn)
            cnt_pos = np.sum(null >= obs[None, :], 0)
            n_pos = np.sum(null >= 0, 0)
            cnt_neg = np.sum(null <= obs[None, :], 0)
            n_neg = np.sum(null < 0, 0)
            p = np.where(obs >= 0, (cnt_pos + 1) / (n_pos + 1), (cnt_neg + 1) / (n_neg + 1))
            return (nes, np.minimum(p, 1.0))

        panel_res, null_sig = ({}, {})
        for panel in plan["panels"]:
            d = common_rm.load(panel)
            Z = common_rm.background(d)
            genes = list(Z.index)
            gi = {g: i for i, g in enumerate(genes)}
            names, ptr, idx = ([], [0], [])
            for nm, members in catalog.items():
                h = sorted({gi[g] for g in members if g in gi})
                if 15 <= len(h) <= 500:
                    names.append(nm)
                    idx += h
                    ptr.append(len(idx))
            ptr = np.array(ptr, np.int64)
            idx = np.array(idx, np.int64)
            C = d["C"]
            Linv = linalg.inv(np.linalg.cholesky(C))
            x = d["mass"]
            Zt = Z.to_numpy().T
            t_obs = gene_t(Zt, x, Linv)
            es_obs = es_for(t_obs, ptr, idx)
            null = np.empty((N, len(names)))
            for b in range(N):
                sim = rng.multivariate_normal(np.zeros(len(x)), C)
                sim = (sim - sim.mean()) / sim.std() * x.std() + x.mean()
                null[b] = es_for(gene_t(Zt, sim, Linv), ptr, idx)
            nes, p = nes_p(es_obs, null)
            q = multipletests(p, method="fdr_bh")[1]
            rk = np.empty_like(null)
            for s in range(null.shape[1]):
                col = null[:, s]
                pos = col >= 0
                sp = np.sort(col[pos])
                sn = np.sort(col[~pos])
                r = np.empty(N)
                r[pos] = (len(sp) - np.searchsorted(sp, col[pos], "left") + 0) / max(1, len(sp))
                r[~pos] = np.searchsorted(sn, col[~pos], "right") / max(1, len(sn))
                rk[:, s] = r
            null_sig[panel] = pd.DataFrame(np.sign(null) * (rk < 0.05), columns=names)
            panel_res[panel] = pd.DataFrame(
                {
                    "set": names,
                    "size": np.diff(ptr),
                    "ES": es_obs,
                    "NES": nes,
                    "p": p,
                    "q_BH_panel": q,
                }
            )
            pd.DataFrame({"gene": genes, "PGLS_t": t_obs}).to_csv(
                RUN / f"results/05_gene_rank_{panel}.tsv", sep="\t", index=False
            )
            print(
                panel,
                len(genes),
                "genes",
                len(names),
                "sets",
                "q<0.05:",
                int((q < 0.05).sum()),
                "min q",
                float(q.min()),
            )
        common = sorted(set.intersection(*[set(v.set) for v in panel_res.values()]))
        W = pd.DataFrame({"set": common})
        for panel, r in panel_res.items():
            r = r.set_index("set").loc[common]
            for c in ["size", "NES", "p", "q_BH_panel"]:
                W[f"{c}_{panel}"] = r[c].to_numpy()
        sig = np.stack([np.sign(W[f"NES_{p}"]) * (W[f"p_{p}"] < 0.05) for p in plan["panels"]])
        W["concordant_3"] = np.abs(sig.sum(0)) == 3
        W["concordant_2plus"] = np.maximum((sig > 0).sum(0), (sig < 0).sum(0)) >= 2
        z = np.stack(
            [np.sign(W[f"NES_{p}"]) * stats.norm.isf(W[f"p_{p}"] / 2) for p in plan["panels"]]
        )
        W["stouffer_z"] = z.sum(0) / np.sqrt(3)
        W["stouffer_p"] = 2 * stats.norm.sf(np.abs(W.stouffer_z))
        W["stouffer_q"] = multipletests(W.stouffer_p, method="fdr_bh")[1]
        nsig = np.stack([null_sig[p][common].to_numpy() for p in plan["panels"]])
        null3 = (np.abs(nsig.sum(0)) == 3).sum(1)
        null2 = (np.maximum((nsig > 0).sum(0), (nsig < 0).sum(0)) >= 2).sum(1)
        summ = dict(
            n_common_sets=len(common),
            observed_concordant_3=int(W.concordant_3.sum()),
            null3_median=float(np.median(null3)),
            null3_q95=float(np.quantile(null3, 0.95)),
            p_concordant_3=float((1 + (null3 >= W.concordant_3.sum()).sum()) / (1 + N)),
            observed_concordant_2plus=int(W.concordant_2plus.sum()),
            null2_median=float(np.median(null2)),
            null2_q95=float(np.quantile(null2, 0.95)),
            p_concordant_2plus=float((1 + (null2 >= W.concordant_2plus.sum()).sum()) / (1 + N)),
            stouffer_q_lt_0_05=int((W.stouffer_q < 0.05).sum()),
            per_panel_q_lt_0_05={p: int((r.q_BH_panel < 0.05).sum()) for p, r in panel_res.items()},
        )
        for p, r in panel_res.items():
            r.to_csv(RUN / f"results/05_gsea_{p}.tsv", sep="\t", index=False)
        W.sort_values("stouffer_p").to_csv(
            RUN / "results/05_gsea_cross_panel.tsv", sep="\t", index=False
        )
        (RUN / "results/05_summary.json").write_text(json.dumps(summ, indent=1))
        print(json.dumps(summ, indent=1))
        pd.set_option("display.width", 250)
        pd.set_option("display.max_colwidth", 60)
        cols = (
            ["set"]
            + [f"NES_{p}" for p in plan["panels"]]
            + [f"p_{p}" for p in plan["panels"]]
            + ["stouffer_q"]
        )
        print("concordant in all three:\n", W[W.concordant_3][cols].round(4).to_string())
        print("top Stouffer:\n", W.sort_values("stouffer_p").head(25)[cols].round(4).to_string())
        print("Ruslan sets:\n", W[W.set.str.startswith("Ruslan")][cols].round(4).to_string())
        for p, r in panel_res.items():
            print(p, "top by q:\n", r.sort_values("p").head(10).round(4).to_string())
        print("DONE05")

    return SimpleNamespace(**locals())


def _common_rm(context):
    """Panel loader shared with 01_sets_mass.py (copied verbatim) plus the fixed per-panel background used in step 1.

    Source: analyses/mammal_maps/code/common_rm.py at d230e56ab942."""
    from pathlib import Path
    import numpy as np
    import pandas as pd

    RUN = context.run_dir("maps")
    RES = RUN.parent
    GI = RES / "growth_regulators_2026-10-21_v1/inputs"
    adm = pd.read_csv(
        RES / "nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv",
        sep="\t",
    )

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
        lin_s = (
            by_species(lin).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
        )
        log_s = (
            by_species(lg).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
        )
        ok = (lin_s > 0).all(axis=1) & np.isfinite(log_s).all(axis=1)
        lin_s, log_s = (lin_s[ok], log_s.loc[lin_s.index[ok]])
        mass = np.log10(a.loc[keep, "accepted_mass_g"].to_numpy(float))
        return dict(species=keep, C=C, mass=mass, lin=lin_s, log=log_s)

    def background(d):
        Z = d["log"].sub(d["log"].mean(1), axis=0).div(d["log"].std(1), axis=0)
        return Z[np.isfinite(Z).all(1)]

    return SimpleNamespace(**locals())


def _describe_correlations(context):
    """Requested descriptive interpretation; retains phylogenetic results separately.

    Source: analyses/mammal_maps/code/describe_correlations.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import argparse, json, hashlib
        import numpy as np, pandas as pd
        from scipy import stats

        ap = argparse.ArgumentParser()
        ap.add_argument("--out", type=Path, required=True)
        ap.add_argument("--results-root", type=Path, required=True)
        ap.add_argument("--genesets-root", type=Path, required=True)
        a = ap.parse_args(context.argv)
        O = a.out
        amendment = dict(
            user_request="Наша задача просто выделить млекопитающие связи; если что-то коррелирует с размером, это уже интересно.",
            primary_interpretation="Descriptive Pearson correlation of log10 linear normalized map sum with log10 body mass. Phylogeny is a separate diagnostic, not a gate.",
            timing="After initial full screen, directly requested by user; Pearson outcomes were already calculated and are not new discoveries from refitting.",
            multiple_testing="BH over all eligible maps within each panel and composition mode; replication shown as context, not an admission criterion.",
            minimum_measured_genes=5,
            human_review="not_completed",
            training_enabled=False,
        )
        (O / "config/descriptive_focus_amendment.json").write_text(
            json.dumps(amendment, ensure_ascii=False, indent=2)
        )
        T = pd.read_csv(O / "results/map_tests.tsv.gz", sep="\t")
        cat = pd.read_csv(O / "config/catalog.tsv.gz", sep="\t")
        hits = T[T.q_Pearson_panel.lt(0.05)].copy().sort_values(["mode", "panel", "p_Pearson"])
        hits.to_csv(O / "results/Pearson_BH_hits.tsv", sep="\t", index=False)
        ids = set(hits.set_id)
        score_rows = []
        diag = []
        contrib = []
        err_r = err_p = 0.0
        checks = 0
        pairrows = []

        def bh_manual(ps):
            p = np.asarray(ps)
            order = np.argsort(p)
            q = np.empty(len(p))
            q[order] = np.minimum(
                1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
            )
            return q

        for (mode, panel), z in T.groupby(["mode", "panel"]):
            design = np.load(O / "results" / f"{panel}_design.npz")
            sc = np.load(O / "results" / f"{mode}_{panel}_scores.npz")
            lookup = {cat.iloc[k].set_id: i for i, k in enumerate(sc["catalog_index"])}
            x = design["mass_log10"]
            sp = list(design["species"])
            assert np.max(abs(bh_manual(z.p_Pearson) - z.q_Pearson_panel)) < 1e-10
            chosen = pd.concat(
                [
                    z[z.set_id.isin(ids)],
                    z.nsmallest(10, "p_Pearson"),
                    z.sample(25, random_state=2307),
                ]
            ).drop_duplicates("set_id")
            src = a.results_root / (
                "nonbird_regulatory_followup_2026-10-21_v1/results/ma2016_admitted_M.tsv.gz"
                if panel == "ma2016"
                else f"nonbird_common_state_2026-10-21_v1/results/{panel}_admission_before_gene_filter_M.tsv.gz"
            )
            M = pd.read_csv(src, sep="\t", index_col=0)[sp]
            for row in chosen.itertuples():
                y = sc["log10_sum"][lookup[row.set_id]]
                res = stats.pearsonr(x, y)
                err_r = max(err_r, abs(res.statistic - row.r_Pearson))
                err_p = max(err_p, abs(res.pvalue - row.p_Pearson))
                checks += 1
                if row.set_id not in ids:
                    continue
                rr = []
                pp = []
                bb = []
                for i, s in enumerate(sp):
                    keep = np.arange(len(sp)) != i
                    fit = stats.linregress(x[keep], y[keep])
                    rr.append(fit.rvalue)
                    pp.append(fit.pvalue)
                    bb.append(fit.slope)
                gm = M.loc[row.genes.split(";")]
                share = gm.div(gm.sum(0), axis=1).mean(1).sort_values(ascending=False)
                diag.append(
                    dict(
                        mode=mode,
                        panel=panel,
                        set_id=row.set_id,
                        n_genes=len(gm),
                        loo_min_r=min(rr),
                        loo_max_r=max(rr),
                        loo_min_abs_r=min(map(abs, rr)),
                        loo_sign_fraction=float(np.mean(np.sign(rr) == np.sign(row.r_Pearson))),
                        loo_max_p=max(pp),
                        top_gene=share.index[0],
                        top_gene_mean_share=float(share.iloc[0]),
                        top3_mean_share=float(share.iloc[:3].sum()),
                    )
                )
                for i, s in enumerate(sp):
                    score_rows.append(
                        dict(
                            mode=mode,
                            panel=panel,
                            set_id=row.set_id,
                            species=s,
                            log10_mass=x[i],
                            log10_sum=y[i],
                            linear_sum=10 ** y[i],
                        )
                    )
                for gene in share.index:
                    r, p = stats.pearsonr(x, np.log10(gm.loc[gene].to_numpy() + 1))
                    contrib.append(
                        dict(
                            mode=mode,
                            panel=panel,
                            set_id=row.set_id,
                            gene=gene,
                            mean_linear_share=float(share[gene]),
                            individual_log1p_r=float(r),
                            individual_log1p_p=float(p),
                        )
                    )
        pd.DataFrame(score_rows).to_csv(
            O / "results/Pearson_hit_species_scores.tsv", sep="\t", index=False
        )
        pd.DataFrame(diag).to_csv(O / "results/Pearson_hit_diagnostics.tsv", sep="\t", index=False)
        pd.DataFrame(contrib).to_csv(O / "results/Pearson_hit_genes.tsv.gz", sep="\t", index=False)
        for mode, z in T.groupby("mode"):
            pv = z.pivot(index="set_id", columns="panel", values="p_Pearson")
            rv = z.pivot(index="set_id", columns="panel", values="r_Pearson")
            same = rv.gt(0).all(1) | rv.lt(0).all(1)
            complete = pv.notna().all(1)
            pmax = pv.max(1).where(same, 1.0)
            q = bh_manual(pmax[complete])
            f = pd.DataFrame(
                dict(
                    set_id=pv.index[complete],
                    mode=mode,
                    p_max_same_direction=pmax[complete],
                    q_all3=q,
                )
            )
            pairrows.append(f)
        pd.concat(pairrows).to_csv(O / "results/Pearson_all3_context.tsv.gz", sep="\t", index=False)
        summary = dict(
            question="Какие суммы экспрессии полных измеренных частей KEGG/GO/Reactome коррелируют с массой млекопитающих?",
            primary_metric=amendment["primary_interpretation"],
            panels=[],
            common_composition_hits=hits[hits["mode"].eq("all3_common")][
                ["panel", "set_id", "n_measured", "r_Pearson", "p_Pearson", "q_Pearson_panel"]
            ].to_dict("records"),
            significant_main_terms=hits[hits["mode"].eq("panel_available")][
                [
                    "panel",
                    "set_id",
                    "n_measured",
                    "n_canonical",
                    "r_Pearson",
                    "beta_OLS",
                    "p_Pearson",
                    "q_Pearson_panel",
                ]
            ].to_dict("records"),
            independent_checks=dict(
                Pearson_fits=checks,
                max_r_error=err_r,
                max_p_error=err_p,
                BH_every_row_verified=True,
            ),
            interpretation="18 whole measured-map associations in GSE105160; contextual differences across panels are visible and not used to dismiss those descriptive associations. Three additional Ma2016 associations refer to the smaller common composition, not the same complete sums.",
            human_review="not_completed",
            training_enabled=False,
        )
        for (mode, panel), z in T.groupby(["mode", "panel"]):
            summary["panels"].append(
                dict(
                    mode=mode,
                    panel=panel,
                    n_species=int(z.n_species.iloc[0]),
                    n_maps=len(z),
                    nominal=int(z.p_Pearson.lt(0.05).sum()),
                    BH05=int(z.q_Pearson_panel.lt(0.05).sum()),
                    min_q=float(z.q_Pearson_panel.min()),
                )
            )
        assert err_r < 1e-09 and err_p < 1e-09
        (O / "results/descriptive_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2)
        )
        (O / "audit/Pearson_verification.json").write_text(
            json.dumps(summary["independent_checks"], indent=2)
        )
        src = a.genesets_root / "mammal_catalogs/msigdb.2026.1.Hs.C5.rds"
        (O / "audit/GO_source.json").write_text(
            json.dumps(
                dict(
                    reference="mammal_catalogs/msigdb.2026.1.Hs.C5.rds",
                    sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                    bytes=src.stat().st_size,
                    conversion="code/export_go.R",
                ),
                indent=2,
            )
        )
        print(
            json.dumps(
                dict(panels=summary["panels"], checks=summary["independent_checks"]),
                ensure_ascii=False,
                indent=2,
            )
        )

    return SimpleNamespace(**locals())


def _run_all_maps(context):
    """Whole measured-map linear sums, no trait-driven membership selection.

    Source: analyses/mammal_maps/code/run_all_maps.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import argparse, json, hashlib, sys, platform, itertools
        import numpy as np, pandas as pd, scipy
        from scipy import stats, sparse

        ap = argparse.ArgumentParser()
        ap.add_argument("--results-root", type=Path, required=True)
        ap.add_argument("--genesets-root", type=Path, required=True)
        ap.add_argument("--out", type=Path, required=True)
        a = ap.parse_args(context.argv)
        R, GS, O = (a.results_root, a.genesets_root, a.out)
        for d in ["config", "results", "audit"]:
            (O / d).mkdir(parents=True, exist_ok=True)
        manifest = []

        def read(p, root, label, **kw):
            b = p.read_bytes()
            manifest.append(
                dict(
                    origin=label,
                    reference=str(p.relative_to(root)),
                    sha256=hashlib.sha256(b).hexdigest(),
                    bytes=len(b),
                )
            )
            return pd.read_csv(p, sep="\t", **kw)

        def bh(p):
            p = np.asarray(p, float)
            q = np.ones(len(p))
            ok = np.flatnonzero(np.isfinite(p))
            q[ok] = stats.false_discovery_control(p[ok])
            return q

        plan = dict(
            primary="log10(sum of all measured fixed-map genes in linear normalized RNA), PGLS lambda=1",
            catalogs="all KEGG and Reactome stored maps; all GO BP/MF/CC in stored MSigDB2026.1 C5; no Hallmark/HPO",
            admission="at least 5 measured genes and positive map sum in every admitted species; measured finite nonnegative genes, true zeros retained; no upper map-size limit",
            coverage="canonical gene count and measured fraction reported, no zero imputation; coverage>=50% secondary stratum",
            composition_modes=["panel_available", "all3_common"],
            multiple_testing="BH within panel and across all panel-map tests separately per composition mode; same-sign max(two-sided p) conjunction tests across all3 or any panel pair, BH over all eligible terms/pairs; no independence assumption",
            sensitivity="ordinary unweighted Pearson/OLS on log10 scales; fixed diagonal covariance lambda0; leave-one-species for globally BH significant or conjunction nominal terms",
            null="1000 Brownian covariance simulated traits per panel; common simulated trait shared by every map within panel; descriptive hit-count calibration only; no cross-panel independence assumption",
            seed=2026102307,
            B=1000,
            trait="corrected species admission; Microtus pennsylvanicus 49g",
            human_review="not_completed",
            training_enabled=False,
        )
        (O / "config/plan.json").write_text(json.dumps(plan, indent=2))
        cats = []
        for source, file in [
            ("KEGG", "kegg_hsa_all_pathways.tsv"),
            ("Reactome", "mammal_catalogs/reactome.tsv"),
        ]:
            t = read(GS / file, GS, "gene_sets")
            for row in t.itertuples():
                name = str(row.pathway) + " " + row.name if source == "KEGG" else row.gs_name
                cats.append(
                    dict(source=source, term=name, genes=sorted(set(str(row.genes).split(";"))))
                )
        t = read(O / "config/go_all.tsv", O, "current_run")
        for row in t.itertuples():
            cats.append(
                dict(
                    source=row.collection.replace(":", "_"),
                    term=row.term,
                    genes=sorted(set(str(row.genes).split(";"))),
                )
            )
        catalog = pd.DataFrame(cats)
        catalog["set_id"] = catalog.source + ":" + catalog.term
        assert catalog.set_id.is_unique
        catalog["n_canonical"] = catalog.genes.map(len)
        catalog.assign(genes=catalog.genes.map(";".join)).to_csv(
            O / "config/catalog.tsv.gz", sep="\t", index=False
        )
        traits = read(
            R / "nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv",
            R,
            "results",
        )
        data = {}
        panels = ["ma2016", "gse105160", "wagner"]
        speciesrows = []
        for panel in panels:
            p = R / (
                "nonbird_regulatory_followup_2026-10-21_v1/results/ma2016_admitted_M.tsv.gz"
                if panel == "ma2016"
                else f"nonbird_common_state_2026-10-21_v1/results/{panel}_admission_before_gene_filter_M.tsv.gz"
            )
            M = read(p, R, "results", index_col=0)
            p = R / f"growth_regulators_2026-10-21_v1/inputs/{panel}_full.cov.npz"
            b = p.read_bytes()
            manifest.append(
                dict(
                    origin="results",
                    reference=str(p.relative_to(R)),
                    sha256=hashlib.sha256(b).hexdigest(),
                    bytes=len(b),
                )
            )
            cv = np.load(p, allow_pickle=True)
            sp = list(map(str, cv["species"]))
            tr = traits[traits.panel.eq(panel) & traits.status.eq("admitted")].set_index("label")
            keep = [s for s in sp if s in tr.index]
            ix = [sp.index(s) for s in keep]
            C = cv["C"][np.ix_(ix, ix)]
            C /= np.diag(C).mean()
            M = M[keep]
            good = np.isfinite(M).all(1) & M.ge(0).all(1)
            M = M.loc[good]
            assert M.index.is_unique
            x = np.log10(tr.loc[keep, "accepted_mass_g"].to_numpy(float))
            data[panel] = dict(M=M, C=C, x=x, species=keep)
            np.savez_compressed(
                O / "results" / f"{panel}_design.npz", C=C, mass_log10=x, species=keep
            )
            for i, s in enumerate(keep):
                speciesrows.append(
                    dict(
                        panel=panel,
                        species=s,
                        scientific_name=tr.loc[s, "scientific_name"],
                        mass_g=10 ** x[i],
                        sum_available=float(M[s].sum()),
                        n_measured=len(M),
                    )
                )
        pd.DataFrame(speciesrows).to_csv(O / "results/species.tsv", sep="\t", index=False)
        common = set.intersection(*(set(d["M"].index) for d in data.values()))

        def fit(Y, x, C):
            Q = np.linalg.inv(C)
            one = np.ones(len(x))
            P = Q - np.outer(Q @ one, Q @ one) / (one @ Q @ one)
            xx = x @ P @ x
            xy = Y @ P @ x
            yy = np.einsum("ij,ij->i", Y @ P, Y)
            r = np.clip(xy / np.sqrt(np.maximum(xx * yy, 1e-300)), -1 + 1e-14, 1 - 1e-14)
            p = 2 * stats.t.sf(np.abs(r) * np.sqrt((len(x) - 2) / (1 - r * r)), len(x) - 2)
            return (xy / xx, p, r, P)

        allrows = []
        scoreblocks = {}
        cover = []
        nullrows = []
        rng = np.random.default_rng(plan["seed"])
        for mode in plan["composition_modes"]:
            for panel, d in data.items():
                M = d["M"] if mode == "panel_available" else d["M"].loc[sorted(common)]
                gindex = {g: i for i, g in enumerate(M.index)}
                members = [sorted(set(g) & set(gindex)) for g in catalog.genes]
                nn = np.array(list(map(len, members)))
                ids = np.flatnonzero(nn >= 5)
                rr = []
                cc = []
                for j in ids:
                    rr.extend([j] * nn[j])
                    cc.extend((gindex[g] for g in members[j]))
                A = sparse.csr_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(catalog), len(M)))
                sums = np.asarray(A @ M.to_numpy())
                positive = np.isfinite(sums).all(1) & (sums > 0).all(1)
                ids = ids[positive[ids]]
                Y = np.log10(sums[ids])
                base = catalog.iloc[ids].drop(columns="genes").copy()
                base["mode"] = mode
                base["panel"] = panel
                base["n_species"] = len(d["x"])
                base["n_measured"] = nn[ids]
                base["coverage"] = nn[ids] / base.n_canonical
                base["genes"] = [";".join(members[i]) for i in ids]
                beta, p, r, P = fit(Y, d["x"], d["C"])
                base["beta"] = beta
                base["p"] = p
                base["r_phylo"] = r
                base["q_panel"] = bh(p)
                b0, p0, r0, _ = fit(Y, d["x"], np.eye(len(d["x"])))
                base["beta_OLS"] = b0
                base["p_Pearson"] = p0
                base["r_Pearson"] = r0
                base["q_Pearson_panel"] = bh(p0)
                _, pdg, _, _ = fit(Y, d["x"], np.diag(np.diag(d["C"])))
                base["p_diagonal"] = pdg
                allrows.append(base)
                testedset = set(ids)
                cover.extend(
                    (
                        dict(
                            mode=mode,
                            panel=panel,
                            set_id=catalog.set_id[i],
                            n_canonical=int(catalog.n_canonical[i]),
                            n_measured=int(nn[i]),
                            status=(
                                "tested" if i in testedset else "fewer_than_5_measured_or_zero_sum"
                            ),
                        )
                        for i in range(len(catalog))
                    )
                )
                scoreblocks[mode, panel] = (ids, Y, members)
                np.savez_compressed(
                    O / "results" / f"{mode}_{panel}_scores.npz", catalog_index=ids, log10_sum=Y
                )
                if mode == "panel_available":
                    xs = rng.multivariate_normal(np.zeros(len(d["x"])), d["C"], size=plan["B"])
                    xx = np.einsum("ij,ij->i", xs @ P, xs)
                    yy = np.einsum("ij,ij->i", Y @ P, Y)
                    nr = np.clip(Y @ P @ xs.T / np.sqrt(yy[:, None] * xx), -1 + 1e-14, 1 - 1e-14)
                    ps = 2 * stats.t.sf(
                        np.abs(nr) * np.sqrt((len(d["x"]) - 2) / (1 - nr * nr)), len(d["x"]) - 2
                    )
                    nom = (ps < 0.05).sum(0)
                    nbh = np.array([(bh(ps[:, i]) < 0.05).sum() for i in range(plan["B"])])
                    np.savez_compressed(
                        O / "results" / f"{panel}_null.npz", nominal_hits=nom, BH_hits=nbh
                    )
                    nullrows.append(
                        dict(
                            panel=panel,
                            observed_nominal=int((p < 0.05).sum()),
                            null_nominal_median=float(np.median(nom)),
                            nominal_count_p=float(
                                (1 + (nom >= (p < 0.05).sum()).sum()) / (1 + plan["B"])
                            ),
                            observed_BH=int((bh(p) < 0.05).sum()),
                            null_BH_median=float(np.median(nbh)),
                            BH_count_p=float(
                                (1 + (nbh >= (bh(p) < 0.05).sum()).sum()) / (1 + plan["B"])
                            ),
                        )
                    )
                print(
                    mode,
                    panel,
                    "genes",
                    len(M),
                    "maps",
                    len(ids),
                    "PGLS_BH",
                    int((base.q_panel < 0.05).sum()),
                    "Pearson_BH",
                    int((base.q_Pearson_panel < 0.05).sum()),
                    flush=True,
                )
        out = pd.concat(allrows, ignore_index=True)
        out["q_all_tests"] = 1.0
        for mode, z in out.groupby("mode"):
            out.loc[z.index, "q_all_tests"] = bh(z.p)
        out.to_csv(O / "results/map_tests.tsv.gz", sep="\t", index=False)
        pd.DataFrame(cover).to_csv(O / "results/coverage.tsv.gz", sep="\t", index=False)
        pd.DataFrame(nullrows).to_csv(O / "results/null_summary.tsv", sep="\t", index=False)
        reps = []
        for mode, z in out.groupby("mode"):
            for pp in [panels] + list(map(list, itertools.combinations(panels, 2))):
                pvt = (
                    z[z.panel.isin(pp)].pivot(index="set_id", columns="panel", values="p").dropna()
                )
                bet = (
                    z[z.panel.isin(pp)]
                    .pivot(index="set_id", columns="panel", values="beta")
                    .loc[pvt.index, pp]
                )
                same = bet.gt(0).all(1) | bet.lt(0).all(1)
                pmax = pvt.max(1).where(same, 1.0)
                r = pd.DataFrame(
                    dict(
                        set_id=pvt.index,
                        mode=mode,
                        panels="+".join(pp),
                        same_sign=same.to_numpy(),
                        p_rep=pmax.to_numpy(),
                    )
                )
                r["q_rep"] = bh(r.p_rep)
                for panel in pp:
                    r[panel + "_p"] = pvt[panel].to_numpy()
                    r[panel + "_beta"] = bet[panel].to_numpy()
                reps.append(r)
        rep = pd.concat(reps, ignore_index=True)
        rep["q_all_pairs"] = np.nan
        for mode, z in rep[rep.panels.str.count("\\+").eq(1)].groupby("mode"):
            rep.loc[z.index, "q_all_pairs"] = bh(z.p_rep)
        rep.to_csv(O / "results/replication.tsv.gz", sep="\t", index=False)
        diagids = set(out.loc[out.q_all_tests.lt(0.05), "set_id"]) | set(
            rep.loc[rep.p_rep.lt(0.05), "set_id"]
        )
        diagnostics = []
        driver = []
        for mode in plan["composition_modes"]:
            for panel, d in data.items():
                ids, Y, members = scoreblocks[mode, panel]
                lookup = {catalog.set_id[i]: j for j, i in enumerate(ids)}
                for sid in sorted(diagids & set(lookup)):
                    j = lookup[sid]
                    i = ids[j]
                    vals = []
                    ps = []
                    for k, s in enumerate(d["species"]):
                        mask = np.arange(len(d["x"])) != k
                        b, p, _, _ = fit(
                            Y[j : j + 1, mask], d["x"][mask], d["C"][np.ix_(mask, mask)]
                        )
                        vals.append(float(b[0]))
                        ps.append(float(p[0]))
                    original = out[
                        (out["mode"] == mode) & out.panel.eq(panel) & out.set_id.eq(sid)
                    ].iloc[0]
                    diagnostics.append(
                        dict(
                            mode=mode,
                            panel=panel,
                            set_id=sid,
                            loo_sign_fraction=float(
                                np.mean(np.sign(vals) == np.sign(original.beta))
                            ),
                            loo_max_p=max(ps),
                            loo_min_beta=min(vals),
                            loo_max_beta=max(vals),
                        )
                    )
                    gs = members[i]
                    gm = d["M"].loc[gs]
                    share = gm.div(gm.sum(0), axis=1).mean(1).sort_values(ascending=False)
                    driver.append(
                        dict(
                            mode=mode,
                            panel=panel,
                            set_id=sid,
                            top_gene=share.index[0],
                            top_gene_mean_share=float(share.iloc[0]),
                            top3_mean_share=float(share.iloc[:3].sum()),
                        )
                    )
        pd.DataFrame(diagnostics).to_csv(
            O / "results/shortlist_leave_one_species.tsv.gz", sep="\t", index=False
        )
        pd.DataFrame(driver).to_csv(
            O / "results/shortlist_gene_dominance.tsv.gz", sep="\t", index=False
        )
        summary = dict(
            catalog_counts=catalog.source.value_counts().to_dict(),
            catalog_total=len(catalog),
            common_genes=len(common),
            panels=[],
            replication=[],
            limits=[
                "Relative normalized RNA sums, not absolute RNA per cell or flux.",
                "Human-reference catalog symbols use inherited mapping; not a new native orthology audit.",
                "Partial measured map composition differs between panels in panel_available; all3_common fixes composition.",
                "Small species counts, different source contexts and pathway overlap remain.",
                "Post hoc whole-catalog screen, not independent validation of previously selected own5.",
            ],
            human_review="not_completed",
            training_enabled=False,
        )
        for (mode, panel), z in out.groupby(["mode", "panel"]):
            summary["panels"].append(
                dict(
                    mode=mode,
                    panel=panel,
                    n_tests=len(z),
                    nominal=int(z.p.lt(0.05).sum()),
                    BH_panel=int(z.q_panel.lt(0.05).sum()),
                    BH_all_tests=int(z.q_all_tests.lt(0.05).sum()),
                    q_min=float(z.q_panel.min()),
                    Pearson_BH=int(z.q_Pearson_panel.lt(0.05).sum()),
                    coverage_median=float(z.coverage.median()),
                    maps_ge50pct=int(z.coverage.ge(0.5).sum()),
                )
            )
        for (mode, pp), z in rep.groupby(["mode", "panels"]):
            summary["replication"].append(
                dict(
                    mode=mode,
                    panels=pp,
                    n_terms=len(z),
                    nominal=int(z.p_rep.lt(0.05).sum()),
                    BH=int(z.q_rep.lt(0.05).sum()),
                    q_min=float(z.q_rep.min()),
                )
            )
        (O / "results/summary.json").write_text(json.dumps(summary, indent=2))
        (O / "audit/input_manifest.json").write_text(json.dumps(manifest, indent=2))
        (O / "audit/environment.json").write_text(
            json.dumps(
                dict(
                    python=sys.version,
                    numpy=np.__version__,
                    pandas=pd.__version__,
                    scipy=scipy.__version__,
                    platform=platform.platform(),
                ),
                indent=2,
            )
        )
        print(json.dumps(summary, indent=2), flush=True)

    return SimpleNamespace(**locals())
