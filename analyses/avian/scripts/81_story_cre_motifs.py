#!/usr/bin/env python3
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
"""
from __future__ import annotations

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

sys.path.insert(0, str(Path(__file__).resolve().parent))
s74 = importlib.import_module("74_nineteen_via_cage")
s75 = importlib.import_module("75_honest_gene_mass_scan")
from oxphos_thread_lib import DATA, FIG, FIMO, RES, SPECIES, SPECIES_ORDER  # noqa: E402

OUT = RES / "story_cre_motifs"
MEME = RES / "module_consistency" / "fimo_chicken" / "JASPAR2024_vertebrates.meme"
PROX_BP = 2000
DISTAL_BP = 100_000
MIN_WIN = 400
FIMO_THRESH = "1e-4"
SETS = [
    ("oxphos_16", "OXPHOS 16"),
    ("folate_3", "folate 3"),
    ("ribo_5", "ribo 5"),
    ("green_hubs_7", "green hubs 7"),
    ("green_sig_51", "green 51"),
]
# Pre-specified names: his PDF + OXPHOS drivers we already used in cis.
WATCH = [
    "SNAI3", "PHOX2B", "ZEB1", "SOX14", "MEF2B", "HEY2",
    "NRF1", "GABPA", "ESRRA", "YY1",
]


def load_py(stem: str, fname: str):
    spec = importlib.util.spec_from_file_location(stem, Path(__file__).resolve().parent / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def genes_of(stem: str) -> list[str]:
    return [ln.strip() for ln in (DATA / "genesets" / f"story_{stem}.txt").read_text().splitlines() if ln.strip()]


def assign_peaks(peaks: pd.DataFrame, loci: pd.DataFrame, max_bp: int) -> pd.DataFrame:
    if loci.empty or peaks.empty:
        return pd.DataFrame()
    index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
        g = g.sort_values("five")
        index[(seq, st)] = (g["five"].to_numpy(float), g["symbol"].tolist())
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
        if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-6:
            continue
        rows.append(dict(
            gene=syms[jbest],
            peak_id=str(r.id),
            seq=str(r.seq),
            start=int(r.start),
            end=int(r.end),
            strand=str(r.strand),
            fib=float(r.fib),
            epi=float(r.epi),
            dist=float(d),
        ))
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
            lo, hi = int(r.start), int(r.end)
            if hi < lo:
                lo, hi = hi, lo
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
    cmd = [str(FIMO), "--text", "--no-qvalue", "--verbosity", "1",
           "--thresh", FIMO_THRESH, str(MEME), str(fasta)]
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
    x, y, names = x[ok], y[ok], [n for n, k in zip(names, ok) if k]
    if len(x) < 5:
        ax.set_visible(False)
        return {"r": np.nan, "p": np.nan, "n": int(len(x))}
    fit = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, fit.intercept + fit.slope * xs, color=s40.LINE, lw=1.6, zorder=2)
    ax.scatter(x, y, s=22, c=s40.POINT, zorder=3, linewidths=0)
    texts = [ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
             for name, xi, yi in zip(names, x, y)]
    s40.adjust_text(
        texts, ax=ax, x=x, y=y, expand=(1.15, 1.25),
        force_text=(0.25, 0.4), force_static=(0.15, 0.2),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    s40.theme_light(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
    ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(len(names)),
            "slope": float(fit.slope), "intercept": float(fit.intercept)}


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
            print(f"  пиков proximal {n_pr}  distal {n_ds}  генов {a.gene.nunique()}", flush=True)
        peaks = pd.concat(chunks, ignore_index=True)
        peaks.to_csv(peak_path, sep="\t", index=False)

    # proximal sums must match the locked fibroblast matrix
    old = pd.read_csv(RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv", sep="\t", index_col=0)
    prox = peaks[peaks.klass == "proximal"]
    rels = []
    for g in all_story:
        for sp in SPECIES_ORDER:
            got = float(prox.loc[(prox.gene == g) & (prox.species == sp), "fib"].sum())
            exp = float(old.loc[g, sp]) if g in old.index else np.nan
            if exp and np.isfinite(exp) and exp > 0:
                rels.append(abs(got - exp) / exp)
    print(f"сверка proximal vs M_fib: max отн. {max(rels) if rels else np.nan:.4g}  n={len(rels)}")

    # CAGE sums proximal / distal, fib / epi
    recs = []
    for stem, title in SETS:
        genes = story[stem]
        for klass_i in ("proximal", "distal"):
            sub = peaks[(peaks.klass == klass_i) & peaks.gene.isin(genes)]
            for tissue, col in (("fib", "fib"), ("epi", "epi")):
                sums, npeaks = [], []
                for sp in SPECIES_ORDER:
                    sl = sub[sub.species == sp]
                    sums.append(float(sl[col].sum()))
                    npeaks.append(int(sl.peak_id.nunique()))
                y = np.array([np.log10(s) if s > 0 else np.nan for s in sums], float)
                ok = np.isfinite(y)
                if ok.sum() >= 5:
                    fit = stats.linregress(logm[ok], y[ok])
                    rr, pp, sl, ice = fit.rvalue, fit.pvalue, fit.slope, fit.intercept
                else:
                    rr = pp = sl = ice = np.nan
                recs.append(dict(
                    set=stem, title=title, klass=klass_i, tissue=tissue,
                    n_genes=len(genes), n_species=int(ok.sum()),
                    r=rr, p=pp, slope=sl, intercept=ice,
                    **{f"sum_{sp.lower()}": sums[i] for i, sp in enumerate(SPECIES_ORDER)},
                    **{f"npeaks_{sp.lower()}": npeaks[i] for i, sp in enumerate(SPECIES_ORDER)},
                ))
                print(f"{stem:16s} {klass_i:9s} {tissue}  r={rr:+.3f} p={pp:.3f}  "
                      f"peaks={npeaks} sums={[round(s,1) for s in sums]}", flush=True)
    cage = pd.DataFrame(recs)
    cage.to_csv(OUT / "cage_prox_distal_fits.tsv", sep="\t", index=False)

    # figures: distal fib (new object) and two-up prox|distal fib
    distal_panels = []
    fig, axes = plt.subplots(5, 2, figsize=(8.4, 16.2))
    for i, (stem, title) in enumerate(SETS):
        for j, klass_i in enumerate(("proximal", "distal")):
            row = cage[(cage.set == stem) & (cage.klass == klass_i) & (cage.tissue == "fib")].iloc[0]
            sums = [float(row[f"sum_{sp.lower()}"]) for sp in SPECIES_ORDER]
            y = np.array([np.log10(s) if s > 0 else np.nan for s in sums])
            ylab = "log10_fib_proximal" if klass_i == "proximal" else "log10_fib_distal"
            fit = scatter(axes[i, j], names, logm, y, "log10_body_mass", ylab, s40)
            ptxt = "NA" if not np.isfinite(fit["p"]) else (
                f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}")
            rtxt = "NA" if not np.isfinite(fit["r"]) else f"{fit['r']:.3f}"
            axes[i, j].set_title(f"{title} {klass_i}   cor. {rtxt}; p = {ptxt}", fontsize=9, pad=6)
            if klass_i == "distal":
                distal_panels.append((title, names, logm, y, fit))
    fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.08, right=0.99, top=0.98, bottom=0.03)
    save(fig, "fig_story_prox_distal_fib")
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
    axes = axes.ravel()
    for i, (title, nm, x, y, fit) in enumerate(distal_panels):
        scatter(axes[i], nm, x, y, "log10_body_mass", "log10_fib_distal", s40)
        ptxt = "NA" if not np.isfinite(fit["p"]) else (
            f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}")
        rtxt = "NA" if not np.isfinite(fit["r"]) else f"{fit['r']:.3f}"
        axes[i].set_title(f"{title}   cor. {rtxt}; p = {ptxt}", fontsize=10, pad=7)
    axes[5].set_visible(False)
    fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.06, right=0.99, top=0.94, bottom=0.08)
    save(fig, "fig_story_five_distal_fib")
    plt.close(fig)

    # FIMO on proximal peaks
    if not MEME.exists() or not FIMO.exists():
        raise SystemExit(f"need FIMO={FIMO.exists()} MEME={MEME.exists()}")
    hit_rows = []
    fish_rows = []
    watch_rows = []
    pending = [
        sp for sp in SPECIES_ORDER
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
        sp_peaks["sid"] = (sp_peaks.gene.astype(str) + "|" + sp_peaks.peak_id.astype(str)
                           + f"|{sp}|proximal")
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
                rec.update(set=stem, title=title, species=sp, tf=tf,
                           n_pos=len(pos), n_neg=len(neg))
                fish_rows.append(rec)
                ps.append(rec["p"])
            # BH later per block
            print(f"  {stem} pos={len(pos)} neg={len(neg)} motifs={len(by_tf)}", flush=True)
        # watch TFs
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
                rec.update(set=stem, title=title, species=sp, tf=w,
                           n_pos=len(pos), n_neg=len(neg))
                watch_rows.append(rec)

    fish = pd.DataFrame(fish_rows)
    if not fish.empty:
        fish["bh"] = np.nan
        for (stem, sp), idx in fish.groupby(["set", "species"]).groups.items():
            p = fish.loc[idx, "p"].to_numpy(float)
            fish.loc[idx, "bh"] = stats.false_discovery_control(np.nan_to_num(p, nan=1.0), method="bh")
        fish.to_csv(OUT / "motif_fisher_proximal.tsv", sep="\t", index=False)
        sig = fish[fish.bh < 0.05]
        print("\nBH<0.05 proximal:")
        print(sig.groupby(["set", "species"]).size().to_string() if len(sig) else "  нет")
        # motifs significant in ≥3 species, per set
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
        # pre-specified: -log10 p vs mass, no top-N
        wsum = []
        for (stem, tf), g in watch.groupby(["set", "tf"]):
            g = g.set_index("species").reindex(SPECIES_ORDER)
            y = -np.log10(np.clip(g["p"].to_numpy(float), 1e-12, 1))
            if np.isfinite(y).sum() >= 5:
                fit = stats.linregress(logm[np.isfinite(y)], y[np.isfinite(y)])
                wsum.append(dict(set=stem, tf=tf, r=float(fit.rvalue), p=float(fit.pvalue),
                                 **{f"p_{sp.lower()}": float(g.loc[sp, "p"]) for sp in SPECIES_ORDER},
                                 **{f"odds_{sp.lower()}": float(g.loc[sp, "odds"]) for sp in SPECIES_ORDER}))
        pd.DataFrame(wsum).to_csv(OUT / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False)
        print("\nwatch TF, −log10 p vs mass:")
        ws = pd.DataFrame(wsum).sort_values("p")
        print(ws[["set", "tf", "r", "p"]].head(20).to_string(index=False))

    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
