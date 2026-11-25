"""Avian CAGE measurement, assignments and expression matrices.

Calculation bodies are retained from the source commit; configuration and execution scope are explicit.
"""

from __future__ import annotations
from types import SimpleNamespace


def _three_checks(context):
    """The three checks the recovered result needs before it can be stated.

    A. Matched null. Is the level decline of the recovered core stronger than that of
       random gene sets matched on post-recovery detection breadth and on chicken
       abundance decile? Without this, -0.880 has nothing to be compared against.

    B. Symmetric measurement. In script 51 the three Augustus species got their genes
       measured by our window rule while the four native species kept the old
       annotation-based rule. Same gene, two procedures, and the procedures differ
       exactly along the axis under test. Here every gene in every species is measured
       by ONE rule: take the gene's 5' end (from the species' own annotation where it
       exists, from the miniprot locus where it does not), collect peaks within
       1000 bp on the gene's strand, threshold each peak at 1 TPM per library, sum,
       average over the tissue's libraries.

    C. The promoter-class filter. Recovered genes bypassed it because that class was
       assigned by the annotation being corrected. Reported both ways: any class, and
       promoter class only, under the symmetric rule.

    Source: analyses/avian/scripts/53_three_checks.py at d230e56ab942."""
    import re
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    RES = context.unit("avian.oxphos_thread_lib").RES
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    load_geneset = context.unit("avian.oxphos_thread_lib").load_geneset
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    EXP = Path(
        str(
            context.resolve_path(
                "/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab"
            )
        )
    )
    GFFDIR = RES / "gene_recovery"
    PROTSRC = Path(
        str(
            context.resolve_path(
                "/path/to/storage/kleiber_birds_gemcat/ref/chicken_prot/query_proteins_all.faa"
            )
        )
    )
    OUT = RES / "three_checks"
    AUG = ["Turkey", "Peacock", "Emu"]
    LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))
    WINDOW = 1000
    THR = 1.0
    N_NULL = 2000
    RNG = np.random.default_rng(20260920)
    LIST_NAME = "его OXPHOS 133"

    def fit(y) -> tuple[float, float]:
        y = np.asarray(y, float)
        ok = np.isfinite(y)
        if ok.sum() < 3:
            return (np.nan, np.nan)
        r = stats.linregress(LOGM[ok], y[ok])
        return (float(r.rvalue), float(r.pvalue))

    def five_prime_table(
        species: str, peaks: pd.DataFrame, use_miniprot: bool = True
    ) -> pd.DataFrame:
        """symbol -> (seqname, 5' position, strand), one row per symbol.

        Annotation first: every peak row carries the coordinates of the model it was
        assigned to. Then the miniprot loci, for the genes the annotation lacks.
        """
        rows = []
        d = peaks.dropna(subset=["geneStart", "geneEnd"]).copy()
        d = d[d["symbol"].ne("")]
        d["strand_g"] = (
            d["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(d["strand"])
        )
        for sym, g in d.groupby("symbol", observed=True):
            r = g.iloc[0]
            st = r["strand_g"]
            rows.append(
                dict(
                    symbol=sym,
                    seqname=r["seqnames"],
                    strand=st,
                    five=float(r["geneStart"] if st == "+" else r["geneEnd"]),
                    src="annotation",
                )
            )
        have = {r["symbol"] for r in rows}
        gff = GFFDIR / f"{species}_miniprot.gff"
        if use_miniprot and gff.exists():
            pid2sym = dict(
                re.findall("^>(\\S+)\\s+(\\S+)\\s+\\[organism=", PROTSRC.read_text(), re.M)
            )
            best: dict[str, tuple[float, dict]] = {}
            for line in gff.read_text().splitlines():
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 9 or f[2] != "mRNA":
                    continue
                at = dict((kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv))
                sym = pid2sym.get(at.get("Target", "").split()[0])
                if not sym or sym in have:
                    continue
                seq = f[0].split()[0]
                if seq.startswith("ENA|"):
                    seq = seq.split("|")[-1]
                sc = float(f[5]) if f[5] != "." else 0.0
                st = f[6]
                cand = dict(
                    symbol=sym,
                    seqname=seq,
                    strand=st,
                    five=float(f[3] if st == "+" else f[4]),
                    src="miniprot",
                )
                if sym not in best or sc > best[sym][0]:
                    best[sym] = (sc, cand)
            rows.extend((v[1] for v in best.values()))
        return pd.DataFrame(rows)

    def measure(
        species: str, promoter_only: bool, use_miniprot: bool = True
    ) -> dict[str, pd.Series]:
        """One rule for every species: window around the 5' end."""
        d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        fp = five_prime_table(species, d, use_miniprot)
        cols = {t: lib_cols(d, t) for t in ("fib", "epi")}
        pool = d[d["class"].astype(str).str.lower().eq("promoter")] if promoter_only else d
        pool = pool.assign(mid=(pool["start"] + pool["end"]) / 2.0)
        by_seq = {k: v for k, v in pool.groupby("seqnames", observed=True)}
        out = {t: {} for t in ("fib", "epi")}
        for _, g in fp.iterrows():
            sub = by_seq.get(g["seqname"])
            if sub is None:
                continue
            near = sub[(sub["strand"] == g["strand"]) & ((sub["mid"] - g["five"]).abs() <= WINDOW)]
            if near.empty:
                continue
            for t in ("fib", "epi"):
                c = cols[t]
                if not c:
                    continue
                m = near[c].apply(pd.to_numeric, errors="coerce").fillna(0.0)
                m = m.where(m >= THR, 0.0)
                v = float(m.sum(axis=0).mean())
                if v > 0:
                    out[t][g["symbol"]] = v
        return {t: pd.Series(out[t], dtype=float) for t in ("fib", "epi")}

    def matched_null_p(wide: pd.DataFrame, core: list[str], obs_r: float) -> tuple[float, float]:
        det = wide.notna() & (wide > 0)
        breadth = det.sum(axis=1)
        chick = wide["Chicken"]
        pool = wide.index[breadth > 0]
        dec = pd.qcut(chick.reindex(pool).rank(method="first"), 10, labels=False, duplicates="drop")
        key = pd.DataFrame({"b": breadth.reindex(pool), "d": dec}, index=pool)
        hit = [g for g in core if g in key.index]
        if len(hit) < 5 or not np.isfinite(obs_r):
            return (np.nan, np.nan)
        want = key.loc[hit].groupby(["b", "d"]).size()
        buckets = {
            k: key.index[(key["b"] == k[0]) & (key["d"] == k[1])].to_numpy() for k in want.index
        }
        null = np.empty(N_NULL)
        for i in range(N_NULL):
            pick: list[str] = []
            for k, need in want.items():
                av = buckets[k]
                if len(av):
                    pick.extend(RNG.choice(av, size=min(need, len(av)), replace=False))
            null[i] = fit(np.log10(wide.loc[pick].mean(axis=0).to_numpy(float)))[0]
        null = null[np.isfinite(null)]
        return (float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null)))

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        genes = sorted({normalize_symbol(g) for g in load_geneset("ruslan_oxphos_133")} - {""})
        rows = []
        for promoter_only in (False, True):
            tag = "только promoter" if promoter_only else "любой класс"
            mats = {sp: measure(sp, promoter_only) for sp in SPECIES_ORDER}
            wide = {
                t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")
            }
            idx = sorted(set(wide["fib"].index) | set(wide["epi"].index))
            ok = None
            for t in ("fib", "epi"):
                m = wide[t].reindex(idx)
                g = m.notna() & (m > 0)
                ok = g if ok is None else ok & g
            breadth = ok.sum(axis=1)
            for k in (7, 6, 5, 4, 3):
                core = [g for g in genes if breadth.get(g, 0) >= k]
                for t in ("fib", "epi"):
                    w = wide[t]
                    r, p = fit(np.log10(w.reindex(core).mean(axis=0).to_numpy(float)))
                    pn, nm = (np.nan, np.nan)
                    if k in (7, 5):
                        pn, nm = matched_null_p(w, core, r)
                    rows.append(
                        dict(
                            правило=tag,
                            k=k,
                            n_core=len(core),
                            tissue=t,
                            r=round(r, 3),
                            p=round(p, 4),
                            p_matched_null=None if np.isnan(pn) else round(pn, 4),
                            null_median_r=None if np.isnan(nm) else round(nm, 3),
                            n_detected=";".join(
                                (str(int(v)) for v in w.reindex(core).notna().sum(axis=0))
                            ),
                        )
                    )
            print(
                f"[{tag}] генов измерено: "
                + ", ".join((f"{sp} {int((wide['fib'][sp] > 0).sum())}" for sp in SPECIES_ORDER))
            )
        tab = pd.DataFrame(rows)
        tab.to_csv(OUT / "three_checks.tsv", sep="\t", index=False)
        pd.set_option("display.width", 260)
        for tag in tab["правило"].unique():
            print("\n" + "=" * 118)
            print(f"ЕДИНОЕ ПРАВИЛО ДЛЯ ВСЕХ СЕМИ ВИДОВ — {tag}")
            print("=" * 118)
            print(tab[tab["правило"] == tag].drop(columns=["правило"]).to_string(index=False))
        print(f"\nwrote {OUT / 'three_checks.tsv'}")

    return SimpleNamespace(**locals())


def _working_sets_round2(context):
    """Working sets and fits after BOTH recovery rounds.

    Gene position priority: the species' own annotation, then the round-1 miniprot
    locus (chicken protein, script 50), then the round-2 locus (four-bird proteins,
    script 56). Measurement rule is the same for all seven species: peaks within
    WINDOW bp of the 5' end on the gene's strand, 1 TPM per library, summed, averaged
    over the tissue's libraries.

    Working set = genes of the KEGG map detected in every one of the seven birds.

    Source: analyses/avian/scripts/57_working_sets_round2.py at d230e56ab942."""
    import argparse
    import importlib
    import re
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd

    s53 = context.unit("avian.53_three_checks")
    RES = context.unit("avian.oxphos_thread_lib").RES
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    load_geneset = context.unit("avian.oxphos_thread_lib").load_geneset
    load_geo_rnaseq = context.unit("avian.oxphos_thread_lib").load_geo_rnaseq
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    EXP = Path(
        str(
            context.resolve_path(
                "/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab"
            )
        )
    )
    R1 = RES / "gene_recovery"
    R2 = RES / "recovery_round2"
    PROTSRC = Path(
        str(
            context.resolve_path(
                "/path/to/storage/kleiber_birds_gemcat/ref/chicken_prot/query_proteins_all.faa"
            )
        )
    )
    OUT = RES / "working_sets_round2"
    SETS = {"hsa00190 OXPHOS": "kegg_00190", "hsa04714 Thermogenesis": "kegg_04714"}
    THR = 1.0

    def loci_from_gff(gff: Path, pid2sym: dict[str, str] | None) -> dict[str, dict]:
        best: dict[str, tuple[float, dict]] = {}
        if not gff.exists():
            return {}
        for line in gff.read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            at = dict((kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv))
            tgt = at.get("Target", "").split()[0]
            sym = tgt.split("##")[1] if "##" in tgt else (pid2sym or {}).get(tgt)
            if not sym:
                continue
            seq = f[0].split()[0]
            if seq.startswith("ENA|"):
                seq = seq.split("|")[-1]
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            cand = dict(seqname=seq, strand=st, five=float(f[3] if st == "+" else f[4]))
            if sym not in best or sc > best[sym][0]:
                best[sym] = (sc, cand)
        return {k: v[1] for k, v in best.items()}

    def measure(species: str, window: int) -> dict[str, pd.Series]:
        d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        fp: dict[str, dict] = {}
        ann = d.dropna(subset=["geneStart", "geneEnd"])
        ann = ann[ann["symbol"].ne("")].copy()
        ann["sg"] = (
            ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
        )
        for sym, g in ann.groupby("symbol", observed=True):
            r = g.iloc[0]
            st = r["sg"]
            fp[sym] = dict(
                seqname=r["seqnames"],
                strand=st,
                five=float(r["geneStart"] if st == "+" else r["geneEnd"]),
            )
        pid2sym = dict(re.findall("^>(\\S+)\\s+(\\S+)\\s+\\[organism=", PROTSRC.read_text(), re.M))
        REANNO = Path(
            str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies"))
        )
        for gff, p2s in (
            (REANNO / species / "reanno.gff", None),
            (R1 / f"{species}_miniprot.gff", pid2sym),
            (R2 / f"{species}_q2.gff", None),
        ):
            for sym, c in loci_from_gff(gff, p2s).items():
                fp.setdefault(sym, c)
        cols = {t: lib_cols(d, t) for t in ("fib", "epi")}
        pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
        by = {k: v for k, v in pool.groupby("seqnames", observed=True)}
        out = {t: {} for t in ("fib", "epi")}
        for sym, c in fp.items():
            sub = by.get(c["seqname"])
            if sub is None:
                continue
            near = sub[(sub["strand"] == c["strand"]) & ((sub["mid"] - c["five"]).abs() <= window)]
            if near.empty:
                continue
            for t in ("fib", "epi"):
                cc = cols[t]
                if not cc:
                    continue
                m = near[cc].apply(pd.to_numeric, errors="coerce").fillna(0.0)
                m = m.where(m >= THR, 0.0)
                v = float(m.sum(axis=0).mean())
                if v > 0:
                    out[t][sym] = v
        return {t: pd.Series(out[t], dtype=float) for t in ("fib", "epi")}

    def main() -> None:
        ap = argparse.ArgumentParser()
        ap.add_argument("--window", type=int, default=1000)
        a = ap.parse_args(context.argv)
        OUT.mkdir(parents=True, exist_ok=True)
        mats = {sp: measure(sp, a.window) for sp in SPECIES_ORDER}
        cage = {t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")}
        idx = sorted(set(cage["fib"].index) | set(cage["epi"].index))
        ok = None
        for t in ("fib", "epi"):
            m = cage[t].reindex(idx)
            g = m.notna() & (m > 0)
            ok = g if ok is None else ok & g
        rna = pd.DataFrame(
            {
                sp: (
                    lambda g: g["gene_mat"].div(g["lib_total"], axis=1).mul(1000000.0).mean(axis=1)
                )(load_geo_rnaseq(sp))
                for sp in SPECIES_ORDER
            }
        )
        rows = []
        for name, key in SETS.items():
            genes = sorted({normalize_symbol(g) for g in load_geneset(key)} - {""})
            det = ok.reindex(genes).fillna(False).infer_objects(copy=False)
            per = det.sum(axis=0)
            core = sorted(det.index[det.all(axis=1)])
            (OUT / f"working_{key}_w{a.window}.txt").write_text("\n".join(core) + "\n")
            print(
                f"\n{name}: карта {len(genes)}, найдено по видам "
                + ", ".join((f"{sp} {int(per[sp])}" for sp in SPECIES_ORDER))
                + f"; РАБОЧИЙ НАБОР {len(core)}"
            )
            for layer, w in (
                ("CAGE фибробласты", cage["fib"]),
                ("CAGE эпибласт", cage["epi"]),
                ("RNA-seq фибробласты", rna),
            ):
                present = [g for g in core if g in w.index]
                r, p = s53.fit(np.log10(w.reindex(present).mean(axis=0).to_numpy(float)))
                pn, nm = s53.matched_null_p(w, present, r)
                rows.append(
                    dict(
                        набор=name,
                        окно=a.window,
                        n_рабочий=len(core),
                        слой=layer,
                        r=round(r, 3),
                        p=round(p, 4),
                        p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                        медиана_r_нуля=None if np.isnan(nm) else round(nm, 3),
                    )
                )
                print(
                    f"  {layer:22s} r={r:+.3f} p={p:.4f}  против нуля {('—' if np.isnan(pn) else round(pn, 4))}"
                )
        t = pd.DataFrame(rows)
        f = OUT / f"fits_w{a.window}.tsv"
        t.to_csv(f, sep="\t", index=False)
        print(f"\nwrote {f}")

    return SimpleNamespace(**locals())


def _build_library_matrices(context):
    """Gene x library matrices for WGCNA, plus the widened gene universe.

    Everything in this project so far collapsed libraries to a per-tissue mean. WGCNA
    needs the libraries themselves, so the same measurement rule is applied per
    library here:

        gene 5' end  <- shipped annotation, then reanno.gff (whole-proteome
                        re-annotation, script 58), then the two targeted recovery
                        rounds (scripts 50 and 56)
        value(gene, library) = sum over peaks within WINDOW bp of that 5' end on the
                        gene's strand, each peak first thresholded at 1 TPM
        stored as log10(value + 1)   -- a log, not a VST: species are compared here

    Universe, widened as planned: the old 2828 required non-zero signal in BOTH
    tissues in all seven species. Here a gene enters the fibroblast universe if it is
    non-zero in at least one fibroblast library in every species, and likewise for
    epiblast; the combined universe is the union.

    Networks are built per tissue across species, which is what removes the
    tissue-contrast bias without any correction:
        M_fib  42 libraries   M_epi  23 libraries   M_all  65 libraries
    plus one matrix per species (7-11 libraries) as a sensitivity.

    Source: analyses/avian/scripts/62_build_library_matrices.py at d230e56ab942."""
    import argparse
    import importlib
    import re
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd

    s57 = context.unit("avian.57_working_sets_round2")
    RES = context.unit("avian.oxphos_thread_lib").RES
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    EXP = Path(
        str(
            context.resolve_path(
                "/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab"
            )
        )
    )
    OUT = RES / "wgcna_input"
    THR = 1.0

    def five_prime(species: str, d: pd.DataFrame) -> dict[str, dict]:
        """Same locus priority as script 57, exposed here so both stay in step."""
        fp: dict[str, dict] = {}
        ann = d.dropna(subset=["geneStart", "geneEnd"])
        ann = ann[ann["symbol"].ne("")].copy()
        ann["sg"] = (
            ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
        )
        for sym, g in ann.groupby("symbol", observed=True):
            r = g.iloc[0]
            st = r["sg"]
            fp[sym] = dict(
                seqname=r["seqnames"],
                strand=st,
                five=float(r["geneStart"] if st == "+" else r["geneEnd"]),
            )
        pid2sym = dict(
            re.findall("^>(\\S+)\\s+(\\S+)\\s+\\[organism=", s57.PROTSRC.read_text(), re.M)
        )
        reanno = Path(
            str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies"))
        )
        for gff, p2s in (
            (reanno / species / "reanno.gff", None),
            (s57.R1 / f"{species}_miniprot.gff", pid2sym),
            (s57.R2 / f"{species}_q2.gff", None),
        ):
            for sym, c in s57.loci_from_gff(gff, p2s).items():
                fp.setdefault(sym, c)
        return fp

    def per_library(species: str, window: int) -> pd.DataFrame:
        d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        fp = five_prime(species, d)
        cols = lib_cols(d, "fib") + lib_cols(d, "epi")
        pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
        num = pool[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        num = num.where(num >= THR, 0.0)
        pool = pool.assign(**{c: num[c] for c in cols})
        by = {k: v for k, v in pool.groupby("seqnames", observed=True)}
        rows: dict[str, np.ndarray] = {}
        for sym, c in fp.items():
            sub = by.get(c["seqname"])
            if sub is None:
                continue
            near = sub[(sub["strand"] == c["strand"]) & ((sub["mid"] - c["five"]).abs() <= window)]
            if near.empty:
                continue
            v = near[cols].to_numpy(float).sum(axis=0)
            if v.sum() <= 0:
                continue
            rows[sym] = v
        m = pd.DataFrame.from_dict(rows, orient="index", columns=cols)
        m.columns = [
            f"{species}|{('fib' if c.startswith('norm.C') else 'epi')}|{c.split('.')[1]}"
            for c in m.columns
        ]
        return np.log10(m + 1.0)

    def main() -> None:
        ap = argparse.ArgumentParser()
        ap.add_argument("--window", type=int, default=1000)
        a = ap.parse_args(context.argv)
        OUT.mkdir(parents=True, exist_ok=True)
        per: dict[str, pd.DataFrame] = {}
        for sp in SPECIES_ORDER:
            per[sp] = per_library(sp, a.window)
            f = OUT / f"M_{sp}_w{a.window}.tsv"
            per[sp].to_csv(f, sep="\t")
            nf = sum((1 for c in per[sp].columns if "|fib|" in c))
            ne = sum((1 for c in per[sp].columns if "|epi|" in c))
            print(
                f"[{sp:8s}] генов {len(per[sp]):>6}  библиотек {per[sp].shape[1]} (фибр {nf}, эпи {ne})"
            )

        def universe(tissue: str) -> list[str]:
            sets = []
            for sp in SPECIES_ORDER:
                c = [x for x in per[sp].columns if f"|{tissue}|" in x]
                m = per[sp][c]
                sets.append(set(m.index[(m > 0).any(axis=1)]))
            return sorted(set.intersection(*sets))

        uf, ue = (universe("fib"), universe("epi"))
        ua = sorted(set(uf) | set(ue))
        print(f"\nвселенная: фибробласты {len(uf)}, эпибласт {len(ue)}, объединение {len(ua)}")
        for nm, u in (("fib", uf), ("epi", ue), ("all", ua)):
            (OUT / f"universe_{nm}_w{a.window}.txt").write_text("\n".join(u) + "\n")

        def combine(tissue: str | None, genes: list[str]) -> pd.DataFrame:
            parts = []
            for sp in SPECIES_ORDER:
                c = [x for x in per[sp].columns if tissue is None or f"|{tissue}|" in x]
                parts.append(per[sp].reindex(genes)[c])
            return pd.concat(parts, axis=1)

        mats = {"fib": combine("fib", uf), "epi": combine("epi", ue), "all": combine(None, ua)}
        for nm, m in mats.items():
            m = m.dropna(how="all")
            m.to_csv(OUT / f"M_{nm}_w{a.window}.tsv", sep="\t")
            print(
                f"M_{nm}: {m.shape[0]} генов x {m.shape[1]} библиотек, пропусков {int(m.isna().sum().sum())}"
            )
        meta = []
        for c in mats["all"].columns:
            sp, tis, lib = c.split("|")
            meta.append(
                dict(
                    sample=c,
                    species=sp,
                    tissue=tis,
                    library=lib,
                    mass_g=SPECIES[sp]["mass_g"],
                    log10_mass=np.log10(SPECIES[sp]["mass_g"]),
                )
            )
        pd.DataFrame(meta).to_csv(OUT / f"samples_w{a.window}.tsv", sep="\t", index=False)
        print(f"\nwrote {OUT}/M_*_w{a.window}.tsv, universe_*, samples_w{a.window}.tsv")

    return SimpleNamespace(**locals())


def _rebuild_matrices_asymmetric(context):
    """Gene x library matrices rebuilt on a symmetric footing.

    Two things are changed against script 62, and both were forced by measurement.

    Window. Script 62 took peaks within 1000 bp of a gene's 5' end in either
    direction. Among genes detected in exactly four of the seven species, the three
    missing are Peacock, Turkey and Emu in 1036, 1016 and 1030 cases out of 1073 —
    the loss runs along the axis under test. The cause is that their annotations mark
    where translation starts while CAGE marks where transcription starts: for genes
    whose nearest peak falls outside the old window the median offset is -1766 bp in
    Peacock, -2066 in Turkey and -2126 in Emu, against +1399 in Chicken with only 41%
    upstream. So the window is opened upstream — for all seven species alike, which
    is what check B of script 53 requires. Chicken has 627 genes outside the old
    window and gains little; Peacock has 2805 and gains most.

    A window reaching 5 kb upstream can cover a neighbouring gene's promoter, so a
    peak is no longer counted for every locus that covers it. Two rules, in order:

      1. If the peak table already names a gene that we have a locus for, the peak
         stays with that gene, on any contig. Distance to the miniprot CDS start is
         not consulted: RGS22's promoter sits 75 kb from its own CDS, and a 20 kb
         cap would refuse the name and also refuse to give the peak to COX6C, which
         is how emu COX6C went to zero. A peak named as a different gene is never
         given to a neighbour (COX6C/RGS22, SDHD/C24H11ORF57, NDUFB1/TRIP11).
      2. Unnamed peaks go to the nearest locus inside the window, provided no other
         gene's 5' lies between them.

    Loci. Script 58 projected a proteome onto three genomes out of seven and gave
    each a different query, which is why locus counts ran 10 953 in Chicken against
    18 618 in Emu. Script 70 replaces that with one combined query — one protein per
    gene from chicken, duck, quail, zebra, turkey and emu, 79 570 proteins over
    16 779 named symbols — projected onto all seven genomes. Symbols matter: the
    first pass kept LOC-prefixed placeholders, which do not agree between species and
    planted some 21 600 phantom loci per genome; with peaks going to the nearest
    locus those phantoms stole real promoters, and ATP5MF and UQCRC2 fell from
    detected in all seven species to detected in one. With LOC dropped the layer
    comes out even: 16 054 to 16 523 symbols per species, under 3% spread.

    Two ways of using it are offered and the choice is made on evidence, not taste:

      annot_first  the species' own annotation first, reanno_common filling the rest.
                   More accurate per species (real TSS where the annotation is good),
                   but the accuracy itself is unequal across species.
      common_only  reanno_common for everyone. Less accurate — every locus is a
                   coding start — but identical in kind for all seven.

    Universe is unchanged: a gene enters if non-zero in at least one library of that
    tissue in EVERY species. That requirement is not an artefact; it is the only
    threshold at which detection rate does not track body mass (r = 0.000 at k = 7
    against -0.670 at k >= 6 and -0.790 at k >= 4 under the old window).

    Source: analyses/avian/scripts/68_rebuild_matrices_asymmetric.py at d230e56ab942."""
    import argparse
    import importlib
    import re
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    s57 = context.unit("avian.57_working_sets_round2")
    s62 = context.unit("avian.62_build_library_matrices")
    RES = context.unit("avian.oxphos_thread_lib").RES
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    ASM = Path(
        str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies"))
    )
    OUT = RES / "wgcna_input_v2"
    THR = 1.0

    def loci_from_reanno(gff: Path) -> dict[str, dict]:
        """Best-scoring model per NORMALISED symbol.

        Normalising is not cosmetic here. The six source annotations disagree on the
        ATP synthase nomenclature, so the query carries ATP5C1 next to ATP5F1C,
        ATP5H next to ATP5PD and so on. Left as they are, one gene gets two loci at
        the same place which then compete for the same peak, and the gene drops out
        of most species: ATP5F1E and ATP5PF survived in one species each before this
        was applied. PREFERRED_ALIAS is the table the rest of the project already
        uses (oxphos_thread_lib.normalize_symbol).
        """
        best: dict[str, tuple[float, dict]] = {}
        if not gff.exists():
            return {}
        with gff.open() as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t")
                if len(f) < 9 or f[2] != "mRNA":
                    continue
                m = re.search("Target=(\\S+)", f[8])
                if not m or "##" not in m.group(1):
                    continue
                sym = normalize_symbol(m.group(1).split("##")[1])
                if not sym:
                    continue
                seq = f[0].split()[0]
                if seq.startswith("ENA|"):
                    seq = seq.split("|")[-1]
                sc = float(f[5]) if f[5] != "." else 0.0
                st = f[6]
                if sym not in best or sc > best[sym][0]:
                    best[sym] = (
                        sc,
                        dict(seqname=seq, strand=st, five=float(f[3] if st == "+" else f[4])),
                    )
        return {k: v[1] for k, v in best.items()}

    def loci_of(species: str, d: pd.DataFrame, mode: str) -> tuple[dict[str, dict], dict[str, int]]:
        """symbol -> 5' locus, plus how many each layer contributed."""
        fp: dict[str, dict] = {}
        src: dict[str, int] = {}
        if mode == "annot_first":
            ann = d.dropna(subset=["geneStart", "geneEnd"])
            ann = ann[ann["symbol"].ne("")].copy()
            ann["sg"] = (
                ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
            )
            for sym, g in ann.groupby("symbol", observed=True):
                r = g.iloc[0]
                st = r["sg"]
                fp[sym] = dict(
                    seqname=r["seqnames"],
                    strand=st,
                    five=float(r["geneStart"] if st == "+" else r["geneEnd"]),
                )
            src["аннотация"] = len(fp)
        common = loci_from_reanno(ASM / species / "reanno_named.gff")
        added = 0
        for sym, c in common.items():
            if sym not in fp:
                fp[sym] = c
                added += 1
        src["reanno_named"] = added
        if mode == "common_then_annot":
            ann = d.dropna(subset=["geneStart", "geneEnd"])
            ann = ann[ann["symbol"].ne("")].copy()
            ann["sg"] = (
                ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
            )
            back = 0
            for sym, g in ann.groupby("symbol", observed=True):
                if sym in fp:
                    continue
                r = g.iloc[0]
                st = r["sg"]
                fp[sym] = dict(
                    seqname=r["seqnames"],
                    strand=st,
                    five=float(r["geneStart"] if st == "+" else r["geneEnd"]),
                )
                back += 1
            src["аннотация добором"] = back
        return (fp, src)

    def _oriented_offset(mid: float, five: float, strand: str) -> float:
        d = mid - five
        return -d if strand == "-" else d

    def _table_fives(peaks: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
        """Every named gene's 5' from the peak table, used only as a jump barrier."""
        if peaks.empty or "geneStart" not in peaks.columns:
            return {}
        sg = peaks["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"})
        if "strand" in peaks.columns:
            sg = sg.fillna(peaks["strand"])
        rows: dict[tuple[str, str], list[float]] = {}
        gs = pd.to_numeric(peaks["geneStart"], errors="coerce").to_numpy(float)
        ge = (
            pd.to_numeric(peaks["geneEnd"], errors="coerce").to_numpy(float)
            if "geneEnd" in peaks.columns
            else np.full(len(peaks), np.nan)
        )
        syms = (
            peaks["symbol"].to_numpy(object)
            if "symbol" in peaks.columns
            else np.array([""] * len(peaks), dtype=object)
        )
        for seq, st, gsv, gev, sym in zip(
            peaks["seqnames"].to_numpy(object), sg.to_numpy(object), gs, ge, syms
        ):
            if not isinstance(sym, str) or not sym or st not in ("+", "-"):
                continue
            if not np.isfinite(gsv):
                continue
            five = float(gsv if st == "+" else gev if np.isfinite(gev) else gsv)
            rows.setdefault((str(seq), str(st)), []).append(five)
        return {k: np.unique(np.asarray(v, float)) for k, v in rows.items()}

    def _has_between(mid: float, five: float, others: np.ndarray) -> bool:
        if others.size == 0:
            return False
        lo, hi = (mid, five) if mid <= five else (five, mid)
        i0 = int(np.searchsorted(others, lo, side="right"))
        i1 = int(np.searchsorted(others, hi, side="left"))
        return i1 > i0

    def assign_peaks(
        fp: dict[str, dict], peaks: pd.DataFrame, up: int, down: int
    ) -> tuple[pd.Series, dict[str, int]]:
        """Named peak stays with its gene, regardless of the miniprot CDS start.

        A peak named as some other gene is never given to a neighbour. Remaining
        unnamed peaks go to the nearest locus in [-up, +down] that does not skip
        over another gene's 5'.
        """
        owner = pd.Series(index=peaks.index, dtype=object)
        locked = np.zeros(len(peaks), dtype=bool)
        n_named = n_locked = 0
        named_i: list = []
        named_s: list[str] = []
        if "symbol" in peaks.columns and (not peaks.empty):
            for k, (i, seq, strand, mid, raw) in enumerate(
                zip(
                    peaks.index.to_numpy(),
                    peaks["seqnames"].to_numpy(object),
                    peaks["strand"].to_numpy(object),
                    peaks["mid"].to_numpy(float),
                    peaks["symbol"].to_numpy(object),
                )
            ):
                if not isinstance(raw, str) or not raw:
                    continue
                sym = normalize_symbol(raw)
                if not sym:
                    continue
                c = fp.get(sym)
                if c is not None:
                    named_i.append(i)
                    named_s.append(sym)
                    locked[k] = True
                    n_named += 1
                else:
                    locked[k] = True
                    n_locked += 1
        if named_i:
            owner.loc[named_i] = named_s
        free = peaks.loc[~locked]
        n_near = 0
        if not free.empty:
            loci = pd.DataFrame(
                [
                    dict(symbol=s, seqname=c["seqname"], strand=c["strand"], five=c["five"])
                    for s, c in fp.items()
                ]
            )
            extra = _table_fives(peaks)
            by_peak = {k: v for k, v in free.groupby(["seqnames", "strand"], observed=True)}
            lab: list = []
            won: list[str] = []
            for (seq, strand), lg in loci.groupby(["seqname", "strand"], observed=True):
                pg = by_peak.get((seq, strand))
                if pg is None or pg.empty:
                    continue
                lg = lg.sort_values("five")
                five = lg["five"].to_numpy(float)
                syms = lg["symbol"].to_numpy(object)
                mid = pg["mid"].to_numpy(float)
                idx = pg.index.to_numpy()
                block = extra.get((str(seq), str(strand)), np.array([], float))
                if five.size:
                    block = np.unique(np.concatenate([block, five])) if block.size else five
                lo = mid - (down if strand == "+" else up)
                hi = mid + (up if strand == "+" else down)
                i0 = np.searchsorted(five, lo, side="left")
                i1 = np.searchsorted(five, hi, side="right")
                for k in np.flatnonzero(i1 > i0):
                    a, b = (int(i0[k]), int(i1[k]))
                    dd = mid[k] - five[a:b]
                    if strand == "-":
                        dd = -dd
                    order = np.argsort(np.abs(dd))
                    chosen = None
                    for j in order:
                        cand_five = five[a:b][j]
                        others = block[block != cand_five]
                        if _has_between(float(mid[k]), float(cand_five), others):
                            continue
                        chosen = syms[a + int(j)]
                        break
                    if chosen is None:
                        continue
                    lab.append(idx[k])
                    won.append(chosen)
            if lab:
                owner.loc[lab] = won
                n_near = len(lab)
        stats_n = dict(якорь=n_named, чужое_имя=n_locked, nearest=n_near)
        return (owner, stats_n)

    def fib_species_means(M: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
        sp = pd.Series([c.split("|")[0] for c in M.columns], index=M.columns)
        per = M.T.groupby(sp.to_numpy()).mean().T
        hit = [g for g in genes if g in per.index]
        y = np.log10(per.reindex(hit)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
        return (y, hit)

    def report_oxphos(
        M: pd.DataFrame, core: list[str], uf: list[str], logm: np.ndarray, tag: str
    ) -> None:
        """Compare the working set on this matrix vs the previous nearest-only run."""
        watch = ["COX6C", "SDHD", "NDUFB1"]
        old_tag = tag.replace("_anchored", "")
        old_path = OUT / f"M_fib_{old_tag}.tsv"
        w1000 = RES / "wgcna_input" / "M_fib_w1000.tsv"

        def one(label: str, mat: pd.DataFrame, genes: list[str]) -> None:
            y, hit = fib_species_means(mat, genes)
            r = stats.linregress(logm, y)
            means = " ".join((f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y)))
            print(f"  {label}: n={len(hit)} r={r.rvalue:+.3f} p={r.pvalue:.4f}")
            print(f"    {means}")

        print("\nOXPHOS, log10(среднее по генам), фибробласты:")
        one(f"якорь {tag}, 55∩вселенная", M, [g for g in core if g in uf])
        one(f"якорь {tag}, все 55 (дыры как есть)", M, core)
        if old_path.exists():
            old = pd.read_csv(old_path, sep="\t", index_col=0)
            u_old = set(old.index)
            g47 = [g for g in core if g in u_old]
            print(f"\n  те же {len(g47)} генов, что были на nearest-only:")
            one("nearest-only", old, g47)
            one("якорь", M, g47)
            print("\n  три гена, которые украли соседний пик:")
            for g in watch:
                print(f"    {g}")
                for lab, mat in (("nearest", old), ("якорь", M)):
                    if g not in mat.index:
                        print(f"      {lab}: нет в матрице")
                        continue
                    sp = pd.Series([c.split("|")[0] for c in mat.columns])
                    y = mat.loc[g].groupby(sp.to_numpy()).mean().reindex(SPECIES_ORDER)
                    print(
                        "      "
                        + lab
                        + ": "
                        + " ".join(
                            (f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y.to_numpy(float)))
                        )
                    )
        if w1000.exists():
            w = pd.read_csv(w1000, sep="\t", index_col=0)
            print("\n  старое окно ±1000 (скрипт 62):")
            one("55 генов", w, core)
            g47 = (
                [g for g in core if g in set(pd.read_csv(old_path, sep="\t", index_col=0).index)]
                if old_path.exists()
                else [g for g in core if g in uf]
            )
            one("те же 47", w, g47)

    def per_library(species: str, up: int, down: int, mode: str) -> pd.DataFrame:
        d = pd.read_csv(s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        fp, src = loci_of(species, d, mode)
        cols = lib_cols(d, "fib") + lib_cols(d, "epi")
        num = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        num = num.where(num >= THR, 0.0)
        pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
        pool = pool.assign(**{c: num[c] for c in cols})
        pool = pool[num.to_numpy().sum(axis=1) > 0].reset_index(drop=True)
        owner, how = assign_peaks(fp, pool, up, down)
        pool = pool.assign(owner=owner).dropna(subset=["owner"])
        m = pool.groupby("owner", observed=True)[cols].sum()
        m = m[m.to_numpy().sum(axis=1) > 0]
        m.columns = [
            f"{species}|{('fib' if c.startswith('norm.C') else 'epi')}|{c.split('.')[1]}"
            for c in cols
        ]
        m.index.name = None
        layers = ", ".join((f"{k} {v}" for k, v in src.items()))
        print(
            f"[{species:8s}] локусов {len(fp):>6} ({layers}); якорь {how['якорь']} / чужое имя {how['чужое_имя']} / nearest {how['nearest']}; генов с сигналом {len(m):>6}",
            flush=True,
        )
        return np.log10(m + 1.0)

    def main() -> None:
        ap = argparse.ArgumentParser()
        ap.add_argument("--up", type=int, default=5000)
        ap.add_argument("--down", type=int, default=1000)
        ap.add_argument(
            "--loci",
            choices=["annot_first", "common_only", "common_then_annot"],
            default="common_then_annot",
        )
        ap.add_argument(
            "--tag-suffix",
            default="anchored",
            help="appended to the output tag so the previous nearest-only matrices are kept",
        )
        a = ap.parse_args(context.argv)
        OUT.mkdir(parents=True, exist_ok=True)
        tag = f"{a.loci}_u{a.up}d{a.down}"
        if a.tag_suffix:
            tag = f"{tag}_{a.tag_suffix}"
        per = {sp: per_library(sp, a.up, a.down, a.loci) for sp in SPECIES_ORDER}
        for sp, m in per.items():
            m.to_csv(OUT / f"M_{sp}_{tag}.tsv", sep="\t")

        def universe(tissue: str) -> list[str]:
            sets = []
            for sp in SPECIES_ORDER:
                c = [x for x in per[sp].columns if f"|{tissue}|" in x]
                m = per[sp][c]
                sets.append(set(m.index[(m > 0).any(axis=1)]))
            return sorted(set.intersection(*sets))

        uf, ue = (universe("fib"), universe("epi"))
        ua = sorted(set(uf) | set(ue))
        print(f"\nвселенная: фибробласты {len(uf)}, эпибласт {len(ue)}, объединение {len(ua)}")
        for nm, u in (("fib", uf), ("epi", ue), ("all", ua)):
            (OUT / f"universe_{nm}_{tag}.txt").write_text("\n".join(u) + "\n")
        logm = np.log10([SPECIES[s]["mass_g"] for s in SPECIES_ORDER])
        nz = [
            int((per[sp][[x for x in per[sp].columns if "|fib|" in x]] > 0).any(axis=1).sum())
            for sp in SPECIES_ORDER
        ]
        r = stats.linregress(logm, np.array(nz, float))
        print(
            "обнаружено в фибробластах: "
            + ", ".join((f"{s} {n}" for s, n in zip(SPECIES_ORDER, nz)))
        )
        print(f"перекос обнаружения по массе: r={r.rvalue:+.3f} p={r.pvalue:.3f}")

        def combine(tissue: str | None, genes: list[str]) -> pd.DataFrame:
            parts = []
            for sp in SPECIES_ORDER:
                c = [x for x in per[sp].columns if tissue is None or f"|{tissue}|" in x]
                parts.append(per[sp].reindex(genes)[c])
            return pd.concat(parts, axis=1)

        mats = {"fib": combine("fib", uf), "epi": combine("epi", ue), "all": combine(None, ua)}
        for nm, m in mats.items():
            m = m.dropna(how="all")
            m.to_csv(OUT / f"M_{nm}_{tag}.tsv", sep="\t")
            print(
                f"M_{nm}: {m.shape[0]} генов x {m.shape[1]} библиотек, пропусков {int(m.isna().sum().sum())}"
            )
        meta = [
            dict(
                sample=c,
                species=c.split("|")[0],
                tissue=c.split("|")[1],
                library=c.split("|")[2],
                mass_g=SPECIES[c.split("|")[0]]["mass_g"],
                log10_mass=np.log10(SPECIES[c.split("|")[0]]["mass_g"]),
            )
            for c in mats["all"].columns
        ]
        pd.DataFrame(meta).to_csv(OUT / f"samples_{tag}.tsv", sep="\t", index=False)
        ox = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
        if ox.exists():
            core = [l.strip() for l in ox.read_text().splitlines() if l.strip()]
            n_in = sum((1 for g in core if g in uf))
            print(f"\nрабочий набор OXPHOS: {n_in} из {len(core)} в новой фибробластной вселенной")
            report_oxphos(mats["fib"], core, uf, logm, tag)
        print(f"\nwrote {OUT}/*_{tag}.tsv")

    return SimpleNamespace(**locals())


def _modules_v2(context):
    """Modules vs mass on the rebuilt panel, plus the centred-within-species network.

    Three things, written next to the old run rather than over it:

    1. The standard module-trait table of script 65, both library-wise (optimistic)
       and species-wise, plus matched-null p with BH across modules — that last
       correction was missing from the first write-up of green and yellow.

    2. The centred-within-species network of script 67, same defaults.

    3. A sweep of mergeCutHeight / deepSplit / minModuleSize so the centred network
       is not judged on the four-module degeneration that made yellow's BH of 0.028
       cheaper than green's 0.054. The question is whether a growth-like module
       survives when 15-20 modules are asked of the same data.

    The OXPHOS working set is the intersection of the old 55 with the new fibroblast
    universe (47 genes). The eight that dropped out are listed, not forced back in.

    Source: analyses/avian/scripts/72_modules_v2.py at d230e56ab942."""
    import importlib
    import subprocess
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    s53 = context.unit("avian.53_three_checks")
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    IN = RES / "wgcna_input_v2"
    WG = RES / "wgcna_out_v2"
    OUT = RES / "modules_v2"
    TAG = "common_only_u5000d1000"
    OXFILE = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
    OLD_YELLOW = RES / "modules_within_species" / "gene_modules_fib_centered.tsv"
    RSCRIPT = str(context.resolve_path("/path/to/user-home/miniconda3/envs/wgcna_r/bin/Rscript"))

    def bh(p) -> np.ndarray:
        return stats.false_discovery_control(
            np.nan_to_num(np.asarray(p, float), nan=1.0), method="bh"
        )

    R_CODE = '\nsuppressMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads(16)\na <- commandArgs(TRUE)\ninf <- a[1]; out <- a[2]\nmch <- as.numeric(a[3]); ds <- as.integer(a[4]); minsz <- as.integer(a[5])\nm <- as.matrix(read.table(inf, sep="\\t", header=TRUE, row.names=1, check.names=FALSE))\nd <- t(m)\ngsg <- goodSamplesGenes(d, verbose=0); if(!gsg$allOK) d <- d[gsg$goodSamples, gsg$goodGenes, drop=FALSE]\nsft <- pickSoftThreshold(d, powerVector=c(1:10, seq(12,30,2)), networkType="signed hybrid", verbose=0)\nok <- sft$fitIndices$SFT.R.sq >= 0.8\np <- if(any(ok)) sft$fitIndices$Power[which(ok)[1]] else 14\ncat(sprintf("softPower=%d  mergeCutHeight=%.2f  deepSplit=%d  minModuleSize=%d  (R2>=0.8: %s)\\n",\n            p, mch, ds, minsz, any(ok)))\nnet <- blockwiseModules(d, power=p, networkType="signed hybrid", TOMType="signed",\n                        minModuleSize=minsz, deepSplit=ds, mergeCutHeight=mch,\n                        numericLabels=TRUE, pamRespectsDendro=FALSE, maxBlockSize=20000, verbose=0)\ncol <- labels2colors(net$colors)\nwrite.table(data.frame(gene=colnames(d), module=col), out, sep="\\t", row.names=FALSE, quote=FALSE)\ncat(sprintf("модулей %d (без grey %d), grey %d\\n", length(unique(col)),\n            length(setdiff(unique(col),"grey")), sum(col=="grey")))\n'

    def trait_table(net: str, samples: pd.DataFrame, ox: list[str]) -> pd.DataFrame:
        gm = pd.read_csv(WG / f"gene_modules_{net}.tsv", sep="\t")
        me = pd.read_csv(WG / f"module_eigengenes_{net}.tsv", sep="\t").set_index("sample")
        s = samples.loc[me.index]
        mass_lib = s["log10_mass"].to_numpy(float)
        rows = []
        for mod in me.columns:
            v = me[mod].to_numpy(float)
            r1, p1 = stats.pearsonr(mass_lib, v)
            per = pd.DataFrame({"sp": s["species"].to_numpy(), "v": v}).groupby("sp")["v"].mean()
            per = per.reindex(SPECIES_ORDER)
            r2, p2 = s53.fit(per.to_numpy(float))
            name = mod.replace("ME", "")
            n_in = int((gm["module"] == name).sum())
            n_ox = int(((gm["module"] == name) & gm["gene"].isin(ox)).sum())
            rows.append(
                dict(
                    сеть=net,
                    модуль=name,
                    n_генов=n_in,
                    n_OXPHOS=n_ox,
                    r_по_библиотекам=round(r1, 3),
                    p_по_библиотекам=p1,
                    r_по_видам=round(r2, 3),
                    p_по_видам=p2,
                )
            )
        t = pd.DataFrame(rows)
        t["bh_по_библиотекам"] = bh(t["p_по_библиотекам"]).round(4)
        t["bh_по_видам"] = bh(t["p_по_видам"]).round(4)
        t["p_по_библиотекам"] = t["p_по_библиотекам"].round(5)
        t["p_по_видам"] = t["p_по_видам"].round(4)
        return t.sort_values("r_по_видам")

    def centred_and_test(
        M: pd.DataFrame,
        meta: pd.DataFrame,
        ox: list[str],
        mch: float,
        ds: int,
        minsz: int,
        tag: str,
        old_yellow: set[str],
    ) -> pd.DataFrame:
        C = M.copy()
        for sp in SPECIES_ORDER:
            cols = [c for c in M.columns if c.startswith(f"{sp}|")]
            C[cols] = M[cols].sub(M[cols].mean(axis=1), axis=0)
        C = C[C.std(axis=1) > 0]
        fmat = OUT / "M_fib_centered.tsv"
        if not fmat.exists():
            C.to_csv(fmat, sep="\t")
        rf = OUT / "_run.R"
        rf.write_text(R_CODE)
        gmf = OUT / f"gene_modules_centered_{tag}.tsv"
        r = subprocess.run(
            [RSCRIPT, str(rf), str(fmat), str(gmf), str(mch), str(ds), str(minsz)],
            capture_output=True,
            text=True,
        )
        print(r.stdout.strip())
        if r.returncode != 0:
            print("R упал:", r.stderr.strip()[-400:])
            return pd.DataFrame()
        gm = pd.read_csv(gmf, sep="\t")
        permean = M.T.groupby(meta["species"].to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
        dist = gm[gm["gene"].isin(ox)]["module"].value_counts()
        rows = []
        for mod in sorted(gm["module"].unique()):
            if mod == "grey":
                continue
            members = [g for g in gm[gm["module"] == mod]["gene"] if g in permean.index]
            if len(members) < 10:
                continue
            r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
            pn, _ = s53.matched_null_p(permean, members, r_set)
            overlap = len(set(members) & old_yellow)
            rows.append(
                dict(
                    настройка=tag,
                    merge=mch,
                    deepSplit=ds,
                    minSize=minsz,
                    модуль=mod,
                    n_генов=len(members),
                    n_OXPHOS=int(dist.get(mod, 0)),
                    overlap_old_yellow=overlap,
                    r_по_видам=round(r_set, 3),
                    p=round(p_set, 4),
                    p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                )
            )
        t = pd.DataFrame(rows)
        if t.empty:
            return t
        t["bh_регрессии"] = bh(t["p"]).round(4)
        t["bh_нуля"] = bh(t["p_против_нуля"].fillna(1.0)).round(4)
        t = t.sort_values("r_по_видам")
        n_mod = gm["module"].nunique() - int("grey" in set(gm["module"]))
        print(
            f"  настроек {tag}: модулей без grey {n_mod}; BH нуля < 0.05 у {int((t['bh_нуля'] < 0.05).sum())}"
        )
        return t

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        samples = pd.read_csv(IN / f"samples_{TAG}.tsv", sep="\t").set_index("sample")
        ox_all = [l.strip() for l in OXFILE.read_text().splitlines() if l.strip()]
        uf = set(pd.read_csv(IN / f"universe_fib_{TAG}.txt", header=None)[0])
        ox = [g for g in ox_all if g in uf]
        lost = [g for g in ox_all if g not in uf]
        print(
            f"рабочий набор OXPHOS: {len(ox)} из {len(ox_all)} в новой вселенной; нет: {', '.join(lost)}\n"
        )
        (OUT / "working_oxphos_in_universe.txt").write_text("\n".join(ox) + "\n")
        old_yellow: set[str] = set()
        if OLD_YELLOW.exists():
            oy = pd.read_csv(OLD_YELLOW, sep="\t")
            old_yellow = set(oy.loc[oy["module"] == "yellow", "gene"])
            print(f"старый yellow центрированной сети: {len(old_yellow)} генов\n")
        for net in ("fib", "epi", "all_tc"):
            t = trait_table(net, samples, ox)
            t.to_csv(OUT / f"module_trait_{net}.tsv", sep="\t", index=False)
            print("=" * 100)
            print(f"сеть {net}: модулей {len(t) - int('grey' in set(t['модуль']))} без grey")
            print(
                f"  BH<0.05 по библиотекам: {int(((t['модуль'] != 'grey') & (t['bh_по_библиотекам'] < 0.05)).sum())}"
            )
            print(
                f"  BH<0.05 по видам:       {int(((t['модуль'] != 'grey') & (t['bh_по_видам'] < 0.05)).sum())}"
            )
            hit = t[t["n_OXPHOS"] > 0]
            if len(hit):
                print(
                    hit[
                        ["модуль", "n_генов", "n_OXPHOS", "r_по_видам", "p_по_видам", "bh_по_видам"]
                    ].to_string(index=False)
                )
        gm = pd.read_csv(WG / "gene_modules_fib.tsv", sep="\t")
        M = pd.read_csv(IN / f"M_fib_{TAG}.tsv", sep="\t", index_col=0)
        meta = samples.loc[M.columns]
        permean = M.T.groupby(meta["species"].to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
        inmod = gm[gm["gene"].isin(ox)]
        print("\nкуда попали гены OXPHOS (сеть фибробластов v2):")
        print(inmod["module"].value_counts().to_string())
        det = []
        for mod, k in inmod["module"].value_counts().items():
            members = gm[gm["module"] == mod]["gene"].tolist()
            r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
            pn, _ = s53.matched_null_p(permean, [g for g in members if g in permean.index], r_set)
            det.append(
                dict(
                    модуль=mod,
                    n_всего=len(members),
                    n_OXPHOS=int(k),
                    r_модуля_по_видам=round(r_set, 3),
                    p=round(p_set, 4),
                    p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                )
            )
            print(
                f"  {mod}: n={len(members)} OXPHOS={k} r={r_set:+.3f} нуль {('—' if np.isnan(pn) else round(pn, 4))}"
            )
        dett = pd.DataFrame(det)
        dett["bh_нуля"] = bh(dett["p_против_нуля"].fillna(1.0)).round(4)
        dett.to_csv(OUT / "oxphos_modules_detail.tsv", sep="\t", index=False)
        print("BH нуля по модулям с OXPHOS:")
        print(dett.sort_values("p_против_нуля").to_string(index=False))
        print("\n=== центрированная сеть, сетка настроек ===")
        grid = [
            (0.25, 2, 30, "default"),
            (0.15, 2, 30, "merge015"),
            (0.1, 2, 20, "merge010_sz20"),
            (0.1, 3, 20, "merge010_ds3"),
            (0.05, 4, 15, "split"),
        ]
        parts = []
        for mch, ds, minsz, tag in grid:
            print(f"\n--- {tag} ---")
            t = centred_and_test(M, meta, ox, mch, ds, minsz, tag, old_yellow)
            if len(t):
                parts.append(t)
                print(t.to_string(index=False))
        if parts:
            allc = pd.concat(parts, ignore_index=True)
            allc.to_csv(OUT / "centered_sweep.tsv", sep="\t", index=False)
            print(f"\nwrote {OUT}/centered_sweep.tsv")

    return SimpleNamespace(**locals())


def _nineteen_via_cage(context):
    """The 19 genes that carried the old ±1000 decline, measured via CAGE.

    The old window called them detected in all seven because a peak sat near a
    miniprot 5'. For Augustus species that 5' was often a neighbour (RGS22, ADCK2,
    GATC, CABP7). This script does not change the set: it asks, for each of those
    19, where the gene actually is and what CAGE belongs to it.

    Two independent placements, then CAGE at that place:

      1. CDS from the CHICKEN protein only (not the highest-scoring turkey/emu
         isoform — those can extend 5' onto a neighbour, which is how emu COX6C
         ate RGS22). Promoter = CAGE on the transcribed strand in the 5' catchment
         up to the next gene, window [-5000, +300] bp of that CDS 5'. A peak already
         named as a different gene is kept only if it overlaps the lifted chicken
         promoter (the name was wrong, the place is right).

      2. The chicken CAGE peak that already carries this gene's name, sequence
         lifted with minimap2 onto each genome. CAGE overlapping the lift (±200 bp)
         is the homologous promoter.

    A species counts as having the gene if the chicken-protein CDS is there AND
    there is CAGE from (1) or (2). Native-annotation species almost always have
    both; Augustus are the test.

    Source: analyses/avian/scripts/74_nineteen_via_cage.py at d230e56ab942."""
    import importlib
    import re
    import subprocess
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    s62 = context.unit("avian.62_build_library_matrices")
    RES = context.unit("avian.oxphos_thread_lib").RES
    RUSLAN = context.unit("avian.oxphos_thread_lib").RUSLAN
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    ASM = Path(
        str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies"))
    )
    OUT = RES / "nineteen_via_cage"
    NINETEEN = [
        "ATP5F1E",
        "ATP6V0E1",
        "ATP6V0E2",
        "COX6A1",
        "COX6C",
        "COX7A2",
        "NDUFA5",
        "NDUFA7",
        "NDUFAB1",
        "NDUFB1",
        "NDUFB2",
        "NDUFB5",
        "NDUFB6",
        "NDUFB8",
        "NDUFC2",
        "NDUFS1",
        "SDHD",
        "UQCR10",
        "UQCRC1",
    ]
    WANT = set(NINETEEN)
    LOGM = np.log10([SPECIES[s]["mass_g"] for s in SPECIES_ORDER])
    THR = 1.0
    UP, DOWN = (5000, 300)
    LIFT_PAD = 200
    GENOME = {
        "Chicken": ASM / "Chicken" / "genome.fa",
        "Duck": ASM / "Duck" / "genome.fa",
        "Quail": ASM / "Quail" / "genome.fa",
        "Zebra": ASM / "Zebra" / "genome.fa",
        "Peacock": ASM / "Peacock" / "genome_localnames.fa",
        "Turkey": ASM / "Turkey" / "genome.fa",
        "Emu": ASM / "Emu" / "genome.fa",
    }
    PRODUCT = {
        "ATP5F1E": "ATP synthase F1 ε, ~50 aa",
        "ATP6V0E1": "V-ATPase c″ subunit",
        "ATP6V0E2": "V-ATPase c″ isoform 2",
        "COX6A1": "Complex IV subunit 6A liver isoform",
        "COX6C": "Complex IV subunit 6C",
        "COX7A2": "Complex IV subunit 7A2",
        "NDUFA5": "Complex I supernumerary A5",
        "NDUFA7": "Complex I supernumerary A7",
        "NDUFAB1": "Complex I ACP / SDAP",
        "NDUFB1": "Complex I supernumerary B1 (MNLL)",
        "NDUFB2": "Complex I supernumerary B2 (AGGG)",
        "NDUFB5": "Complex I supernumerary B5",
        "NDUFB6": "Complex I supernumerary B6",
        "NDUFB8": "Complex I supernumerary B8",
        "NDUFC2": "Complex I supernumerary C2 (B14.5b)",
        "NDUFS1": "Complex I core Fe-S 75 kDa",
        "SDHD": "Complex II cytochrome b small subunit",
        "UQCR10": "Complex III subunit 10 (QCR9)",
        "UQCRC1": "Complex III core subunit 1",
    }

    def contig(s: str) -> str:
        s = str(s).split()[0]
        return s.split("|")[-1]

    def complex_of(g: str) -> str:
        if g.startswith("NDUF"):
            return "CI"
        if g.startswith("SDH"):
            return "CII"
        if g.startswith("UQCR"):
            return "CIII"
        if g.startswith("COX"):
            return "CIV"
        if g.startswith("ATP5"):
            return "CV"
        if g.startswith("ATP6V"):
            return "V-ATPase"
        return "other"

    def sh(cmd: list[str], **kw):
        return subprocess.run(cmd, check=True, **kw)

    def chicken_genes(gff: Path) -> pd.DataFrame:
        rows = []
        with gff.open() as fh:
            for line in fh:
                if line[0] == "#" or "\tgene\t" not in line:
                    continue
                f = line.split("\t")
                if len(f) < 9:
                    continue
                if (
                    "gene_biotype=protein_coding" not in f[8]
                    and "biotype=protein_coding" not in f[8]
                ):
                    if "gene_biotype=" in f[8] and "protein_coding" not in f[8]:
                        continue
                m = re.search("(?:gene=|Name=)([^;]+)", f[8])
                if not m:
                    continue
                sym = normalize_symbol(m.group(1))
                if not sym:
                    continue
                st = f[6]
                five = int(f[3] if st == "+" else f[4])
                rows.append(
                    dict(
                        symbol=sym,
                        seqname=contig(f[0]),
                        strand=st,
                        start=int(f[3]),
                        end=int(f[4]),
                        five=five,
                        product=re.search("description=([^;]+)", f[8]),
                    )
                )
        g = pd.DataFrame(rows).drop_duplicates("symbol", keep="first")
        g["product"] = g["product"].map(lambda m: m.group(1) if m else "")
        return g

    def flanks_of(genes: pd.DataFrame, sym: str) -> tuple[str, str]:
        row = genes[genes["symbol"] == sym]
        if row.empty:
            return ("", "")
        r = row.iloc[0]
        same = genes[genes["seqname"] == r["seqname"]].sort_values("start").reset_index(drop=True)
        hits = same.index[same["symbol"] == sym].tolist()
        if not hits:
            return ("", "")
        i = int(hits[0])
        left = same.iloc[i - 1]["symbol"] if i > 0 else ""
        right = same.iloc[i + 1]["symbol"] if i < len(same) - 1 else ""
        return (str(left), str(right))

    def chicken_protein_loci(gff: Path, symbols: set[str]) -> dict[str, dict]:
        """One locus per symbol: the chicken_* query with the highest identity."""
        best: dict[str, tuple[float, dict]] = {}
        if not gff.exists():
            return {}
        with gff.open() as fh:
            for line in fh:
                if line[0] == "#" or "\tmRNA\t" not in line:
                    continue
                f = line.split("\t")
                if len(f) < 9:
                    continue
                m = re.search("Target=(\\S+)", f[8])
                if not m or "##" not in m.group(1):
                    continue
                src, raw = m.group(1).split("##", 1)
                if not src.lower().startswith("chicken"):
                    continue
                sym = normalize_symbol(raw)
                if sym not in symbols:
                    continue
                ident = re.search("Identity=([0-9.]+)", f[8])
                idv = float(ident.group(1)) if ident else 0.0
                sc = float(f[5]) if f[5] != "." else 0.0
                st = f[6]
                seq = contig(f[0])
                rec = dict(
                    seqname=seq,
                    strand=st,
                    start=int(f[3]),
                    end=int(f[4]),
                    five=float(f[3] if st == "+" else f[4]),
                    ident=idv,
                    score=sc,
                    query=src,
                )
                key = (idv, sc)
                prev = best.get(sym)
                if prev is None or key > (prev[0], prev[1]["score"]):
                    best[sym] = (idv, rec)
        return {k: v[1] for k, v in best.items()}

    def all_fives(gff: Path) -> dict[str, list[tuple[str, float, str]]]:
        """seqname -> list of (symbol, five, strand) for catchment barriers."""
        best: dict[str, tuple[float, tuple]] = {}
        with gff.open() as fh:
            for line in fh:
                if line[0] == "#" or "\tmRNA\t" not in line:
                    continue
                f = line.split("\t")
                if len(f) < 9:
                    continue
                m = re.search("Target=(\\S+)", f[8])
                if not m or "##" not in m.group(1):
                    continue
                src = m.group(1)
                if not src.lower().startswith("chicken"):
                    continue
                sym = normalize_symbol(src.split("##", 1)[1])
                if not sym:
                    continue
                sc = float(f[5]) if f[5] != "." else 0.0
                ident = re.search("Identity=([0-9.]+)", f[8])
                idv = float(ident.group(1)) if ident else 0.0
                st = f[6]
                seq = contig(f[0])
                five = float(f[3] if st == "+" else f[4])
                rec = (idv, sc, (sym, five, st, seq))
                if sym not in best or (idv, sc) > (best[sym][0], best[sym][1]):
                    best[sym] = rec
        by: dict[str, list[tuple[str, float, str]]] = {}
        for _, _, (sym, five, st, seq) in best.values():
            by.setdefault(seq, []).append((sym, five, st))
        for seq in by:
            by[seq].sort(key=lambda t: t[1])
        return by

    def oriented(mid: float, five: float, strand: str) -> float:
        d = mid - five
        return -d if strand == "-" else d

    def catchment_lohi(
        five: float, strand: str, barriers: list[tuple[str, float, str]], self_sym: str
    ) -> tuple[float, float]:
        """Window [-UP, +DOWN] of CDS 5', clipped by the next gene 5' upstream/down."""
        if strand == "+":
            up_b = [b[1] for b in barriers if b[0] != self_sym and b[1] < five]
            down_b = [b[1] for b in barriers if b[0] != self_sym and b[1] > five]
            lo = max(five - UP, max(up_b) if up_b else five - UP)
            hi = min(five + DOWN, min(down_b) if down_b else five + DOWN)
        else:
            up_b = [b[1] for b in barriers if b[0] != self_sym and b[1] > five]
            down_b = [b[1] for b in barriers if b[0] != self_sym and b[1] < five]
            lo = max(five - DOWN, max(down_b) if down_b else five - DOWN)
            hi = min(five + UP, min(up_b) if up_b else five + UP)
        return (lo, hi)

    def load_peaks(species: str) -> pd.DataFrame:
        d = pd.read_csv(s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
        cols = lib_cols(d, "fib")
        num = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        num = num.where(num >= THR, 0.0)
        d = d.assign(
            mid=(d["start"] + d["end"]) / 2.0, fib=num.mean(axis=1), seq=d["seqnames"].map(contig)
        )
        return d[d["fib"] > 0].copy()

    def ensure_faidx(fa: Path) -> None:
        if not Path(str(fa) + ".fai").exists():
            sh(["samtools", "faidx", str(fa)])

    def chicken_promoter_fasta(peaks: pd.DataFrame, fa: Path, out_fa: Path) -> pd.DataFrame:
        """One strongest named CAGE peak per gene on chicken."""
        ensure_faidx(fa)
        rows = []
        bed = OUT / "chicken_promoters.bed"
        with bed.open("w") as fh:
            for g in NINETEEN:
                hit = peaks[peaks["symbol"] == g].sort_values("fib", ascending=False)
                if hit.empty:
                    continue
                r = hit.iloc[0]
                start = max(0, int(r["start"]) - 50)
                end = int(r["end"]) + 50
                name = f"{g}|{r['id']}"
                fh.write(f"{r['seqnames']}\t{start}\t{end}\t{name}\t0\t{r['strand']}\n")
                rows.append(
                    dict(
                        gene=g,
                        peak_id=r["id"],
                        seqname=r["seqnames"],
                        start=int(r["start"]),
                        end=int(r["end"]),
                        strand=r["strand"],
                        fib=float(r["fib"]),
                    )
                )
        sh(
            [
                "bedtools",
                "getfasta",
                "-fi",
                str(fa),
                "-bed",
                str(bed),
                "-s",
                "-nameOnly",
                "-fo",
                str(out_fa),
            ]
        )
        return pd.DataFrame(rows)

    def lift_one(genome: Path, query: Path, paf: Path) -> pd.DataFrame:
        ensure_faidx(genome)
        with paf.open("w") as fh:
            sh(
                [
                    "minimap2",
                    "-t",
                    "8",
                    "-x",
                    "asm20",
                    "-k",
                    "11",
                    "-w",
                    "5",
                    "--secondary=no",
                    "-N",
                    "3",
                    str(genome),
                    str(query),
                ],
                stdout=fh,
                stderr=subprocess.DEVNULL,
            )
        cols = [
            "q",
            "qlen",
            "qs",
            "qe",
            "strand",
            "t",
            "tlen",
            "ts",
            "te",
            "nmatch",
            "alen",
            "mapq",
        ]
        if paf.stat().st_size == 0:
            return pd.DataFrame(columns=cols + ["gene"])
        raw = pd.read_csv(paf, sep="\t", header=None, usecols=range(12), names=cols)
        raw["gene"] = raw["q"].str.split("|").str[0]
        raw["t"] = raw["t"].map(contig)
        raw["ident"] = raw["nmatch"] / raw["alen"].clip(lower=1)
        raw = raw.sort_values(["gene", "mapq", "ident"], ascending=[True, False, False])
        return raw.drop_duplicates("gene", keep="first")

    def score_peaks(near: pd.DataFrame) -> float:
        if near is None or near.empty:
            return float("nan")
        v = float(near["fib"].sum())
        return v if v > 0 else float("nan")

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        print("курица: гены NCBI + CAGE промоторы", flush=True)
        ch_genes = chicken_genes(ASM / "Chicken" / "genomic.gff")
        peaks = {sp: load_peaks(sp) for sp in SPECIES_ORDER}
        ch_prom = chicken_promoter_fasta(
            peaks["Chicken"], GENOME["Chicken"], OUT / "chicken_promoters.fa"
        )
        ch_prom.to_csv(OUT / "chicken_promoters.tsv", sep="\t", index=False)
        cards = []
        for g in NINETEEN:
            left, right = flanks_of(ch_genes, g)
            prow = ch_prom[ch_prom["gene"] == g]
            grow = ch_genes[ch_genes["symbol"] == g]
            cards.append(
                dict(
                    gene=g,
                    complex=complex_of(g),
                    product=PRODUCT.get(g, ""),
                    chicken_chr="" if grow.empty else grow.iloc[0]["seqname"],
                    chicken_tss="" if grow.empty else int(grow.iloc[0]["five"]),
                    chicken_strand="" if grow.empty else grow.iloc[0]["strand"],
                    flank_left=left,
                    flank_right=right,
                    cage_tpm=None if prow.empty else round(float(prow.iloc[0]["fib"]), 1),
                    cage_peak="" if prow.empty else prow.iloc[0]["peak_id"],
                )
            )
        pd.DataFrame(cards).to_csv(OUT / "identity_cards.tsv", sep="\t", index=False)
        print("карточки:", OUT / "identity_cards.tsv", flush=True)
        lifts = {}
        for sp in SPECIES_ORDER:
            if sp == "Chicken":
                continue
            print(f"lift → {sp}", flush=True)
            lifts[sp] = lift_one(GENOME[sp], OUT / "chicken_promoters.fa", OUT / f"lift_{sp}.paf")
            print(
                f"  mapped {(lifts[sp]['gene'].nunique() if len(lifts[sp]) else 0)} / 19",
                flush=True,
            )
        loci = {
            sp: chicken_protein_loci(ASM / sp / "reanno_named.gff", WANT) for sp in SPECIES_ORDER
        }
        barriers = {sp: all_fives(ASM / sp / "reanno_named.gff") for sp in SPECIES_ORDER}
        lifts["Chicken"] = pd.DataFrame(
            [
                dict(
                    gene=r.gene,
                    t=r.seqname,
                    ts=r.start - 50,
                    te=r.end + 50,
                    strand=r.strand,
                    mapq=60,
                    ident=1.0,
                )
                for r in ch_prom.itertuples()
            ]
        )
        rows = []
        for g in NINETEEN:
            for sp in SPECIES_ORDER:
                rec = dict(gene=g, species=sp, complex=complex_of(g))
                loc = loci[sp].get(g)
                rec["cds"] = False
                if loc:
                    rec.update(
                        cds=True,
                        cds_ident=round(loc["ident"], 3),
                        cds_seq=loc["seqname"],
                        cds_strand=loc["strand"],
                        cds_five=int(loc["five"]),
                    )
                pk = peaks[sp]
                catch_val = np.nan
                catch_n = 0
                catch_names = ""
                if loc:
                    barr = barriers[sp].get(loc["seqname"], [])
                    lo, hi = catchment_lohi(loc["five"], loc["strand"], barr, g)
                    rec["catch_lo"] = int(lo)
                    rec["catch_hi"] = int(hi)
                    sub = pk[
                        (pk["seq"] == loc["seqname"])
                        & (pk["strand"] == loc["strand"])
                        & (pk["mid"] >= lo)
                        & (pk["mid"] <= hi)
                    ].copy()
                    lf = lifts.get(sp, pd.DataFrame())
                    lift_hit = lf[lf["gene"] == g] if len(lf) else lf
                    lift_iv = None
                    if lift_hit is not None and len(lift_hit):
                        lr = lift_hit.iloc[0]
                        lift_iv = (
                            contig(lr["t"]),
                            int(lr["ts"]) - LIFT_PAD,
                            int(lr["te"]) + LIFT_PAD,
                        )
                    keep_idx = []
                    for i, r in sub.iterrows():
                        nom = r["symbol"] if isinstance(r["symbol"], str) else ""
                        if nom == g or nom == "":
                            keep_idx.append(i)
                            continue
                        if (
                            lift_iv
                            and r["seq"] == lift_iv[0]
                            and (r["end"] >= lift_iv[1])
                            and (r["start"] <= lift_iv[2])
                        ):
                            keep_idx.append(i)
                            continue
                    kept = sub.loc[keep_idx] if keep_idx else sub.iloc[0:0]
                    catch_val = score_peaks(kept)
                    catch_n = int(len(kept))
                    catch_names = ",".join(
                        sorted({str(x) for x in kept["symbol"] if isinstance(x, str) and x})[:5]
                    )
                rec["cage_catchment"] = catch_val
                rec["catch_n"] = catch_n
                rec["catch_names"] = catch_names
                lf = lifts.get(sp, pd.DataFrame())
                rec["lifted"] = False
                lift_val = np.nan
                if lf is not None and len(lf) and (g in set(lf["gene"])):
                    lr = lf[lf["gene"] == g].iloc[0]
                    rec.update(
                        lifted=True,
                        lift_seq=contig(lr["t"]),
                        lift_start=int(lr["ts"]),
                        lift_end=int(lr["te"]),
                        lift_strand=lr["strand"],
                        lift_mapq=int(lr["mapq"]),
                        lift_ident=round(float(lr["ident"]), 3),
                    )
                    sub = pk[
                        (pk["seq"] == contig(lr["t"]))
                        & (pk["end"] >= int(lr["ts"]) - LIFT_PAD)
                        & (pk["start"] <= int(lr["te"]) + LIFT_PAD)
                    ]
                    st = str(lr["strand"])
                    if st in "+-":
                        sub = sub[sub["strand"] == st]
                    lift_val = score_peaks(sub)
                    rec["lift_peak_names"] = ",".join(
                        sorted({str(x) for x in sub["symbol"] if isinstance(x, str) and x})[:5]
                    )
                rec["cage_lift"] = lift_val
                chosen = np.nan
                how = "нет"
                if loc and rec.get("lifted") and (rec.get("lift_seq") == loc["seqname"]):
                    dist = abs((int(rec["lift_start"]) + int(rec["lift_end"])) / 2 - loc["five"])
                    rec["lift_vs_cds_bp"] = int(dist)
                    if dist <= 20000 and np.isfinite(lift_val):
                        chosen, how = (lift_val, "lift+CDS")
                if not np.isfinite(chosen) and np.isfinite(catch_val):
                    chosen, how = (catch_val, "catchment")
                if not np.isfinite(chosen) and loc and rec.get("lifted") and np.isfinite(lift_val):
                    chosen, how = (lift_val, "lift_only")
                rec["tpm"] = chosen
                rec["how"] = how if np.isfinite(chosen) else "CDS_без_CAGE" if loc else "нет_локуса"
                rows.append(rec)
        t = pd.DataFrame(rows)
        t.to_csv(OUT / "per_gene_species.tsv", sep="\t", index=False)
        M = t.pivot(index="gene", columns="species", values="tpm").reindex(
            index=NINETEEN, columns=SPECIES_ORDER
        )
        M.to_csv(OUT / "M_nineteen_tpm.tsv", sep="\t")
        det = M.notna() & (M > 0)
        print("\nобнаружено (CAGE при курином CDS или lift):")
        print(" ", " ".join((f"{sp} {int(det[sp].sum())}" for sp in SPECIES_ORDER)))
        n7 = [g for g in NINETEEN if det.loc[g].all()]
        print(f" во всех семи: {len(n7)} — {(', '.join(n7) if n7 else '—')}")

        def fit_of(genes: list[str], label: str) -> None:
            if not genes:
                print(f"\n{label}: пусто")
                return
            y = np.log10(M.reindex(genes)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
            r = stats.linregress(LOGM, y)
            print(f"\n{label}: n={len(genes)} r={r.rvalue:+.3f} p={r.pvalue:.4f}")
            print(" ", " ".join((f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y))))

        fit_of(n7, "19, только гены с CAGE во всех семи")
        y_found = []
        ns = []
        for sp in SPECIES_ORDER:
            v = M[sp].to_numpy(float)
            ok = np.isfinite(v) & (v > 0)
            ns.append(int(ok.sum()))
            y_found.append(np.log10(v[ok].mean()) if ok.any() else np.nan)
        r = stats.linregress(LOGM, np.array(y_found, float))
        print(
            f"\nсреднее только найденных (дыры завышают): n={ns} r={r.rvalue:+.3f} p={r.pvalue:.4f}"
        )
        print(" ", " ".join((f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y_found))))
        old = pd.read_csv(RES / "wgcna_input" / "M_fib_w1000.tsv", sep="\t", index_col=0)
        spcol = pd.Series([c.split("|")[0] for c in old.columns], index=old.columns)
        per = old.T.groupby(spcol.to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
        y = np.log10(per.reindex(NINETEEN).mean(axis=0).to_numpy(float))
        r = stats.linregress(LOGM, y)
        print(f"\nте же 19 на СТАРОМ окне ±1000: r={r.rvalue:+.3f} p={r.pvalue:.4f}")
        print(" ", " ".join((f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y))))
        print("\nкак измерен каждый ген × вид (how):")
        print(
            t.groupby(["species", "how"])
            .size()
            .unstack(fill_value=0)
            .reindex(SPECIES_ORDER)
            .to_string()
        )
        print(f"\nwrote {OUT}")

    return SimpleNamespace(**locals())


def _honest_gene_mass_scan(context):
    """Per-gene correlation with adult mass under the homologous-CAGE rule.

    Same assignment as the 19-gene rescue (nearest CDS ≤ 2 kb, chicken protein
    only, peak goes to the closer 5' if two genes compete). A gene enters the
    scan only if that rule gives a fibroblast CAGE number in all seven species.
    Empty cells are not filled with zero and are not averaged over.

    Null: exact permutation of the seven mass labels (7! = 5040). BH over genes.

    Source: analyses/avian/scripts/75_honest_gene_mass_scan.py at d230e56ab942."""
    import importlib
    import itertools
    import re
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    s74 = context.unit("avian.74_nineteen_via_cage")
    DATA = context.unit("avian.oxphos_thread_lib").DATA
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    load_geneset = context.unit("avian.oxphos_thread_lib").load_geneset
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    ASM = Path(
        str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies"))
    )
    OUT = RES / "honest_gene_scan"
    MAX_BP = 2000
    LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
    NINETEEN = list(s74.NINETEEN)
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

    def contig(s: str) -> str:
        return s74.contig(s)

    def all_chicken_loci(gff: Path) -> pd.DataFrame:
        """One locus per symbol: chicken_* query with highest Identity, then score."""
        best: dict[str, tuple[float, float, dict]] = {}
        with gff.open() as fh:
            for line in fh:
                if "\tmRNA\t" not in line:
                    continue
                f = line.split("\t")
                if len(f) < 9:
                    continue
                m = re.search("Target=(\\S+)", f[8])
                if not m or "##" not in m.group(1):
                    continue
                src, raw = m.group(1).split("##", 1)
                if not src.lower().startswith("chicken"):
                    continue
                sym = normalize_symbol(raw)
                if not sym:
                    continue
                ident = re.search("Identity=([0-9.]+)", f[8])
                idv = float(ident.group(1)) if ident else 0.0
                sc = float(f[5]) if f[5] != "." else 0.0
                st = f[6]
                rec = dict(
                    symbol=sym,
                    seqname=contig(f[0]),
                    strand=st,
                    five=float(f[3] if st == "+" else f[4]),
                    ident=idv,
                    score=sc,
                )
                prev = best.get(sym)
                if prev is None or (idv, sc) > (prev[0], prev[1]):
                    best[sym] = (idv, sc, rec)
        if not best:
            return pd.DataFrame(columns=["symbol", "seqname", "strand", "five", "ident", "score"])
        return pd.DataFrame([v[2] for v in best.values()])

    def assign_species(peaks: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
        """Sum fibroblast TPM of peaks whose midpoint is closer to this CDS 5'
        than to any other, and within MAX_BP."""
        if loci.empty or peaks.empty:
            return pd.DataFrame(
                columns=["gene", "tpm", "n_peaks", "dist", "peak_names", "frac_foreign"]
            )
        index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
        for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
            g = g.sort_values("five")
            index[seq, st] = (g["five"].to_numpy(float), g["symbol"].tolist())
        acc: dict[str, list] = {}
        for r in peaks.itertuples():
            key = (r.seq, r.strand)
            packed = index.get(key)
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
            if d > MAX_BP:
                continue
            if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-06:
                continue
            gene = syms[jbest]
            nom = r.symbol if isinstance(r.symbol, str) else ""
            acc.setdefault(gene, []).append((float(r.fib), d, nom))
        rows = []
        for gene, hits in acc.items():
            tpm = float(sum((h[0] for h in hits)))
            if tpm <= 0:
                continue
            names = sorted({h[2] for h in hits if h[2]})
            foreign = sum((h[0] for h in hits if h[2] and h[2] != gene))
            rows.append(
                dict(
                    gene=gene,
                    tpm=tpm,
                    n_peaks=len(hits),
                    dist=min((h[1] for h in hits)),
                    peak_names=",".join(names[:6]),
                    frac_foreign=foreign / tpm,
                )
            )
        return pd.DataFrame(rows)

    def r_matrix(X: np.ndarray, y: np.ndarray) -> np.ndarray:
        yc = y - y.mean()
        Xc = X - X.mean(axis=1, keepdims=True)
        num = Xc @ yc
        den = np.sqrt((Xc**2).sum(axis=1) * (yc**2).sum())
        return num / np.where(den == 0, np.nan, den)

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        ox = {normalize_symbol(g) for g in load_geneset("kegg_00190")} - {""}
        print("локусы куриного белка + CAGE ≤ 2 т.п.н.", flush=True)
        wide = {}
        meta_rows = []
        ident_wide = {}
        for sp in SPECIES_ORDER:
            print(f"  {sp}", flush=True)
            loci = all_chicken_loci(ASM / sp / "reanno_named.gff")
            peaks = s74.load_peaks(sp)
            assigned = assign_species(peaks, loci)
            wide[sp] = (
                assigned.set_index("gene")["tpm"] if len(assigned) else pd.Series(dtype=float)
            )
            ident_wide[sp] = (
                loci.set_index("symbol")["ident"] if len(loci) else pd.Series(dtype=float)
            )
            if len(assigned):
                a = assigned.copy()
                a["species"] = sp
                meta_rows.append(a)
            print(f"    локусов {len(loci)}  генов с CAGE {len(assigned)}", flush=True)
        M = pd.DataFrame(wide).reindex(columns=SPECIES_ORDER)
        I = pd.DataFrame(ident_wide).reindex(index=M.index, columns=SPECIES_ORDER)
        M.to_csv(OUT / "M_fib_nearest2kb.tsv", sep="\t")
        pd.concat(meta_rows, ignore_index=True).to_csv(
            OUT / "peak_assignment.tsv", sep="\t", index=False
        )
        det = M.notna() & (M > 0)
        n7 = det.all(axis=1)
        genes = M.index[n7].tolist()
        print(
            f"\nво всех семи с CAGE: {len(genes)}  (локусов хотя бы у одного: {int(M.notna().any(axis=1).sum())})"
        )
        prev = RES / "nineteen_via_cage" / "M_nineteen_nearest_cds.tsv"
        if prev.exists():
            P = pd.read_csv(prev, sep="\t", index_col=0)
            print("сверка с 19 nearest_cds:")
            for g in NINETEEN:
                if g not in M.index:
                    print(f"  {g}: нет в скане")
                    continue
                a = M.loc[g].to_numpy(float)
                b = (
                    P.loc[g].reindex(SPECIES_ORDER).to_numpy(float)
                    if g in P.index
                    else np.full(7, np.nan)
                )
                ok = np.isfinite(a) & np.isfinite(b) & (b > 0)
                if ok.any():
                    rel = np.nanmax(np.abs(a[ok] - b[ok]) / np.clip(b[ok], 1e-06, None))
                else:
                    rel = np.nan
                print(
                    f"  {g}: max отн. разница {rel:.4f}  скан={np.array2string(np.nan_to_num(a), precision=0)}"
                )
        X = np.log10(M.loc[genes, SPECIES_ORDER].to_numpy(float))
        r_obs = r_matrix(X, LOGM)
        p_obs = 2 * stats.t.sf(
            np.abs(r_obs) * np.sqrt(5 / np.clip(1 - r_obs**2, 1e-12, None)), df=5
        )
        bh = stats.false_discovery_control(np.nan_to_num(p_obs, nan=1.0), method="bh")
        orders = list(itertools.permutations(range(7)))
        n_ge = np.zeros(len(genes), dtype=int)
        thr = (0.7, 0.75, 0.8, 0.85, 0.9)
        null_counts = {t: [] for t in thr}
        for perm in orders:
            rp = r_matrix(X, LOGM[list(perm)])
            n_ge += (np.abs(rp) >= np.abs(r_obs)).astype(int)
            for t in thr:
                null_counts[t].append(int((np.abs(rp) >= t).sum()))
        p_perm = n_ge / len(orders)
        min_ident = I.loc[genes].min(axis=1)
        foreign = (
            pd.concat(meta_rows, ignore_index=True)
            .groupby("gene")["frac_foreign"]
            .mean()
            .reindex(genes)
        )
        fit = pd.DataFrame(
            {
                "gene": genes,
                "r": r_obs,
                "p": p_obs,
                "bh": bh,
                "p_perm": p_perm,
                "min_ident": min_ident.to_numpy(float),
                "frac_foreign": foreign.fillna(0).to_numpy(float),
                "in_hsa00190": [g in ox for g in genes],
                "in_16": [g in SIXTEEN for g in genes],
                "in_19": [g in NINETEEN for g in genes],
                **{sp: M.loc[genes, sp].to_numpy(float) for sp in SPECIES_ORDER},
            }
        ).sort_values("r")
        fit["rank"] = np.arange(1, len(fit) + 1)
        fit.to_csv(OUT / "per_gene_mass_fit.tsv", sep="\t", index=False)
        rows = []
        for t in thr:
            obs = int((np.abs(r_obs) >= t).sum())
            nc = np.array(null_counts[t], float)
            pn = float((nc >= obs).mean())
            rows.append(
                dict(
                    порог=t,
                    наблюдено=obs,
                    медиана_нуля=float(np.median(nc)),
                    перц95_нуля=float(np.percentile(nc, 95)),
                    p_перестановки=pn,
                )
            )
        pd.DataFrame(rows).to_csv(OUT / "threshold_vs_null.tsv", sep="\t", index=False)
        print(f"\nгенов {len(genes)}  p<0.05 {(p_obs < 0.05).sum()}  BH<0.05 {(bh < 0.05).sum()}")
        print("перестановка |r|:")
        for rec in rows:
            print(
                f"  |r|≥{rec['порог']}  наблюдано {rec['наблюдено']}  нуль {rec['медиана_нуля']:.0f} (95% {rec['перц95_нуля']:.0f})  p={rec['p_перестановки']:.3f}"
            )
        print("\n20 самых отрицательных r:")
        show = fit.head(20)
        for r in show.itertuples():
            flag = []
            if r.in_16:
                flag.append("16")
            elif r.in_hsa00190:
                flag.append("00190")
            if r.frac_foreign >= 0.5:
                flag.append("чужое_имя")
            print(
                f"  {int(r.rank):4d} {r.gene:12s} r={r.r:+.3f} p={r.p:.4f} bh={r.bh:.3f}  {(','.join(flag) if flag else '')}"
            )
        print("\nгде 16:")
        sub = fit[fit.in_16]
        for r in sub.itertuples():
            print(
                f"  rank {int(r.rank):4d} {r.gene:10s} r={r.r:+.3f} p={r.p:.4f} bh={r.bh:.3f}  чужое_имя {r.frac_foreign:.2f}"
            )
        y16 = np.log10(M.reindex(SIXTEEN)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
        rr = stats.linregress(LOGM, y16)
        print(f"\n16 как набор: r={rr.rvalue:+.3f} p={rr.pvalue:.4f}")
        print(f"wrote {OUT}")

    return SimpleNamespace(**locals())


def _epi_own_promoters(context):
    """Epiblast promoters with the same 2 kb rule, without the fibroblast filter.

    Script 80 kept a peak only when fibroblast CAGE was above the threshold and
    then read epiblast on that peak. Here a peak enters if its own epiblast mean
    is above the threshold. The catchment is script 75: chicken-protein CDS, same
    strand, closer 5' end, within 2 kb. Missing gene is dropped. y = log10(sum).

    Source: analyses/avian/scripts/98_epi_own_promoters.py at d230e56ab942."""
    import sys
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats
    import importlib

    s74 = context.unit("avian.74_nineteen_via_cage")
    s75 = context.unit("avian.75_honest_gene_mass_scan")
    DATA = context.unit("avian.oxphos_thread_lib").DATA
    RES = context.unit("avian.oxphos_thread_lib").RES
    SPECIES = context.unit("avian.oxphos_thread_lib").SPECIES
    SPECIES_ORDER = context.unit("avian.oxphos_thread_lib").SPECIES_ORDER
    annotate_table = context.unit("avian.oxphos_thread_lib").annotate_table
    lib_cols = context.unit("avian.oxphos_thread_lib").lib_cols
    load_geneset = context.unit("avian.oxphos_thread_lib").load_geneset
    normalize_symbol = context.unit("avian.oxphos_thread_lib").normalize_symbol
    peak_to_symbol_de = context.unit("avian.oxphos_thread_lib").peak_to_symbol_de
    OUT = RES / "epi_own_promoters"
    ORDER = SPECIES_ORDER
    LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in ORDER], float))
    SIXTEEN = s75.SIXTEEN

    def load_peaks_epi(species: str) -> pd.DataFrame:
        d = pd.read_csv(s74.s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
        d = annotate_table(
            d,
            peak_to_symbol_de(
                Path(str(context.resolve_path("/path/to/storage/deviatiiarov/chicken/bird_sizes")))
                / SPECIES[species]["de"]
            ),
        )
        ecols = lib_cols(d, "epi")
        fcols = lib_cols(d, "fib")
        enum = d[ecols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        enum = enum.where(enum >= s74.THR, 0.0)
        fnum = d[fcols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        fnum = fnum.where(fnum >= s74.THR, 0.0)
        d = d.assign(
            mid=(d["start"] + d["end"]) / 2.0,
            fib=enum.mean(axis=1),
            had_fib=fnum.mean(axis=1) > 0,
            seq=d["seqnames"].map(s74.contig),
        )
        return d[d["fib"] > 0].copy()

    def fit_sum(wide: pd.DataFrame, genes: list[str]) -> tuple[int, float, float, list[str]]:
        hit = []
        for g in genes:
            if g not in wide.index:
                continue
            y = wide.loc[g, ORDER].to_numpy(float)
            if np.all(np.isfinite(y)) and np.all(y > 0):
                hit.append(g)
        if len(hit) < 3:
            return (len(hit), float("nan"), float("nan"), hit)
        total = np.log10(wide.loc[hit, ORDER].sum(axis=0).to_numpy(float))
        r, p = stats.pearsonr(LOGM, total)
        return (len(hit), float(r), float(p), hit)

    def main() -> None:
        OUT.mkdir(parents=True, exist_ok=True)
        ox = [normalize_symbol(g) for g in load_geneset("kegg_00190")]
        ox = [g for g in dict.fromkeys(ox) if g]
        wide = {}
        only_epi_peaks = {}
        for sp in ORDER:
            print(f"{sp}", flush=True)
            loci = s75.all_chicken_loci(s74.ASM / sp / "reanno_named.gff")
            peaks = load_peaks_epi(sp)
            only_epi_peaks[sp] = int((~peaks["had_fib"]).sum())
            print(
                f"  epi peaks {len(peaks)}  of which no fibroblast {only_epi_peaks[sp]}  loci {len(loci)}",
                flush=True,
            )
            assigned = s75.assign_species(peaks, loci)
            wide[sp] = (
                assigned.set_index("gene")["tpm"] if len(assigned) else pd.Series(dtype=float)
            )
            print(f"  genes {len(assigned)}", flush=True)
        M = pd.DataFrame(wide).reindex(columns=ORDER)
        M.to_csv(OUT / "M_epi_own_nearest2kb.tsv", sep="\t")
        old = pd.read_csv(RES / "story_epi_traits" / "M_epi_nearest2kb.tsv", sep="\t", index_col=0)
        old.index = old.index.astype(str).str.upper()
        M.index = M.index.astype(str).str.upper()
        both = [g for g in ox if g in M.index and g in old.index]
        same = 0
        changed = []
        for g in both:
            a = M.loc[g, ORDER].to_numpy(float)
            b = old.loc[g, ORDER].to_numpy(float)
            if np.allclose(np.nan_to_num(a), np.nan_to_num(b), rtol=1e-06, atol=0.0001):
                same += 1
            else:
                changed.append(g)
        new_genes = sorted(set(M.index) - set(old.index))
        new_ox = [g for g in ox if g in M.index and g not in old.index]
        print(f"\nOXPHOS in both matrices {len(both)} identical {same} changed {len(changed)}")
        print("changed", " ".join(changed[:40]))
        print(
            f"genes only in the epiblast call {len(new_genes)}  of them OXPHOS {len(new_ox)} {' '.join(new_ox)}"
        )
        print("epi-only peaks", only_epi_peaks)
        for label, genes in (("16", SIXTEEN), ("00190", ox)):
            n, r, p, hit = fit_sum(M, genes)
            n0, r0, p0, hit0 = fit_sum(old, genes)
            print(
                f"SUM {label:8s} own-epi n={n:3d} r={r:+.3f} p={p:.3g}   old-epi n={n0:3d} r={r0:+.3f} p={p0:.3g}"
            )
            gained = sorted(set(hit) - set(hit0))
            lost = sorted(set(hit0) - set(hit))
            if gained or lost:
                print(f"  membership gained {gained} lost {lost}")
        n, r, p, hit = fit_sum(M, ox)
        rows = []
        for g in hit:
            y = np.log10(M.loc[g, ORDER].to_numpy(float))
            rr, pp = stats.pearsonr(LOGM, y)
            rows.append((g, rr, pp, g in SIXTEEN))
        tab = pd.DataFrame(rows, columns=["gene", "r", "p", "in16"]).sort_values("p")
        ps = tab["p"].to_numpy(float)
        order = np.argsort(ps)
        q = np.empty(len(ps))
        prev = 1.0
        m = len(ps)
        for rank, i in enumerate(order[::-1]):
            k = m - rank
            prev = min(prev, ps[i] * m / k)
            q[i] = prev
        tab["bh"] = q
        tab.to_csv(OUT / "oxphos_epi_own_per_gene.tsv", sep="\t", index=False)
        nom = tab[tab.p < 0.05]
        print(
            f"nominal p<0.05 among {len(tab)} all-7 genes: {len(nom)}  BH<0.05 {(tab.bh < 0.05).sum()}"
        )
        print(nom.to_string(index=False))
        print(f"wrote {OUT}")

    return SimpleNamespace(**locals())


def _oxphos_thread_lib(context):
    """Shared helpers for OXPHOS species-thread pipeline.

    Source: analyses/avian/scripts/oxphos_thread_lib.py at d230e56ab942."""
    import re
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from scipy import stats

    ROOT = Path(str(context.resolve_path("/path/to/user-home/chicken/kleiber_birds_gemcat")))
    DATA = ROOT / "data"
    EXP = Path(
        str(
            context.resolve_path(
                "/path/to/storage/kleiber_birds_gemcat/dropbox_20260808/unpacked/exp_tab"
            )
        )
    )
    RUSLAN = Path(str(context.resolve_path("/path/to/storage/deviatiiarov/chicken/bird_sizes")))
    RES = ROOT / "results" / "oxphos_species_thread"
    FIG = ROOT / "figures"
    GENOME_ROOT = Path(
        str(context.resolve_path("/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt"))
    )
    FIMO = Path(str(context.resolve_path("/path/to/user-home/miniconda3/envs/meme/bin/fimo")))
    SALMON = ROOT / "results" / "salmon_all7_chicken" / "gene_tpm_by_species.tsv"
    MTDNA = {
        "ND1",
        "ND2",
        "ND3",
        "ND4",
        "ND4L",
        "ND5",
        "ND6",
        "CYTB",
        "COX1",
        "COX2",
        "COX3",
        "ATP6",
        "ATP8",
    }
    PREFERRED_ALIAS = {
        "ATP5A1": "ATP5F1A",
        "ATP5B": "ATP5F1B",
        "ATP5C1": "ATP5F1C",
        "ATP5D": "ATP5F1D",
        "ATP5E": "ATP5F1E",
        "ATP5F1": "ATP5PB",
        "ATP5G1": "ATP5MC1",
        "ATP5G3": "ATP5MC3",
        "ATP5H": "ATP5PD",
        "ATP5J": "ATP5PF",
        "ATP5J2": "ATP5MF",
        "ATP5L": "ATP5MG",
        "ATP5O": "ATP5PO",
        "COX4": "COX4I1",
        "PGC1A": "PPARGC1A",
        "NRF2": "NFE2L2",
        "ERRA": "ESRRA",
    }
    SPECIES = {
        "Zebra": {
            "file": "Zebra_fib_epi_exp.txt",
            "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 11.93,
            "salmon": "ZebraFinch",
            "geo_rna": "GSE213253_Zebrafinch_RNAseq.xlsx",
        },
        "Quail": {
            "file": "Quail_fib_epi_exp.txt",
            "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 96.6,
            "salmon": "Quail",
            "geo_rna": "GSE213253_Quail_RNAseq.xlsx",
        },
        "Chicken": {
            "file": "Chicken_fib_epi_exp.txt",
            "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 912.5,
            "salmon": "Chicken",
            "geo_rna": "GSE213253_Chicken_RNAseq.xlsx",
        },
        "Duck": {
            "file": "Duck_fib_epi_exp.txt",
            "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 1121.0,
            "salmon": "Duck",
            "geo_rna": "GSE213253_Duck_RNAseq.xlsx",
        },
        "Peacock": {
            "file": "Peacock_fib_epi_exp.txt",
            "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 4093.75,
            "salmon": "Peacock",
            "geo_rna": "GSE213253_Peacock_RNAseq.xlsx",
        },
        "Turkey": {
            "file": "Turkey_fib_epi_exp.txt",
            "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 5811.0,
            "salmon": "Turkey",
            "geo_rna": "GSE213253_Turkey_RNAseq.xlsx",
        },
        "Emu": {
            "file": "Emu_fib_epi_exp.txt",
            "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx",
            "mass_g": 36200.0,
            "salmon": "Emu",
            "geo_rna": "GSE213253_Emu_RNAseq.xlsx",
        },
    }
    SPECIES_ORDER = list(SPECIES.keys())
    MATCHED_HIGH_CODES = {
        "Chicken": ["C10", "C12"],
        "Quail": ["C26", "C28"],
        "Duck": ["C32", "C34"],
        "Peacock": ["C06", "C04"],
        "Turkey": ["C16", "C18"],
        "Emu": ["C22", "C24"],
        "Zebra": ["C39", "C40"],
    }
    DRIVERS = ["NRF1", "GABPA", "ESRRA", "YY1"]
    DRIVERS_EXTRA = ["TFAM", "POLRMT", "SIRT1", "PPRC1"]
    CHICKEN_NC_TO_CHR = {
        "NC_006088.5": "chr1",
        "NC_006089.5": "chr2",
        "NC_006090.5": "chr3",
        "NC_006091.5": "chr4",
        "NC_006092.5": "chr5",
        "NC_006093.5": "chr6",
        "NC_006094.5": "chr7",
        "NC_006095.5": "chr8",
        "NC_006096.5": "chr9",
        "NC_006097.5": "chr10",
        "NC_006098.5": "chr11",
        "NC_006099.5": "chr12",
        "NC_006100.5": "chr13",
        "NC_006101.5": "chr14",
        "NC_006102.5": "chr15",
        "NC_006103.5": "chr16",
        "NC_006104.5": "chr17",
        "NC_006105.5": "chr18",
        "NC_006106.5": "chr19",
        "NC_006107.5": "chr20",
        "NC_006108.5": "chr21",
        "NC_006109.5": "chr22",
        "NC_006110.5": "chr23",
        "NC_006111.5": "chr24",
        "NC_006112.5": "chr25",
        "NC_006113.5": "chr26",
        "NC_006114.5": "chr27",
        "NC_006115.5": "chr28",
        "NC_006119.5": "chr32",
        "NC_006126.5": "chr33",
        "NC_006127.5": "chrW",
        "NC_006088": "chr1",
        "NC_001323.1": "chrM",
    }

    def load_geneset(name: str) -> set[str]:
        path = DATA / "genesets" / f"{name}.txt"
        return {ln.strip().upper() for ln in path.read_text().splitlines() if ln.strip()}

    def nuclear_oxphos() -> set[str]:
        return {g for g in load_geneset("OXPHOS") if g not in MTDNA}

    def is_augustus_id(s: str) -> bool:
        u = str(s).upper().strip()
        return bool(re.match("^G\\d+(\\.T\\d+)?$", u)) or u in {"NAN", "NA", "", "NONE", "."}

    def normalize_symbol(sym: str) -> str:
        s = str(sym).upper().strip()
        if not s or s in {"NAN", "NA", "NONE", "."} or is_augustus_id(s):
            return ""
        return PREFERRED_ALIAS.get(s, s)

    def peak_to_symbol_de(de_path: Path) -> pd.Series:
        de = pd.read_excel(de_path, usecols=["id", "genes", "geneId"])
        for col in ("genes", "geneId"):
            de[col] = de[col].astype(str).str.upper()
        de["sym"] = de["genes"].where(~de["genes"].isin(["NAN", "NA", ""]), de["geneId"])
        de["sym"] = de["sym"].map(normalize_symbol)
        de = de[de["sym"] != ""].drop_duplicates("id")
        return de.set_index("id")["sym"]

    def annotate_table(df: pd.DataFrame, id2sym: pd.Series) -> pd.DataFrame:
        out = df.copy()
        gid = out["geneId"].astype(str)
        cg = (
            out["chicken_genes"].astype(str)
            if "chicken_genes" in out.columns
            else pd.Series([""] * len(out))
        )
        sym = []
        for g, c in zip(gid, cg):
            cu = str(c).upper().strip()
            gu = str(g).upper().strip()
            pick = ""
            if cu and (not is_augustus_id(cu)):
                pick = cu
            elif not is_augustus_id(gu):
                pick = gu
            sym.append(normalize_symbol(pick))
        out["symbol"] = sym
        fill = out["symbol"].eq("")
        out.loc[fill, "symbol"] = out.loc[fill, "id"].map(id2sym).fillna("").map(normalize_symbol)
        return out

    def lib_cols(df: pd.DataFrame, tissue: str) -> list[str]:
        pref = "norm.C" if tissue == "fib" else "norm.E"
        return [c for c in df.columns if str(c).startswith(pref)]

    def filter_matched_cols(cols: list[str], codes: list[str]) -> list[str]:
        keep = []
        for c in cols:
            for code in codes:
                if re.search(f"(?:^|[._]){code}(?:[._]|$)", str(c), flags=re.I) or code in str(c):
                    keep.append(c)
                    break
        return keep

    def gene_by_lib_matrix(
        df: pd.DataFrame,
        genes: set[str],
        tissue: str,
        thr: float = 1.0,
        cols: list[str] | None = None,
    ) -> pd.DataFrame:
        """Rows=genes, cols=libs; total same-strand promoter TPM per lib."""
        cols = cols or lib_cols(df, tissue)
        if not cols:
            return pd.DataFrame()
        sub = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].isin(genes)]
        if sub.empty:
            return pd.DataFrame(columns=cols)
        sub = sub.copy()
        mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sub = sub.assign(**{c: mat[c] for c in cols})
        sub["_tot"] = mat.sum(axis=1)
        strand_tot = sub.groupby(["symbol", "strand"], observed=True)["_tot"].sum().reset_index()
        maj = strand_tot.sort_values("_tot", ascending=False).drop_duplicates("symbol")
        maj_map = maj.set_index("symbol")["strand"]
        sub = sub[sub["strand"].values == sub["symbol"].map(maj_map).values]
        for c in cols:
            sub[c] = sub[c].where(sub[c] >= thr, 0.0)
        out = sub.groupby("symbol", observed=True)[cols].sum()
        return out

    def median_pairwise_spearman(mat: pd.DataFrame) -> float:
        """mat: genes x samples. Return median pairwise Spearman among genes with variance."""
        m = mat.loc[mat.std(axis=1) > 1e-12]
        if len(m) < 3:
            return float("nan")
        ranks = m.rank(axis=1, method="average").to_numpy(dtype=float)
        c = np.corrcoef(ranks)
        if c.ndim < 2:
            return float("nan")
        iu = np.triu_indices(c.shape[0], k=1)
        vals = c[iu]
        vals = vals[np.isfinite(vals)]
        return float(np.median(vals)) if len(vals) else float("nan")

    def ols_logmass(y: pd.Series, masses: pd.Series) -> dict:
        x = np.log10(masses.astype(float).values)
        yy = y.astype(float).values
        mask = np.isfinite(x) & np.isfinite(yy)
        if mask.sum() < 4:
            return {
                "n": int(mask.sum()),
                "slope": np.nan,
                "r": np.nan,
                "p": np.nan,
                "se": np.nan,
                "intercept": np.nan,
            }
        sl, ic, r, p, se = stats.linregress(x[mask], yy[mask])
        return {
            "n": int(mask.sum()),
            "slope": float(sl),
            "r": float(r),
            "p": float(p),
            "se": float(se),
            "intercept": float(ic),
        }

    def residual_after_covariate(y: pd.Series, cov: pd.Series) -> pd.Series:
        x = cov.astype(float).values
        yy = y.astype(float).values
        mask = np.isfinite(x) & np.isfinite(yy)
        out = pd.Series(np.nan, index=y.index)
        if mask.sum() < 3:
            return out
        sl, ic, *_ = stats.linregress(x[mask], yy[mask])
        pred = ic + sl * cov.astype(float)
        out = y.astype(float) - pred
        return out

    def spearman_logmass(y: pd.Series, masses: pd.Series) -> dict:
        x = np.log10(masses.astype(float).values)
        yy = y.astype(float).values
        mask = np.isfinite(x) & np.isfinite(yy)
        if mask.sum() < 4:
            return {"n": int(mask.sum()), "rho": np.nan, "p": np.nan}
        rho, p = stats.spearmanr(x[mask], yy[mask])
        return {"n": int(mask.sum()), "rho": float(rho), "p": float(p)}

    def de_augustus_to_symbol(de_path: Path) -> dict[str, str]:
        """Map Augustus transcript/gene IDs (g9190.t1, g9190) → symbol via Ruslan DE."""
        de = pd.read_excel(de_path, usecols=["geneId", "transcriptId", "genes"])
        out: dict[str, str] = {}
        for _, r in de.iterrows():
            sym = normalize_symbol(r["genes"]) or normalize_symbol(r["geneId"])
            if not sym:
                continue
            tid = str(r["transcriptId"]).strip()
            gid = str(r["geneId"]).strip()
            keys = {tid, gid}
            if "." in tid:
                keys.add(tid.split(".", 1)[0])
            for key in keys:
                if key and key.lower() not in {"nan", "none", ""}:
                    out[key.upper()] = sym
        return out

    def load_geo_rnaseq(species: str) -> dict:
        """Official GSE213253 RNA-seq: symbol-indexed gene×lib norm matrix + library totals.

        Peacock/Turkey/Emu use Augustus IDs; mapped through DE transcriptId.
        GEO 'norm' totals are NOT a shared 1e6 TPM scale across species — use
        within-table fractions, not raw interspecific means, as the primary metric.
        """
        meta = SPECIES[species]
        df = pd.read_excel(DATA / meta["geo_rna"])
        norms = [c for c in df.columns if str(c).startswith("norm.")]
        raw_id = df["gene"].astype(str)
        mat = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        lib_total = mat.sum(axis=0)
        if raw_id.map(is_augustus_id).mean() < 0.5:
            sym = raw_id.map(normalize_symbol)
            mapped_how = "native_symbol"
        else:
            id2 = de_augustus_to_symbol(RUSLAN / meta["de"])
            keys = raw_id.str.upper()
            keys2 = keys.str.replace("\\.T\\d+$", "", regex=True)
            mapped = keys.map(id2)
            miss = mapped.isna() | (mapped == "")
            mapped = mapped.where(~miss, keys2.map(id2))
            sym = mapped.fillna("")
            mapped_how = "de_transcriptId"
        keep = sym != ""
        tmp = mat.loc[keep].copy()
        tmp["_sym"] = sym[keep].values
        gene_mat = tmp.groupby("_sym", observed=True)[norms].mean()
        return {
            "gene_mat": gene_mat,
            "lib_total": lib_total,
            "mapped_how": mapped_how,
            "n_raw": int(len(df)),
            "n_mapped_rows": int(keep.sum()),
            "n_symbols": int(gene_mat.shape[0]),
            "norm_cols": norms,
        }

    return SimpleNamespace(**locals())
