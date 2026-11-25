"""Existing pySCENIC preparation, inference, pruning and sensitivity.

Calculation bodies are retained from the source commit; configuration and execution scope are explicit.
"""

from __future__ import annotations
from types import SimpleNamespace

_GRN_STATE = {}


def _infer_target(i):
    return _GRN_STATE["infer_partial_network"](
        regressor_type="GBM",
        regressor_kwargs=_GRN_STATE["SGBM_KWARGS"],
        tf_matrix=_GRN_STATE["TFX"],
        tf_matrix_gene_names=_GRN_STATE["TFN"],
        target_gene_name=_GRN_STATE["GN"][i],
        target_gene_expression=_GRN_STATE["X"][:, i],
        include_meta=False,
        early_stop_window_length=_GRN_STATE["EARLY_STOP_WINDOW_LENGTH"],
        seed=_GRN_STATE["SEED"],
    )


def _prepare(context):
    """Existing analysis stage; calculation statements retained from the source.

    Source: analyses/mammal_scenic/code/01_prepare.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import os, json, hashlib, sys, re, gzip, math, importlib.metadata as md

        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        import numpy as np, pandas as pd

        R = context.run_dir("scenic")
        P = R.parent
        H = Path(str(context.resolve_path("/path/to/user-home/chicken/kleiber_birds_gemcat")))
        inputs = []

        def sha(p):
            h = hashlib.sha256()
            with Path(p).open("rb") as f:
                for b in iter(lambda: f.read(1 << 22), b""):
                    h.update(b)
            return h.hexdigest()

        def use(p):
            p = Path(p)
            inputs.append(dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size))
            return p

        for d in ["inputs", "config", "reference", "results", "audit", "logs", "figures"]:
            (R / d).mkdir(exist_ok=True)
        plan = dict(
            run="ma2016_pyscenic_2026-10-23_v1",
            authorization="User: Ну так в чем проблема это сделать?",
            experiment="Ma2016 primary mammalian fibroblasts, same source as frozen WGCNA",
            technical_pairs=[
                ["RedSquirrel.1A", "RedSquirrel.1B"],
                ["DeerMouse.2A", "DeerMouse.2B"],
            ],
            normalization="sum raw counts within explicit technical pairs; CPM from all source gene counts; log2(CPM+1)",
            filter="CPM >=1 in >=20% of26 collapsed sample-label groups, nonconstant; no body mass or WGCNA selection",
            orthology="unique one-to-one mouse-to-human mappings, unique human symbols",
            seeds=[20261023, 20261024, 20261025],
            primary_seed=20261023,
            workers=24,
            regressor="Official arboreto GRNBoost2 defaults",
            cisTarget="Official custom database builder, cluster-buster, JASPAR2024 one-per-symbol816 motifs; standard downloads and official mirror refused connection",
            reference_windows={"500up100down": [-500, 100], "10kb_each_side": [-10000, 10000]},
            reference="Human GRCh38.99 longest exonic transcript on primary chromosomes; protein-coding genes plus expression-input symbols; unique gene symbols",
            pruning=dict(
                rank_threshold=1500,
                auc_threshold=0.05,
                nes_threshold=3.0,
                annotation="direct JASPAR TF label, heterodimer components expanded, HGNC unambiguous aliases normalized",
            ),
            overlap="all inferred primary regulons x4 fixed WGCNA modules; measured/ranked WGCNA-network gene background; one-sided hypergeometric, BH across all tested module-regulon pairs",
            activity="AUCell continuous per donor-labelled group; species means for mass correlation; BH across primary regulons; exploratory same-study analysis",
            stability="two additional predeclared full GRN/prune runs, target Jaccard and TF recovery; not independent data",
            wgcna_rebuilt=False,
            slides_edited=False,
        )
        (R / "config/plan.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
        raw = pd.read_excel(
            use(H / "data/file_panels/ma2016/elife19130_supp1.xlsx"),
            sheet_name="A. raw counts",
            index_col=0,
        )
        meta0 = pd.read_excel(
            use(H / "data/file_panels/ma2016/elife19130_fig1.xlsx"),
            sheet_name="B. RNAseq & Mapping to ortholog",
        )
        assert raw.index.is_unique and raw.columns.is_unique and (raw.shape == (16816, 28))
        map0 = pd.read_csv(use(H / "data/ensembl_orthologs/mouse_to_human.tsv"), sep="\t").dropna()
        map0 = map0[map0["Human homology type"].eq("ortholog_one2one")]
        map0 = map0[["Gene stable ID", "Human gene name"]].drop_duplicates()
        map0 = map0[
            ~map0["Gene stable ID"].duplicated(keep=False)
            & ~map0["Human gene name"].duplicated(keep=False)
        ]
        map0["Human gene name"] = map0["Human gene name"].str.upper()
        hgnc = pd.read_csv(
            use(P / "sizes5_execution_2026-10-16_v1/sources/hgnc_complete_set.txt"),
            sep="\t",
            dtype=str,
        ).fillna("")
        official = set(hgnc.symbol.str.upper())
        aliases = {}
        for row in hgnc.itertuples():
            for col in ["alias_symbol", "prev_symbol"]:
                for a in str(getattr(row, col, "")).upper().split("|"):
                    if a:
                        aliases.setdefault(a, set()).add(row.symbol.upper())
        rename = {a: next(iter(s)) for a, s in aliases.items() if len(s) == 1 and a not in official}
        map0["Human gene name"] = map0["Human gene name"].replace(rename)
        map0 = map0[~map0["Human gene name"].duplicated(keep=False)]
        meta = meta0.set_index("Sample ID")
        cols = {s: s for s in raw.columns}
        for aa, bb in plan["technical_pairs"]:
            cols[bb] = aa[:-1]
            cols[aa] = aa[:-1]
        X = raw.T.groupby(pd.Index([cols[s] for s in raw.columns]), sort=False).sum().T
        lib = X.sum(axis=0)
        CPM = X.div(lib, axis=1) * 1000000.0
        m = map0.set_index("Gene stable ID")["Human gene name"]
        ids = [g for g in X.index if g in m]
        CPM = CPM.loc[ids]
        CPM.index = m.loc[ids].tolist()
        assert CPM.index.is_unique
        keep = (CPM.ge(1).sum(axis=1) >= math.ceil(0.2 * CPM.shape[1])) & (
            np.log2(CPM + 1).std(axis=1) > 1e-08
        )
        CPM = CPM.loc[keep]
        E = np.log2(CPM + 1).T
        assert len(E) == 26 and E.columns.is_unique
        E.to_pickle(R / "inputs/expression.pkl")
        E.to_csv(R / "inputs/expression.tsv.gz", sep="\t")
        CPM.to_csv(R / "inputs/CPM.tsv.gz", sep="\t")
        map0.to_csv(R / "inputs/orthology_used.tsv", sep="\t", index=False)
        tfs = sorted(set(use(R / "reference/hs_hgnc_tfs.txt").read_text().splitlines()))
        tfs = sorted({rename.get(t.upper(), t.upper()) for t in tfs} & set(E.columns))
        (R / "inputs/tf_names.txt").write_text("\n".join(tfs) + "\n")
        rows = []
        label = {
            "Shrew": "shrew",
            "BigBrownBat": "big brown bat",
            "LittleBrownBat": "little brown bat",
            "GuineaPig": "guinea pig",
            "Porcupine": "porcupine",
            "Chinchilla": "chinchilla",
            "Chipmunk": "chipmunk",
            "RedSquirrel": "red squirrel",
            "FoxSquirrel": "fox squirrel",
            "Beaver": "beaver",
            "AfricanGrassRat": "grass rat",
            "MeadowVole": "meadow vole",
            "CottonRat": "cotton rat",
            "WhiteFootedMouse": "white-footed mouse",
            "DeerMouse": "deer mouse",
        }
        for s in E.index:
            original = [k for k, v in cols.items() if v == s]
            sp = s.split(".")[0]
            rows.append(
                dict(
                    sample=s,
                    species=label[sp],
                    source_columns=";".join(original),
                    technical_libraries=len(original),
                    source_library_sum=int(lib[s]),
                )
            )
        T = pd.DataFrame(rows)
        T.to_csv(R / "inputs/samples.tsv", sep="\t", index=False)
        assert T.species.nunique() == 15
        memb = pd.read_csv(
            use(P / "growth_regulators_2026-10-21_v1/results/wgcna/ma2016_full__raw/members.tsv"),
            sep="\t",
        )
        memb.to_csv(R / "inputs/frozen_WGCNA_members.tsv", sep="\t", index=False)
        gtf = use(P / "sizes5_execution_2026-10-16_v1/sources/Homo_sapiens.GRCh38.99.gtf.gz")
        genome = use(
            Path(
                str(context.resolve_path("/path/to/storage/human_cs6_embryo/track2/ref/GRCh38.fa"))
            )
        )
        use(Path(str(genome) + ".fai"))
        genes = {}
        trans = {}
        exlen = {}
        biotype = {}
        with gzip.open(gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip().split("\t")
                if len(c) < 9 or c[2] not in ["gene", "transcript", "exon"]:
                    continue
                d = dict(re.findall('(\\w+) "([^"]+)"', c[8]))
                name = d.get("gene_name", "").upper()
                name = rename.get(name, name)
                gid = d.get("gene_id", "")
                tid = d.get("transcript_id")
                if not name:
                    continue
                genes.setdefault(name, set()).add(gid)
                if c[2] == "gene":
                    biotype[name] = d.get("gene_biotype", "")
                if c[2] == "transcript" and tid:
                    trans[tid] = (
                        name,
                        c[0] if c[0].startswith("chr") else "chr" + c[0],
                        int(c[3]),
                        int(c[4]),
                        c[6],
                    )
                if c[2] == "exon" and tid:
                    exlen[tid] = exlen.get(tid, 0) + int(c[4]) - int(c[3]) + 1
        best = {}
        for tid, (name, chrom, start, end, strand) in trans.items():
            if len(genes[name]) != 1 or (
                biotype.get(name) != "protein_coding" and name not in E.columns
            ):
                continue
            if not re.match("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chrom):
                continue
            length = exlen.get(tid, 0)
            v = (length, tid, chrom, start, end, strand)
            if name not in best or (length, tid) > (best[name][0], best[name][1]):
                best[name] = v
        fai = {
            c[0]: list(map(int, c[1:5]))
            for c in (
                line.split("\t") for line in Path(str(genome) + ".fai").read_text().splitlines()
            )
        }
        comp = str.maketrans("ACGTN", "TGCAN")
        genomef = genome.open("rb")
        regions = []
        for tag, (up, down) in plan["reference_windows"].items():
            n = 0
            with (R / "reference" / f"human_{tag}.fa").open("w") as out:
                for name, (length, tid, chrom, start, end, strand) in sorted(best.items()):
                    lo, hi = (
                        (start - 1 + up, start - 1 + down)
                        if strand == "+"
                        else (end - down, end - up)
                    )
                    if chrom not in fai or lo < 0 or hi > fai[chrom][0]:
                        continue
                    size, offset, bases, width = fai[chrom]
                    p0 = offset + lo // bases * width + lo % bases
                    p1 = offset + (hi - 1) // bases * width + (hi - 1) % bases + 1
                    genomef.seek(p0)
                    seq = (
                        genomef.read(p1 - p0)
                        .replace(b"\n", b"")
                        .replace(b"\r", b"")
                        .decode()
                        .upper()
                    )
                    assert len(seq) == down - up
                    if seq.count("N") > 0.2 * len(seq):
                        continue
                    if strand == "-":
                        seq = seq.translate(comp)[::-1]
                    out.write(">" + name + "\n" + seq + "\n")
                    n += 1
                    regions.append(
                        dict(
                            window=tag,
                            gene=name,
                            transcript=tid,
                            chrom=chrom,
                            start0=lo,
                            end=hi,
                            strand=strand,
                            N_fraction=seq.count("N") / len(seq),
                        )
                    )
            print("FASTA", tag, n, flush=True)
        genomef.close()
        pd.DataFrame(regions).to_csv(R / "reference/promoter_regions.tsv.gz", sep="\t", index=False)
        meme = use(
            P / "size_program_regulons_2026-10-21_v1/reference/JASPAR2024_one_per_symbol.meme"
        )
        lines = meme.read_text().splitlines()
        motifs = []
        annot = []
        (R / "reference/motifs_cb").mkdir(exist_ok=True)
        i = 0
        while i < len(lines):
            if lines[i].startswith("MOTIF "):
                mid, tf = lines[i].split()[1:3]
                i += 1
                while i < len(lines) and "letter-probability matrix:" not in lines[i]:
                    i += 1
                length = int(re.search("w=\\s*(\\d+)", lines[i]).group(1))
                mat = np.array(
                    [[float(x) for x in lines[i + j + 1].split()] for j in range(length)]
                )
                counts = np.rint(mat * 10000).astype(int)
                (R / "reference/motifs_cb" / f"{mid}.cb").write_text(
                    ">"
                    + mid
                    + "\n"
                    + "\n".join(("\t".join(map(str, row)) for row in counts))
                    + "\n"
                )
                motifs.append(mid)
                for gene in tf.upper().split("::"):
                    gene = rename.get(gene, gene)
                    annot.append(
                        {
                            "#motif_id": mid,
                            "gene_name": gene,
                            "motif_similarity_qvalue": 0.0,
                            "orthologous_identity": 1.0,
                            "description": "JASPAR2024 direct TF-label assignment; heterodimer component if applicable; no native-mammal occupancy measurement",
                        }
                    )
            i += 1
        assert len(motifs) == 816
        pd.DataFrame(annot).to_csv(R / "reference/motif_annotations.tsv", sep="\t", index=False)
        (R / "reference/motifs.txt").write_text("\n".join(motifs) + "\n")
        out = dict(
            n_profiles=E.shape[0],
            n_species=T.species.nunique(),
            n_genes=E.shape[1],
            n_TFs=len(tfs),
            filter_min_samples=math.ceil(0.2 * E.shape[0]),
            E2F2_present="E2F2" in E.columns,
            E2F4_present="E2F4" in E.columns,
            motif_count=len(motifs),
            selected_WGCNA_coverage={
                mod: {
                    "total": int(memb.module.eq(mod).sum()),
                    "RNA_input": int(memb[memb.module.eq(mod)].feature.isin(E.columns).sum()),
                }
                for mod in ["orange", "violet", "cyan", "skyblue3"]
            },
            code_sha256=sha(Path(__file__)),
            inputs=inputs,
            versions={
                n: md.version(n)
                for n in [
                    "numpy",
                    "pandas",
                    "scipy",
                    "pyscenic",
                    "arboreto",
                    "ctxcore",
                    "scikit-learn",
                    "pyarrow",
                ]
            },
        )
        (R / "audit/preparation.json").write_text(json.dumps(out, indent=2))
        print(
            json.dumps({k: v for k, v in out.items() if k not in ["inputs"]}, indent=2), flush=True
        )

    return SimpleNamespace(**locals())


def _grn(context):
    """Existing functions and settings retained in their original local scope.

    Source: analyses/mammal_scenic/code/02_grn.py at d230e56ab942."""
    import pathlib, json, time, hashlib, inspect, multiprocessing as mp, os, argparse

    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np

    np.object = object
    import pandas as pd
    from arboreto.core import (
        infer_partial_network,
        SGBM_KWARGS,
        EARLY_STOP_WINDOW_LENGTH,
        to_tf_matrix,
    )
    from arboreto.algo import _prepare_input
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)
    R = context.run_dir("scenic")

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--seed", type=int, required=True)
        args = ap.parse_args(context.argv)
        SEED = args.seed
        out = R / "results" / str(SEED)
        out.mkdir(exist_ok=True)
        start = time.time()
        expr = pd.read_pickle(R / "inputs/expression.pkl")
        tfs = (R / "inputs/tf_names.txt").read_text().splitlines()
        X, GN, tfs = _prepare_input(expr, None, tfs)
        TFX, TFN = to_tf_matrix(X, GN, tfs)
        _GRN_STATE.update(
            dict(
                infer_partial_network=infer_partial_network,
                SGBM_KWARGS=SGBM_KWARGS,
                EARLY_STOP_WINDOW_LENGTH=EARLY_STOP_WINDOW_LENGTH,
                X=X,
                TFX=TFX,
                TFN=TFN,
                GN=GN,
                SEED=SEED,
            )
        )
        info = dict(
            seed=SEED,
            workers=24,
            n_profiles=X.shape[0],
            genes=X.shape[1],
            TFs=len(TFN),
            regressor_kwargs=SGBM_KWARGS,
            early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
            code_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
            input_sha256=hashlib.sha256((R / "inputs/expression.pkl").read_bytes()).hexdigest(),
            official_core="arboreto.core.infer_partial_network; GRNBoost2 defaults; fork-shared dispatcher",
        )
        (out / "GRN_before.json").write_text(json.dumps(info, indent=2))
        (R / "audit/arboreto_infer_partial_network.py").write_text(
            inspect.getsource(infer_partial_network)
        )
        print(json.dumps(info), flush=True)
        parts = []
        with mp.get_context("fork").Pool(24) as pool:
            for i, a in enumerate(
                pool.imap_unordered(_infer_target, range(len(GN)), chunksize=2), 1
            ):
                parts.append(a)
                if i % 1000 == 0 or i == len(GN):
                    print(
                        json.dumps(
                            dict(
                                finished_targets=i,
                                total=len(GN),
                                elapsed_seconds=time.time() - start,
                            )
                        ),
                        flush=True,
                    )
        adj = pd.concat(parts).sort_values(
            ["importance", "TF", "target"], ascending=[False, True, True], kind="stable"
        )
        assert not adj.duplicated(["TF", "target"]).any()
        adj.to_csv(out / "adjacencies.tsv.gz", sep="\t", index=False)
        (out / "GRN_summary.json").write_text(
            json.dumps(
                dict(
                    positive_edges=len(adj),
                    seconds=time.time() - start,
                    seed=SEED,
                    official_core_GRNBoost2=True,
                ),
                indent=2,
            )
        )
        print("GRN_DONE", SEED, len(adj), time.time() - start, flush=True)

    return SimpleNamespace(**locals())


def _build_db(context):
    """Existing analysis stage; calculation statements retained from the source.

    Source: analyses/mammal_scenic/code/03_build_db.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import subprocess, json, hashlib, time

        R = context.run_dir("scenic")
        P = R.parent / "planaria_pyscenic_2026-10-21_v1/reference"
        tool = P / "create_cistarget_tools/create_cistarget_motif_databases.py"
        py = str(
            context.resolve_path(
                "/path/to/user-home/miniconda3/envs/create_cistarget_databases/bin/python"
            )
        )
        sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
        for tag in ["500up100down", "10kb_each_side"]:
            t = time.time()
            prefix = R / "reference" / f"human_{tag}"
            fa = Path(str(prefix) + ".fa")
            args = [
                py,
                str(tool),
                "-f",
                str(fa),
                "-M",
                str(R / "reference/motifs_cb"),
                "-m",
                str(R / "reference/motifs.txt"),
                "-o",
                str(prefix),
                "-c",
                str(P / "cbust"),
                "-t",
                "24",
                "-s",
                "20261023",
            ]
            (R / "audit" / f"build_{tag}.json").write_text(
                json.dumps(
                    dict(
                        command=args,
                        tool_sha256=sha(tool),
                        cbust_sha256=sha(P / "cbust"),
                        fasta_sha256=sha(fa),
                        code_sha256=sha(Path(__file__)),
                    ),
                    indent=2,
                )
            )
            with (R / "logs" / f"build_{tag}.log").open("w") as f:
                subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)
            print("DATABASE_DONE", tag, time.time() - t, flush=True)

    return SimpleNamespace(**locals())


def _prune(context):
    """Existing functions and settings retained in their original local scope.

    Source: analyses/mammal_scenic/code/04_prune.py at d230e56ab942."""
    from pathlib import Path
    import os, json, pickle, time, hashlib, argparse

    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np

    np.object = object
    import pandas as pd
    from ctxcore.rnkdb import FeatherRankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)
    R = context.run_dir("scenic")

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--seed", type=int, required=True)
        args = ap.parse_args(context.argv)
        seed = args.seed
        O = R / "results" / str(seed)
        expr = pd.read_pickle(R / "inputs/expression.pkl")
        adj = pd.read_csv(O / "adjacencies.tsv.gz", sep="\t")
        t = time.time()
        modules = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
        pd.DataFrame(
            [
                dict(
                    module=i,
                    TF=m.transcription_factor,
                    genes=len(m),
                    context=";".join(sorted(m.context)),
                )
                for i, m in enumerate(modules)
            ]
        ).to_csv(O / "coexpression_modules.tsv", sep="\t", index=False)
        with (O / "coexpression_modules.pkl").open("wb") as f:
            pickle.dump(modules, f)
        print("MODULES", seed, len(modules), flush=True)
        dfs = []
        windows = {}
        summ = []
        for tag in ["500up100down", "10kb_each_side"]:
            p = R / "reference" / f"human_{tag}.regions_vs_motifs.rankings.feather"
            db = FeatherRankingDatabase(str(p), name="human_GRCh38_JASPAR2024_" + tag)
            info = dict(
                seed=seed,
                window=tag,
                ranking_genes=db.total_genes,
                ranking_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                rank_threshold=1500,
                auc_threshold=0.05,
                NES_threshold=3.0,
                num_workers=16,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            )
            (O / f"prune_{tag}_before.json").write_text(json.dumps(info, indent=2))
            df = prune2df(
                [db],
                modules,
                str(R / "reference/motif_annotations.tsv"),
                rank_threshold=1500,
                auc_threshold=0.05,
                nes_threshold=3.0,
                num_workers=16,
                client_or_address="custom_multiprocessing",
                module_chunksize=100,
                filter_for_annotation=True,
            )
            df.to_pickle(O / f"motif_enrichment_{tag}.pkl")
            df.to_csv(O / f"motif_enrichment_{tag}.csv")
            regs = df2regulons(df) if len(df) else []
            with (O / f"regulons_{tag}.pkl").open("wb") as f:
                pickle.dump(regs, f)
            windows[tag] = {r.name: set(r.gene2weight) - {r.transcription_factor} for r in regs}
            dfs.append(df)
            summ.append(dict(window=tag, regulons=len(regs), motif_enrichment_rows=len(df)))
            print("PRUNED", seed, tag, len(regs), flush=True)
        df = pd.concat(dfs)
        df.to_pickle(O / "motif_enrichment_union.pkl")
        regs = df2regulons(df) if len(df) else []
        with (O / "regulons.pkl").open("wb") as f:
            pickle.dump(regs, f)
        edges = []
        inv = []
        for r in regs:
            tf = r.transcription_factor
            targets = set(r.gene2weight) - {tf}
            inv.append(
                dict(
                    regulon=r.name,
                    TF=tf,
                    targets_nonself=len(targets),
                    windows=sum((r.name in v for v in windows.values())),
                    targets_both_windows=len(
                        set.intersection(*[v.get(r.name, set()) for v in windows.values()])
                    ),
                )
            )
            for g in sorted(targets):
                edges.append(
                    dict(
                        regulon=r.name,
                        TF=tf,
                        target=g,
                        importance=r.gene2weight[g],
                        window_count=sum((g in v.get(r.name, set()) for v in windows.values())),
                    )
                )
        pd.DataFrame(
            inv, columns=["regulon", "TF", "targets_nonself", "windows", "targets_both_windows"]
        ).to_csv(O / "regulon_inventory.tsv", sep="\t", index=False)
        pd.DataFrame(
            edges, columns=["regulon", "TF", "target", "importance", "window_count"]
        ).to_csv(O / "regulon_edges.tsv.gz", sep="\t", index=False)
        if regs:
            auc = aucell(expr, regs, auc_threshold=0.05, num_workers=8, seed=seed)
            auc.to_csv(O / "AUCell.tsv", sep="\t")
        (O / "SCENIC_summary.json").write_text(
            json.dumps(
                dict(
                    seed=seed,
                    coexpression_modules=len(modules),
                    windows=summ,
                    union_regulons=len(regs),
                    nonself_edges=len(edges),
                    seconds=time.time() - t,
                    expression="26 bulk donor-labelled groups; not single cells",
                    reference="custom human GRCh38+JASPAR2024; not native species occupancy",
                ),
                indent=2,
            )
        )
        print("SCENIC_DONE", seed, len(regs), len(edges), time.time() - t, flush=True)

    return SimpleNamespace(**locals())


def _prune_existing(context):
    """Existing functions and settings retained in their original local scope.

    Source: analyses/mammal_scenic/code/04_prune_existing.py at d230e56ab942."""
    from pathlib import Path
    import os, json, pickle, time, hashlib, argparse

    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    import numpy as np

    np.object = object
    import pandas as pd
    from ctxcore.rnkdb import FeatherRankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)
    R = context.run_dir("scenic")

    def sha(p):
        h = hashlib.sha256()
        with Path(p).open("rb") as f:
            for b in iter(lambda: f.read(1 << 22), b""):
                h.update(b)
        return h.hexdigest()

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--seed", type=int, required=True)
        args = ap.parse_args(context.argv)
        seed = args.seed
        O = R / "results" / str(seed)
        expr = pd.read_pickle(R / "inputs/expression.pkl")
        adj = pd.read_csv(O / "adjacencies.tsv.gz", sep="\t")
        t = time.time()
        cfg = json.loads((R / "config/existing_reference.json").read_text())
        modules = list(modules_from_adjacencies(adj, expr, rho_mask_dropouts=False))
        pd.DataFrame(
            [
                dict(
                    module=i,
                    TF=m.transcription_factor,
                    genes=len(m),
                    context=";".join(sorted(m.context)),
                )
                for i, m in enumerate(modules)
            ]
        ).to_csv(O / "coexpression_modules.tsv", sep="\t", index=False)
        with (O / "coexpression_modules.pkl").open("wb") as f:
            pickle.dump(modules, f)
        print("MODULES", seed, len(modules), flush=True)
        p = Path(cfg["database"])
        ann = Path(cfg["annotation"])
        assert sha(p) == cfg["database_sha256"]
        assert sha(ann) == cfg["annotation_sha256"]
        db = FeatherRankingDatabase(str(p), name="human_hg38_v10_clustered_500up100down_existing")
        info = dict(
            seed=seed,
            ranking_genes=db.total_genes,
            ranking_sha256=cfg["database_sha256"],
            annotation_sha256=cfg["annotation_sha256"],
            rank_threshold=1500,
            auc_threshold=0.05,
            NES_threshold=3.0,
            num_workers=16,
            code_sha256=sha(Path(__file__)),
            reference="Existing official hg38 full-tx 500bpUp100down v10-clustered; all5876motifs match HGNC snapshot m1e-5",
        )
        (O / "prune_before.json").write_text(json.dumps(info, indent=2))
        df = prune2df(
            [db],
            modules,
            str(ann),
            rank_threshold=1500,
            auc_threshold=0.05,
            nes_threshold=3.0,
            num_workers=16,
            client_or_address="custom_multiprocessing",
            module_chunksize=100,
            filter_for_annotation=True,
        )
        df.to_pickle(O / "motif_enrichment.pkl")
        df.to_csv(O / "motif_enrichment.csv")
        regs = df2regulons(df) if len(df) else []
        with (O / "regulons.pkl").open("wb") as f:
            pickle.dump(regs, f)
        edges = []
        inv = []
        for r in regs:
            tf = r.transcription_factor
            targets = set(r.gene2weight) - {tf}
            inv.append(dict(regulon=r.name, TF=tf, targets_nonself=len(targets)))
            for g in sorted(targets):
                edges.append(dict(regulon=r.name, TF=tf, target=g, importance=r.gene2weight[g]))
        pd.DataFrame(inv, columns=["regulon", "TF", "targets_nonself"]).to_csv(
            O / "regulon_inventory.tsv", sep="\t", index=False
        )
        pd.DataFrame(edges, columns=["regulon", "TF", "target", "importance"]).to_csv(
            O / "regulon_edges.tsv.gz", sep="\t", index=False
        )
        if regs:
            auc = aucell(expr, regs, auc_threshold=0.05, num_workers=8, seed=seed)
            auc.to_csv(O / "AUCell.tsv", sep="\t")
        (O / "SCENIC_summary.json").write_text(
            json.dumps(
                dict(
                    seed=seed,
                    coexpression_modules=len(modules),
                    motif_enrichment_rows=len(df),
                    union_regulons=len(regs),
                    nonself_edges=len(edges),
                    seconds=time.time() - t,
                    expression="26 bulk donor-labelled groups; not single cells",
                    reference="existing human hg38 v10-clustered short-promoter DB; not native species occupancy",
                    two_window_sensitivity=False,
                ),
                indent=2,
            )
        )
        print("SCENIC_DONE", seed, len(regs), len(edges), time.time() - t, flush=True)

    return SimpleNamespace(**locals())


def _analyze(context):
    """Existing analysis stage; calculation statements retained from the source.

    Source: analyses/mammal_scenic/code/05_analyze.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import os, json, pickle, hashlib

        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        import numpy as np, pandas as pd

        np.object = object
        from scipy.stats import hypergeom, pearsonr, t as tdist
        from ctxcore.rnkdb import FeatherRankingDatabase
        from pyscenic.aucell import aucell

        R = context.run_dir("scenic")
        P = R.parent
        O = R / "results/analysis"
        O.mkdir(exist_ok=True)
        mods = ["orange", "violet", "cyan", "skyblue3"]
        seeds = [20261023, 20261024, 20261025]

        def bh(x):
            x = np.asarray(x, float)
            assert np.isfinite(x).all()
            ix = np.argsort(x, kind="stable")
            q = np.minimum.accumulate((x[ix] * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
            out = np.empty_like(q)
            out[ix] = np.minimum(q, 1)
            return out

        cfg = json.loads((R / "config/existing_reference.json").read_text())
        db = FeatherRankingDatabase(cfg["database"], name="existing")
        E = pd.read_pickle(R / "inputs/expression.pkl")
        meta = pd.read_csv(R / "inputs/samples.tsv", sep="\t").set_index("sample")
        M = pd.read_csv(R / "inputs/frozen_WGCNA_members.tsv", sep="\t")
        all_selected = set(M.loc[M.module.isin(mods), "feature"])
        groups = {mod: set(M.loc[M.module.eq(mod), "feature"]) for mod in mods}
        U = set(M.feature) & set(E.columns) & set(db.genes)
        (O / "background_genes.txt").write_text("\n".join(sorted(U)) + "\n")
        traits = pd.read_csv(
            P / "growth_regulators_2026-10-21_v1/inputs/ma2016_full.samples.tsv", sep="\t"
        ).set_index("sample")
        species = traits.index.tolist()
        mass = np.log10(traits.loc[species, "mass_g"].to_numpy())
        assert set(species) == set(meta.species)
        eig = (
            pd.read_csv(
                P / "growth_regulators_2026-10-21_v1/results/wgcna/ma2016_full__raw/eigengenes.tsv",
                sep="\t",
            )
            .pivot(index="sample", columns="module", values="eigengene")
            .loc[species]
        )
        edges = {}
        sets = {}
        tests = {}
        for seed in seeds:
            a = pd.read_csv(R / "results" / str(seed) / "regulon_edges.tsv.gz", sep="\t")
            edges[seed] = a
            sets[seed] = {reg: set(t.target) for reg, t in a.groupby("regulon")}
            tfmap = a.groupby("regulon").TF.first().to_dict()
            rows = []
            for reg, gg in sets[seed].items():
                tf = tfmap[reg]
                u = U - {tf}
                tg = gg & u
                for mod in mods:
                    group = groups[mod] & u
                    hits = group & tg
                    N = len(u)
                    K = len(tg)
                    n = len(group)
                    k = len(hits)
                    p = float(hypergeom.sf(k - 1, N, K, n))
                    rows.append(
                        dict(
                            seed=seed,
                            module=mod,
                            regulon=reg,
                            TF=tf,
                            N_background=N,
                            K_targets_background=K,
                            n_module_covered=n,
                            k_shared=k,
                            n_module_total=len(groups[mod]),
                            p=p,
                            fold_enrichment=k / (n * K / N) if n * K else np.nan,
                            genes=";".join(sorted(hits)),
                        )
                    )
            t = pd.DataFrame(rows)
            t["q"] = bh(t.p)
            t.to_csv(O / f"overlap_{seed}.tsv", sep="\t", index=False)
            tests[seed] = t
        T = tests[seeds[0]].copy()
        eligible = T.K_targets_background.ge(5)
        T["q_ge5_diagnostic"] = np.nan
        T.loc[eligible, "q_ge5_diagnostic"] = bh(T.loc[eligible, "p"])
        st = []
        for row in T.itertuples():
            base = sets[seeds[0]][row.regulon]
            ss = [sets[s].get(row.regulon, set()) for s in seeds]
            core = set.intersection(*ss)
            target_j = [
                len(base & ss[i]) / len(base | ss[i]) if base | ss[i] else np.nan for i in [1, 2]
            ]
            st.append(
                dict(
                    regulon=row.regulon,
                    module=row.module,
                    TF=row.TF,
                    seeds_recovered=sum((bool(g) for g in ss)),
                    Jaccard_seed2=target_j[0],
                    Jaccard_seed3=target_j[1],
                    targets_all3=len(core),
                    module_genes_all3=len(core & groups[row.module]),
                    module_genes_all3_names=";".join(sorted(core & groups[row.module])),
                    overlap_q_seed2=float(
                        tests[seeds[1]]
                        .set_index(["regulon", "module"])
                        .q.get((row.regulon, row.module), 1.0)
                    ),
                    overlap_q_seed3=float(
                        tests[seeds[2]]
                        .set_index(["regulon", "module"])
                        .q.get((row.regulon, row.module), 1.0)
                    ),
                )
            )
        S = pd.DataFrame(st)
        S.to_csv(O / "seed_stability.tsv", sep="\t", index=False)
        T = T.merge(S, on=["regulon", "module", "TF"])
        T = T.sort_values(["q", "p", "module", "TF"])
        T.to_csv(O / "module_regulon_results.tsv", sep="\t", index=False)
        a = edges[seeds[0]].copy()
        a = a[a.target.isin(all_selected)].merge(
            M[["feature", "module", "kME", "weight"]].rename(columns={"feature": "target"}),
            on="target",
        )
        a = a[a.module.isin(mods)]
        a["recovered_seed2"] = [
            g in sets[seeds[1]].get(r, set()) for r, g in zip(a.regulon, a.target)
        ]
        a["recovered_seed3"] = [
            g in sets[seeds[2]].get(r, set()) for r, g in zip(a.regulon, a.target)
        ]
        a["binding_measured"] = False
        a.to_csv(O / "module_regulon_genes.tsv", sep="\t", index=False)
        auc = pd.read_csv(R / "results" / str(seeds[0]) / "AUCell.tsv", sep="\t", index_col=0)
        auc_species = auc.groupby(meta.loc[auc.index, "species"]).mean().loc[species]
        auc_species.to_csv(O / "AUCell_species.tsv", sep="\t")
        exprsp = E.groupby(meta.loc[E.index, "species"]).mean().loc[species]
        activity = []
        for reg in auc_species.columns:
            y = auc_species[reg].to_numpy()
            r, p = pearsonr(mass, y) if np.ptp(y) > 1e-12 else (np.nan, 1.0)
            tf = edges[seeds[0]].loc[edges[seeds[0]].regulon.eq(reg), "TF"].iloc[0]
            activity.append(
                dict(
                    regulon=reg,
                    TF=tf,
                    n_species=len(species),
                    mass_r=r,
                    mass_p=p,
                    own_RNA_mass_r=(
                        pearsonr(mass, exprsp[tf])[0]
                        if tf in exprsp.columns and np.ptp(exprsp[tf]) > 1e-12
                        else np.nan
                    ),
                )
            )
        A = pd.DataFrame(activity)
        A["mass_q"] = bh(A.mass_p)
        A.to_csv(O / "regulon_mass_associations.tsv", sep="\t", index=False)
        with (R / "results" / str(seeds[0]) / "regulons.pkl").open("rb") as f:
            regs = pickle.load(f)
        proxy = []
        proxyrows = []
        for reg in regs:
            keep = set(reg.gene2weight) - all_selected - {reg.transcription_factor}
            if len(keep) >= 5:
                proxy.append(reg.copy(gene2weight={g: reg.gene2weight[g] for g in keep}))
        if proxy:
            pa = aucell(E, proxy, auc_threshold=0.05, num_workers=8, seed=20261023)
            ps = pa.groupby(meta.loc[pa.index, "species"]).mean().loc[species]
            ps.to_csv(O / "disjoint_AUCell_species.tsv", sep="\t")
            for row in T.itertuples():
                if row.regulon not in ps.columns:
                    continue
                y = ps[row.regulon].to_numpy()
                z = eig[row.module].to_numpy()
                r, p = pearsonr(y, z) if np.ptp(y) > 1e-12 else (np.nan, 1.0)
                design = np.column_stack([np.ones(len(mass)), mass])
                yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
                zr = z - design @ np.linalg.lstsq(design, z, rcond=None)[0]
                pr = (
                    float(pearsonr(yr, zr)[0])
                    if np.std(yr) > 1e-12 and np.std(zr) > 1e-12
                    else np.nan
                )
                pp = (
                    float(
                        2
                        * tdist.sf(
                            abs(pr) * np.sqrt((len(mass) - 3) / max(1e-16, 1 - pr * pr)),
                            len(mass) - 3,
                        )
                    )
                    if np.isfinite(pr)
                    else 1.0
                )
                proxyrows.append(
                    dict(
                        regulon=row.regulon,
                        TF=row.TF,
                        module=row.module,
                        n_nonmodule_targets=len(
                            sets[seeds[0]][row.regulon] - all_selected - {row.TF}
                        ),
                        module_r=r,
                        module_p=p,
                        partial_mass_r=pr,
                        partial_mass_p=pp,
                    )
                )
            Q = pd.DataFrame(proxyrows)
            Q["module_q"] = bh(Q.module_p)
            Q["partial_mass_q"] = bh(Q.partial_mass_p)
            Q.to_csv(O / "disjoint_activity_support.tsv", sep="\t", index=False)
        project = []
        for mod in mods:
            members = M[M.module.eq(mod)].set_index("feature")
            gg = sorted(set(members.index) & set(exprsp.columns))
            x = exprsp[gg]
            Z = (x - x.mean()) / x.std(ddof=1)
            score = Z.to_numpy() @ members.loc[gg, "weight"].to_numpy()
            project.append(
                dict(
                    module=mod,
                    n_genes=len(gg),
                    r_with_frozen_eigengene=float(pearsonr(score, eig[mod])[0]),
                    r_mass=float(pearsonr(score, mass)[0]),
                    scope="Diagnostic projection on new normalization, frozen genes/weights; WGCNA not refit",
                )
            )
        pd.DataFrame(project).to_csv(
            O / "normalization_projection_diagnostic.tsv", sep="\t", index=False
        )
        summary = dict(
            background=len(U),
            module_coverage={m: len(groups[m] & U) for m in mods},
            seeds={
                s: json.loads((R / "results" / str(s) / "SCENIC_summary.json").read_text())
                for s in seeds
            },
            primary_overlap_tests=len(T),
            primary_overlap_hits=int(T.q.lt(0.05).sum()),
            TFs_with_primary_module_hits=sorted(T.loc[T.q.lt(0.05), "TF"].unique()),
            mass_association_tests=len(A),
            mass_association_hits=int(A.mass_q.lt(0.05).sum()),
            top_module_results=T.head(15).replace({np.nan: None}).to_dict("records"),
            normalization_diagnostic=project,
            E2F_results=T[T.TF.str.startswith("E2F")].replace({np.nan: None}).to_dict("records"),
            causal_regulator_established=False,
            method="Full pySCENIC on26 Ma2016 bulk donor-labelled groups; preexisting human v10 short-promoter ranking; exploratory within-study results",
            postprocessing_revision="v2: evaluate all primary regulons x4 as declared; >=5-background-target family retained only as sensitivity; constant AUCell profiles assigned p=1 before BH after independent R check exposed NaNs in v1",
            code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        )
        (O / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False))
        print(
            json.dumps(
                {
                    k: v
                    for k, v in summary.items()
                    if k not in ["top_module_results", "E2F_results"]
                },
                indent=2,
            ),
            flush=True,
        )

    return SimpleNamespace(**locals())


def _stability(context):
    """Existing analysis stage; calculation statements retained from the source.

    Source: analyses/mammal_scenic/code/08_stability.py at d230e56ab942."""

    def main():
        from pathlib import Path
        import os, json, hashlib

        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        import pandas as pd, numpy as np
        from scipy.stats import pearsonr

        R = context.run_dir("scenic")
        O = R / "results/analysis"
        seeds = [20261023, 20261024, 20261025]
        meta = pd.read_csv(R / "inputs/samples.tsv", sep="\t").set_index("sample")
        traits = pd.read_csv(R / "inputs/species_traits.tsv", sep="\t").set_index("sample")
        sp = traits.index.tolist()
        x = np.log10(traits.mass_g.to_numpy())

        def bh(a):
            a = np.asarray(a)
            ix = np.argsort(a)
            q = np.minimum.accumulate((a[ix] * len(a) / np.arange(1, len(a) + 1))[::-1])[::-1]
            r = np.empty_like(q)
            r[ix] = np.minimum(q, 1)
            return r

        sets = {}
        rows = []
        for seed in seeds:
            e = pd.read_csv(R / "results" / str(seed) / "regulon_edges.tsv.gz", sep="\t")
            sets[seed] = {r: set(t.target) for r, t in e.groupby("regulon")}
            auc = pd.read_csv(R / "results" / str(seed) / "AUCell.tsv", sep="\t", index_col=0)
            a = auc.groupby(meta.loc[auc.index, "species"]).mean().loc[sp]
            temp = []
            for reg in a:
                y = a[reg].to_numpy()
                r, p = pearsonr(x, y) if np.ptp(y) > 1e-12 else (np.nan, 1.0)
                temp.append(dict(seed=seed, regulon=reg, mass_r=r, mass_p=p))
            t = pd.DataFrame(temp)
            t["mass_q"] = bh(t.mass_p)
            rows.append(t)
        A = pd.concat(rows)
        A.to_csv(O / "mass_association_all_seeds.tsv", sep="\t", index=False)
        base = sets[seeds[0]]
        j = []
        for reg, ss in base.items():
            gg = [sets[s].get(reg, set()) for s in seeds]
            r = dict(
                regulon=reg,
                TF=reg.split("(")[0],
                recovered=sum((bool(z) for z in gg)),
                primary_targets=len(ss),
                targets_all3=len(set.intersection(*gg)),
            )
            for i, s in enumerate(seeds[1:], 2):
                r[f"Jaccard_seed{i}"] = (
                    len(ss & gg[i - 1]) / len(ss | gg[i - 1]) if ss | gg[i - 1] else 0.0
                )
            j.append(r)
        J = pd.DataFrame(j)
        J.to_csv(O / "regulon_stability.tsv", sep="\t", index=False)
        T = pd.read_csv(O / "module_regulon_results.tsv", sep="\t")
        T["q_positive_runs"] = (
            T.q.lt(0.05).astype(int)
            + T.overlap_q_seed2.lt(0.05).astype(int)
            + T.overlap_q_seed3.lt(0.05).astype(int)
        )
        T.to_csv(O / "module_regulon_results.tsv", sep="\t", index=False)
        s = {
            "primary_regulons": len(base),
            "recovered_in_all3": int(J.recovered.eq(3).sum()),
            "median_target_Jaccard_seed2": float(J.Jaccard_seed2.median()),
            "median_target_Jaccard_seed3": float(J.Jaccard_seed3.median()),
            "primary_significant_module_pairs": T[T.q.lt(0.05)]
            .replace({np.nan: None})
            .to_dict("records"),
            "module_pairs_q_positive_all3": int(T.q_positive_runs.eq(3).sum()),
            "mass_ERG_FOXP4_all_seeds": A[A.regulon.isin(["ERG(+)", "FOXP4(+)"])]
            .replace({np.nan: None})
            .to_dict("records"),
            "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
        (O / "stability_summary.json").write_text(json.dumps(s, indent=2, allow_nan=False))
        print(json.dumps(s, indent=2))

    return SimpleNamespace(**locals())
