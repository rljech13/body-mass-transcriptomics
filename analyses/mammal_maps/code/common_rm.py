"""Panel loader shared with 01_sets_mass.py (copied verbatim) plus the fixed per-panel background used in step 1."""
from pathlib import Path
import numpy as np
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
RES = RUN.parent
GI = RES / "growth_regulators_2026-10-21_v1/inputs"
adm = pd.read_csv(RES / "nonbird_common_state_2026-10-21_v1/config/species_trait_admission.corrected.tsv", sep="\t")


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
    lin_s = by_species(lin).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
    log_s = by_species(lg).loc[lambda d: d.index.isin(fs.index)].rename(index=fs.symbols.to_dict())
    # identical composition: measured with positive linear expression in every admitted species
    ok = (lin_s > 0).all(axis=1) & np.isfinite(log_s).all(axis=1)
    lin_s, log_s = lin_s[ok], log_s.loc[lin_s.index[ok]]
    mass = np.log10(a.loc[keep, "accepted_mass_g"].to_numpy(float))
    return dict(species=keep, C=C, mass=mass, lin=lin_s, log=log_s)


def background(d):
    Z = d["log"].sub(d["log"].mean(1), axis=0).div(d["log"].std(1), axis=0)
    return Z[np.isfinite(Z).all(1)]
