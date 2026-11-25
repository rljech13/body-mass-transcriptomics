#!/usr/bin/env python3
"""Distance of proximal ELA to the promoter PLA of the same gene.

The ELA in question are the near-TSS fragments already assigned to the
rising set (green 51) and the falling set (OXPHOS 16). Distance is the
gap between CREate cores on the same chromosome. A gap of 0 means the
cores overlap or touch. The full CREate region is scored the same way.

Missing loci are omitted. They are not filled with zero.
"""

from __future__ import annotations

from pathlib import Path

BASE = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
PEAKS = (
    BASE
    / "results/oxphos_species_thread/create_to_genes/peaks_nearest_tss.tsv"
)
OUT = BASE / "results/oxphos_species_thread/create_to_genes/ela_to_pla"
ORDER = ["Zebra", "Quail", "Chicken", "Duck", "Peacock", "Turkey", "Emu"]
BINS = [
    ("overlap", 0, 0),
    ("1-50", 1, 50),
    ("51-200", 51, 200),
    ("201-500", 201, 500),
    ("501-1000", 501, 1000),
    ("1001-2000", 1001, 2000),
    (">2000", 2001, 10**12),
]


def gap(a0: int, a1: int, b0: int, b1: int) -> int:
    if a1 < a0:
        a0, a1 = a1, a0
    if b1 < b0:
        b0, b1 = b1, b0
    if a0 < b1 and b0 < a1:
        return 0
    if a1 <= b0:
        return b0 - a1
    return a0 - b1


def load_sets() -> dict[str, set[str]]:
    out = {}
    for key, name in (
        ("green_51", "story_green_sig_51.txt"),
        ("oxphos_16", "story_oxphos_16.txt"),
    ):
        text = (BASE / "data/genesets" / name).read_text().splitlines()
        out[key] = {line.strip().upper() for line in text if line.strip()}
    return out


def main() -> None:
    sets = load_sets()
    pla: dict[tuple[str, str], list[dict]] = {}
    ela_rows: list[dict] = []
    header = PEAKS.read_text().splitlines()[0].split("\t")
    for line in PEAKS.read_text().splitlines()[1:]:
        rec = dict(zip(header, line.split("\t")))
        gene = rec["gene"].upper()
        rec["gene"] = gene
        rec["start"] = int(rec["start"])
        rec["end"] = int(rec["end"])
        rec["core_start"] = int(rec["core_start"])
        rec["core_end"] = int(rec["core_end"])
        rec["dist_bp"] = int(rec["dist_bp"])
        rec["tpm"] = float(rec["tpm"])
        membership = []
        for key, genes in sets.items():
            if gene in genes:
                membership.append(key)
        if not membership:
            continue
        rec["membership"] = membership
        if rec["cre_class"] == "PLA":
            pla.setdefault((rec["species"], gene), []).append(rec)
        elif rec["cre_class"] == "ELA" and rec["shelf"] == "proximal":
            ela_rows.append(rec)

    OUT.mkdir(parents=True, exist_ok=True)
    detail = OUT / "proximal_ela_to_pla.tsv"
    with detail.open("w") as handle:
        handle.write(
            "set\tspecies\tgene\ttpm\tdist_tss\t"
            "core_gap\tregion_gap\tnearest_pla_shelf\tnearest_pla_tpm\t"
            "n_pla_same_gene\tchrom\tela_start\tela_end\t"
            "ela_core_start\tela_core_end\n"
        )
        scored: dict[str, list[dict]] = {key: [] for key in sets}
        for rec in ela_rows:
            mates = [
                p
                for p in pla.get((rec["species"], rec["gene"]), [])
                if p["chrom"] == rec["chrom"]
            ]
            if mates:
                best = min(
                    mates,
                    key=lambda p: (
                        gap(
                            rec["core_start"],
                            rec["core_end"],
                            p["core_start"],
                            p["core_end"],
                        ),
                        gap(rec["start"], rec["end"], p["start"], p["end"]),
                    ),
                )
                core_gap = gap(
                    rec["core_start"],
                    rec["core_end"],
                    best["core_start"],
                    best["core_end"],
                )
                region_gap = gap(
                    rec["start"], rec["end"], best["start"], best["end"]
                )
                pla_shelf = best["shelf"]
                pla_tpm = best["tpm"]
                n_pla = len(mates)
            else:
                core_gap = region_gap = -1
                pla_shelf = ""
                pla_tpm = 0.0
                n_pla = 0
            row = {
                **rec,
                "core_gap": core_gap,
                "region_gap": region_gap,
                "pla_shelf": pla_shelf,
                "pla_tpm": pla_tpm,
                "n_pla": n_pla,
            }
            for key in rec["membership"]:
                scored[key].append(row)
                handle.write(
                    f"{key}\t{rec['species']}\t{rec['gene']}\t{rec['tpm']:.6g}\t"
                    f"{rec['dist_bp']}\t{core_gap}\t{region_gap}\t{pla_shelf}\t"
                    f"{pla_tpm:.6g}\t{n_pla}\t{rec['chrom']}\t{rec['start']}\t"
                    f"{rec['end']}\t{rec['core_start']}\t{rec['core_end']}\n"
                )

    lines = [
        "Proximal ELA (≤2 kb from the assigned TSS) versus the nearest",
        "PLA of the same gene on the same chromosome.",
        "core_gap 0 = cores overlap or the intervals interleave.",
        "A gene with no PLA is counted as no_pla, not as a large distance.",
        "",
    ]
    for key in ("green_51", "oxphos_16"):
        rows = scored[key]
        lines.append(f"## {key}  proximal ELA n={len(rows)}")
        lines.extend(block(rows))
        lines.append("")
        lines.append(f"### {key} by species")
        for sp in ORDER:
            sub = [r for r in rows if r["species"] == sp]
            lines.append(f"{sp} n={len(sub)}")
            lines.extend(block(sub, indent="  "))
        lines.append("")
    text = "\n".join(lines) + "\n"
    (OUT / "summary.txt").write_text(text)
    print(text, end="")


def block(rows: list[dict], indent: str = "") -> list[str]:
    if not rows:
        return [indent + "(none)"]
    have = [r for r in rows if r["core_gap"] >= 0]
    none = [r for r in rows if r["core_gap"] < 0]
    tpm_all = sum(r["tpm"] for r in rows)
    tpm_have = sum(r["tpm"] for r in have)
    out = [
        f"{indent}with a PLA on the same gene: {len(have)}/{len(rows)} peaks, "
        f"TPM {tpm_have:.1f}/{tpm_all:.1f}",
        f"{indent}no PLA on that gene: {len(none)} peaks, "
        f"TPM {sum(r['tpm'] for r in none):.1f}",
    ]
    if not have:
        return out
    gaps = sorted(r["core_gap"] for r in have)
    rgaps = sorted(r["region_gap"] for r in have)
    tss = sorted(r["dist_bp"] for r in have)

    def pct(xs: list[int], q: float) -> int:
        return xs[min(len(xs) - 1, int(q * (len(xs) - 1)))]

    out.append(
        f"{indent}core gap bp  p10={pct(gaps, 0.1)}  p50={pct(gaps, 0.5)}  "
        f"p90={pct(gaps, 0.9)}  max={gaps[-1]}"
    )
    out.append(
        f"{indent}region gap bp  p50={pct(rgaps, 0.5)}  p90={pct(rgaps, 0.9)}"
    )
    out.append(
        f"{indent}ELA core to TSS bp  p50={pct(tss, 0.5)}  p90={pct(tss, 0.9)}"
    )
    out.append(f"{indent}core gap, peaks and TPM share among those with a PLA:")
    for name, lo, hi in BINS:
        sub = [r for r in have if lo <= r["core_gap"] <= hi]
        share = 100.0 * len(sub) / len(have)
        tshare = 100.0 * sum(r["tpm"] for r in sub) / tpm_have
        out.append(
            f"{indent}  {name:10} peaks {len(sub):4} ({share:5.1f}%)  "
            f"TPM {sum(r['tpm'] for r in sub):7.1f} ({tshare:5.1f}%)"
        )
    inside = [r for r in have if r["region_gap"] == 0]
    out.append(
        f"{indent}ELA region touches the PLA region: {len(inside)}/{len(have)} "
        f"({100.0 * len(inside) / len(have):.1f}% of peaks, "
        f"{100.0 * sum(r['tpm'] for r in inside) / tpm_have:.1f}% of TPM)"
    )
    prox_pla = [r for r in have if r["pla_shelf"] == "proximal"]
    out.append(
        f"{indent}nearest PLA is itself proximal: {len(prox_pla)}/{len(have)}"
    )
    return out


if __name__ == "__main__":
    main()
