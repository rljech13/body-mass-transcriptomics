#!/usr/bin/env python3
"""Report for the mammal size scan. Does not edit the frozen bird nine."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES, normalize_symbol  # noqa: E402

OUT = RES / "mammal_size_sets"
CORPUS = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/corpus_kleiber_2026-08-06/02_база_знаний")
NINE = (
    ("oxphos_16", "OXPHOS 16"),
    ("folate_12", "фолат 12"),
    ("ribo_5", "ribo 5"),
    ("green_hubs_7", "green hubs 7"),
    ("green_51_rest", "green, остальные 41"),
    ("carbon_7", "углеродный метаболизм 7"),
    ("amino_5", "биосинтез аминокислот 5"),
    ("nafld_4", "NAFLD 4"),
    ("epi_lipid_3", "липиды и атеросклероз 3"),
)


def load_nine() -> dict[str, set[str]]:
    frames = [
        pd.read_csv(CORPUS / "working_sets_disjoint_2026-09-22.tsv", sep="\t"),
        pd.read_csv(CORPUS / "working_sets_overlap_ok_2026-09-22.tsv", sep="\t"),
    ]
    raw = pd.concat(frames, ignore_index=True)
    out = {}
    for sid, _label in NINE:
        genes = {normalize_symbol(g) for g in raw.loc[raw.set_id == sid, "gene"]}
        genes.discard("")
        out[sid] = genes
    return out


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    u = len(a | b)
    return len(a & b) / u if u else float("nan")


def fmt_r(r: float) -> str:
    return f"{r:+.3f}"


def main() -> None:
    sets = pd.read_csv(OUT / "sets_vs_mass.tsv", sep="\t")
    anchors = pd.read_csv(OUT / "anchors_frozen_vs_honest.tsv", sep="\t")
    cam = pd.read_csv(OUT / "camera.tsv", sep="\t")
    gsea = pd.read_csv(OUT / "gsea.tsv", sep="\t")
    phy = pd.read_csv(OUT / "pgls.tsv", sep="\t")
    nine = load_nine()
    union = set().union(*nine.values())

    passed = sets[sets.bh < 0.05].copy()
    key = ["panel", "catalog", "id"]
    wide = passed.merge(cam[key + ["rho", "p_camera", "bh_camera"]], on=key, how="left")
    wide = wide.merge(gsea[key + ["es", "p_gsea", "bh_gsea"]], on=key, how="left")
    if len(phy):
        wide = wide.merge(
            phy[key + ["pgls_tt_p", "pgls_gr_p", "lambda_pegel", "p_brownian", "pic_p", "note"]],
            on=key, how="left",
        )
    else:
        wide["note"] = ""

    # cross-panel: same catalog and id, BH in at least two panels
    if len(passed):
        cnt = passed.groupby(["catalog", "id", "name"]).panel.nunique().reset_index(name="n_panels")
        shared = cnt[cnt.n_panels >= 2]
    else:
        shared = pd.DataFrame(columns=["catalog", "id", "name", "n_panels"])
    shared.to_csv(OUT / "shared_names.tsv", sep="\t", index=False)

    jac_rows = []
    for rec in passed.itertuples(index=False):
        genes = {normalize_symbol(g) for g in str(rec.genes).split(";") if g}
        genes.discard("")
        cover = len(genes & union) / len(genes) if genes else float("nan")
        for sid, label in NINE:
            jac_rows.append(dict(
                panel=rec.panel, catalog=rec.catalog, id=rec.id, name=rec.name,
                nine=sid, nine_label=label, jaccard=jaccard(genes, nine[sid]),
                fraction_in_nine=cover,
                same_observation=bool(cover >= 0.8),
            ))
    jac = pd.DataFrame(jac_rows)
    jac.to_csv(OUT / "jaccard_nine.tsv", sep="\t", index=False)

    lines = []
    lines.append("# Группы генов и масса у млекопитающих")
    lines.append("")
    lines.append("22 сентября 2026. Три панели посчитаны отдельно. В одну регрессию они не входили.")
    lines.append("Девять птичьих наборов не пополнялись.")
    lines.append("")
    lines.append("Счёт набора — log10 суммы видовых значений против log10 массы.")
    lines.append("Ген остаётся, только если он больше нуля у всех видов панели. Пропавший ген вынут, не обнулён.")
    lines.append("Видовой счёт — среднее по библиотекам. Библиотека не единица регрессии.")
    lines.append("Нуль суммы — 2000 наборов той же величины и той же децили медианной обильности.")
    lines.append("Поправка BH — внутри каталога и внутри панели.")
    lines.append("")
    lines.append("## Вселенные")
    lines.append("")
    for panel, label in (
        ("ma2016", "Ma 2016, фибробласты брюшной кожи"),
        ("wagner", "Wagner, кожные фибробласты"),
        ("gse105160", "GSE105160, mock дермальных фибробластов"),
    ):
        m = pd.read_csv(RES / "mammal_honest" / f"{panel}_M.tsv", sep="\t", index_col=0)
        mass = pd.read_csv(RES / "mammal_honest" / f"{panel}_mass.tsv", sep="\t")
        lines.append(f"- {label}: {m.shape[0]} генов, {m.shape[1]} видов, масса от {mass.mass_g.min():.0f} до {mass.mass_g.max():.0f} г.")
    lines.append("")
    lines.append("У Wagner нет третьей библиотеки лошади (`Ec_SF_3`), в матрице остаются две.")
    lines.append("Девять из десяти видов GSE105160 — приматы.")
    lines.append("")
    lines.append("## Якоря: старый счёт и честная вселенная")
    lines.append("")
    lines.append("Старые r скриптов 42 и 43 не перезаписаны. Рядом — те же две карты на честной вселенной.")
    lines.append("")
    lines.append("| панель | карта | r старый | p старый | n честных | r честный | p OLS | p нуля | BH |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for rec in anchors.sort_values(["panel", "id"]).itertuples(index=False):
        lines.append(
            f"| {rec.panel} | {rec.id} | {rec.r_frozen:+.3f} | {rec.p_frozen:.3f} | {int(rec.n_honest)} | "
            f"{rec.r_honest:+.3f} | {rec.p_honest:.3f} | {rec.p_null_honest:.3f} | {rec.bh_honest:.3f} |"
        )
    lines.append("")
    lines.append("## Сколько наборов испытано")
    lines.append("")
    lines.append("| панель | каталог | испытано | BH < 0.05 |")
    lines.append("|---|---|---:|---:|")
    for (panel, cat), g in sets.groupby(["panel", "catalog"]):
        lines.append(f"| {panel} | {cat} | {len(g)} | {int((g.bh < 0.05).sum())} |")
    lines.append("")
    lines.append("Каталоги: KEGG `hsa` из `kegg_hsa_all_pathways.tsv`; Hallmark, Reactome (`CP:REACTOME`) и GO Biological Process — MSigDB 2026.1.Hs (Zenodo 18968178).")
    lines.append("GO взят в окне 8–50 генов каталога, затем оставлены наборы с не меньше чем 8 генами внутри вселенной панели.")
    lines.append("CORUM — человеческие комплексы, снимок Zenodo 17419058 (21 января 2025): живой сервер Helmholtz на этом ходе не ответил.")
    lines.append("Иммунные, позиционные и микроРНК-коллекции MSigDB не брались.")
    lines.append("p нуля 0 значит, что ни один из 2000 подобранных наборов не достиг такого модуля корреляции. Разрешение этого нуля — 0.0005.")
    lines.append("")
    if passed.empty:
        lines.append("## Что прошло поправку суммы")
        lines.append("")
        lines.append("Ни один набор не прошёл BH < 0.05 по нулю обильности. Пустой FDR — ответ этого прогона, порог не ослаблялся.")
        lines.append("")
    else:
        lines.append("## Что прошло поправку суммы")
        lines.append("")
        lines.append("Прохождение — BH нуля обильности. CAMERA, GSEA и PGLS эту строку не стирают.")
        lines.append("")
        lines.append("| панель | каталог | id | n | r | p нуля | BH | BH CAMERA | BH GSEA | PGLS | Грейфен | пометка |")
        lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        show = wide.sort_values(["panel", "bh", "id"])
        for rec in show.itertuples(index=False):
            pgls = "" if pd.isna(rec.pgls_tt_p) else f"{rec.pgls_tt_p:.4f}"
            graf = "" if pd.isna(rec.pgls_gr_p) else f"{rec.pgls_gr_p:.4f}"
            note = "" if not isinstance(getattr(rec, "note", ""), str) else rec.note
            lines.append(
                f"| {rec.panel} | {rec.catalog} | {rec.id} | {int(rec.n)} | {rec.r:+.3f} | {rec.p_null:.4f} | {rec.bh:.3f} | "
                f"{rec.bh_camera:.3f} | {rec.bh_gsea:.3f} | {pgls} | {graf} | {note} |"
            )
        lines.append("")
        if (show.bh_camera >= 0.05).all() and (show.bh_gsea >= 0.05).all():
            lines.append("Ни один из этих наборов не проходит BH в CAMERA и не проходит BH в GSEA.")
        if len(jac) and float(jac.jaccard.max()) == 0.0:
            lines.append("Пересечение каждого из них с девяткой пустое: Жаккар по всем девяти равен нулю (`jaccard_nine.tsv`).")
        lines.append("")

    lines.append("## Имена в более чем одной панели")
    lines.append("")
    if shared.empty:
        lines.append("Нет id, который прошёл BH суммы в двух панелях. Совпадение знака без поправки сюда не входит.")
    else:
        lines.append("Общее имя — один и тот же каталог и id, BH < 0.05 в каждой из панелей, где оно перечислено.")
        lines.append("")
        for rec in shared.itertuples(index=False):
            lines.append(f"- {rec.catalog} `{rec.id}` ({rec.name}): {int(rec.n_panels)} панели.")
    lines.append("")
    lines.append("## Вторые чтения")
    lines.append("")
    lines.append("CAMERA: видовые средние, log2, ковариата log10 массы, корреляция между генами оценена внутри набора и обрезана снизу нулём.")
    lines.append("Пакет limma на этом R не установлен, статистика посчитана по Wu и Smyth (2012) и сжатию дисперсий Smyth (2004).")
    lines.append(f"CAMERA, BH < 0.05: {int((cam.bh_camera < 0.05).sum())} из {len(cam)}.")
    lines.append("")
    lines.append("Ранговый GSEA: вес гена — корреляция log-экспрессии с log10 массы, 10 000 перестановок вектора масс.")
    lines.append("Перестановка состава набора основным нулём не была.")
    lines.append(f"GSEA, BH < 0.05: {int((gsea.bh_gsea < 0.05).sum())} из {len(gsea)}.")
    lines.append("")
    if phy.empty:
        lines.append("PGLS считался только для наборов с BH суммы < 0.05. Таких нет, таблица `pgls.tsv` пустая. OLS других наборов из-за этого не стирался: его просто не водили на дерево.")
    else:
        n_clade = int((phy.note.fillna("") == "клад, не размер").sum())
        lines.append(f"PGLS на приблизительном ультраметрическом дереве панели: {len(phy)} наборов. Пометка «клад, не размер» у {n_clade}: датированный PGLS даёт p ≥ 0.05, OLS-строка при этом остаётся.")
        disagree = phy[(phy.pgls_tt_p >= 0.05) != (phy.pgls_gr_p >= 0.05)]
        if len(disagree):
            bits = [
                f"`{r.id}` датированный {r.pgls_tt_p:.4f}, Грейфен {r.pgls_gr_p:.4f}"
                for r in disagree.itertuples(index=False)
            ]
            lines.append("Датированное дерево и длины Грейфена расходятся на пороге 0.05 у: " + "; ".join(bits) + ".")
        lines.append("Броуновские симуляции лежат в `pgls.tsv`. Дерево помечено как приблизительное (`trees_note.tsv`).")
    lines.append("")
    lines.append("## Чего этот прогон не делал")
    lines.append("")
    lines.append("Сеть внутри вида и спуск Meinshausen не строились: у Ma и Wagner по две–три библиотеки на вид.")
    lines.append("Дерево по видовым средним не строилось.")
    lines.append("Список генов, значимых по одному, не собирался.")
    lines.append("Панели, птицы, FANTOM MSC и одиночные библиотеки не склеивались.")
    lines.append("Официальные `hsa00190` и `hsa04714` не заменялись урезанными ядрами.")
    lines.append("Членство девятки не заморожено заново и не расширено. Приёмка корпуса человеком не закрыта.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    dest = OUT / "ОТЧЁТ_группы_и_масса_маммалианы_2026-09-22.md"
    dest.write_text(text)
    print(dest, flush=True)


if __name__ == "__main__":
    main()
