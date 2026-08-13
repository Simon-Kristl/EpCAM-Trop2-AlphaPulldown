#!/usr/bin/env python3
"""
Analiza pogostosti kontaktnih ostankov glede na iLIS razrede.

KAJ SKRIPTA NAREDI
------------------
Za vsak ostanek na strani vabe prešteje:
1) koliko vseh modelov ima ta ostanek v stolpcu cLIR_indices_i;
2) koliko modelov z iLIS >= 0.223 ima ta ostanek v kontaktu;
3) v katerih iLIS razredih se ostanek pojavlja:
   - iLIS < 0.223
   - 0.223 <= iLIS < 0.4
   - 0.4 <= iLIS < 0.6
   - iLIS >= 0.6
4) vsoto iLIS vrednosti za modele, kjer je ostanek v kontaktu;
5) povprečni iLIS za modele, kjer je ostanek v kontaktu.

POMEMBNO
--------
Analiza šteje MODELE, ne unikatnih tarčnih proteinov.
Če ima isti protein več rangiranih modelov, lahko prispeva večkrat.

KAKO POGNATI ZA EpCAM
---------------------
EpCAM AFM-LIS cLIR_indices_i uporablja modelsko oštevilčenje 1-242.
Ker model 1-242 ustreza UniProt 24-265, uporabimo offset 23.

Primer:
python residue_contacts_by_iLIS_analysis_with_usage.py \
  --csv EpCAM_iLIS.csv \
  --out_prefix EpCAM_contacts_by_iLIS \
  --residue_min 1 \
  --residue_max 242 \
  --offset 23 \
  --threshold 0.223

KAKO POGNATI ZA Trop2
---------------------
Trop2 AFM-LIS cLIR_indices_i uporablja modelsko oštevilčenje 1-248.
Ker model 1-248 ustreza UniProt 27-274, uporabimo offset 26.

Primer:
python residue_contacts_by_iLIS_analysis_with_usage.py \
  --csv "Trop2_iLIS(2).csv" \
  --out_prefix Trop2_contacts_by_iLIS \
  --residue_min 1 \
  --residue_max 248 \
  --offset 26 \
  --threshold 0.223

IZHODNE DATOTEKE
----------------
Skripta ustvari:
1) *_residue_iLIS_summary.csv
   - glavna tabela za vse ostanke

2) *_top_residues_by_iLIS_metrics.csv
   - top ostanki po različnih kriterijih

3) *_counts_by_iLIS_bin.csv
   - število kontaktov po iLIS razredih

4) *_heatmap_counts_by_iLIS_bin.png
   - heatmap kontaktov po iLIS razredih

5) *_heatmap_sum_iLIS_per_residue.png
   - heatmap, kjer so kontakti uteženi z vsoto iLIS

6) *_barplot_contacts_iLIS_ge_threshold.png
   - barplot kontaktov samo za modele z iLIS >= threshold
"""

import argparse
from pathlib import Path
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


TOKEN_RE = re.compile(r"\d+\s*(?:-\s*\d+)?")


def parse_indices(cell):
    """
    Pretvori celico, kot je:
        [120,122,124,162,164,167,169,225,227-230]
    v množico ostankov:
        {120, 122, 124, 162, 164, 167, 169, 225, 227, 228, 229, 230}
    """
    if pd.isna(cell):
        return set()

    text = str(cell).strip()

    if not text or text == "[]":
        return set()

    residues = set()

    for match in TOKEN_RE.finditer(text):
        token = match.group(0).replace(" ", "")

        if "-" in token:
            start, end = map(int, token.split("-", 1))
            if start > end:
                start, end = end, start
            residues.update(range(start, end + 1))
        else:
            residues.add(int(token))

    return residues


def partner_from_name(name):
    """
    Iz imena modela izlušči UniProt kodo partnerja.

    Primer:
        P16422_and_Q9Y5Y6 -> Q9Y5Y6
        P09758_and_Q9Y5Y6 -> Q9Y5Y6
    """
    text = str(name)

    if "_and_" in text:
        return text.split("_and_", 1)[1]

    return text


def analyze(csv_path, out_prefix, residue_min, residue_max, offset=0, threshold=0.223):
    """
    Glavna funkcija za analizo.

    Parametri:
    ----------
    csv_path : str
        Pot do AFM-LIS CSV datoteke.

    out_prefix : str
        Predpona za izhodne datoteke.

    residue_min, residue_max : int
        Razpon ostankov v MODELSKEM oštevilčenju, kot se uporablja v cLIR_indices_i.

    offset : int
        Korekcija za pretvorbo modelskega oštevilčenja v UniProt oštevilčenje.
        EpCAM: offset = 23
        Trop2: offset = 26

    threshold : float
        Prag iLIS za pomembne napovedane interakcije.
    """
    csv_path = Path(csv_path)
    out_prefix = Path(out_prefix)

    df = pd.read_csv(csv_path, low_memory=False)

    if "iLIS" not in df.columns or "cLIR_indices_i" not in df.columns:
        raise ValueError("CSV mora vsebovati stolpca 'iLIS' in 'cLIR_indices_i'.")

    if "name" not in df.columns:
        raise ValueError("CSV mora vsebovati stolpec 'name'.")

    df["iLIS"] = pd.to_numeric(df["iLIS"], errors="coerce").fillna(0)
    df["partner"] = df["name"].apply(partner_from_name)
    df["contact_residues_model"] = df["cLIR_indices_i"].apply(parse_indices)
    df["has_contacts"] = df["contact_residues_model"].apply(len) > 0

    # iLIS razredi.
    # right=False pomeni, da je spodnja meja vključena:
    # 0.223 spada v razred 0.223-0.4.
    bins = [-np.inf, threshold, 0.4, 0.6, np.inf]
    labels = [f"<{threshold}", f"{threshold}-0.4", "0.4-0.6", ">=0.6"]
    df["iLIS_bin"] = pd.cut(df["iLIS"], bins=bins, labels=labels, right=False)

    residues_model = list(range(residue_min, residue_max + 1))
    rows = []

    for residue_model in residues_model:
        residue_display = residue_model + offset

        # contains je logični vektor: True za modele, kjer je ta ostanek v cLIR_indices_i.
        contains = df["contact_residues_model"].apply(lambda residue_set: residue_model in residue_set)
        subset = df[contains]
        high_subset = subset[subset["iLIS"] >= threshold]

        row = {
            "residue_model": residue_model,
            "residue_display": residue_display,

            # Osnovna frekvenca.
            "contact_count_all_models": int(len(subset)),

            # Frekvenca samo pri modelih z iLIS nad pragom.
            "contact_count_iLIS_ge_threshold": int(len(high_subset)),

            # Število unikatnih partnerjev; dodatna informacija, ne glavna metrika.
            "unique_partners_all_models": int(subset["partner"].nunique()),
            "unique_partners_iLIS_ge_threshold": int(high_subset["partner"].nunique()),

            # iLIS utežene metrike.
            "sum_iLIS_for_contacting_models": float(subset["iLIS"].sum()),
            "mean_iLIS_for_contacting_models": float(subset["iLIS"].mean()) if len(subset) else 0.0,
            "max_iLIS_for_contacting_models": float(subset["iLIS"].max()) if len(subset) else 0.0,

            # Delež modelov s tem ostankom, ki so nad pragom.
            "fraction_contacting_models_above_threshold": float(len(high_subset) / len(subset)) if len(subset) else 0.0,
        }

        # Frekvence po iLIS razredih.
        for label in labels:
            row[f"count_iLIS_bin_{label}"] = int(((df["iLIS_bin"] == label) & contains).sum())

        rows.append(row)

    summary = pd.DataFrame(rows)

    # Glavna tabela.
    summary_csv = Path(f"{out_prefix}_residue_iLIS_summary.csv")
    summary.to_csv(summary_csv, index=False)

    # Top ostanki po več kriterijih.
    top_by_count = summary.sort_values("contact_count_all_models", ascending=False).head(30)
    top_by_high = summary.sort_values("contact_count_iLIS_ge_threshold", ascending=False).head(30)
    top_by_sum = summary.sort_values("sum_iLIS_for_contacting_models", ascending=False).head(30)

    # Mean iLIS je lahko zavajajoč, če je ostanek prisoten samo v 1-2 modelih.
    # Zato zahtevamo vsaj 20 kontaktnih modelov.
    top_by_mean = (
        summary[summary["contact_count_all_models"] >= 20]
        .sort_values("mean_iLIS_for_contacting_models", ascending=False)
        .head(30)
    )

    top_csv = Path(f"{out_prefix}_top_residues_by_iLIS_metrics.csv")
    with open(top_csv, "w", encoding="utf-8", newline="") as f:
        f.write("# Top residues by all-model contact count\n")
        top_by_count.to_csv(f, index=False)

        f.write("\n# Top residues by count in models with iLIS >= threshold\n")
        top_by_high.to_csv(f, index=False)

        f.write("\n# Top residues by sum of iLIS values\n")
        top_by_sum.to_csv(f, index=False)

        f.write("\n# Top residues by mean iLIS, requiring at least 20 contacting models\n")
        top_by_mean.to_csv(f, index=False)

    # Tabela po iLIS razredih.
    bin_cols = [f"count_iLIS_bin_{label}" for label in labels]
    bin_matrix = summary[["residue_model", "residue_display"] + bin_cols].copy()

    bin_csv = Path(f"{out_prefix}_counts_by_iLIS_bin.csv")
    bin_matrix.to_csv(bin_csv, index=False)

    # Graf 1: heatmap po iLIS razredih.
    heatmap_bin_png = Path(f"{out_prefix}_heatmap_counts_by_iLIS_bin.png")
    mat = bin_matrix[bin_cols].T.values

    plt.figure(figsize=(16, 4.5))
    plt.imshow(mat, aspect="auto", interpolation="nearest")
    plt.yticks(range(len(labels)), labels)

    tick_step = max(1, len(residues_model) // 25)
    plt.xticks(
        ticks=range(0, len(residues_model), tick_step),
        labels=summary["residue_display"].iloc[::tick_step],
        rotation=90,
    )

    plt.xlabel("Residue position")
    plt.ylabel("iLIS bin")
    plt.title("Contact residue frequency by iLIS bin")
    plt.colorbar(label="Number of models")
    plt.tight_layout()
    plt.savefig(heatmap_bin_png, dpi=300)
    plt.close()

    # Graf 2: kontakti uteženi z vsoto iLIS.
    heatmap_sum_png = Path(f"{out_prefix}_heatmap_sum_iLIS_per_residue.png")

    plt.figure(figsize=(16, 2.5))
    plt.imshow(
        summary["sum_iLIS_for_contacting_models"].values.reshape(1, -1),
        aspect="auto",
        interpolation="nearest",
    )
    plt.yticks([0], ["sum(iLIS)"])

    plt.xticks(
        ticks=range(0, len(residues_model), tick_step),
        labels=summary["residue_display"].iloc[::tick_step],
        rotation=90,
    )

    plt.xlabel("Residue position")
    plt.title("Residue contact frequency weighted by iLIS")
    plt.colorbar(label="Sum of iLIS over contacting models")
    plt.tight_layout()
    plt.savefig(heatmap_sum_png, dpi=300)
    plt.close()

    # Graf 3: samo modeli nad pragom.
    bar_high_png = Path(f"{out_prefix}_barplot_contacts_iLIS_ge_threshold.png")

    plt.figure(figsize=(16, 4.5))
    plt.bar(summary["residue_display"], summary["contact_count_iLIS_ge_threshold"], width=1)
    plt.xlabel("Residue position")
    plt.ylabel(f"Number of models with iLIS >= {threshold}")
    plt.title("Contact frequency restricted to significant iLIS models")
    plt.tight_layout()
    plt.savefig(bar_high_png, dpi=300)
    plt.close()

    print(f"Input: {csv_path}")
    print(f"Rows: {len(df):,}")
    print(f"Rows with contacts: {df['has_contacts'].sum():,}")
    print(f"Rows with iLIS >= {threshold}: {(df['iLIS'] >= threshold).sum():,}")
    print()
    print("Top residues by count in models with iLIS >= threshold:")
    print(
        top_by_high[
            [
                "residue_display",
                "residue_model",
                "contact_count_iLIS_ge_threshold",
                "contact_count_all_models",
                "sum_iLIS_for_contacting_models",
                "mean_iLIS_for_contacting_models",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    print()
    print("Created files:")
    for path in [summary_csv, top_csv, bin_csv, heatmap_bin_png, heatmap_sum_png, bar_high_png]:
        print(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count bait contact residues by iLIS bins from AFM-LIS CSV output."
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Input AFM-LIS CSV file.",
    )

    parser.add_argument(
        "--out_prefix",
        required=True,
        help="Output prefix for CSV and PNG files.",
    )

    parser.add_argument(
        "--residue_min",
        type=int,
        required=True,
        help="First bait residue in model numbering.",
    )

    parser.add_argument(
        "--residue_max",
        type=int,
        required=True,
        help="Last bait residue in model numbering.",
    )

    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Offset added to model numbering for displayed/UniProt numbering.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.223,
        help="iLIS threshold. Default: 0.223.",
    )

    args = parser.parse_args()

    analyze(
        csv_path=args.csv,
        out_prefix=args.out_prefix,
        residue_min=args.residue_min,
        residue_max=args.residue_max,
        offset=args.offset,
        threshold=args.threshold,
    )
