#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


TOKEN_RE = re.compile(r"\d+\s*(?:-\s*\d+)?")


def parse_indices(cell):
    """
    Parse strings like:
    [120,122,124,162,164,167,169,225,227-230,232-236,238]
    into a sorted set of residue numbers.
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
    Extract partner from names like:
    P16422_and_Q9Y5Y6
    P09758_and_Q9Y5Y6
    """
    text = str(name)
    if "_and_" in text:
        return text.split("_and_", 1)[1]
    return text


def main():
    parser = argparse.ArgumentParser(
        description="Create residue contact heatmaps from cLIR_indices_i."
    )

    parser.add_argument("--csv", required=True, help="Input AFM-LIS CSV file.")
    parser.add_argument("--out_prefix", required=True, help="Output file prefix.")
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Offset added to model residue numbers for UniProt numbering. EpCAM: 23, Trop2: 26.",
    )
    parser.add_argument(
        "--min_ilis",
        type=float,
        default=None,
        help="Optional minimum iLIS threshold, for example 0.223.",
    )
    parser.add_argument(
        "--top_model_per_partner",
        action="store_true",
        help="Keep only the highest-iLIS model for each partner.",
    )

    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    if "cLIR_indices_i" not in df.columns:
        raise ValueError("Input CSV must contain column: cLIR_indices_i")

    if "iLIS" not in df.columns:
        raise ValueError("Input CSV must contain column: iLIS")

    if "name" not in df.columns:
        raise ValueError("Input CSV must contain column: name")

    df["iLIS"] = pd.to_numeric(df["iLIS"], errors="coerce").fillna(0)
    df["partner"] = df["name"].apply(partner_from_name)

    if args.min_ilis is not None:
        df = df[df["iLIS"] >= args.min_ilis].copy()

    if args.top_model_per_partner:
        df = (
            df.sort_values("iLIS", ascending=False)
              .drop_duplicates("partner", keep="first")
              .copy()
        )

    df["contact_residues"] = df["cLIR_indices_i"].apply(parse_indices)

    # Apply UniProt offset
    df["contact_residues_offset"] = df["contact_residues"].apply(
        lambda residues: {r + args.offset for r in residues}
    )

    all_residues = sorted(set().union(*df["contact_residues_offset"]))

    if not all_residues:
        raise ValueError("No contact residues found after filtering.")

    # Build binary matrix: rows = partners, columns = residues
    matrix = pd.DataFrame(0, index=df["partner"], columns=all_residues)

    for partner, residues in zip(df["partner"], df["contact_residues_offset"]):
        for residue in residues:
            matrix.loc[partner, residue] = 1

    # Save binary matrix
    matrix.to_csv(f"{args.out_prefix}_binary_contact_matrix.csv")

    # Frequency table
    frequency = matrix.sum(axis=0).reset_index()
    frequency.columns = ["residue", "contact_count"]
    frequency["contact_fraction"] = frequency["contact_count"] / len(matrix)
    frequency.to_csv(f"{args.out_prefix}_contact_frequency.csv", index=False)

    # Plot heatmap protein × residue
    plt.figure(figsize=(16, max(4, 0.25 * len(matrix))))
    plt.imshow(matrix.values, aspect="auto", interpolation="nearest")
    plt.yticks(range(len(matrix.index)), matrix.index)
    plt.xticks(
        ticks=range(0, len(matrix.columns), max(1, len(matrix.columns) // 20)),
        labels=matrix.columns[::max(1, len(matrix.columns) // 20)],
        rotation=90,
    )
    plt.xlabel("Residue position on bait protein")
    plt.ylabel("Target protein")
    plt.title("cLIR_indices_i contact heatmap")
    plt.colorbar(label="Contact present")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}_protein_by_residue_heatmap.png", dpi=300)
    plt.close()

    # Plot frequency heatmap as one-row heatmap
    freq_values = frequency["contact_count"].values.reshape(1, -1)

    plt.figure(figsize=(16, 2.5))
    plt.imshow(freq_values, aspect="auto", interpolation="nearest")
    plt.yticks([0], ["Contact count"])
    plt.xticks(
        ticks=range(0, len(frequency), max(1, len(frequency) // 25)),
        labels=frequency["residue"][::max(1, len(frequency) // 25)],
        rotation=90,
    )
    plt.xlabel("Residue position on bait protein")
    plt.title("Frequency of cLIR_indices_i contacts")
    plt.colorbar(label="Number of models")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}_residue_frequency_heatmap.png", dpi=300)
    plt.close()

    print(f"Rows used: {len(df)}")
    print(f"Residues plotted: {len(all_residues)}")
    print(f"Saved: {args.out_prefix}_protein_by_residue_heatmap.png")
    print(f"Saved: {args.out_prefix}_residue_frequency_heatmap.png")
    print(f"Saved: {args.out_prefix}_contact_frequency.csv")


if __name__ == "__main__":
    main()