#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert residue contact frequency CSV to ChimeraX .defattr file."
    )
    parser.add_argument("--frequency_csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model", default="#1")
    parser.add_argument("--chain", default="A")
    parser.add_argument(
        "--residue_column",
        default="residue",
        help="Column with residue numbers."
    )
    parser.add_argument(
        "--value_column",
        default="contact_count",
        help="Column with values used for coloring."
    )
    parser.add_argument(
        "--attribute_name",
        default="contact_count",
        help="Name of the ChimeraX residue attribute."
    )

    args = parser.parse_args()

    df = pd.read_csv(args.frequency_csv)
    df = df[[args.residue_column, args.value_column]].dropna()

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(f"attribute: {args.attribute_name}\n")
        f.write("recipient: residues\n")
        f.write("match mode: any\n\n")

        for _, row in df.iterrows():
            residue = int(row[args.residue_column])
            value = float(row[args.value_column])
            f.write(f"\t{args.model}/{args.chain}:{residue}\t{value}\n")


if __name__ == "__main__":
    main()