#!/usr/bin/env python3

"""
Run like this

python filter_csv_by_tsv_uniprot.py EpCAM_iLIS_imena.csv 2026-07-02-18-30.tsv -o filtered_output.csv
"""
import argparse
import csv
import re
from pathlib import Path

import pandas as pd


# UniProt accession pattern, including optional isoform suffixes like P16422-2.
UNIPROT_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)(?:-\d+)?\b",
    re.IGNORECASE,
)


def extract_uniprot_accessions(text):
    """
    Extract UniProt accessions from any text.

    Examples:
        uniprotkb:P16422 -> P16422
        P16422-2 -> P16422-2 and P16422
    """
    if text is None:
        return set()

    found = set()
    for match in UNIPROT_RE.finditer(str(text)):
        acc = match.group(0).upper().strip()
        if acc:
            found.add(acc)
            # Also add base accession so P12345-2 can match P12345.
            if "-" in acc:
                found.add(acc.split("-", 1)[0])

    return found


def clean_cell(value):
    """
    Keep values as text and prevent hidden newlines from making a CSV appear broken.
    """
    if value is None:
        return ""
    value = str(value)
    value = value.replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    return value.strip()


def clean_output_columns(df, keep_empty_unnamed=False):
    """
    Make the output easier to open in Excel/LibreOffice.

    The input file has an empty Unnamed: 3 column. Keeping a fully empty column can
    make rows look shifted in spreadsheet viewers, especially when the description
    column is sometimes blank. By default this drops only fully empty Unnamed columns.
    """
    df = df.copy()

    # Rename the populated description column.
    if "Unnamed: 2" in df.columns and "description" not in df.columns:
        df = df.rename(columns={"Unnamed: 2": "description"})

    if not keep_empty_unnamed:
        empty_unnamed_cols = []
        for col in df.columns:
            col_name = str(col)
            if col_name.startswith("Unnamed"):
                is_empty = df[col].map(clean_cell).eq("").all()
                if is_empty:
                    empty_unnamed_cols.append(col)
        if empty_unnamed_cols:
            df = df.drop(columns=empty_unnamed_cols)

    # Clean all cell contents but preserve all values as text.
    for col in df.columns:
        df[col] = df[col].map(clean_cell)

    return df


def validate_written_csv(path):
    """
    Confirm every written row has exactly the same number of fields as the header.
    This catches true column-shift/misalignment problems.
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        lengths = []
        for row_number, row in enumerate(reader, start=1):
            lengths.append((row_number, len(row)))

    if not lengths:
        raise ValueError(f"Output file is empty: {path}")

    expected = lengths[0][1]
    bad_rows = [(row_number, count) for row_number, count in lengths if count != expected]

    if bad_rows:
        preview = ", ".join(f"line {r}: {c} fields" for r, c in bad_rows[:10])
        raise ValueError(
            f"Output CSV validation failed. Expected {expected} fields per row. "
            f"Mismatched rows: {preview}"
        )

    return expected, len(lengths)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Keep CSV rows whose UniProt accession is present in a TSV, "
            "and write a spreadsheet-safe CSV with no shifted rows."
        )
    )

    parser.add_argument("csv_input", help="Input CSV file")
    parser.add_argument("tsv_input", help="Input TSV file")
    parser.add_argument(
        "-c",
        "--csv-uniprot-column",
        default="protein",
        help="CSV column containing UniProt accessions. Default: protein",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="rows_in_both_no_shift.csv",
        help="Output CSV file. Default: rows_in_both_no_shift.csv",
    )
    parser.add_argument(
        "--keep-empty-unnamed-columns",
        action="store_true",
        help=(
            "Keep fully empty Unnamed columns instead of dropping them. "
            "By default they are removed to avoid apparent left-shift in spreadsheets."
        ),
    )

    args = parser.parse_args()

    csv_df = pd.read_csv(args.csv_input, dtype=str, keep_default_na=False)
    tsv_df = pd.read_csv(args.tsv_input, sep="\t", dtype=str, keep_default_na=False)

    if args.csv_uniprot_column not in csv_df.columns:
        available = ", ".join(map(str, csv_df.columns))
        raise ValueError(
            f"CSV column {args.csv_uniprot_column!r} not found. "
            f"Available columns: {available}"
        )

    # Extract all UniProt accessions from the TSV.
    tsv_uniprots = set()
    for col in tsv_df.columns:
        for value in tsv_df[col]:
            tsv_uniprots.update(extract_uniprot_accessions(value))

    # Match CSV rows by UniProt accession in the selected column.
    def csv_cell_matches_tsv(value):
        accessions = extract_uniprot_accessions(value)
        if not accessions:
            accessions = {clean_cell(value).upper()}
        return bool(accessions & tsv_uniprots)

    match_mask = csv_df[args.csv_uniprot_column].map(csv_cell_matches_tsv)
    filtered_df = csv_df.loc[match_mask].copy()

    filtered_df = clean_output_columns(
        filtered_df,
        keep_empty_unnamed=args.keep_empty_unnamed_columns,
    )

    # The important part:
    # QUOTE_ALL keeps comma-containing range cells like [1,2,5-9] inside one cell.
    # utf-8-sig makes the file easier for Excel to detect as UTF-8.
    filtered_df.to_csv(
        args.output,
        index=False,
        encoding="utf-8-sig",
        quoting=csv.QUOTE_ALL,
        doublequote=True,
        lineterminator="\n",
    )

    field_count, physical_rows = validate_written_csv(args.output)

    print(f"CSV rows read: {len(csv_df)}")
    print(f"TSV rows read: {len(tsv_df)}")
    print(f"Unique UniProt accessions found in TSV: {len(tsv_uniprots)}")
    print(f"Matching CSV rows written: {len(filtered_df)}")
    print(f"Output columns written: {field_count}")
    print(f"CSV validation: passed, all {physical_rows} lines have {field_count} fields")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
