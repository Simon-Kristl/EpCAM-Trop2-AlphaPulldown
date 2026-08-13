#!/usr/bin/env python3
"""
Run like this

python remove_dimer_surface_rows.py EpCAM_iLIS_imena.csv \
  --column cLIR_indices_i \
  --indices-file dimer_indices.txt \
  -o EpCAM_no_dimer_surface.csv \
  --keep-removed EpCAM_removed_dimer_surface.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Set

import pandas as pd

DEFAULT_BAIT_COLUMN_CANDIDATES = [
    "cLIR_indices_i",   # correct spelling
    "cLIR_indicies_i",  # tolerated misspelling
]


def sniff_separator(path: str | Path, encoding: str) -> str:
    with open(path, "r", encoding=encoding, errors="replace", newline="") as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters=[",", ";", "\t"]).delimiter
    except csv.Error:
        return ","


def read_csv_flexible(path: str | Path) -> pd.DataFrame:
    encodings = ["utf-8", "utf-8-sig", "cp1250", "iso-8859-2", "windows-1252", "latin1"]
    last_error: Exception | None = None
    for enc in encodings:
        try:
            sep = sniff_separator(path, enc)
            return pd.read_csv(path, encoding=enc, sep=sep, dtype=str, keep_default_na=False)
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Could not read {path} with common encodings/separators") from last_error


def resolve_column(df: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested in df.columns:
            return requested
        raise SystemExit("Column not found: " + requested + "\nAvailable columns:\n" + "\n".join(df.columns))
    for candidate in DEFAULT_BAIT_COLUMN_CANDIDATES:
        if candidate in df.columns:
            return candidate
    raise SystemExit(
        "Could not find bait-side cLIR column. Expected one of: "
        + ", ".join(DEFAULT_BAIT_COLUMN_CANDIDATES)
        + "\nAvailable columns:\n"
        + "\n".join(df.columns)
    )


def parse_index_cell(value) -> Set[int]:
    text = str(value).strip()
    if not text or text.lower() in {"[]", "nan", "none"}:
        return set()
    text = text.strip("[](){}")
    found: Set[int] = set()
    # Split on commas, semicolons, or whitespace. Ranges like 67-71 are expanded.
    for part in re.split(r"[,;\s]+", text):
        if not part:
            continue
        match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", part)
        if match:
            a, b = int(match.group(1)), int(match.group(2))
            if a <= b:
                found.update(range(a, b + 1))
            else:
                found.update(range(b, a + 1))
        elif re.fullmatch(r"\d+", part):
            found.add(int(part))
    return found


def read_indices_from_file(path: str | Path) -> list[int]:
    text = Path(path).read_text(encoding="utf-8")
    return [int(x) for x in re.findall(r"\d+", text)]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove rows where EpCAM/bait-side cLIR_indices_i overlaps supplied residues."
    )
    parser.add_argument("input", help="Input CSV")
    parser.add_argument("-o", "--output", required=True, help="Filtered output CSV")
    parser.add_argument("--column", default=None, help="Column to inspect; default auto-detects cLIR_indices_i")
    parser.add_argument("--indices", nargs="*", type=int, default=[], help="Residues to remove")
    parser.add_argument("--indices-file", help="Text file with residues to remove")
    parser.add_argument("--keep-removed", help="Optional file containing removed rows")
    parser.add_argument(
        "--sep",
        default=",",
        choices=[",", ";", "tab"],
        help="Output separator. Default comma. Use --sep ';' for European spreadsheet settings or --sep tab for TSV.",
    )
    args = parser.parse_args()

    indices = list(args.indices)
    if args.indices_file:
        indices.extend(read_indices_from_file(args.indices_file))
    blocked = set(indices)
    if not blocked:
        raise SystemExit("No indices supplied. Use --indices or --indices-file.")

    df = read_csv_flexible(args.input)
    column = resolve_column(df, args.column)

    remove_mask = df[column].apply(lambda cell: bool(parse_index_cell(cell) & blocked))
    removed = df[remove_mask].copy()
    kept = df[~remove_mask].copy()

    sep = "\t" if args.sep == "tab" else args.sep
    kept.to_csv(args.output, index=False, encoding="utf-8-sig", sep=sep, quoting=csv.QUOTE_ALL, lineterminator="\n")
    if args.keep_removed:
        removed.to_csv(args.keep_removed, index=False, encoding="utf-8-sig", sep=sep, quoting=csv.QUOTE_ALL, lineterminator="\n")

    print(f"Checked column: {column}")
    print(f"Input rows:     {len(df)}")
    print(f"Removed rows:   {len(removed)}")
    print(f"Kept rows:      {len(kept)}")
    print(f"Output:         {args.output}")
    print("Note: output uses UTF-8 with BOM and quotes every field to prevent spreadsheet column-shifting.")
    if args.keep_removed:
        print(f"Removed file:   {args.keep_removed}")


if __name__ == "__main__":
    main()
