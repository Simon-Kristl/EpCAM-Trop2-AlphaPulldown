import pandas as pd
import re
from pathlib import Path

# =========================
# NASTAVITVE
# =========================

input_file = "EpCAM_iLIS_imena.csv"

# Neželene regije na chain_i.
# Lahko jih je več.
BAD_I = [(225, 242)]


# =========================
# FUNKCIJE
# =========================

def parse_indices(value):
    if pd.isna(value):
        return set()

    value = str(value).strip()
    value = value.strip("[]")
    value = value.replace(" ", "")

    if value == "" or value.lower() == "nan":
        return set()

    indices = set()
    parts = re.split(r"[;,]", value)

    for part in parts:
        if not part:
            continue

        if "-" in part:
            start, end = part.split("-")
            start = int(start)
            end = int(end)
            indices.update(range(start, end + 1))
        else:
            indices.add(int(part))

    return indices


def expand_regions(regions):
    positions = set()
    for start, end in regions:
        positions.update(range(start, end + 1))
    return positions


def classify_region(bad_frac):
    if bad_frac == 0:
        return "clean"
    elif bad_frac <= 0.25:
        return "minor_bad_region"
    elif bad_frac <= 0.50:
        return "mixed"
    elif bad_frac < 1.00:
        return "bad_region_dominated"
    else:
        return "only_bad_region"


def add_region_stats(row, bad_i):
    clir_i = parse_indices(row["cLIR_indices_i"])
    total_clir = len(clir_i)

    if total_clir == 0:
        bad_hits = 0
        bad_frac = 0.0
    else:
        bad_hits = len(clir_i & bad_i)
        bad_frac = bad_hits / total_clir

    regional_cLIR_score = 1 - bad_frac

    return pd.Series({
        "total_cLIR_indices_i": total_clir,
        "bad_hits_i": bad_hits,
        "bad_frac_i": bad_frac,
        "regional_cLIR_score": regional_cLIR_score,
        "region_class": classify_region(bad_frac)
    })


# =========================
# GLAVNI DEL
# =========================

input_path = Path(input_file)

if not input_path.exists():
    raise FileNotFoundError(
        f"Datoteka {input_file} ni najdena. "
        "Preveri, da sta CSV in ta skripta v isti mapi."
    )

df = pd.read_csv(input_path, sep=None, engine="python")

required_columns = ["cLIR_indices_i", "name"]

missing = [col for col in required_columns if col not in df.columns]

if missing:
    raise ValueError(f"Manjkajo stolpci v CSV: {missing}")

bad_i = expand_regions(BAD_I)

df = df.join(
    df.apply(
        lambda row: add_region_stats(row, bad_i),
        axis=1
    )
)

# Sortiranje: najprej regijska čistost, nato originalne metrike
sort_columns = ["regional_cLIR_score"]

for col in ["iLIS", "cLIS", "ipTM"]:
    if col in df.columns:
        sort_columns.append(col)

df_sorted = df.sort_values(sort_columns, ascending=False)

# Najboljši model za vsak proteinski par
best_per_pair = (
    df_sorted
    .groupby("name", as_index=False)
    .head(1)
)

df_sorted.to_csv("all_models_with_regional_cLIR_score.csv", index=False)
best_per_pair.to_csv("best_per_pair_by_regional_cLIR_score.csv", index=False)

print("Končano.")
print(f"Vseh vrstic: {len(df)}")
print(f"Najboljši model na par: {len(best_per_pair)}")
print()
print("Ustvarjeni datoteki:")
print("all_models_with_regional_cLIR_score.csv")
print("best_per_pair_by_regional_cLIR_score.csv")
