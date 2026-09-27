#!/usr/bin/env python3

"""
Reconstruct taxonomic distribution and association statistics for
conventional HDC and HDC2 across the final 497 Photobacterium genomes.

Inputs
------
results/reproducibility/Photobacterium_primary497_HDC_HDC2_reconstructed.tsv
results/phylogeny/itol_taxon_groups_colorstrip.txt

Analyses
--------
1. Recover the archived 15-group taxonomic classification.
2. Construct the 15 x 4 taxon-by-HDC-system contingency table.
3. Calculate Pearson chi-square statistic.
4. Calculate Cramer's V.
5. Perform 100,000 Monte Carlo permutations preserving:
      - taxon group sizes
      - global HDC-system state totals
6. Calculate adjusted standardized residuals.

Random seed
-----------
20260919
"""

from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

CALLS_FILE = (
    ROOT
    / "results/reproducibility/"
      "Photobacterium_primary497_HDC_HDC2_reconstructed.tsv"
)

TAXON_FILE = (
    ROOT
    / "results/phylogeny/itol_taxon_groups_colorstrip.txt"
)

OUTDIR = ROOT / "results/reproducibility"
OUTDIR.mkdir(parents=True, exist_ok=True)

SUMMARY_OUT = (
    OUTDIR
    / "Photobacterium_primary497_taxon_HDC_summary_reconstructed.tsv"
)

PREVALENCE_OUT = (
    OUTDIR
    / "Photobacterium_primary497_taxon_HDC_prevalence_reconstructed.tsv"
)

RESIDUAL_OUT = (
    OUTDIR
    / "Photobacterium_primary497_taxon_HDC_residuals_reconstructed.tsv"
)

TEST_OUT = (
    OUTDIR
    / "Photobacterium_primary497_taxon_HDC_permutation_test_reconstructed.txt"
)

B = 100000
SEED = 20260919

STATE_ORDER = [
    "HDC_ONLY",
    "HDC2_ONLY",
    "BOTH",
    "NEITHER",
]


def load_taxon_mapping(path):
    rows = []
    data_started = False

    for line in path.read_text().splitlines():

        if line.strip() == "DATA":
            data_started = True
            continue

        if not data_started or not line.strip():
            continue

        parts = line.split("\t")

        if (
            len(parts) >= 3
            and (
                parts[0].startswith("GCF_")
                or parts[0].startswith("GCA_")
            )
        ):
            rows.append(
                {
                    "Assembly Accession": parts[0].strip(),
                    "Taxon_Group": parts[2].strip(),
                    "Taxon_Color": parts[1].strip(),
                }
            )

    return pd.DataFrame(rows)


def chi_square_stat(table):
    observed = table.astype(float)

    row_totals = observed.sum(axis=1, keepdims=True)
    col_totals = observed.sum(axis=0, keepdims=True)
    total = observed.sum()

    expected = row_totals @ col_totals / total

    statistic = np.sum(
        (observed - expected) ** 2 / expected
    )

    return statistic, expected


print("=" * 76)
print("TAXONOMIC DISTRIBUTION AND HDC-SYSTEM ASSOCIATION")
print("=" * 76)

calls = pd.read_csv(
    CALLS_FILE,
    sep="\t",
    dtype=str
)

taxa = load_taxon_mapping(TAXON_FILE)

print()
print("Input genomes:", len(calls))
print("Taxon mappings:", len(taxa))

assert len(calls) == 497
assert len(taxa) == 497
assert taxa["Assembly Accession"].nunique() == 497

df = calls.merge(
    taxa,
    on="Assembly Accession",
    how="left",
    validate="one_to_one",
)

assert len(df) == 497
assert df["Taxon_Group"].notna().all()

print("Taxonomic groups:", df["Taxon_Group"].nunique())
assert df["Taxon_Group"].nunique() == 15

# ------------------------------------------------------------
# Contingency table
# ------------------------------------------------------------

table = pd.crosstab(
    df["Taxon_Group"],
    df["HDC_System"]
)

table = table.reindex(
    columns=STATE_ORDER,
    fill_value=0
)

table["Total"] = table.sum(axis=1)

print()
print("=== TAXON x HDC-SYSTEM COUNTS ===")
print(table.to_string())

# Locked taxonomic totals
expected_taxon_totals = {
    "Other Photobacterium": 53,
    "P. angustum": 28,
    "P. aquimaris": 9,
    "P. carnosum": 21,
    "P. damselae": 56,
    "P. damselae subsp. damselae": 52,
    "P. damselae subsp. piscicida": 14,
    "P. iliopiscarium": 18,
    "P. kishitanii": 24,
    "P. leiognathi": 10,
    "P. leiognathi subsp. mandapamensis": 108,
    "P. phosphoreum": 60,
    "P. piscicola": 14,
    "P. swingsii": 8,
    "Photobacterium sp.": 22,
}

for group, expected_n in expected_taxon_totals.items():
    observed_n = int(table.loc[group, "Total"])
    assert observed_n == expected_n, (
        group,
        observed_n,
        expected_n,
    )

# Locked system totals
assert int(table["HDC_ONLY"].sum()) == 95
assert int(table["HDC2_ONLY"].sum()) == 33
assert int(table["BOTH"].sum()) == 3
assert int(table["NEITHER"].sum()) == 366
assert int(table["Total"].sum()) == 497

# ------------------------------------------------------------
# Save counts + prevalence
# ------------------------------------------------------------

summary = table.reset_index()
summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False,
)

prevalence = summary.copy()

for state in STATE_ORDER:
    prevalence[state + "_pct"] = (
        100
        * prevalence[state]
        / prevalence["Total"]
    )

prevalence.to_csv(
    PREVALENCE_OUT,
    sep="\t",
    index=False,
)

# ------------------------------------------------------------
# Pearson chi-square
# ------------------------------------------------------------

observed = table[STATE_ORDER].to_numpy(dtype=float)

chi2, expected = chi_square_stat(observed)

n = observed.sum()
r, c = observed.shape

cramers_v = np.sqrt(
    chi2
    / (
        n
        * min(r - 1, c - 1)
    )
)

print()
print("=== OBSERVED ASSOCIATION ===")
print(f"Chi-square:  {chi2:.4f}")
print(f"df:          {(r - 1) * (c - 1)}")
print(f"Cramer's V: {cramers_v:.4f}")

assert abs(chi2 - 509.0521) < 0.001
assert abs(cramers_v - 0.5843) < 0.0001

# Expected-count diagnostics
below5 = int((expected < 5).sum())
below1 = int((expected < 1).sum())
cells = expected.size

print()
print("Expected-count diagnostics")
print(
    f"  cells <5: {below5}/{cells} "
    f"({100 * below5 / cells:.2f}%)"
)
print(
    f"  cells <1: {below1}/{cells} "
    f"({100 * below1 / cells:.2f}%)"
)

assert below5 == 38
assert below1 == 20

# ------------------------------------------------------------
# Monte Carlo permutation test
# ------------------------------------------------------------

print()
print(
    f"Running {B:,} permutations "
    f"(seed={SEED})..."
)

rng = np.random.default_rng(SEED)

taxon_codes, taxon_levels = pd.factorize(
    df["Taxon_Group"],
    sort=True
)

states = df["HDC_System"].to_numpy()

state_to_int = {
    state: i
    for i, state in enumerate(STATE_ORDER)
}

state_codes = np.array(
    [state_to_int[x] for x in states],
    dtype=np.int16,
)

n_taxa = len(taxon_levels)
n_states = len(STATE_ORDER)

ge_observed = 0

for _ in range(B):

    permuted = rng.permutation(state_codes)

    perm_table = np.zeros(
        (n_taxa, n_states),
        dtype=np.int32,
    )

    np.add.at(
        perm_table,
        (taxon_codes, permuted),
        1,
    )

    perm_chi2, _ = chi_square_stat(
        perm_table
    )

    if perm_chi2 >= chi2:
        ge_observed += 1

mc_p = (
    ge_observed + 1
) / (
    B + 1
)

print()
print("=== MONTE CARLO PERMUTATION TEST ===")
print("Permutations:", B)
print("Permuted >= observed:", ge_observed)
print(f"Monte Carlo p-value: {mc_p:.7g}")

assert ge_observed == 0
assert abs(mc_p - (1 / 100001)) < 1e-15

# ------------------------------------------------------------
# Adjusted standardized residuals
# ------------------------------------------------------------

row_props = (
    observed.sum(axis=1, keepdims=True)
    / n
)

col_props = (
    observed.sum(axis=0, keepdims=True)
    / n
)

adjusted_residuals = (
    observed - expected
) / np.sqrt(
    expected
    * (1 - row_props)
    * (1 - col_props)
)

residual_rows = []

for i, taxon in enumerate(table.index):
    for j, state in enumerate(STATE_ORDER):

        residual_rows.append(
            {
                "Taxon_Group": taxon,
                "HDC_System": state,
                "Observed": int(observed[i, j]),
                "Expected": expected[i, j],
                "Adjusted_Standardized_Residual":
                    adjusted_residuals[i, j],
            }
        )

residual_df = pd.DataFrame(residual_rows)

residual_df["Abs_Residual"] = (
    residual_df[
        "Adjusted_Standardized_Residual"
    ].abs()
)

residual_df = residual_df.sort_values(
    "Abs_Residual",
    ascending=False,
)

residual_df.to_csv(
    RESIDUAL_OUT,
    sep="\t",
    index=False,
)

print()
print("=== LARGEST ADJUSTED STANDARDIZED RESIDUALS ===")

print(
    residual_df[
        [
            "Taxon_Group",
            "HDC_System",
            "Observed",
            "Expected",
            "Adjusted_Standardized_Residual",
        ]
    ]
    .head(12)
    .to_string(
        index=False,
        formatters={
            "Expected": "{:.2f}".format,
            "Adjusted_Standardized_Residual":
                "{:.2f}".format,
        },
    )
)

# ------------------------------------------------------------
# Save statistical report
# ------------------------------------------------------------

with open(TEST_OUT, "w") as out:

    out.write(
        "Photobacterium taxon x HDC-system association\n"
    )
    out.write("=" * 50 + "\n\n")

    out.write(f"Genomes: {int(n)}\n")
    out.write(f"Taxon groups: {r}\n")
    out.write(f"HDC-system states: {c}\n\n")

    out.write(
        f"Pearson chi-square: {chi2:.6f}\n"
    )
    out.write(
        f"Degrees of freedom: {(r - 1) * (c - 1)}\n"
    )
    out.write(
        f"Cramer's V: {cramers_v:.6f}\n\n"
    )

    out.write(
        f"Expected cells <5: "
        f"{below5}/{cells} "
        f"({100 * below5 / cells:.2f}%)\n"
    )

    out.write(
        f"Expected cells <1: "
        f"{below1}/{cells} "
        f"({100 * below1 / cells:.2f}%)\n\n"
    )

    out.write(
        f"Monte Carlo permutations: {B}\n"
    )
    out.write(
        f"Random seed: {SEED}\n"
    )
    out.write(
        f"Permuted statistics >= observed: "
        f"{ge_observed}\n"
    )
    out.write(
        f"Monte Carlo p-value: {mc_p:.10g}\n"
    )

print()
print("Saved:")
print(" ", SUMMARY_OUT.relative_to(ROOT))
print(" ", PREVALENCE_OUT.relative_to(ROOT))
print(" ", RESIDUAL_OUT.relative_to(ROOT))
print(" ", TEST_OUT.relative_to(ROOT))

print()
print("STATUS: PASS")
