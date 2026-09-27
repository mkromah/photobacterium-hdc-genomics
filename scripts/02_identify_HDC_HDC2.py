#!/usr/bin/env python3

"""
Reconstruct conventional HDC and HDC2 genome-level calls.

Inputs
------
metadata/photobacterium_primary_cohort_2026-09-12.tsv
results/conventional_HDC_blastp_all.tsv
results/HDC2_AK3_blastp_all.tsv

Operational candidate rule reconstructed from archived outputs:
    amino-acid identity >= 80%
    query coverage >= 80%

Expected final 497-genome cohort:
    Conventional HDC positive = 98
    HDC2 positive             = 36
    BOTH                      = 3
"""

from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

COHORT = ROOT / "metadata/photobacterium_primary_cohort_2026-09-12.tsv"
EXCLUSIONS = ROOT / "results/ubcg_primary_exclusions_2026-09-13.tsv"
HDC_BLAST = ROOT / "results/conventional_HDC_blastp_all.tsv"
HDC2_BLAST = ROOT / "results/HDC2_AK3_blastp_all.tsv"

OUTDIR = ROOT / "results/reproducibility"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "Photobacterium_primary497_HDC_HDC2_reconstructed.tsv"

PIDENT_MIN = 80.0
QCOV_MIN = 80.0

EXPECTED_HDC_QUERIES = {
    "KP728803.1",
    "KP728804.1",
    "KP728805.1",
}

EXPECTED_HDC2_QUERY = {"WP_065193859.1"}

BLAST_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "evalue",
    "bitscore",
    "qcovs",
    "stitle",
]


def extract_accession(value):
    m = re.search(r"GCF_\d+\.\d+", str(value))
    return m.group(0) if m else None


def find_column(df, candidates):
    lower = {c.lower(): c for c in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    return None


def load_blast(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=BLAST_COLS,
        dtype={"qseqid": str, "sseqid": str},
    )

    for col in [
        "pident",
        "length",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
        "qcovs",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["Assembly Accession"] = df["sseqid"].map(extract_accession)

    return df


print("=" * 72)
print("RECONSTRUCTING HDC / HDC2 GENOME-LEVEL CALLS")
print("=" * 72)

# ---------------------------------------------------------------------
# 1. Load final primary cohort
# ---------------------------------------------------------------------

cohort = pd.read_csv(COHORT, sep="\t", dtype=str)

acc_col = find_column(
    cohort,
    [
        "Assembly Accession",
        "assembly_accession",
        "Accession",
    ],
)

if acc_col is None:
    raise SystemExit(
        "ERROR: could not identify assembly-accession column "
        f"in {COHORT}"
    )

cohort[acc_col] = cohort[acc_col].astype(str)

if cohort[acc_col].duplicated().any():
    duplicates = cohort.loc[
        cohort[acc_col].duplicated(keep=False),
        acc_col,
    ].tolist()

    raise SystemExit(
        "ERROR: duplicate assembly accessions in primary cohort:\n"
        + "\n".join(duplicates[:20])
    )

initial_cohort_set = set(cohort[acc_col])

print("\nInitial primary cohort:", len(initial_cohort_set), "genomes")

if len(initial_cohort_set) != 500:
    raise SystemExit(
        f"ERROR: expected 500 pre-UBCG primary genomes, "
        f"found {len(initial_cohort_set)}"
    )

# Reconstruct the final 497-genome cohort by applying recorded UBCG exclusions.
exclusions = pd.read_csv(EXCLUSIONS, sep=chr(9), dtype=str)

if "Assembly_Accession" not in exclusions.columns:
    raise SystemExit(
        "ERROR: Assembly_Accession column not found in UBCG exclusion table"
    )

excluded_set = set(exclusions["Assembly_Accession"].dropna())

if len(excluded_set) != 3:
    raise SystemExit(
        f"ERROR: expected 3 UBCG exclusions, found {len(excluded_set)}"
    )

missing_exclusions = excluded_set - initial_cohort_set

if missing_exclusions:
    raise SystemExit(
        "ERROR: UBCG exclusion accession(s) absent from 500-genome cohort: "
        + ", ".join(sorted(missing_exclusions))
    )

cohort = cohort.loc[
    ~cohort[acc_col].isin(excluded_set)
].copy()

cohort_set = set(cohort[acc_col])

print("UBCG exclusions:      ", len(excluded_set))
for acc in sorted(excluded_set):
    reason = exclusions.loc[
        exclusions["Assembly_Accession"] == acc,
        "Reason"
    ].iloc[0]
    print(f"  {acc}: {reason}")

print("Final primary cohort: ", len(cohort_set), "genomes")

if len(cohort_set) != 497:
    raise SystemExit(
        f"ERROR: expected 497 final primary genomes, found {len(cohort_set)}"
    )

# ---------------------------------------------------------------------
# 2. Conventional HDC
# ---------------------------------------------------------------------

hdc = load_blast(HDC_BLAST)

observed_hdc_queries = set(hdc["qseqid"].dropna())

if observed_hdc_queries != EXPECTED_HDC_QUERIES:
    raise SystemExit(
        "ERROR: unexpected conventional-HDC query IDs.\n"
        f"Observed: {sorted(observed_hdc_queries)}"
    )

hdc_pass = hdc[
    (hdc["pident"] >= PIDENT_MIN)
    & (hdc["qcovs"] >= QCOV_MIN)
].copy()

hdc_candidates = set(
    hdc_pass["Assembly Accession"].dropna()
)

hdc_final = hdc_candidates & cohort_set

print("\nConventional HDC")
print("  query proteins:       ", ", ".join(sorted(EXPECTED_HDC_QUERIES)))
print("  raw BLAST rows:       ", len(hdc))
print("  candidates at 80/80:  ", len(hdc_candidates))
print("  final cohort positive:", len(hdc_final))
print("  outside final cohort: ", len(hdc_candidates - cohort_set))

if hdc_candidates - cohort_set:
    print(
        "  excluded accession(s):",
        ", ".join(sorted(hdc_candidates - cohort_set))
    )

# ---------------------------------------------------------------------
# 3. HDC2
# ---------------------------------------------------------------------

hdc2 = load_blast(HDC2_BLAST)

observed_hdc2_queries = set(hdc2["qseqid"].dropna())

if observed_hdc2_queries != EXPECTED_HDC2_QUERY:
    raise SystemExit(
        "ERROR: unexpected HDC2 query ID.\n"
        f"Observed: {sorted(observed_hdc2_queries)}"
    )

hdc2_pass = hdc2[
    (hdc2["pident"] >= PIDENT_MIN)
    & (hdc2["qcovs"] >= QCOV_MIN)
].copy()

hdc2_candidates = set(
    hdc2_pass["Assembly Accession"].dropna()
)

hdc2_final = hdc2_candidates & cohort_set

print("\nHDC2")
print("  query protein:        ", next(iter(EXPECTED_HDC2_QUERY)))
print("  raw BLAST rows:       ", len(hdc2))
print("  candidates at 80/80:  ", len(hdc2_candidates))
print("  final cohort positive:", len(hdc2_final))
print("  outside final cohort: ", len(hdc2_candidates - cohort_set))

# ---------------------------------------------------------------------
# 4. Build genome-level classification
# ---------------------------------------------------------------------

out = cohort.copy()

out["Conventional_HDC"] = out[acc_col].isin(hdc_final)
out["HDC2"] = out[acc_col].isin(hdc2_final)

def classify(row):
    hdc_pos = bool(row["Conventional_HDC"])
    hdc2_pos = bool(row["HDC2"])

    if hdc_pos and hdc2_pos:
        return "BOTH"
    if hdc_pos:
        return "HDC_ONLY"
    if hdc2_pos:
        return "HDC2_ONLY"
    return "NEITHER"

out["HDC_System"] = out.apply(classify, axis=1)

counts = out["HDC_System"].value_counts()

print("\nGenome-level states")
for state in ["HDC_ONLY", "HDC2_ONLY", "BOTH", "NEITHER"]:
    print(f"  {state:10s}: {int(counts.get(state, 0))}")

print("\nFinal totals")
print("  Conventional HDC:", int(out["Conventional_HDC"].sum()))
print("  HDC2:            ", int(out["HDC2"].sum()))
print("  BOTH:            ", int((out["HDC_System"] == "BOTH").sum()))

# ---------------------------------------------------------------------
# 5. Hard validation against locked project totals
# ---------------------------------------------------------------------

assert len(out) == 497
assert int(out["Conventional_HDC"].sum()) == 98
assert int(out["HDC2"].sum()) == 36
assert int((out["HDC_System"] == "BOTH").sum()) == 3

# Expected state totals implied by the final project:
assert int((out["HDC_System"] == "HDC_ONLY").sum()) == 95
assert int((out["HDC_System"] == "HDC2_ONLY").sum()) == 33
assert int((out["HDC_System"] == "NEITHER").sum()) == 366

out.to_csv(OUTFILE, sep="\t", index=False)

print("\nSaved:")
print(OUTFILE.relative_to(ROOT))
print("\nSTATUS: PASS")
