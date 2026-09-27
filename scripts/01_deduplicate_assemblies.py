import pandas as pd
from pathlib import Path

INPUT = Path("metadata/photobacterium_candidate_manifest_2026-09-12.tsv")
OUTPUT = Path("metadata/photobacterium_unique_assemblies_2026-09-12.tsv")
DUPLICATES = Path("metadata/photobacterium_removed_paired_duplicates_2026-09-12.tsv")

df = pd.read_csv(INPUT, sep="\t", dtype=str).fillna("")

print(f"Input assembly records: {len(df)}")

# Keep current assemblies only
current = df[
    df["Assembly Status"].str.lower().eq("current")
].copy()

print(f"Current assembly records: {len(current)}")

# A GCA and its paired GCF represent the same underlying assembly.
# Build one common identifier for each pair.
def make_pair_key(row):
    accession = row["Assembly Accession"].strip()
    paired = row["Assembly Paired Assembly Accession"].strip()

    if paired:
        return "|".join(sorted([accession, paired]))

    return accession

current["Assembly Pair Key"] = current.apply(make_pair_key, axis=1)

# Prefer RefSeq over GenBank when both representations exist.
priority = {
    "SOURCE_DATABASE_REFSEQ": 0,
    "SOURCE_DATABASE_GENBANK": 1
}

current["Database Priority"] = (
    current["Source Database"]
    .map(priority)
    .fillna(2)
)

current = current.sort_values(
    ["Assembly Pair Key", "Database Priority", "Assembly Accession"]
)

# Keep one preferred assembly from every GCA/GCF pair
unique = current.drop_duplicates(
    subset="Assembly Pair Key",
    keep="first"
).copy()

# Save records removed only because they were paired duplicates
removed = current.loc[~current.index.isin(unique.index)].copy()

unique.drop(columns=["Database Priority"]).to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

removed.drop(columns=["Database Priority"]).to_csv(
    DUPLICATES,
    sep="\t",
    index=False
)

print(f"Unique biological assemblies retained: {len(unique)}")
print(f"Paired GCA/GCF duplicates removed: {len(removed)}")

print("\nRetained source databases:")
print(unique["Source Database"].value_counts())

print("\nRetained assembly levels:")
print(unique["Assembly Level"].value_counts())

print("\nTop 20 organisms after deduplication:")
print(unique["Organism Name"].value_counts().head(20))

print("\nCheckM completeness available:")
print(unique["CheckM completeness"].replace("", pd.NA).notna().sum())

print("\nCheckM contamination available:")
print(unique["CheckM contamination"].replace("", pd.NA).notna().sum())
