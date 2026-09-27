#!/usr/bin/env python3

"""
Reconstruct protein FASTA inputs used for conventional-HDC and HDC2
phylogenetic analyses directly from the archived NCBI protein files.

Expected datasets
-----------------
Conventional HDC: 98 proteins
HDC2:             36 proteins

The reconstructed sequences are also compared with the archived
gene-phylogeny FASTA files to verify exact sequence identity.
"""

from pathlib import Path
import re
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]

RAW_ROOT = (
    ROOT
    / "data/raw/photobacterium_primary_500_ncbi/ncbi_dataset/data"
)

HDC_CANDIDATES = (
    ROOT
    / "results/conventional_HDC_candidate_neighborhoods_2026-09-13.tsv"
)

HDC2_CANDIDATES = (
    ROOT
    / "results/HDC2_candidate_neighborhoods_2026-09-13.tsv"
)

EXCLUSIONS = (
    ROOT
    / "results/ubcg_primary_exclusions_2026-09-13.tsv"
)

ARCHIVED_HDC = (
    ROOT
    / "results/gene_phylogeny/"
      "Photobacterium_conventional_HDC_final497.faa"
)

ARCHIVED_HDC2 = (
    ROOT
    / "results/gene_phylogeny/"
      "Photobacterium_HDC2_final497.faa"
)

OUTDIR = ROOT / "results/reproducibility"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_HDC = (
    OUTDIR
    / "Photobacterium_conventional_HDC_final497_reconstructed.faa"
)

OUT_HDC2 = (
    OUTDIR
    / "Photobacterium_HDC2_final497_reconstructed.faa"
)


def clean_label(text):
    text = "" if pd.isna(text) else str(text)
    return text.replace(" ", "_").replace("|", "_")


def load_proteome(accession):
    path = RAW_ROOT / accession / "protein.faa"

    if not path.exists():
        raise FileNotFoundError(
            f"Missing protein.faa for {accession}"
        )

    records = {}

    for record in SeqIO.parse(path, "fasta"):
        records[record.id] = record

    return records


def reconstruct(table, excluded, outfile, expected_n, system):

    rows = []

    for _, r in table.iterrows():

        accession = str(r["Assembly Accession"])

        if accession in excluded:
            continue

        protein = str(r["Protein"])
        organism = str(r.get("Organism Name", ""))
        strain = str(r.get("Strain", ""))

        rows.append(
            {
                "Assembly": accession,
                "Protein": protein,
                "Organism": organism,
                "Strain": strain,
            }
        )

    # Candidate tables are genome-level.
    df = pd.DataFrame(rows).drop_duplicates(
        subset=["Assembly"]
    )

    if len(df) != expected_n:
        raise RuntimeError(
            f"{system}: expected {expected_n} genomes, "
            f"found {len(df)}"
        )

    output_records = []
    missing = []

    for _, row in df.iterrows():

        accession = row["Assembly"]
        protein = row["Protein"]

        proteome = load_proteome(accession)

        if protein not in proteome:
            missing.append((accession, protein))
            continue

        source = proteome[protein]

        organism_label = clean_label(row["Organism"])
        strain_label = clean_label(row["Strain"])

        fields = [
            accession,
            protein,
            organism_label,
        ]

        if strain_label and strain_label.lower() != "nan":
            fields.append(strain_label)

        source.id = "|".join(fields)
        source.name = source.id
        source.description = ""

        output_records.append(source)

    if missing:
        print("\nMissing target proteins:")
        for item in missing:
            print(" ", item)

        raise RuntimeError(
            f"{system}: {len(missing)} protein(s) missing"
        )

    if len(output_records) != expected_n:
        raise RuntimeError(
            f"{system}: expected {expected_n} sequences, "
            f"reconstructed {len(output_records)}"
        )

    SeqIO.write(
        output_records,
        outfile,
        "fasta"
    )

    return output_records


def protein_id_from_header(record):
    """
    Recover WP_/NP_/YP_ etc. protein accession from archived
    or reconstructed FASTA labels.
    """
    m = re.search(
        r"(?:WP|NP|YP|XP|AP|KP)_?\d+(?:\.\d+)?",
        record.id
    )

    if m:
        return m.group(0)

    # Most project labels use accession|protein|...
    parts = record.id.split("|")

    if len(parts) >= 2:
        return parts[1]

    raise RuntimeError(
        f"Cannot recover protein ID from {record.id}"
    )


def sequence_map(path):
    out = {}

    for record in SeqIO.parse(path, "fasta"):
        pid = protein_id_from_header(record)

        seq = str(record.seq)

        if pid in out and out[pid] != seq:
            raise RuntimeError(
                f"Protein ID {pid} occurs with different "
                f"sequences in {path}"
            )

        out[pid] = seq

    return out


def compare_fastas(reconstructed, archived, name):

    new = sequence_map(reconstructed)
    old = sequence_map(archived)

    print()
    print(f"=== {name} ARCHIVED COMPARISON ===")
    print("Reconstructed unique proteins:", len(new))
    print("Archived unique proteins:     ", len(old))

    missing = sorted(set(old) - set(new))
    extra = sorted(set(new) - set(old))

    mismatched = sorted(
        pid
        for pid in set(new) & set(old)
        if new[pid] != old[pid]
    )

    print("Missing from reconstruction:", len(missing))
    print("Extra in reconstruction:    ", len(extra))
    print("Sequence mismatches:         ", len(mismatched))

    if missing:
        print("Missing:", missing)

    if extra:
        print("Extra:", extra)

    if mismatched:
        print("Mismatched:", mismatched)

    if missing or extra or mismatched:
        raise RuntimeError(
            f"{name}: archived sequence comparison failed"
        )

    print("Exact sequence agreement: YES")


print("=" * 76)
print("RECONSTRUCTING GENE-PHYLOGENY PROTEIN INPUTS")
print("=" * 76)

hdc = pd.read_csv(
    HDC_CANDIDATES,
    sep="\t",
    dtype=str
)

hdc2 = pd.read_csv(
    HDC2_CANDIDATES,
    sep="\t",
    dtype=str
)

ex = pd.read_csv(
    EXCLUSIONS,
    sep="\t",
    dtype=str
)

excluded = set(
    ex["Assembly_Accession"].dropna()
)

print()
print("Candidate genomes")
print("  Conventional HDC:", len(hdc))
print("  HDC2:            ", len(hdc2))
print("  UBCG exclusions: ", len(excluded))

assert len(hdc) == 99
assert len(hdc2) == 36
assert len(excluded) == 3

hdc_records = reconstruct(
    hdc,
    excluded,
    OUT_HDC,
    98,
    "Conventional HDC"
)

hdc2_records = reconstruct(
    hdc2,
    excluded,
    OUT_HDC2,
    36,
    "HDC2"
)

print()
print("Reconstructed datasets")
print("  Conventional HDC:", len(hdc_records))
print("  HDC2:            ", len(hdc2_records))

hdc_lengths = [
    len(x.seq)
    for x in hdc_records
]

hdc2_lengths = [
    len(x.seq)
    for x in hdc2_records
]

print()
print("Sequence lengths")
print(
    "  Conventional HDC:",
    min(hdc_lengths),
    "-",
    max(hdc_lengths),
    "aa"
)
print(
    "  HDC2:",
    min(hdc2_lengths),
    "-",
    max(hdc2_lengths),
    "aa"
)

assert len(hdc_records) == 98
assert len(hdc2_records) == 36

assert min(hdc_lengths) == 374
assert max(hdc_lengths) == 380

assert min(hdc2_lengths) == 680
assert max(hdc2_lengths) == 771

compare_fastas(
    OUT_HDC,
    ARCHIVED_HDC,
    "CONVENTIONAL HDC"
)

compare_fastas(
    OUT_HDC2,
    ARCHIVED_HDC2,
    "HDC2"
)

print()
print("Saved:")
print(" ", OUT_HDC.relative_to(ROOT))
print(" ", OUT_HDC2.relative_to(ROOT))

print()
print("STATUS: PASS")
