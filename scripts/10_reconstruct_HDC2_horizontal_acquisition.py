#!/usr/bin/env python3

"""
Reconstruct evidence for recent horizontal acquisition of an HDC2-associated
genomic region shared by:

  Photobacterium phosphoreum FS3_1  GCF_001676275.1
  Photobacterium toruni S2MW2       GCF_036040135.1

Analyses:
  1. HDC2 protein identity
  2. Ranking among cross-species HDC2 comparisons
  3. Reciprocal whole-genome ANI
  4. 16-kb HDC2-locus BLASTN comparison
  5. Full HDC2-bearing contig comparison
  6. Distribution of the FS3_1 cassette across all 36 HDC2-positive genomes
  7. Per-position prevalence across the 16-kb cassette

Interpretation is deliberately conservative: the analysis can support recent
horizontal acquisition/shared acquisition, but cannot determine transfer
direction or exclude an unsampled donor.
"""

from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]

RAW = (
    ROOT
    / "data/raw/photobacterium_primary_500_ncbi/ncbi_dataset/data"
)

HDC2_CANDIDATES = (
    ROOT
    / "results/HDC2_candidate_neighborhoods_2026-09-13.tsv"
)

ANNOTATION = (
    ROOT
    / "results/phylogeny/"
      "photobacterium_primary497_HDC_HDC2_annotation.tsv"
)

ALIGNMENT = (
    ROOT
    / "results/reproducibility/"
      "Photobacterium_HDC2_final497_reconstructed_aligned.faa"
)

BASE = (
    ROOT
    / "results/gene_phylogeny/HDC2_FS3_S2MW2_DNA"
)

FS3_WINDOW = (
    BASE
    / "P_phosphoreum_FS3_1_HDC2_window.fna"
)

S2_WINDOW = (
    BASE
    / "P_toruni_S2MW2_HDC2_window.fna"
)

FS3_CONTIG = (
    BASE
    / "full_contigs/FS3_1_HDC2_contig_normalized.fna"
)

S2_CONTIG = (
    BASE
    / "full_contigs/S2MW2_HDC2_contig_normalized.fna"
)

ARCHIVED_DISTRIBUTION = (
    BASE
    / "FS3_1_cassette_vs_all_HDC2_genomes.tsv"
)

OUT = (
    ROOT
    / "results/reproducibility/HDC2_horizontal_acquisition"
)

BLAST_OUT = OUT / "cassette_blastn"

OUT.mkdir(parents=True, exist_ok=True)
BLAST_OUT.mkdir(parents=True, exist_ok=True)

FS3 = "GCF_001676275.1"
S2 = "GCF_036040135.1"

ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


def find_accession_column(df):
    best = None
    best_n = -1

    for col in df.columns:
        n = (
            df[col]
            .astype(str)
            .str.contains(
                r"GC[AF]_\d+\.\d+",
                regex=True,
                na=False
            )
            .sum()
        )

        if n > best_n:
            best = col
            best_n = n

    if best is None or best_n == 0:
        raise RuntimeError(
            "Could not identify assembly-accession column"
        )

    return best


def find_organism_column(df):
    candidates = []

    for col in df.columns:
        n = (
            df[col]
            .astype(str)
            .str.contains(
                "Photobacterium",
                case=False,
                regex=False,
                na=False
            )
            .sum()
        )

        candidates.append((n, col))

    n, col = max(candidates)

    if n == 0:
        raise RuntimeError(
            "Could not identify organism-name column"
        )

    return col


def accession_from_text(text):
    m = ACC_RE.search(str(text))
    return m.group(0) if m else None


def genome_fasta(accession):
    folder = RAW / accession
    files = sorted(folder.glob("*_genomic.fna"))

    if len(files) != 1:
        raise RuntimeError(
            f"{accession}: expected one genomic FASTA; "
            f"found {len(files)}"
        )

    return files[0]


def fasta_length(path):
    return sum(
        len(record.seq)
        for record in SeqIO.parse(path, "fasta")
    )


def aligned_identity(a, b):
    if len(a) != len(b):
        raise RuntimeError(
            "Alignment sequences have unequal lengths"
        )

    identical = 0
    comparable = 0

    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue

        comparable += 1

        if x == y:
            identical += 1

    if comparable == 0:
        return np.nan, 0, 0

    return (
        100.0 * identical / comparable,
        identical,
        comparable,
    )


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def run_blastn(query, subject, outfile):
    subprocess.run(
        [
            "blastn",
            "-query", str(query),
            "-subject", str(subject),
            "-outfmt",
            (
                "6 qseqid sseqid pident length mismatch "
                "gapopen qstart qend sstart send evalue "
                "bitscore qlen slen"
            ),
            "-out", str(outfile),
        ],
        check=True,
    )


def read_blast(path):
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    return pd.read_csv(
        path,
        sep="\t",
        names=BLAST_COLUMNS,
    )


def merged_length(intervals):
    if not intervals:
        return 0

    intervals = sorted(
        (
            min(int(a), int(b)),
            max(int(a), int(b)),
        )
        for a, b in intervals
    )

    start, end = intervals[0]
    total = 0

    for a, b in intervals[1:]:
        if a <= end + 1:
            end = max(end, b)
        else:
            total += end - start + 1
            start, end = a, b

    total += end - start + 1
    return total


def blast_summary(df, min_identity):
    selected = df.loc[
        df["pident"] >= min_identity
    ].copy()

    if selected.empty:
        return 0, 0.0, np.nan

    qlen = int(selected["qlen"].iloc[0])

    covered = merged_length(
        list(
            zip(
                selected["qstart"],
                selected["qend"],
            )
        )
    )

    coverage = 100.0 * covered / qlen

    weighted_identity = (
        (
            selected["pident"]
            * selected["length"]
        ).sum()
        /
        selected["length"].sum()
    )

    return covered, coverage, weighted_identity


def run_fastani(query, ref, outfile):
    subprocess.run(
        [
            "fastANI",
            "--query", str(query),
            "--ref", str(ref),
            "--output", str(outfile),
        ],
        check=True,
    )

    fields = (
        outfile
        .read_text()
        .strip()
        .split("\t")
    )

    if len(fields) != 5:
        raise RuntimeError(
            "Unexpected fastANI output"
        )

    return {
        "ANI": float(fields[2]),
        "Matched": int(fields[3]),
        "Total": int(fields[4]),
    }


print("=" * 78)
print("HDC2 HORIZONTAL-ACQUISITION RECONSTRUCTION")
print("=" * 78)

for executable in ["blastn", "fastANI"]:
    location = shutil.which(executable)

    if location is None:
        raise RuntimeError(
            f"{executable} not found"
        )

    print(f"{executable}: {location}")

required = [
    HDC2_CANDIDATES,
    ANNOTATION,
    ALIGNMENT,
    FS3_WINDOW,
    S2_WINDOW,
    FS3_CONTIG,
    S2_CONTIG,
]

for path in required:
    if not path.exists():
        raise FileNotFoundError(path)


candidate_df = pd.read_csv(
    HDC2_CANDIDATES,
    sep="\t",
    dtype=str,
)

candidate_acc_col = find_accession_column(
    candidate_df
)

hdc2_accessions = sorted(
    {
        accession_from_text(x)
        for x in candidate_df[candidate_acc_col]
        if accession_from_text(x)
    }
)

assert len(hdc2_accessions) == 36

annotation = pd.read_csv(
    ANNOTATION,
    sep="\t",
    dtype=str,
)

ann_acc_col = find_accession_column(
    annotation
)

ann_org_col = find_organism_column(
    annotation
)

organism_map = {}

for _, row in annotation.iterrows():
    acc = accession_from_text(
        row[ann_acc_col]
    )

    if acc:
        organism_map[acc] = str(
            row[ann_org_col]
        )

missing_org = [
    x
    for x in hdc2_accessions
    if x not in organism_map
]

if missing_org:
    raise RuntimeError(
        "Missing organism mapping: "
        + ", ".join(missing_org)
    )


print()
print("=== 1. HDC2 PROTEIN IDENTITY ===")

records = list(
    SeqIO.parse(
        ALIGNMENT,
        "fasta"
    )
)

seq_by_accession = {}

for record in records:
    acc = accession_from_text(
        record.description
    )

    if acc:
        if acc in seq_by_accession:
            raise RuntimeError(
                f"Duplicate aligned HDC2 sequence: {acc}"
            )

        seq_by_accession[acc] = str(
            record.seq
        )

assert len(seq_by_accession) == 36

identity, identical, comparable = (
    aligned_identity(
        seq_by_accession[FS3],
        seq_by_accession[S2],
    )
)

print(
    f"FS3_1 vs S2MW2 HDC2 identity: "
    f"{identity:.3f}%"
)

print(
    f"Identical/comparable residues: "
    f"{identical}/{comparable}"
)

assert abs(identity - 97.058824) < 0.01


print()
print("=== 2. CROSS-SPECIES HDC2 COMPARISONS ===")

pairwise = []

accessions = sorted(seq_by_accession)

for i in range(len(accessions)):
    a = accessions[i]

    for j in range(i + 1, len(accessions)):
        b = accessions[j]

        if organism_map[a] == organism_map[b]:
            continue

        pid, _, ncomp = aligned_identity(
            seq_by_accession[a],
            seq_by_accession[b],
        )

        pairwise.append(
            {
                "Assembly_1": a,
                "Organism_1": organism_map[a],
                "Assembly_2": b,
                "Organism_2": organism_map[b],
                "Identity_pct": pid,
                "Comparable_sites": ncomp,
            }
        )

pairwise = pd.DataFrame(pairwise)

pairwise.to_csv(
    OUT / "HDC2_cross_species_pairwise_identity.tsv",
    sep="\t",
    index=False,
)

print(
    "Cross-species comparisons:",
    len(pairwise)
)

assert len(pairwise) == 162

target_pair = pairwise.loc[
    (
        (
            pairwise["Assembly_1"] == FS3
        )
        &
        (
            pairwise["Assembly_2"] == S2
        )
    )
    |
    (
        (
            pairwise["Assembly_1"] == S2
        )
        &
        (
            pairwise["Assembly_2"] == FS3
        )
    )
]

assert len(target_pair) == 1

target_identity = float(
    target_pair["Identity_pct"].iloc[0]
)

higher = int(
    (
        pairwise["Identity_pct"]
        > target_identity + 1e-9
    ).sum()
)

print(
    f"FS3_1-S2MW2: "
    f"{target_identity:.3f}%"
)

rank = higher + 1

print(
    "Cross-species comparisons "
    "with higher global identity:",
    higher
)

print(
    f"FS3_1-S2MW2 global-identity rank: "
    f"{rank}/{len(pairwise)}"
)


print()
print("=== 3. WHOLE-GENOME ANI ===")

fs3_genome = genome_fasta(FS3)
s2_genome = genome_fasta(S2)

ani_fs3_s2 = run_fastani(
    fs3_genome,
    s2_genome,
    OUT / "FS3_1_vs_S2MW2_fastANI.txt",
)

ani_s2_fs3 = run_fastani(
    s2_genome,
    fs3_genome,
    OUT / "S2MW2_vs_FS3_1_fastANI.txt",
)

mean_ani = (
    ani_fs3_s2["ANI"]
    + ani_s2_fs3["ANI"]
) / 2

print(
    f"FS3_1 -> S2MW2: "
    f"{ani_fs3_s2['ANI']:.4f}%"
)

print(
    f"S2MW2 -> FS3_1: "
    f"{ani_s2_fs3['ANI']:.4f}%"
)

print(
    f"Mean reciprocal ANI: "
    f"{mean_ani:.4f}%"
)

assert 86.2 < mean_ani < 86.5


print()
print("=== 4. 16-kb HDC2 LOCUS ===")

query_length = fasta_length(
    FS3_WINDOW
)

print(
    "FS3_1 window length:",
    query_length,
    "bp"
)

assert query_length == 16000

locus_out = (
    OUT
    / "FS3_1_vs_S2MW2_16kb_blastn.tsv"
)

run_blastn(
    FS3_WINDOW,
    S2_WINDOW,
    locus_out,
)

locus = read_blast(
    locus_out
)

assert len(locus) >= 1

top = (
    locus
    .sort_values(
        [
            "bitscore",
            "length",
        ],
        ascending=False,
    )
    .iloc[0]
)

covered_locus, locus_cov, _ = (
    blast_summary(
        locus,
        90.0
    )
)

print(
    f"Best HSP identity: "
    f"{top['pident']:.3f}%"
)

print(
    f"Alignment length: "
    f"{int(top['length'])} bp"
)

print(
    "Mismatches:",
    int(top["mismatch"])
)

print(
    "Gap openings:",
    int(top["gapopen"])
)

print(
    f"Query coverage: "
    f"{locus_cov:.2f}%"
)

assert abs(
    float(top["pident"])
    - 99.513
) < 0.01

assert int(top["length"]) == 16001
assert int(top["mismatch"]) == 74
assert int(top["gapopen"]) == 4
assert abs(locus_cov - 100.0) < 0.01


print()
print("=== 5. FULL HDC2-BEARING CONTIG ===")

fs3_contig_length = fasta_length(
    FS3_CONTIG
)

print(
    "FS3_1 contig length:",
    fs3_contig_length,
    "bp"
)

assert fs3_contig_length == 23907

contig_out = (
    OUT
    / "FS3_1_vs_S2MW2_fullcontig_blastn.tsv"
)

run_blastn(
    FS3_CONTIG,
    S2_CONTIG,
    contig_out,
)

contig = read_blast(
    contig_out
)

covered_contig, contig_cov, contig_pid = (
    blast_summary(
        contig,
        99.0
    )
)

print(
    "Non-overlapping >=99% identity coverage:",
    covered_contig,
    "bp"
)

print(
    f"Contig coverage: "
    f"{contig_cov:.2f}%"
)

print(
    f"HSP-length-weighted identity: "
    f"{contig_pid:.3f}%"
)

assert abs(covered_contig - 22897) <= 5
assert abs(contig_cov - 95.78) < 0.05
assert abs(contig_pid - 99.590) < 0.03


print()
print(
    "=== 6. CASSETTE DISTRIBUTION ACROSS "
    "36 HDC2 GENOMES ==="
)

distribution_rows = []
coverage_masks = []

for i, accession in enumerate(
    hdc2_accessions,
    1
):
    print(
        f"[{i:02d}/36] {accession}"
    )

    outfile = (
        BLAST_OUT
        / f"{accession}.tsv"
    )

    run_blastn(
        FS3_WINDOW,
        genome_fasta(accession),
        outfile,
    )

    hits = read_blast(
        outfile
    )

    if hits.empty:
        best_identity = np.nan
        best_length = 0
        best_contig = ""
    else:
        best = (
            hits
            .sort_values(
                [
                    "bitscore",
                    "length",
                ],
                ascending=False,
            )
            .iloc[0]
        )

        best_identity = float(
            best["pident"]
        )

        best_length = int(
            best["length"]
        )

        best_contig = str(
            best["sseqid"]
        )

    covered, coverage, weighted = (
        blast_summary(
            hits,
            90.0
        )
    )

    distribution_rows.append(
        {
            "Assembly": accession,
            "Organism": organism_map[accession],
            "Best_HSP_identity": best_identity,
            "Best_HSP_length": best_length,
            "Weighted_identity_ge90": weighted,
            "Query_coverage_ge90": coverage,
            "Best_contig": best_contig,
        }
    )

    mask = np.zeros(
        query_length,
        dtype=np.uint8,
    )

    selected = hits.loc[
        hits["pident"] >= 90.0
    ]

    for _, hit in selected.iterrows():
        start = (
            min(
                int(hit["qstart"]),
                int(hit["qend"]),
            )
            - 1
        )

        end = max(
            int(hit["qstart"]),
            int(hit["qend"]),
        )

        mask[start:end] = 1

    coverage_masks.append(
        mask
    )

distribution = pd.DataFrame(
    distribution_rows
).sort_values(
    "Query_coverage_ge90",
    ascending=False,
)

distribution_file = (
    OUT
    / "FS3_1_cassette_vs_all_HDC2_genomes_reconstructed.tsv"
)

distribution.to_csv(
    distribution_file,
    sep="\t",
    index=False,
)

print()
print(
    distribution[
        [
            "Assembly",
            "Organism",
            "Best_HSP_identity",
            "Weighted_identity_ge90",
            "Query_coverage_ge90",
        ]
    ]
    .head(10)
    .to_string(index=False)
)

near_full = distribution.loc[
    (
        distribution["Query_coverage_ge90"] >= 95
    )
    &
    (
        distribution["Weighted_identity_ge90"] >= 99
    )
]

print()
print(
    "Near-full / near-identical matches:",
    len(near_full)
)

print(
    near_full[
        [
            "Assembly",
            "Organism",
            "Weighted_identity_ge90",
            "Query_coverage_ge90",
        ]
    ].to_string(index=False)
)

assert set(
    near_full["Assembly"]
) == {
    FS3,
    S2,
}


print()
print("=== 7. CASSETTE POSITION PREVALENCE ===")

matrix = np.vstack(
    coverage_masks
)

position_counts = matrix.sum(
    axis=0
)

prevalence = pd.DataFrame(
    {
        "Query_position":
            np.arange(
                1,
                query_length + 1
            ),
        "HDC2_genomes_covered":
            position_counts,
        "Prevalence_pct":
            (
                100.0
                * position_counts
                / 36
            ),
    }
)

prevalence_file = (
    OUT
    / "FS3_1_16kb_per_position_prevalence.tsv"
)

prevalence.to_csv(
    prevalence_file,
    sep="\t",
    index=False,
)

pair_specific = (
    position_counts <= 2
)

blocks = []
start = None

for pos, flag in enumerate(
    pair_specific,
    start=1
):
    if flag and start is None:
        start = pos

    elif (
        not flag
        and start is not None
    ):
        blocks.append(
            (
                start,
                pos - 1,
            )
        )
        start = None

if start is not None:
    blocks.append(
        (
            start,
            query_length,
        )
    )

blocks.sort(
    key=lambda x:
        x[1] - x[0] + 1,
    reverse=True,
)

largest_start, largest_end = (
    blocks[0]
)

largest_length = (
    largest_end
    - largest_start
    + 1
)

print(
    "Largest <=2-genome block:",
    f"{largest_start}-{largest_end}",
)

print(
    "Length:",
    largest_length,
    "bp"
)

assert abs(largest_start - 7915) <= 5
assert abs(largest_end - 15778) <= 5
assert abs(largest_length - 7864) <= 10


print()
print(
    "=== 8. ARCHIVED DISTRIBUTION VALIDATION ==="
)

if ARCHIVED_DISTRIBUTION.exists():
    archived = pd.read_csv(
        ARCHIVED_DISTRIBUTION,
        sep="\t",
        dtype=str,
    )

    archived_acc_col = (
        find_accession_column(
            archived
        )
    )

    archived_set = {
        accession_from_text(x)
        for x in archived[
            archived_acc_col
        ]
        if accession_from_text(x)
    }

    reconstructed_set = set(
        distribution["Assembly"]
    )

    print(
        "Archived genomes:",
        len(archived_set)
    )

    print(
        "Reconstructed genomes:",
        len(reconstructed_set)
    )

    assert archived_set == reconstructed_set

    print(
        "Assembly-set agreement: YES"
    )


summary_file = (
    OUT
    / "HDC2_horizontal_acquisition_reconstructed_summary.txt"
)

with open(
    summary_file,
    "w"
) as handle:

    handle.write(
        "HDC2 horizontal-acquisition evidence reconstruction\n"
    )
    handle.write("=" * 55 + "\n\n")
    handle.write(f"FS3_1: {FS3}\n")
    handle.write(f"S2MW2: {S2}\n\n")
    handle.write(
        f"HDC2 protein identity: "
        f"{identity:.3f}%\n"
    )
    handle.write(
        f"Cross-species HDC2 comparisons: "
        f"{len(pairwise)}\n"
    )
    handle.write(
        f"Cross-species comparisons with higher global identity: "
        f"{higher}\n"
    )

    handle.write(
        f"FS3_1-S2MW2 global-identity rank: "
        f"{rank}/{len(pairwise)}\n\n"
    )
    handle.write(
        f"FS3_1 -> S2MW2 whole-genome ANI: "
        f"{ani_fs3_s2['ANI']:.4f}%\n"
    )
    handle.write(
        f"S2MW2 -> FS3_1 whole-genome ANI: "
        f"{ani_s2_fs3['ANI']:.4f}%\n"
    )
    handle.write(
        f"Mean reciprocal ANI: "
        f"{mean_ani:.4f}%\n\n"
    )
    handle.write(
        f"16-kb locus identity: "
        f"{float(top['pident']):.3f}%\n"
    )
    handle.write(
        f"16-kb query coverage: "
        f"{locus_cov:.2f}%\n"
    )
    handle.write(
        f"16-kb alignment length: "
        f"{int(top['length'])} bp\n"
    )
    handle.write(
        f"Mismatches: "
        f"{int(top['mismatch'])}\n"
    )
    handle.write(
        f"Gap openings: "
        f"{int(top['gapopen'])}\n\n"
    )
    handle.write(
        f"FS3_1 HDC2-bearing contig length: "
        f"{fs3_contig_length} bp\n"
    )
    handle.write(
        f">=99%-identity non-overlapping coverage: "
        f"{covered_contig} bp "
        f"({contig_cov:.2f}%)\n"
    )
    handle.write(
        f"HSP-length-weighted identity: "
        f"{contig_pid:.3f}%\n\n"
    )
    handle.write(
        f"Near-full near-identical cassette matches: "
        f"{len(near_full)}/36\n"
    )
    handle.write(
        f"Largest pair-specific/low-prevalence block: "
        f"{largest_start}-{largest_end} "
        f"({largest_length} bp)\n\n"
    )
    handle.write(
        "Interpretation: the localized near-identity of the HDC2-associated "
        "region despite substantially lower whole-genome ANI is strongly "
        "consistent with recent horizontal acquisition/shared recent "
        "acquisition. Transfer direction is unresolved, and an unsampled "
        "third donor cannot be excluded.\n"
    )

print()
print("Saved:")
print(" ", distribution_file.relative_to(ROOT))
print(" ", prevalence_file.relative_to(ROOT))
print(" ", summary_file.relative_to(ROOT))

print()
print("STATUS: PASS")
