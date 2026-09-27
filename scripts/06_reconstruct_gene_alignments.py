#!/usr/bin/env python3

"""
Reconstruct MAFFT alignments for conventional HDC and HDC2 and perform
alignment-quality checks.

Expected primary alignments
---------------------------
Conventional HDC:
    sequences = 98
    alignment length = 380 aa
    columns containing gaps = 6
    variable columns = 84
    invariant columns = 296

HDC2:
    sequences = 36
    alignment length = 771 aa
    columns containing gaps = 112
    variable columns = 134
    invariant columns = 637

The shorter P. phosphoreum FS3_1 HDC2 sequence is expected to contain
approximately 11.80% gaps in the final HDC2 alignment.
"""

from pathlib import Path
import shutil
import subprocess
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]

INDIR = ROOT / "results/reproducibility"

HDC_IN = (
    INDIR
    / "Photobacterium_conventional_HDC_final497_reconstructed.faa"
)

HDC2_IN = (
    INDIR
    / "Photobacterium_HDC2_final497_reconstructed.faa"
)

HDC_ALN = (
    INDIR
    / "Photobacterium_conventional_HDC_final497_reconstructed_aligned.faa"
)

HDC2_ALN = (
    INDIR
    / "Photobacterium_HDC2_final497_reconstructed_aligned.faa"
)

HDC_LOG = (
    INDIR
    / "Photobacterium_conventional_HDC_MAFFT_reconstructed.log"
)

HDC2_LOG = (
    INDIR
    / "Photobacterium_HDC2_MAFFT_reconstructed.log"
)


def run_mafft(infile, outfile, logfile):

    cmd = [
        "mafft",
        "--auto",
        "--thread",
        "-1",
        str(infile),
    ]

    print()
    print("Running:")
    print(" ".join(cmd))

    with open(outfile, "w") as out, open(logfile, "w") as err:
        subprocess.run(
            cmd,
            stdout=out,
            stderr=err,
            check=True,
        )


def alignment_qc(path, name):

    records = list(SeqIO.parse(path, "fasta"))

    if not records:
        raise RuntimeError(
            f"{name}: alignment contains no sequences"
        )

    lengths = {len(r.seq) for r in records}

    if len(lengths) != 1:
        raise RuntimeError(
            f"{name}: alignment sequences have unequal lengths"
        )

    aln_len = next(iter(lengths))
    nseq = len(records)

    seqs = [str(r.seq).upper() for r in records]

    total_cells = nseq * aln_len
    total_gaps = sum(seq.count("-") for seq in seqs)

    gap_fraction = total_gaps / total_cells

    gap_columns = 0
    variable_columns = 0

    for i in range(aln_len):

        column = [seq[i] for seq in seqs]

        if "-" in column:
            gap_columns += 1

        residues = {
            x for x in column
            if x not in {"-", "?", "X"}
        }

        if len(residues) > 1:
            variable_columns += 1

    invariant_columns = aln_len - variable_columns

    sequence_gap_stats = []

    for record in records:

        seq = str(record.seq)
        gaps = seq.count("-")
        pct = 100 * gaps / aln_len

        sequence_gap_stats.append(
            (record.id, gaps, pct)
        )

    high_gap = [
        x for x in sequence_gap_stats
        if x[2] > 5
    ]

    print()
    print("=" * 72)
    print(name)
    print("=" * 72)
    print("Sequences:             ", nseq)
    print("Alignment length:      ", aln_len, "aa")
    print(
        "Overall gap fraction: ",
        f"{100 * gap_fraction:.2f}%"
    )
    print(
        "Gap-containing columns:",
        f"{gap_columns}/{aln_len}",
        f"({100 * gap_columns / aln_len:.1f}%)"
    )
    print(
        "Variable columns:      ",
        f"{variable_columns}/{aln_len}",
        f"({100 * variable_columns / aln_len:.1f}%)"
    )
    print(
        "Invariant columns:     ",
        f"{invariant_columns}/{aln_len}",
        f"({100 * invariant_columns / aln_len:.1f}%)"
    )

    print()
    print("Sequences with >5% gaps:", len(high_gap))

    for record_id, gaps, pct in high_gap:
        print(
            f"  {record_id}: "
            f"{gaps} gaps ({pct:.2f}%)"
        )

    return {
        "nseq": nseq,
        "aln_len": aln_len,
        "gap_fraction": gap_fraction,
        "gap_columns": gap_columns,
        "variable_columns": variable_columns,
        "invariant_columns": invariant_columns,
        "gap_stats": sequence_gap_stats,
    }


print("=" * 76)
print("RECONSTRUCTING HDC / HDC2 MAFFT ALIGNMENTS")
print("=" * 76)

mafft = shutil.which("mafft")

if mafft is None:
    raise SystemExit(
        "ERROR: mafft not found in current environment"
    )

print()
print("MAFFT:", mafft)

version = subprocess.run(
    ["mafft", "--version"],
    capture_output=True,
    text=True,
)

print(
    "Version:",
    (version.stderr or version.stdout).strip()
)

for f in [HDC_IN, HDC2_IN]:
    if not f.exists():
        raise FileNotFoundError(f)

run_mafft(
    HDC_IN,
    HDC_ALN,
    HDC_LOG
)

run_mafft(
    HDC2_IN,
    HDC2_ALN,
    HDC2_LOG
)

hdc = alignment_qc(
    HDC_ALN,
    "CONVENTIONAL HDC ALIGNMENT"
)

hdc2 = alignment_qc(
    HDC2_ALN,
    "HDC2 ALIGNMENT"
)

# ------------------------------------------------------------
# Locked QC validation
# ------------------------------------------------------------

assert hdc["nseq"] == 98
assert hdc["aln_len"] == 380
assert hdc["gap_columns"] == 6
assert hdc["variable_columns"] == 84
assert hdc["invariant_columns"] == 296

assert hdc2["nseq"] == 36
assert hdc2["aln_len"] == 771
assert hdc2["gap_columns"] == 112
assert hdc2["variable_columns"] == 134
assert hdc2["invariant_columns"] == 637

# FS3_1 sensitivity target
fs3 = [
    x for x in hdc2["gap_stats"]
    if x[0].startswith("GCF_001676275.1|")
]

if len(fs3) != 1:
    raise RuntimeError(
        f"Expected one FS3_1 sequence, found {len(fs3)}"
    )

fs3_id, fs3_gaps, fs3_gap_pct = fs3[0]

print()
print("=== FS3_1 ALIGNMENT CHECK ===")
print("ID:", fs3_id)
print("Gap residues:", fs3_gaps)
print(f"Gap fraction: {fs3_gap_pct:.2f}%")

assert fs3_gaps == 91
assert abs(fs3_gap_pct - 11.80) < 0.02

print()
print("Saved:")
print(" ", HDC_ALN.relative_to(ROOT))
print(" ", HDC2_ALN.relative_to(ROOT))
print(" ", HDC_LOG.relative_to(ROOT))
print(" ", HDC2_LOG.relative_to(ROOT))

print()
print("STATUS: PASS")
