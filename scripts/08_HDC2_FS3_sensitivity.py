#!/usr/bin/env python3

"""
Sensitivity analysis for the shorter P. phosphoreum FS3_1 HDC2 sequence.

Workflow
--------
1. Remove GCF_001676275.1 (FS3_1) from the reconstructed 36-sequence HDC2 FASTA.
2. Realign the remaining 35 sequences independently with MAFFT.
3. Reconstruct an ML tree with IQ-TREE using:
      - ModelFinder
      - 1000 UFBoot
      - 1000 SH-aLRT
4. Compare strongly supported splits between:
      - the primary 36-sequence HDC2 tree, restricted to the 35 common taxa
      - the 35-sequence sensitivity tree

Strong support:
    SH-aLRT >= 80
    UFBoot >= 95
"""

from pathlib import Path
import re
import shutil
import subprocess

from Bio import SeqIO, Phylo

ROOT = Path(__file__).resolve().parents[1]

INDIR = ROOT / "results/reproducibility"

PRIMARY_FASTA = (
    INDIR
    / "Photobacterium_HDC2_final497_reconstructed.faa"
)

PRIMARY_TREE = (
    INDIR
    / "IQTREE/HDC2/"
      "Photobacterium_HDC2_primary36_reconstructed.treefile"
)

OUTDIR = INDIR / "HDC2_sensitivity_no_FS3_1"
OUTDIR.mkdir(parents=True, exist_ok=True)

SENS_FASTA = OUTDIR / "Photobacterium_HDC2_no_FS3_1.faa"
SENS_ALN = OUTDIR / "Photobacterium_HDC2_no_FS3_1_aligned.faa"

PREFIX = OUTDIR / "Photobacterium_HDC2_no_FS3_1"

FS3_ACCESSION = "GCF_001676275.1"
SEED = 20260927


def canonical_split(side, universe):
    side = frozenset(side)
    other = frozenset(universe - side)

    if len(side) < 2 or len(other) < 2:
        return None

    # Canonical representation of the bipartition.
    if len(side) < len(other):
        return side

    if len(other) < len(side):
        return other

    return min(
        side,
        other,
        key=lambda x: tuple(sorted(x))
    )


def get_support(clade):
    """
    IQ-TREE commonly stores support as SH-aLRT/UFBoot.
    Biopython may place this in clade.name.
    """

    candidates = []

    if clade.name:
        candidates.append(str(clade.name))

    if clade.confidence is not None:
        candidates.append(str(clade.confidence))

    for text in candidates:

        m = re.search(
            r"([0-9]+(?:\.[0-9]+)?)/"
            r"([0-9]+(?:\.[0-9]+)?)",
            text,
        )

        if m:
            return float(m.group(1)), float(m.group(2))

    return None


def strong_splits(treefile, universe, remove=None):

    tree = Phylo.read(treefile, "newick")

    results = {}

    for clade in tree.get_nonterminals():

        support = get_support(clade)

        if support is None:
            continue

        sh_alrt, ufboot = support

        if sh_alrt < 80 or ufboot < 95:
            continue

        leaves = {
            x.name
            for x in clade.get_terminals()
        }

        if remove is not None:
            leaves.discard(remove)

        split = canonical_split(
            leaves,
            universe
        )

        if split is None:
            continue

        results[split] = {
            "SH_aLRT": sh_alrt,
            "UFBoot": ufboot,
        }

    return results


print("=" * 76)
print("HDC2 FS3_1 SENSITIVITY ANALYSIS")
print("=" * 76)

# ------------------------------------------------------------
# Remove FS3_1
# ------------------------------------------------------------

records = list(
    SeqIO.parse(PRIMARY_FASTA, "fasta")
)

assert len(records) == 36

removed = [
    r for r in records
    if r.id.startswith(FS3_ACCESSION + "|")
]

retained = [
    r for r in records
    if not r.id.startswith(FS3_ACCESSION + "|")
]

print()
print("Primary sequences:", len(records))
print("FS3_1 removed:   ", len(removed))
print("Retained:        ", len(retained))

assert len(removed) == 1
assert len(retained) == 35

print(
    "Removed ID:",
    removed[0].id
)
print(
    "Removed length:",
    len(removed[0].seq),
    "aa"
)

assert len(removed[0].seq) == 680

SeqIO.write(
    retained,
    SENS_FASTA,
    "fasta"
)

# ------------------------------------------------------------
# MAFFT realignment
# ------------------------------------------------------------

if shutil.which("mafft") is None:
    raise SystemExit("ERROR: mafft not found")

print()
print("Running independent MAFFT alignment...")

with open(SENS_ALN, "w") as out:
    subprocess.run(
        [
            "mafft",
            "--auto",
            "--thread",
            "-1",
            str(SENS_FASTA),
        ],
        stdout=out,
        stderr=subprocess.DEVNULL,
        check=True,
    )

aligned = list(
    SeqIO.parse(SENS_ALN, "fasta")
)

assert len(aligned) == 35
assert len({len(x.seq) for x in aligned}) == 1

print(
    "Sensitivity alignment:",
    len(aligned),
    "sequences x",
    len(aligned[0].seq),
    "aa"
)

# ------------------------------------------------------------
# IQ-TREE sensitivity reconstruction
# ------------------------------------------------------------

if shutil.which("iqtree") is None:
    raise SystemExit("ERROR: iqtree not found")

print()
print("Running IQ-TREE sensitivity phylogeny...")

subprocess.run(
    [
        "iqtree",
        "-s", str(SENS_ALN),
        "-m", "MFP",
        "-B", "1000",
        "--alrt", "1000",
        "-T", "AUTO",
        "-seed", str(SEED),
        "-pre", str(PREFIX),
        "-redo",
    ],
    check=True,
)

SENS_TREE = Path(str(PREFIX) + ".treefile")
SENS_REPORT = Path(str(PREFIX) + ".iqtree")

# ------------------------------------------------------------
# Taxon validation
# ------------------------------------------------------------

common_taxa = {
    r.id for r in retained
}

sens_tree = Phylo.read(
    SENS_TREE,
    "newick"
)

sens_taxa = {
    x.name
    for x in sens_tree.get_terminals()
}

print()
print("Sensitivity FASTA taxa:", len(common_taxa))
print("Sensitivity tree taxa: ", len(sens_taxa))
print(
    "Missing taxa:",
    len(common_taxa - sens_taxa)
)
print(
    "Extra taxa:",
    len(sens_taxa - common_taxa)
)

assert common_taxa == sens_taxa

# ------------------------------------------------------------
# Strong split comparison
# ------------------------------------------------------------

primary_strong = strong_splits(
    PRIMARY_TREE,
    common_taxa,
    remove=removed[0].id,
)

sensitivity_strong = strong_splits(
    SENS_TREE,
    common_taxa,
)

primary_set = set(primary_strong)
sensitivity_set = set(sensitivity_strong)

shared = primary_set & sensitivity_set
primary_only = primary_set - sensitivity_set
sensitivity_only = sensitivity_set - primary_set

print()
print("=== STRONGLY SUPPORTED SPLITS ===")
print(
    "Primary tree on 35 common taxa:",
    len(primary_set)
)
print(
    "Sensitivity tree:             ",
    len(sensitivity_set)
)
print(
    "Shared:                       ",
    len(shared)
)
print(
    "Primary-only:                 ",
    len(primary_only)
)
print(
    "Sensitivity-only:             ",
    len(sensitivity_only)
)

# Previously established sensitivity result.
assert len(sensitivity_set) == 3
assert len(shared) == 3
assert len(primary_only) == 0
assert len(sensitivity_only) == 0

print()
print(
    "All strongly supported sensitivity-tree "
    "splits are retained from the primary analysis: YES"
)

# ------------------------------------------------------------
# Recover selected sensitivity model
# ------------------------------------------------------------

report = SENS_REPORT.read_text()

m = re.search(
    r"Best-fit model according to BIC:\s*(\S+)",
    report
)

model = m.group(1) if m else "UNRESOLVED"

print()
print("Sensitivity-tree model:", model)

print()
print("Saved:")
print(" ", SENS_FASTA.relative_to(ROOT))
print(" ", SENS_ALN.relative_to(ROOT))
print(" ", SENS_TREE.relative_to(ROOT))
print(" ", SENS_REPORT.relative_to(ROOT))

print()
print("Random seed:", SEED)
print("STATUS: PASS")
