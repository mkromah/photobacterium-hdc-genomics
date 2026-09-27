#!/usr/bin/env python3

"""
Reconstruct maximum-likelihood phylogenies for conventional HDC and HDC2.

Method
------
IQ-TREE ModelFinder (-m MFP)
1000 ultrafast bootstrap replicates (-B 1000)
1000 SH-aLRT replicates (--alrt 1000)

An explicit seed is supplied for reproducible reruns.
"""

from pathlib import Path
import re
import shutil
import subprocess
from Bio import SeqIO, Phylo

ROOT = Path(__file__).resolve().parents[1]

INDIR = ROOT / "results/reproducibility"

HDC_ALN = (
    INDIR
    / "Photobacterium_conventional_HDC_final497_reconstructed_aligned.faa"
)

HDC2_ALN = (
    INDIR
    / "Photobacterium_HDC2_final497_reconstructed_aligned.faa"
)

OUTROOT = INDIR / "IQTREE"
HDC_DIR = OUTROOT / "HDC"
HDC2_DIR = OUTROOT / "HDC2"

HDC_DIR.mkdir(parents=True, exist_ok=True)
HDC2_DIR.mkdir(parents=True, exist_ok=True)

HDC_PREFIX = (
    HDC_DIR
    / "Photobacterium_conventional_HDC_primary98_reconstructed"
)

HDC2_PREFIX = (
    HDC2_DIR
    / "Photobacterium_HDC2_primary36_reconstructed"
)

SEED = 20260927


def run_iqtree(alignment, prefix):

    cmd = [
        "iqtree",
        "-s", str(alignment),
        "-m", "MFP",
        "-B", "1000",
        "--alrt", "1000",
        "-T", "AUTO",
        "-seed", str(SEED),
        "-pre", str(prefix),
        "-redo",
    ]

    print()
    print("Running:")
    print(" ".join(cmd))
    print()

    subprocess.run(cmd, check=True)


def get_model(iqtree_file):

    text = iqtree_file.read_text()

    patterns = [
        r"Best-fit model according to BIC:\s*(\S+)",
        r"Model of substitution:\s*(\S+)",
    ]

    for pattern in patterns:
        m = re.search(pattern, text)

        if m:
            return m.group(1)

    raise RuntimeError(
        f"Could not recover model from {iqtree_file}"
    )


def validate(name, alignment, prefix, expected_n, expected_model):

    treefile = Path(str(prefix) + ".treefile")
    iqtree_file = Path(str(prefix) + ".iqtree")

    if not treefile.exists():
        raise FileNotFoundError(treefile)

    if not iqtree_file.exists():
        raise FileNotFoundError(iqtree_file)

    alignment_ids = {
        r.id
        for r in SeqIO.parse(alignment, "fasta")
    }

    tree = Phylo.read(treefile, "newick")

    tree_ids = {
        terminal.name
        for terminal in tree.get_terminals()
    }

    model = get_model(iqtree_file)

    print()
    print("=" * 72)
    print(name)
    print("=" * 72)

    print("Alignment taxa:", len(alignment_ids))
    print("Tree taxa:     ", len(tree_ids))
    print("Model:         ", model)
    print("Expected model:", expected_model)

    missing = sorted(alignment_ids - tree_ids)
    extra = sorted(tree_ids - alignment_ids)

    print("Missing tree taxa:", len(missing))
    print("Extra tree taxa:  ", len(extra))

    if missing:
        print("Missing:", missing)

    if extra:
        print("Extra:", extra)

    assert len(alignment_ids) == expected_n
    assert len(tree_ids) == expected_n
    assert not missing
    assert not extra
    assert model == expected_model

    print("Validation: PASS")


print("=" * 76)
print("RECONSTRUCTING HDC / HDC2 MAXIMUM-LIKELIHOOD PHYLOGENIES")
print("=" * 76)

iqtree = shutil.which("iqtree")

if iqtree is None:
    raise SystemExit(
        "ERROR: iqtree not found in current environment"
    )

print()
print("IQ-TREE executable:", iqtree)

version = subprocess.run(
    ["iqtree", "--version"],
    capture_output=True,
    text=True,
)

version_text = (
    version.stdout.strip()
    or version.stderr.strip()
)

print("IQ-TREE version:")
print(version_text.splitlines()[0] if version_text else "unknown")

for f in [HDC_ALN, HDC2_ALN]:
    if not f.exists():
        raise FileNotFoundError(f)

run_iqtree(
    HDC_ALN,
    HDC_PREFIX
)

run_iqtree(
    HDC2_ALN,
    HDC2_PREFIX
)

validate(
    "CONVENTIONAL HDC TREE",
    HDC_ALN,
    HDC_PREFIX,
    expected_n=98,
    expected_model="WAG+G4",
)

validate(
    "HDC2 TREE",
    HDC2_ALN,
    HDC2_PREFIX,
    expected_n=36,
    expected_model="JTT+I",
)

print()
print("Saved tree files:")
print(
    " ",
    Path(str(HDC_PREFIX) + ".treefile").relative_to(ROOT)
)
print(
    " ",
    Path(str(HDC2_PREFIX) + ".treefile").relative_to(ROOT)
)

print()
print("Random seed:", SEED)
print("STATUS: PASS")
