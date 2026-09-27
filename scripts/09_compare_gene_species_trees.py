#!/usr/bin/env python3

"""
Compare strongly supported HDC/HDC2 gene-tree splits with corresponding
pruned core-genome species trees.

Primary evolutionary results are based on the archived IQ-TREE analyses
used during the original analysis.

Seeded reconstructed IQ-TREE analyses are treated as reproducibility /
robustness checks because bootstrap support values may vary slightly
between independent runs.

Strong support:
    SH-aLRT >= 80
    UFBoot >= 95
"""

from pathlib import Path
import re
import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]

# ------------------------------------------------------------
# Primary archived trees
# ------------------------------------------------------------

PRIMARY_HDC = (
    ROOT
    / "results/gene_phylogeny/HDC_IQTREE/"
      "Photobacterium_conventional_HDC_primary98.treefile"
)

PRIMARY_HDC2 = (
    ROOT
    / "results/gene_phylogeny/HDC2_IQTREE/"
      "Photobacterium_HDC2_primary36.treefile"
)

# ------------------------------------------------------------
# Seeded reproducibility trees
# ------------------------------------------------------------

RECON_HDC = (
    ROOT
    / "results/reproducibility/IQTREE/HDC/"
      "Photobacterium_conventional_HDC_primary98_reconstructed.treefile"
)

RECON_HDC2 = (
    ROOT
    / "results/reproducibility/IQTREE/HDC2/"
      "Photobacterium_HDC2_primary36_reconstructed.treefile"
)

# ------------------------------------------------------------
# Species trees
# ------------------------------------------------------------

SPECIES_HDC = (
    ROOT
    / "results/gene_phylogeny/species_tree_pruned/"
      "Photobacterium_species_tree_HDC98.nwk"
)

SPECIES_HDC2 = (
    ROOT
    / "results/gene_phylogeny/species_tree_pruned/"
      "Photobacterium_species_tree_HDC2_36.nwk"
)

OUTDIR = ROOT / "results/reproducibility"

DETAIL_OUT = (
    OUTDIR
    / "gene_species_tree_split_comparison.tsv"
)

SUMMARY_OUT = (
    OUTDIR
    / "gene_species_tree_comparison_summary.tsv"
)

SH_MIN = 80
UF_MIN = 95


def accession(label):
    m = re.search(
        r"(GC[AF]_\d+\.\d+)",
        str(label)
    )

    if not m:
        raise RuntimeError(
            f"Cannot recover accession from {label}"
        )

    return m.group(1)


def support_pair(clade):

    values = []

    if clade.name is not None:
        values.append(str(clade.name))

    if clade.confidence is not None:
        values.append(str(clade.confidence))

    for text in values:

        m = re.search(
            r"([0-9]+(?:\.[0-9]+)?)/"
            r"([0-9]+(?:\.[0-9]+)?)",
            text
        )

        if m:
            return (
                float(m.group(1)),
                float(m.group(2))
            )

    return None


def canonical_split(side, universe):

    side = frozenset(side)
    other = frozenset(universe - side)

    if len(side) < 2 or len(other) < 2:
        return None

    if len(side) < len(other):
        return side

    if len(other) < len(side):
        return other

    return min(
        side,
        other,
        key=lambda x: tuple(sorted(x))
    )


def taxon_set(tree):

    return {
        accession(x.name)
        for x in tree.get_terminals()
    }


def species_splits(tree, universe):

    out = set()

    for clade in tree.get_nonterminals():

        side = {
            accession(x.name)
            for x in clade.get_terminals()
        }

        s = canonical_split(
            side,
            universe
        )

        if s is not None:
            out.add(s)

    return out


def supported_gene_splits(tree, universe):

    out = {}

    for clade in tree.get_nonterminals():

        pair = support_pair(clade)

        if pair is None:
            continue

        sh, uf = pair

        side = {
            accession(x.name)
            for x in clade.get_terminals()
        }

        split = canonical_split(
            side,
            universe
        )

        if split is None:
            continue

        current = out.get(split)

        row = {
            "SH_aLRT": sh,
            "UFBoot": uf,
            "strong":
                sh >= SH_MIN
                and uf >= UF_MIN,
        }

        if (
            current is None
            or (sh, uf) >
            (
                current["SH_aLRT"],
                current["UFBoot"]
            )
        ):
            out[split] = row

    return out


def compare_to_species(
    system,
    analysis,
    gene_file,
    species_file,
):

    gene = Phylo.read(
        gene_file,
        "newick"
    )

    species = Phylo.read(
        species_file,
        "newick"
    )

    gene_taxa = taxon_set(gene)
    species_taxa = taxon_set(species)

    if gene_taxa != species_taxa:
        raise RuntimeError(
            f"{system} {analysis}: "
            "gene/species taxon sets differ"
        )

    species_bp = species_splits(
        species,
        gene_taxa
    )

    all_gene = supported_gene_splits(
        gene,
        gene_taxa
    )

    strong = {
        split: values
        for split, values in all_gene.items()
        if values["strong"]
    }

    rows = []

    for i, (split, values) in enumerate(
        sorted(
            strong.items(),
            key=lambda x: (
                len(x[0]),
                tuple(sorted(x[0]))
            )
        ),
        1
    ):

        concordant = (
            split in species_bp
        )

        rows.append(
            {
                "System": system,
                "Analysis": analysis,
                "Strong_Split_ID": i,
                "SH_aLRT":
                    values["SH_aLRT"],
                "UFBoot":
                    values["UFBoot"],
                "Species_tree_status":
                    (
                        "CONCORDANT"
                        if concordant
                        else "DISCORDANT"
                    ),
                "Smaller_split_size":
                    len(split),
                "Smaller_split_accessions":
                    ",".join(sorted(split)),
            }
        )

    concordant = sum(
        r["Species_tree_status"]
        == "CONCORDANT"
        for r in rows
    )

    discordant = sum(
        r["Species_tree_status"]
        == "DISCORDANT"
        for r in rows
    )

    return (
        rows,
        len(strong),
        concordant,
        discordant,
        all_gene,
    )


def strong_set(data):

    return {
        split
        for split, x in data.items()
        if x["strong"]
    }


print("=" * 78)
print("GENE-TREE / SPECIES-TREE COMPARISON")
print("=" * 78)

detail_rows = []
summary_rows = []

systems = [
    (
        "Conventional_HDC",
        PRIMARY_HDC,
        RECON_HDC,
        SPECIES_HDC,
    ),
    (
        "HDC2",
        PRIMARY_HDC2,
        RECON_HDC2,
        SPECIES_HDC2,
    ),
]

for (
    system,
    primary_file,
    recon_file,
    species_file,
) in systems:

    print()
    print("=" * 78)
    print(system)
    print("=" * 78)

    (
        primary_rows,
        primary_n,
        primary_con,
        primary_dis,
        primary_all,
    ) = compare_to_species(
        system,
        "PRIMARY_ARCHIVED",
        primary_file,
        species_file,
    )

    (
        recon_rows,
        recon_n,
        recon_con,
        recon_dis,
        recon_all,
    ) = compare_to_species(
        system,
        "SEEDED_RECONSTRUCTION",
        recon_file,
        species_file,
    )

    detail_rows.extend(
        primary_rows
    )

    detail_rows.extend(
        recon_rows
    )

    pstrong = strong_set(
        primary_all
    )

    rstrong = strong_set(
        recon_all
    )

    shared = (
        pstrong & rstrong
    )

    primary_only = (
        pstrong - rstrong
    )

    recon_only = (
        rstrong - pstrong
    )

    print()
    print("PRIMARY ARCHIVED")
    print(
        "  Strong splits:",
        primary_n
    )
    print(
        "  Concordant:",
        primary_con
    )
    print(
        "  Discordant:",
        primary_dis
    )

    print()
    print("SEEDED RECONSTRUCTION")
    print(
        "  Strong splits:",
        recon_n
    )
    print(
        "  Concordant:",
        recon_con
    )
    print(
        "  Discordant:",
        recon_dis
    )

    print()
    print("REPRODUCIBILITY")
    print(
        "  Shared strong splits:",
        len(shared)
    )
    print(
        "  Primary-only strong:",
        len(primary_only)
    )
    print(
        "  Reconstruction-only strong:",
        len(recon_only)
    )

    summary_rows.append(
        {
            "System": system,
            "Primary_strong_splits":
                primary_n,
            "Primary_concordant":
                primary_con,
            "Primary_discordant":
                primary_dis,
            "Reconstructed_strong_splits":
                recon_n,
            "Reconstructed_concordant":
                recon_con,
            "Reconstructed_discordant":
                recon_dis,
            "Shared_strong_splits":
                len(shared),
            "Primary_only_strong":
                len(primary_only),
            "Reconstruction_only_strong":
                len(recon_only),
        }
    )

# ------------------------------------------------------------
# Validate locked primary results
# ------------------------------------------------------------

s = {
    x["System"]: x
    for x in summary_rows
}

assert (
    s["Conventional_HDC"]
    ["Primary_strong_splits"]
    == 4
)

assert (
    s["Conventional_HDC"]
    ["Primary_concordant"]
    == 3
)

assert (
    s["Conventional_HDC"]
    ["Primary_discordant"]
    == 1
)

assert (
    s["HDC2"]
    ["Primary_strong_splits"]
    == 4
)

assert (
    s["HDC2"]
    ["Primary_concordant"]
    == 1
)

assert (
    s["HDC2"]
    ["Primary_discordant"]
    == 3
)

# HDC reproducibility result just established:
assert (
    s["Conventional_HDC"]
    ["Reconstructed_strong_splits"]
    == 3
)

assert (
    s["Conventional_HDC"]
    ["Shared_strong_splits"]
    == 3
)

assert (
    s["Conventional_HDC"]
    ["Primary_only_strong"]
    == 1
)

assert (
    s["Conventional_HDC"]
    ["Reconstruction_only_strong"]
    == 0
)

pd.DataFrame(
    detail_rows
).to_csv(
    DETAIL_OUT,
    sep="\t",
    index=False
)

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)

print()
print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print(
    summary.to_string(
        index=False
    )
)

print()
print("Saved:")
print(
    " ",
    DETAIL_OUT.relative_to(ROOT)
)
print(
    " ",
    SUMMARY_OUT.relative_to(ROOT)
)

print()
print(
    f"Strong-support criterion: "
    f"SH-aLRT >= {SH_MIN}, "
    f"UFBoot >= {UF_MIN}"
)

print()
print("STATUS: PASS")
