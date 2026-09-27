#!/usr/bin/env python3

"""
Reconstruct and validate genomic context of conventional HDC and HDC2 loci.

The script:
  1. Parses NCBI GFF annotations directly.
  2. Locates the target HDC/HDC2 CDS in each candidate genome.
  3. Calculates edge-to-edge distance to the nearest hisS-associated CDS.
  4. Calculates distance to the nearest transporter/permease/antiporter CDS.
  5. Reproduces final cohort/context summary statistics.
  6. Does not assume an unrecovered historical fixed proximity cutoff.

Expected final results
----------------------
Conventional HDC:
    final genomes = 98
    hisS found = 98/98
    transporter/antiporter found = 98/98

HDC2:
    genomes = 36
    hisS found = 36/36
    transporter evaluable/present = 32/36
    COMPLETE_CLUSTER = 32
    BOUNDARY_SUPPORTED = 4
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

RAW_ROOT = (
    ROOT
    / "data/raw/photobacterium_primary_500_ncbi/ncbi_dataset/data"
)

HDC_TABLE = (
    ROOT
    / "results/conventional_HDC_candidate_neighborhoods_2026-09-13.tsv"
)

HDC2_TABLE = (
    ROOT
    / "results/HDC2_candidate_neighborhoods_2026-09-13.tsv"
)

HDC2_CONTEXT = (
    ROOT
    / "results/HDC2_final_context_classification_2026-09-13.tsv"
)

EXCLUSIONS = (
    ROOT
    / "results/ubcg_primary_exclusions_2026-09-13.tsv"
)

OUTDIR = ROOT / "results/reproducibility"
OUTDIR.mkdir(parents=True, exist_ok=True)

HDC_OUT = OUTDIR / "Conventional_HDC_context_reconstructed.tsv"
HDC2_OUT = OUTDIR / "HDC2_context_reconstructed.tsv"
SUMMARY_OUT = OUTDIR / "HDC_HDC2_context_reconstructed_summary.tsv"


def parse_attrs(text):
    out = {}

    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value

    return out


def interval_distance(a1, a2, b1, b2):
    """
    Edge-to-edge distance between two genomic intervals.
    Returns zero if they overlap.
    """
    if a2 < b1:
        return b1 - a2 - 1

    if b2 < a1:
        return a1 - b2 - 1

    return 0


def load_gff_cds(accession):
    folder = RAW_ROOT / accession
    gffs = sorted(folder.glob("*.gff"))

    if not gffs:
        raise FileNotFoundError(
            f"No GFF found for {accession}: {folder}"
        )

    cds = []

    with open(gffs[0]) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9 or parts[2] != "CDS":
                continue

            attrs = parse_attrs(parts[8])

            cds.append(
                {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "protein": attrs.get("protein_id", ""),
                    "gene": attrs.get("gene", ""),
                    "product": attrs.get("product", ""),
                }
            )

    return cds


def is_hiss(feature):
    gene = feature["gene"].lower()
    text = (
        feature["gene"] + " " + feature["product"]
    ).lower()

    return (
        gene == "hiss"
        or "histidine--trna ligase" in text
        or "histidine-trna ligase" in text
    )


def is_transporter(feature):
    text = (
        feature["gene"] + " " + feature["product"]
    ).lower()

    return (
        "permease" in text
        or "transporter" in text
        or "antiporter" in text
    )


def reconstruct(table, system):
    records = []

    for _, row in table.iterrows():
        accession = str(row["Assembly Accession"])
        target_protein = str(row["Protein"])

        cds = load_gff_cds(accession)

        targets = [
            x for x in cds
            if x["protein"] == target_protein
        ]

        if not targets:
            raise RuntimeError(
                f"{accession}: no target CDS found for "
                f"{target_protein}"
            )

        # Some assemblies contain duplicated copies of the same
        # HDC protein accession. Select the copy with the strongest
        # local canonical context: nearest hisS + transporter.
        # Ties are resolved deterministically by contig and start.
        def target_context_score(t):
            same = [
                x for x in cds
                if x["contig"] == t["contig"]
                and not (
                    x["protein"] == target_protein
                    and x["start"] == t["start"]
                    and x["end"] == t["end"]
                )
            ]

            hiss_distances = []
            transporter_distances = []

            for feature in same:
                d = interval_distance(
                    t["start"],
                    t["end"],
                    feature["start"],
                    feature["end"],
                )

                if is_hiss(feature):
                    hiss_distances.append(d)

                if is_transporter(feature):
                    transporter_distances.append(d)

            hiss_min = (
                min(hiss_distances)
                if hiss_distances
                else 10**12
            )

            transporter_min = (
                min(transporter_distances)
                if transporter_distances
                else 10**12
            )

            missing_components = (
                int(not hiss_distances)
                + int(not transporter_distances)
            )

            return (
                missing_components,
                hiss_min + transporter_min,
                hiss_min,
                transporter_min,
                t["contig"],
                t["start"],
            )

        target = min(
            targets,
            key=target_context_score
        )

        if len(targets) > 1:
            print(
                f"  NOTE: {accession} {target_protein}: "
                f"{len(targets)} CDS copies; selected "
                f"{target['contig']}:{target['start']}-{target['end']}"
            )

        same_contig = [
            x for x in cds
            if x["contig"] == target["contig"]
            and x["protein"] != target_protein
        ]

        hiss_hits = []
        transporter_hits = []

        for feature in same_contig:
            distance = interval_distance(
                target["start"],
                target["end"],
                feature["start"],
                feature["end"],
            )

            if is_hiss(feature):
                hiss_hits.append((distance, feature))

            if is_transporter(feature):
                transporter_hits.append((distance, feature))

        hiss_hits.sort(key=lambda x: x[0])
        transporter_hits.sort(key=lambda x: x[0])

        nearest_hiss = hiss_hits[0] if hiss_hits else None
        nearest_transporter = (
            transporter_hits[0]
            if transporter_hits
            else None
        )

        records.append(
            {
                "System": system,
                "Assembly Accession": accession,
                "Organism Name": row.get("Organism Name", ""),
                "Strain": row.get("Strain", ""),
                "Target Protein": target_protein,
                "Target Copy Count": len(targets),
                "Contig": target["contig"],
                "Target Start": target["start"],
                "Target End": target["end"],
                "Target Strand": target["strand"],

                "hisS_found_same_contig":
                    nearest_hiss is not None,

                "nearest_hisS_distance_bp":
                    nearest_hiss[0]
                    if nearest_hiss else pd.NA,

                "nearest_hisS_protein":
                    nearest_hiss[1]["protein"]
                    if nearest_hiss else "",

                "nearest_hisS_product":
                    nearest_hiss[1]["product"]
                    if nearest_hiss else "",

                "transporter_found_same_contig":
                    nearest_transporter is not None,

                "nearest_transporter_distance_bp":
                    nearest_transporter[0]
                    if nearest_transporter else pd.NA,

                "nearest_transporter_protein":
                    nearest_transporter[1]["protein"]
                    if nearest_transporter else "",

                "nearest_transporter_product":
                    nearest_transporter[1]["product"]
                    if nearest_transporter else "",
            }
        )

    return pd.DataFrame(records)


print("=" * 76)
print("RECONSTRUCTING HDC / HDC2 GENOMIC CONTEXT")
print("=" * 76)

hdc = pd.read_csv(HDC_TABLE, sep=chr(9))
hdc2 = pd.read_csv(HDC2_TABLE, sep=chr(9))
context = pd.read_csv(HDC2_CONTEXT, sep=chr(9))
exclusions = pd.read_csv(EXCLUSIONS, sep=chr(9), dtype=str)

excluded = set(
    exclusions["Assembly_Accession"].dropna()
)

print()
print("Archived candidate loci")
print("  Conventional HDC:", len(hdc))
print("  HDC2:            ", len(hdc2))

assert len(hdc) == 99
assert len(hdc2) == 36

# ------------------------------------------------------------
# Reconstruct actual annotated genomic neighborhoods
# ------------------------------------------------------------

hdc_reconstructed = reconstruct(
    hdc,
    "Conventional_HDC"
)

hdc2_reconstructed = reconstruct(
    hdc2,
    "HDC2"
)

# ------------------------------------------------------------
# Apply the recorded UBCG exclusions to conventional HDC
# ------------------------------------------------------------

hdc_final = hdc_reconstructed.loc[
    ~hdc_reconstructed["Assembly Accession"].isin(excluded)
].copy()

print()
print("Final cohort")
print("  Conventional HDC:", len(hdc_final))
print("  HDC2:            ", len(hdc2_reconstructed))

assert len(hdc_final) == 98
assert len(hdc2_reconstructed) == 36

# ------------------------------------------------------------
# Conventional HDC validation
# ------------------------------------------------------------

hdc_hiss_n = int(
    hdc_final["hisS_found_same_contig"].sum()
)

hdc_transporter_n = int(
    hdc_final["transporter_found_same_contig"].sum()
)

hdc_hiss_max = int(
    hdc_final["nearest_hisS_distance_bp"].max()
)

hdc_transporter_max = int(
    hdc_final["nearest_transporter_distance_bp"].max()
)

print()
print("Conventional HDC context")
print(
    f"  hisS on target contig:           "
    f"{hdc_hiss_n}/98"
)
print(
    f"  transporter/antiporter detected: "
    f"{hdc_transporter_n}/98"
)
print(
    f"  maximum nearest hisS distance:   "
    f"{hdc_hiss_max} bp"
)
print(
    f"  maximum nearest transporter distance: "
    f"{hdc_transporter_max} bp"
)

assert hdc_hiss_n == 98
assert hdc_transporter_n == 98

# Locked observed maxima
assert hdc_hiss_max == 183
assert hdc_transporter_max == 402

# ------------------------------------------------------------
# HDC2 validation
# ------------------------------------------------------------

hdc2_hiss_n = int(
    hdc2_reconstructed["hisS_found_same_contig"].sum()
)

hdc2_transporter_n = int(
    hdc2_reconstructed[
        "transporter_found_same_contig"
    ].sum()
)

hdc2_hiss_max = int(
    hdc2_reconstructed[
        "nearest_hisS_distance_bp"
    ].max()
)

hdc2_transporter_max = int(
    hdc2_reconstructed.loc[
        hdc2_reconstructed[
            "transporter_found_same_contig"
        ],
        "nearest_transporter_distance_bp"
    ].max()
)

print()
print("HDC2 context")
print(
    f"  hisS on target contig:           "
    f"{hdc2_hiss_n}/36"
)
print(
    f"  transporter detected:            "
    f"{hdc2_transporter_n}/36"
)
print(
    f"  maximum nearest hisS distance:   "
    f"{hdc2_hiss_max} bp"
)
print(
    f"  maximum nearest transporter distance: "
    f"{hdc2_transporter_max} bp"
)

assert hdc2_hiss_n == 36
assert hdc2_transporter_n == 32
assert hdc2_hiss_max == 376
assert hdc2_transporter_max == 3722

# ------------------------------------------------------------
# Verify archived HDC2 flags
# ------------------------------------------------------------

stored_hiss = hdc2[
    "has_hisS_nearby"
].astype(bool)

stored_transporter = hdc2[
    "has_transporter_nearby"
].astype(bool)

stored_cluster = hdc2[
    "HDC2_cluster_like"
].astype(bool)

assert stored_hiss.sum() == 36
assert stored_transporter.sum() == 32
assert stored_cluster.sum() == 32

assert (
    stored_cluster
    == (stored_hiss & stored_transporter)
).all()

print()
print("Archived HDC2 flag relationship")
print(
    "  HDC2_cluster_like = "
    "has_hisS_nearby AND has_transporter_nearby: YES"
)

# ------------------------------------------------------------
# Final HDC2 context classification
# ------------------------------------------------------------

status_counts = (
    context["Context_Status"]
    .value_counts()
    .to_dict()
)

complete_n = int(
    status_counts.get("COMPLETE_CLUSTER", 0)
)

boundary_n = int(
    status_counts.get("BOUNDARY_SUPPORTED", 0)
)

print()
print("Final HDC2 context classification")
print("  COMPLETE_CLUSTER:   ", complete_n)
print("  BOUNDARY_SUPPORTED: ", boundary_n)

assert complete_n == 32
assert boundary_n == 4

boundary_accessions = sorted(
    context.loc[
        context["Context_Status"]
        == "BOUNDARY_SUPPORTED",
        "Assembly Accession"
    ].tolist()
)

print()
print("Boundary-supported HDC2 assemblies:")
for accession in boundary_accessions:
    print(" ", accession)

# ------------------------------------------------------------
# Save reconstructed outputs
# ------------------------------------------------------------

hdc_final.to_csv(
    HDC_OUT,
    sep="\t",
    index=False
)

hdc2_reconstructed.to_csv(
    HDC2_OUT,
    sep="\t",
    index=False
)

summary = pd.DataFrame(
    [
        {
            "System": "Conventional_HDC",
            "Final_genomes": 98,
            "hisS_found_n": 98,
            "Transporter_found_n": 98,
            "Boundary_supported_n": 0,
            "Max_hisS_distance_bp": 183,
            "Max_transporter_distance_bp": 402,
        },
        {
            "System": "HDC2",
            "Final_genomes": 36,
            "hisS_found_n": 36,
            "Transporter_found_n": 32,
            "Boundary_supported_n": 4,
            "Max_hisS_distance_bp": 376,
            "Max_transporter_distance_bp": 3722,
        },
    ]
)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)

print()
print("Saved:")
print(" ", HDC_OUT.relative_to(ROOT))
print(" ", HDC2_OUT.relative_to(ROOT))
print(" ", SUMMARY_OUT.relative_to(ROOT))

print()
print("STATUS: PASS")
