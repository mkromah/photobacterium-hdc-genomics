# Comparative Genomics and Phylogenomics of Histidine Decarboxylase Systems in *Photobacterium*

This repository contains the reproducible analysis workflow and curated results for a comparative-genomics study of conventional histidine decarboxylase (HDC) and the distinct HDC2 system across *Photobacterium*.

The broader food-safety motivation is to understand the genomic distribution and evolution of histamine-forming potential in *Photobacterium*, a genus relevant to seafood-associated histamine formation.

## Research questions

1. How are conventional HDC and HDC2 distributed across *Photobacterium*?
2. Are HDC/HDC2 states associated with particular taxonomic groups?
3. How conserved are the genomic neighborhoods surrounding the two systems?
4. Do HDC and HDC2 gene phylogenies agree with the core-genome species phylogeny?
5. Is there genomic evidence consistent with horizontal acquisition of HDC2-associated regions?

## Study cohort

- Initial primary-genome cohort: 500 assemblies
- Final analyzed cohort: 497 assemblies
- UBCG exclusions: 3 assemblies
- Conventional HDC-positive genomes: 98
- HDC2-positive genomes: 36
- HDC-only: 95
- HDC2-only: 33
- Both systems: 3
- Neither system: 366

The accession lists and cohort metadata are provided under `metadata/`.

## Main findings

### Taxonomic distribution

HDC-system state was strongly associated with taxonomic group.

- Chi-square: 509.0521
- Degrees of freedom: 42
- Cramer's V: 0.5843
- Monte Carlo permutation test: 100,000 permutations
- Empirical p-value: 9.9999e-06

Because many expected contingency-table counts were small, the Monte Carlo permutation result is emphasized rather than the asymptotic chi-square p-value.

### Genomic context

Conventional HDC:

- 98/98 loci had nearby `hisS`
- 98/98 loci had a nearby histidine/histamine transporter or antiporter

HDC2:

- 36/36 loci had nearby `hisS`
- 32/36 had an evaluable nearby transporter
- 32 loci were classified as `COMPLETE_CLUSTER`
- 4 were `BOUNDARY_SUPPORTED` because assembly boundaries prevented full neighborhood evaluation

### Gene-tree / species-tree comparison

Strong support was defined as:

- SH-aLRT >= 80
- UFBoot >= 95

Primary archived analyses:

- Conventional HDC: 4 strong splits; 3 concordant and 1 discordant with the species tree
- HDC2: 4 strong splits; 1 concordant and 3 discordant with the species tree

Seeded reconstruction reproduced all 4 HDC2 strong splits. For conventional HDC, 3/4 remained above the formal threshold; the fourth retained the same bipartition but shifted from UFBoot 96 to 94.

### HDC2 horizontal-acquisition case study

The strongest localized signal involved:

- *P. phosphoreum* FS3_1 — GCF_001676275.1
- *P. toruni* S2MW2 — GCF_036040135.1

Validated results:

- HDC2 global comparable-site amino-acid identity: 97.059%
- Cross-species global-identity rank: 42/162
- BLASTP local alignment: 655 aa at 100.000% identity
- Whole-genome ANI: 86.3171% and 86.3201% reciprocally
- Mean reciprocal ANI: 86.3186%
- 16-kb HDC2-associated locus identity: 99.513%
- 16-kb query coverage: 100%
- FS3_1 HDC2-bearing contig: 23,907 bp
- Near-identical contig coverage: 22,897 bp (95.78%)
- HSP-length-weighted identity: 99.590%
- Near-full, near-identical cassette matches: 2/36 HDC2-positive genomes
- Largest low-prevalence block: 7,915-15,778 bp (7,864 bp)

These data are strongly consistent with recent horizontal acquisition or shared recent acquisition of an HDC2-associated genomic region. Transfer direction cannot be resolved, and acquisition from an unsampled donor cannot be excluded.

## Reproducibility workflow

The main analysis scripts are in `scripts/`.

1. `01_deduplicate_assemblies.py` — removes redundant assemblies from the starting cohort.
2. `run_ubcg500.sh` — runs UBCG2 core-gene analysis across the primary cohort.
3. `02_identify_HDC_HDC2.py` — reconstructs genome-level HDC/HDC2 calls.
4. `03_reconstruct_genomic_context.py` — reconstructs HDC and HDC2 genomic neighborhoods.
5. `04_taxon_distribution_statistics.py` — reproduces taxon-by-HDC-state summaries, chi-square statistics, Monte Carlo permutations, and adjusted standardized residuals.
6. `05_reconstruct_gene_phylogeny_inputs.py` — reconstructs HDC and HDC2 protein FASTA files from assembly annotations.
7. `06_reconstruct_gene_alignments.py` — recreates MAFFT alignments and alignment QC.
8. `07_reconstruct_gene_phylogenies.py` — rebuilds IQ-TREE gene phylogenies with model selection and branch support.
9. `08_HDC2_FS3_sensitivity.py` — repeats the HDC2 phylogeny after excluding the shortened FS3_1 HDC2 sequence.
10. `09_compare_gene_species_trees.py` — compares strongly supported gene-tree splits with corresponding species-tree splits.
11. `10_reconstruct_HDC2_horizontal_acquisition.py` — reconstructs the FS3_1-S2MW2 HDC2 horizontal-acquisition evidence.

## Software environment

A Conda environment is provided in `environment.yml`.

Key validated software versions:

- Python 3.11.16
- BLAST+ 2.17.0+
- fastANI 1.34
- MAFFT 7.526
- IQ-TREE 3.1.3
- pandas 3.0.5
- NumPy 2.4.6
- Biopython 1.88

Create the environment with:

```bash
conda env create -f environment.yml
conda activate photobacterium
```

## UBCG2

UBCG2 is required for the core-genome phylogeny but is not redistributed in this repository.

The local workflow expects the UBCG2 JAR and associated database files to be available separately. See `scripts/run_ubcg500.sh` for the expected invocation.

## Data availability

Raw NCBI genome files are not included because of repository size.

The repository includes:

- accession and cohort metadata under `metadata/`
- HDC/HDC2 reference/query sequences under `references/`
- curated final and intermediate analysis outputs under `results/`
- publication-quality figures under `figures/`
- reproducibility scripts under `scripts/`

The raw assemblies can be re-obtained from NCBI using the assembly accessions provided in the metadata files.

## Repository structure

```text
.
├── environment.yml
├── figures/
├── metadata/
├── references/
├── results/
│   ├── gene_phylogeny/
│   ├── phylogeny/
│   ├── reproducibility/
│   └── synteny/
└── scripts/
```

## Notes on interpretation

- HDC2 transporter absence was not inferred for loci truncated by contig boundaries.
- No historical fixed genomic-neighborhood distance threshold was invented when it could not be recovered.
- The HDC2 FS3_1-S2MW2 analysis does not establish transfer direction.
- The previously used value of 98.09% HDC2 amino-acid identity was not reproducible and was corrected. The validated global comparable-site identity is 97.059%, while BLASTP identifies a 655-aa local region at 100.000% identity.
- Bootstrap support can vary slightly between independent phylogenetic reruns; the repository preserves both the primary archived analysis and seeded reconstruction checks.

## Citation

Citation information will be added after manuscript publication.
