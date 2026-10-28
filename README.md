# VKaryo

An automated Snakemake workflow that integrates total RNA-seq and whole-exome
sequencing (WES) to detect gene fusions, chromosomal translocations, copy
number variations (CNVs) and single-nucleotide variants (SNVs) in a single
coordinated run, and to support molecular classification and prognostic
stratification of haematological neoplasms.

The detection steps do not depend on the disease under analysis and apply to
cancer and non-cancer samples alike. The classification modules are specific
to haematological neoplasms.

> **Status.** This is release `v0.1.0`, the version accompanying the
> manuscript. The local web interface is not part of this release and is
> planned for a later version.

## What it does

**From RNA-seq** — quality control and adapter trimming (fastp), splice-aware
alignment (STAR, two-pass), fusion detection (Arriba) with confidence
stratification.

**From WES** — trimming (fastp), alignment (BWA-MEM), duplicate marking and
base quality score recalibration (GATK4), somatic SNVs and indels (Mutect2,
tumour-only, with a panel of normals and gnomAD germline filtering), low-VAF
variants for minimal residual disease (VarDict), structural variants from two
orthogonal callers (Manta and DELLY) merged by SURVIVOR and annotated with
AnnotSV, and copy number from CNVkit.

**Copy number, beyond the ratio** — with a flat reference the log2 ratio loses
quantitative meaning as segment coverage approaches zero, so deep losses are
assessed directly by observed versus expected coverage per measurement window.
Three steps build on this: per-bin cohort statistics, focal event detection
(homozygous loss, heterozygous loss, gain), and chromosome-arm imbalance by
amplitude combined with extent. Events falling in immunoglobulin and T-cell
receptor loci are flagged, not silently discarded, since V(D)J recombination
physically excises those segments in any lymphoid sample.

**Classification** — a curated driver gene panel for prioritisation, entity
assignment under WHO-HAEM5/ICC resolved jointly by the defining alteration and
the cell lineage of the sample, and risk stratification where the system
applies.

## Requirements

- Linux (Ubuntu 20.04 or newer)
- Conda or Mamba
- Snakemake >= 7.0
- Python >= 3.9
- 64 GB RAM minimum
- Roughly 500 GB of disk for reference files

Per-rule software is declared in `envs/` and installed by Snakemake itself
with `--use-conda`. You do not need to install the tools by hand.

## Installation

```bash
git clone https://github.com/brunorsci/VKaryo0.1.git
cd VKaryo0.1
conda env create -f envs/vkaryo.yaml
conda activate vkaryo
```

## Reference files

Reference data is **not** shipped with this repository: the genome, its
indexes and the annotation databases are tens of gigabytes, and most groups
already hold their own copies. `docs/references.md` lists every file the
workflow expects, where to obtain it, and how to build the indexes.

## Configuration

Copy the template and edit the paths:

```bash
cp config/config.yaml config/my_config.yaml
```

Every path the workflow uses is resolved from this file. `docs/references.md`
explains each entry.

## Running

```bash
snakemake --configfile config/my_config.yaml \
          --cores 16 \
          --use-conda \
          --conda-frontend conda \
          --rerun-triggers mtime
```

Adjust `--cores` to the machine. Peak memory observed on a 100-sample cohort
was 29 GB per sample.

To inspect the workflow graph before running anything:

```bash
snakemake --configfile config/my_config.yaml --rulegraph | dot -Tsvg > rulegraph.svg
```

## Example

`docs/example.md` describes a small self-contained example: synthetic reads
carrying a planted fusion and a planted homozygous deletion, against a
minimal reference. It runs in minutes and the expected output is stated
explicitly, so it doubles as a regression check.

## Outputs

Per sample, under the output directory given in the config:

| Path | Content |
|---|---|
| `arriba/` | fusion calls with Arriba confidence |
| `mutect2/`, `vardict/` | somatic SNVs and indels, annotated |
| `survivor/`, `annotsv/` | consensus structural variants, annotated |
| `cnvkit/` | copy number segments |
| `cnv_focal/`, `cnv_arm/` | coverage-based focal events and arm imbalance |
| `classification/` | entity assignment and risk stratification |

Filtering and clinical interpretation of the final tables require review by a
qualified bioinformatician. The workflow deliberately preserves calls flagged
by heuristic filters, with the flag recorded, rather than discarding them
before anyone can look.

## Citing

If you use VKaryo, please cite the software through its archived release and
the accompanying article. Metadata for the software citation is in
`CITATION.cff`.

## Licence

MIT. See `LICENSE`.
