# Reference files

VKaryo does not ship reference data. This page lists every file the workflow
expects, where to obtain it, and how to build the indexes. Set the resulting
paths in `config/my_config.yaml`.

Plan for roughly **500 GB** of disk and several hours for index building.

## 1. Genome and annotation

The workflow is built for **GRCh38/hg38**.

| File | Source |
|---|---|
| Genome FASTA | Ensembl: `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`, or GENCODE / UCSC `hg38.fa` |
| Gene annotation GTF | Ensembl release 113 (the release used for the published results) |

```bash
mkdir -p reference && cd reference
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Keep the FASTA, the `.fai` and the `.dict` side by side; GATK requires all
three.

Config entries: `refgenome.ref`, `refgenome.gtf`.

## 2. BWA index (exome branch)

```bash
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Several hours. Writes alongside the FASTA, so no separate config entry.

## 3. STAR index (transcriptome branch)

Needs about 32 GB of RAM and 30 GB of disk.

```bash
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir reference/star_index \
     --genomeFastaFiles reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile reference/Homo_sapiens.GRCh38.113.gtf \
     --sjdbOverhang 150
```

`--sjdbOverhang` is read length minus one; 150 suits 2 x 151 bp reads.

Config entry: `refgenome.star_index`.

## 4. Capture kit target BED

Must match the kit actually used to prepare the libraries. The published
results used the **Agilent SureSelect Human All Exon V6** target set in hg38
coordinates (244,655 intervals, 60.5 Mb). Agilent distributes the BED through
SureDesign, which requires a free account.

If your libraries used a different kit, use that kit's BED. The expected
coverage statistics the copy number steps rely on are kit-specific: feeding
the wrong target set produces wrong expectations, not merely noisier ones.

Config entry: `refgenome.exome_bed`.

## 5. Known variant databases (GATK)

For base quality score recalibration and germline filtering.

| Database | Source |
|---|---|
| dbSNP (v156) | NCBI, or the GATK resource bundle |
| gnomAD (v4.1) | https://gnomad.broadinstitute.org/downloads |
| Panel of normals | GATK resource bundle, `1000g_pon.hg38.vcf.gz` |

The GATK bundle is on Google Cloud:

```bash
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/ reference/gatk_bundle/
```

Without `gsutil`, the same files are downloadable over HTTPS from the
`storage.googleapis.com/genomics-public-data/` prefix.

Config entries: `refgenome.dbsnp`, `refgenome.gnomad`, `mutect2.pon`.

## 6. ANNOVAR databases

ANNOVAR requires separate registration and is **not** installed by conda. Get
it from https://annovar.openbioinformatics.org/ and download the databases
the workflow annotates against:

```bash
cd annovar
for db in refGene clinvar_20240917 cosmic101 gnomad41_genome dbnsfp47a avsnp154; do
  perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar $db humandb/
done
```

`cosmic` requires a COSMIC licence and is prepared separately; see the ANNOVAR
documentation.

Config entry: `annovar.humandb`.

## 7. AnnotSV

Installed by conda, but its annotation bundle is downloaded separately:

```bash
conda activate annotsv
export ANNOTSV=$CONDA_PREFIX/share/AnnotSV
cd $ANNOTSV && make PREFIX=. install-human-annotation
```

Config entry: `annotsv.path`.

## 8. Arriba resources

**Nothing to do.** The blacklist, the known-fusion list and the protein domain
file are installed with the Arriba conda package, and the workflow resolves
them at run time from `$CONDA_PREFIX/var/lib/arriba`. This is why they are not
config entries: the environment prefix differs on every machine.

If you need them outside the workflow, they are in that same directory after
`conda activate arriba`.

## 9. Files shipped with the repository

These are small, curated and versioned here — nothing to download:

| File | Purpose |
|---|---|
| `reference/driver_genes_leukemia.bed` | curated driver gene panel, used for prioritisation and for the low-VAF VarDict scan |
| `reference/linhagem_celular_por_amostra.tsv` | per-sample cell lineage, input to entity assignment |

### Cell lineage input

Entity assignment under WHO-HAEM5/ICC depends jointly on the defining genetic
alteration and on the lineage of the cell, because the same alteration defines
different entities in different lineages. The workflow therefore takes lineage
as a **declared input**, one row per sample:

```
sample_RNA	sample_DNA	label	lineage_stage	lineage_family	source
```

`lineage_stage` is one of `B_precursor`, `B_mature`, `plasmacytic`, `T`, `NK`,
`myeloid`. In clinical use this comes from flow cytometry immunophenotyping,
which is routine on peripheral blood and bone marrow aspirate and is specific
to the clone because it is read with a gate on the blasts.

Lineage is not inferred from expression: a patient sample is a cellular
mixture, so marker expression would report the composition of the aspirate
rather than the neoplastic clone.

If lineage is unknown for a sample, leave the fields empty — the entity is
then reported as assigned by the alteration alone, with a caveat, instead of
being discarded.

## Checklist

Before the first run, confirm each config path resolves:

```bash
snakemake --configfile config/my_config.yaml --rulegraph > /dev/null && echo OK
```

This builds the dependency graph without executing anything. It fails fast if
a reference file is missing or a path is wrong.
