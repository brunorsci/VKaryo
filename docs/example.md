# Worked example

A self-contained example that exercises both branches of the workflow against
a minimal reference, with reads carrying a **planted fusion** and a **planted
homozygous deletion**. Because the answer is known by construction, the
example doubles as a regression check.

It is deliberately small: the full human reference takes hours to index and
tens of gigabytes of disk, which is why nobody runs a "quick example" against
it. Here a single chromosome is used, and the whole thing finishes in minutes
on a laptop.

## What is planted

| Event | Where | Expected call |
|---|---|---|
| Gene fusion | breakpoint joining two genes on the mini reference | reported by Arriba with supporting reads |
| Homozygous deletion | a ~100 kb interval with zero coverage | reported by `cnv_focal_events` as homozygous loss |
| Neutral region | the rest of the chromosome | no event |

## Building the test data

```bash
python scripts/make_test_dataset.py --outdir test/
```

This writes:

```
test/
  reference/   mini FASTA + GTF + target BED, and the indexes
  reads/       sample1 and sample2, paired FASTQ for RNA and WES
  expected/    the calls the run should produce
```

The script needs only Python, `samtools`, `bwa` and `STAR`, all of which the
conda environments already provide. Building the mini STAR index takes a few
minutes and a few hundred MB — not the 32 GB the full genome needs.

## Running

```bash
snakemake --configfile config/config.example.yaml \
          --cores 4 \
          --use-conda \
          --conda-frontend conda
```

## Checking the result

```bash
python scripts/check_test_results.py --outdir results_example/ --expected test/expected/
```

It compares the calls against `test/expected/` and exits non-zero on any
mismatch, so it can be wired into CI.

Reading the output by hand instead:

- `results_example/output_vkaryo/arriba/sample1_fusions.tsv` should contain the
  planted pair, with Arriba confidence `high`.
- `results_example/output_vkaryo/cnv_focal/sample1_focal_events.tsv` should
  contain one `homozygous_loss` event spanning the planted interval, and no
  event over the neutral region.

## Scope

The example proves the workflow is correctly installed, that every reference
path resolves and that both branches run end to end. It does not establish
analytical sensitivity: that comes from validation against a characterised
cohort, reported in the accompanying article.
