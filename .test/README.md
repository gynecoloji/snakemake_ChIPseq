# Catalog test case (`.test/`)

This directory exists so the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
can render the workflow's **tube map** (rule graph) with
[snakevision](https://github.com/snakemake/snakevision). The catalog runs:

```bash
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph
```

`--rulegraph` only resolves the **rule dependency graph** — no rule is executed —
so the reference genomes and FASTQs here are **empty 0-byte placeholders**. They
exist only so the DAG resolves; nothing reads their contents. CI also runs a
`snakemake -d .test -n` dry run against this case.

The sample sheet deliberately exercises every code path: a **control-only** input,
a **narrow** single-replicate IP, and a **broad** 2-replicate IP (so the IDR-,
ratio-track- and both peak-mode rules all appear in the map). `config.yaml` also
defines a differential-binding `contrast`, so the downstream stage (annotation,
motifs, DESeq2, overlap, heatmaps) is covered too.

Contents:

- `config/config.yaml` — ChIP-seq config (paths resolve under `.test/`)
- `config/samples.csv` — 4 samples: 1 input, 1 narrow IP, 2 broad IP replicates
- `data/*.fastq.gz` — empty placeholder paired-end reads
- `ref/*` — empty placeholder reference genomes / annotations / BEDs

**This is not an end-to-end integration test.** Turning it into one (so the
workflow also earns the catalog's "tests" ranking) would require real miniature
reference genomes + reads that actually run through the tools, plus `picard.jar`
and the per-rule conda environments — a larger, separate effort.
