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
a **narrow** 2-replicate IP condition, and a **broad** 2-replicate IP condition.
So both peak modes, both IDR / reproducibility branches (narrow *and* broad),
the input-control + ratio-track rules, and a runnable **2-vs-2 differential-binding
contrast** all appear in the map — i.e. the QC and downstream stages (annotation,
motifs, DESeq2, overlap, heatmaps) are fully covered.

Contents:

- `config/config.yaml` — ChIP-seq config (paths resolve under `.test/`)
- `config/samples.csv` — 5 samples: 1 input, 2 narrow IP reps (H3K4me3), 2 broad IP reps (H3K27me3)
- `data/*.fastq.gz` — empty placeholder paired-end reads
- `ref/*` — empty placeholder reference genomes / annotations / BEDs

**This is not an end-to-end integration test.** Turning it into one (so the
workflow also earns the catalog's "tests" ranking) would require real miniature
reference genomes + reads that actually run through the tools, plus `picard.jar`
and the per-rule conda environments — a larger, separate effort.
