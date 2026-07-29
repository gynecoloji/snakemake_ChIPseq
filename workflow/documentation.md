# Technical documentation

Step-by-step documentation of the `snakemake_ChIPseq` workflow. For installation,
container usage, and the full narrative, see the top-level
[`README.md`](../README.md); for every configuration parameter, see
[`config/README.md`](../config/README.md) and the schema
[`workflow/schemas/config.schema.yaml`](schemas/config.schema.yaml).

The rule graph is rendered as a "tube map" on the workflow's Snakemake Workflow
Catalog page.

## Overview

A single `snakemake -s workflow/Snakefile --use-conda` run builds a **unified
DAG** covering two stages in dependency order:

1. **Primary** (`chipseq_all` target) — alignment → filtering → MACS2 peak calling
   (with a matched input/IgG control, per-sample narrow/broad) → RPGC + IP-over-input
   bigWigs → reproducible consensus peaks + fragment counts.
2. **QC** (`qc_all` target) — deepTools QC, FRiP, IDR, library complexity, TSS
   signal, reads-in-annotation, ENCODE metrics, and an interactive HTML QC report.
3. **Downstream** (`downstream_all` target) — peak annotation + GO, motif
   enrichment, differential binding, peak-set overlap, and signal heatmaps/metagene.

## Inputs

| Input | Location | Notes |
|---|---|---|
| Paired-end reads | `data/<sample_id>_R1_001.fastq.gz`, `_R2_001.fastq.gz` | one pair per sample (IP + control) |
| Sample sheet | `config/samples.csv` | columns `sample_id, condition, replicate, input_control, igg_control, peak_mode, notes` |
| Human genome FASTA | `ref/hg38.fa` | chr-prefixed UCSC |
| Blacklist BED | `ref/hg38_blacklist_regions.bed` | ENCODE, chr-prefixed |
| GTF / 2bit / promoter+enhancer BEDs | `ref/…` | QC references |
| Picard | `ref/picard.jar` | duplicate marking |

Configuration is read from `config/config.yaml` and validated against the schema
at parse time (missing/invalid parameters fail fast). The sample sheet is also
validated (valid peak modes, existing input controls, one peak mode per condition).

## Steps (primary stage)

1. **`fastqc`** — raw-read quality.
2. **`fastp`** — adapter trimming + quality filtering (auto-detects adapters).
3. **`build_bowtie2_index`** — build one Bowtie2 index from the human genome
   (optionally subset to `align_chroms`).
4. **`bowtie2_align`** — align paired-end reads to the human index.
5. **`samtools_sort_filter_index`** — keep uniquely-mapped, properly-paired reads;
   record mitochondrial-% QC; restrict to the analysis chromosomes.
6. **`remove_duplicates`** — Picard MarkDuplicates.
7. **`filter_blacklist`** — fragment-level ENCODE blacklist removal.
8. **`call_peaks_narrow` / `call_peaks_broad`** — MACS2 (BAMPE, `-q`), with
   `-c <control>` when the sample has a resolved control (Input or IgG, per
   `control_type`; default Input); `--broad --broad-cutoff` for broad samples.
   Control-only samples (empty `peak_mode`) are skipped here.
9. **Signal tracks** — `create_bigwig` (RPGC depth-normalized, all samples) and
   `create_ratio_bigwig` (`bamCompare` log2 IP/control, for IP samples that have a
   resolved control).
10. **Consensus peaks** (`relaxed_peaks_*`, `reproducible_idr_*`, `consensus_peaks`,
    `count_fragments_consensus`) — per-condition reproducibility (majority vote for
    ≥3 replicates, IDR for exactly 2), a fixed-width consensus set that handles both
    narrowPeak and broadPeak, and a featureCounts fragment matrix over all samples.

## Steps (QC stage)

deepTools coverage/fragment-size/fingerprint/correlation/PCA/GC/TSS, a numeric
TSS-signal score, FRiP (per IP sample), IDR on relaxed peaks, library complexity
(NRF/PBC1/PBC2), reads-in-annotation and peak summaries, and the **ENCODE** QC
panel — strand cross-correlation (**NSC/RSC**, phantompeakqualtools), fingerprint
**Jensen-Shannon distance** (deepTools `--JSDsample`), and IDR reproducibility
(**self-consistency + rescue ratios** from self/pooled pseudo-replicates) — plus a
FastQC-only MultiQC report and a self-contained interactive HTML QC report
(`results/qc/chipseq_qc_report.html`).

ENCODE QC rules: `cross_correlation` (+ `cross_correlation_summary`),
`fingerprint_jsd` (+ `fingerprint_jsd_summary`), and the pseudo-replicate chain
`repro_pool` → `repro_split` → `repro_relaxed_{narrow,broad}` →
`repro_idr_{narrow,broad}` → `idr_reproducibility_summary` (2-replicate conditions).

## Steps (downstream stage)

`annotate_peaks` (ChIPseeker annotation + clusterProfiler GO, per IP sample),
`motif_enrichment` (HOMER `findMotifsGenome.pl`, per IP sample),
`differential_binding` (DESeq2 per `contrasts:` entry over the consensus matrix,
with ChIPseeker-annotated significant peaks), `peak_jaccard` + `peak_overlap_matrix`
(bedtools Jaccard matrix + heatmap), and `deeptools_peak_heatmap` +
`deeptools_metagene` (peak-centered and gene-body signal). Differential binding is
skipped when `contrasts` is empty.

## Outputs

All outputs are written under `results/` (peaks, bigWigs, consensus matrix, QC
tables and reports); per-rule logs under `logs/`. See the README's "Output Files"
section for the full tree.

## Running the tests

```bash
python -m pytest tests/ -q                               # unit tests
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph   # DAG/tube map
snakemake -s workflow/Snakefile -d .test -n              # dry run
```
