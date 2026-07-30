# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1](https://github.com/gynecoloji/snakemake_ChIPseq/compare/v0.1.0...v0.1.1) (2026-07-29)


### Fixed

* trigger release-please after 1.0.0 bootstrap ([fc70709](https://github.com/gynecoloji/snakemake_ChIPseq/commit/fc70709504dd270134b58ae03f13bf0368511f7a))

## [0.1.0] - 2026-07-14

Initial ChIP-seq workflow, adapted from the `snakemake_ATACseq_spikein` template.

### Added

- Primary ChIP-seq pipeline (`chipseq_all` target): FastQC/fastp → human Bowtie2
  alignment → unique/proper-pair filtering → Picard dedup → fragment-level ENCODE
  blacklist removal → MACS2 peak calling **against a matched input/IgG control**,
  with **per-sample narrow or broad** mode chosen in `config/samples.csv`.
- Per-sample signal tracks: RPGC depth-normalized bigWigs (all samples) and
  **IP-over-input `log2` ratio** bigWigs (`bamCompare`) for IP samples with a
  matched control.
- Reproducible fixed-width **consensus peaks** (Corces-2018 SPM iterative overlap;
  majority-vote / IDR / single reproducibility by replicate count) that handle
  **both narrowPeak and broadPeak**, plus a `featureCounts` fragment matrix.
- QC pipeline (`qc_all` target): fingerprint, fragment size, correlation/PCA, GC
  bias, TSS signal + numeric score, FRiP, library complexity, IDR, reads in
  promoters/enhancers, a FastQC-only MultiQC, and a self-contained **interactive
  HTML QC report**.
- Sample sheet columns `sample_id, condition, replicate, input_control, igg_control,
  peak_mode, notes` with fail-fast validation (valid peak modes, existing controls,
  one peak mode per condition).
- **Choice of control for peak calling**: `control_type` config parameter
  (`input` by default, or `igg`) selects whether each IP uses its `input_control`
  or `igg_control` as the MACS2 `-c` control, with fallback to the other column
  when the chosen one is empty. The IP-over-control `log2` ratio bigWig uses the
  same resolved control.
- **ENCODE QC panel**: strand cross-correlation (**NSC/RSC** + fragment length +
  QualityTag via phantompeakqualtools), fingerprint **Jensen-Shannon distance**
  (deepTools `plotFingerprint --JSDsample`), **IDR reproducibility** self-consistency
  + rescue ratios (self- and pooled-pseudoreplicate IDR for 2-replicate conditions),
  and ENCODE per-mode usable-read-depth flags (narrow ≥ 20M, broad ≥ 45M) in the
  interactive report. Adds a `phantompeakqualtools` conda env (6 envs total).
- **Downstream analysis stage** (`downstream_all` target): **peak annotation + GO**
  (ChIPseeker + clusterProfiler, per IP sample), **motif enrichment** (HOMER, per IP
  sample), **differential binding** (DESeq2 over the consensus count matrix per
  `contrasts:` entry, with ChIPseeker-annotated significant peaks + MA/volcano/PCA),
  **peak-set overlap** (bedtools Jaccard matrix + heatmap), and **signal heatmaps /
  metagene** (deepTools peak-centered + gene-body). Adds `chipseeker` (R/Bioconductor)
  and `homer` conda envs (8 envs total) and a `contrasts` config parameter.
  Differential binding runs only for contrasts with ≥2 replicates per condition;
  single-replicate contrasts are skipped automatically (DESeq2 needs replicates).
- Config schema (`workflow/schemas/config.schema.yaml`) with parameter validation,
  a Docker/Apptainer image, `.test/` catalog case, and CI (unit tests + dry runs).
