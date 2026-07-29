# Configuration

This workflow is configured through two files in this directory:

- `config.yaml` — all workflow parameters (see below)
- `samples.csv` — the sample sheet

plus reference data you download into `ref/` (not tracked in git; see
[Reference data](#reference-data)).

## Sample sheet (`config/samples.csv`)

CSV with one row per sample and these columns:

| column          | description |
|-----------------|-------------|
| `sample_id`     | Sample name. Raw reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`. |
| `condition`     | IP target / biological condition label. **All IP rows that share a `condition` are treated as replicates of one group** (drives consensus/IDR reproducibility). Control samples typically use `Input`. |
| `replicate`     | Replicate index within the condition (1, 2, …). |
| `input_control` | `sample_id` of the matched **Input** control for this IP. |
| `igg_control`   | `sample_id` of the matched **IgG** control for this IP. |
| `peak_mode`     | `narrow`, `broad`, or **empty**. Empty ⇒ a **control-only** sample (aligned + bigWig, usable as a control, but no peaks and not in the consensus). `narrow`/`broad` ⇒ an **IP** sample, peak-called in that mode. |
| `notes`         | Free text. |

The **choice between IgG and Input** as the control is made once for the run with
the `control_type` parameter in `config.yaml` (`input` by default, or `igg`); see
[Choosing the control](#choosing-the-control-igg-vs-input) below.

Example (the shipped sheet — OVCAR3 cJUN/IgG ChIP with matched inputs; each cJUN IP
lists both its Input and its IgG control):

```csv
sample_id,condition,replicate,input_control,igg_control,peak_mode,notes
GSF2801-ChIPseq-OVCAR3-Control-Input_S3,Input,1,,,,Control input (control-only)
GSF2801-ChIPseq-OVCAR3-Control-IP-cJun_S1,Ctrl_cJUN,1,GSF2801-ChIPseq-OVCAR3-Control-Input_S3,GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,narrow,Control cJUN
GSF2801-ChIPseq-OVCAR3-Control-IP-IgG_S2,Ctrl_IgG,1,GSF2801-ChIPseq-OVCAR3-Control-Input_S3,,narrow,Control IgG
GSF2801-ChIPseq-OVCAR3-3D-Input_S6,Input,1,,,,3D input (control-only)
GSF2801-ChIPseq-OVCAR3-3D-IP-cJun_S4,3D_cJUN,1,GSF2801-ChIPseq-OVCAR3-3D-Input_S6,GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,narrow,3D cJUN
GSF2801-ChIPseq-OVCAR3-3D-IP-IgG_S5,3D_IgG,1,GSF2801-ChIPseq-OVCAR3-3D-Input_S6,,narrow,3D IgG
```

**How the columns drive the pipeline:**

- **Peak mode is per sample.** Set `peak_mode` to `broad` for broad marks
  (e.g. H3K27me3, H3K9me3, H3K36me3) and `narrow` for point-source factors and
  sharp marks (e.g. transcription factors, H3K4me3). Every IP row can choose
  independently — MACS2 runs `--broad --broad-cutoff` for broad rows.
- **Both controls are per sample.** List each IP's matched Input in `input_control`
  and its matched IgG in `igg_control` (either may be empty). Which one MACS2
  actually uses as `-c` is selected run-wide by `control_type` (see below).
- **Control-only samples** (empty `peak_mode`, e.g. Input, or an IgG you only use
  as a control) are still aligned, deduplicated and turned into bigWigs, and can be
  named as another sample's `input_control` / `igg_control`, but they are never
  peak-called and never enter the consensus.

## Choosing the control (IgG vs. Input)

`control_type` in `config.yaml` picks which control each IP uses as its MACS2 `-c`,
for the whole run:

```yaml
control_type: "input"   # "input" (default) or "igg"
```

- `input` (**default**) → each IP uses its `input_control`.
- `igg` → each IP uses its `igg_control`.
- **Fallback:** if the selected column is empty for a sample, the other column is
  used; if both are empty, that IP is called treatment-only (no `-c`).

So to compare Input- vs. IgG-based calls, keep both columns filled in `samples.csv`
and just flip `control_type` (or override per run without editing the file:
`snakemake ... --config control_type=igg`). The IP-over-control `log2` ratio bigWig
(`results/ratio_bigwig/`) uses the same resolved control.

**Per-condition reproducibility** is derived automatically from the number of IP
replicates sharing a `condition`:

- **≥ 3 replicates** → majority vote (a peak is kept if it recurs in ≥
  `consensus_min_replicates` replicates).
- **exactly 2 replicates** → IDR (`idr_threshold`).
- **1 replicate** → the sample's own peaks are used as-is.

All replicates of a condition must share one `peak_mode` (validated on load) so
the per-group consensus/IDR is well-defined. If a single ChIP target spans several
biological conditions, give each condition a distinct `condition` name (e.g.
`Ctrl_cJUN` vs `3D_cJUN`, as above) so replicates group correctly.

## Differential-binding contrasts

The downstream stage runs DESeq2 over the consensus count matrix for each contrast
listed under `contrasts:` in `config.yaml`. Each entry names two `condition` values
from the sample sheet (A = test, B = reference; log2FC > 0 means higher in A):

```yaml
contrasts:
  - name: cJUN_3D_vs_Ctrl
    condition_a: "3D_cJUN"
    condition_b: "Ctrl_cJUN"
```

Leave the list empty (`contrasts: []`) to skip differential binding. Replicates per
condition are recommended — a 1-vs-1 contrast still runs but yields fold-changes with
unreliable p-values. The other downstream analyses (peak annotation + GO, motif
enrichment, peak overlap, signal heatmaps) need no configuration.

## Parameters (`config/config.yaml`)

Every parameter — with its type, default, and description — is defined once in the
config schema, [`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml).
That schema is the single source of truth: the workflow validates `config.yaml`
against it on every run (and fills in defaults for anything you omit), and the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
renders it as a parameter table on the workflow page.

To configure a run, edit `config.yaml` directly — it ships with working defaults
and an inline comment on every parameter. At minimum, point the reference-file
paths (`human_fasta`, `blacklist`, `gtf`, `promoter_bed`, `enhancer_bed`) at the
files you provide (see [Reference data](#reference-data)). Peak mode and the input
control are set **per sample in `samples.csv`**, not here; `config.yaml` holds only
the shared MACS2 parameters (`macs2_genome`, `macs2_qvalue`, `broad_cutoff`).

## Reference data

Genomes, indexes and large annotations are **not** shipped in the repo (they are
`.gitignore`d). Download / place them under `ref/` before running, matching the
paths in `config.yaml`:

- `ref/hg38.fa` — chr-prefixed UCSC human genome
- `ref/hg38_blacklist_regions.bed` — ENCODE hg38 blacklist (shipped)
- `ref/gencode.v36.annotation.gtf` — GENCODE annotation (for TSS QC)
- `ref/hg38.2bit` — for `computeGCBias`
- `ref/picard.jar` — Picard (used by MarkDuplicates)

The human Bowtie2 index (`ref/BOWTIE2/`) is built automatically by the
`build_bowtie2_index` rule from `human_fasta`.

See the top-level `README.md` for full setup and run instructions.
