# Shared setup for the primary (chipseq.smk) and QC (qc.smk) rule files:
# imports, config-derived values, output-directory constants, the ChIP-seq
# sample sheet (IP vs. control, per-sample peak mode, matched input control),
# replicate groups and helper functions. Included first by workflow/Snakefile,
# so every name defined here is visible to the rules in the other files.

import pandas as pd
import re
import os
import sys
from itertools import combinations
from snakemake.utils import validate

# ── Config validation ───────────────────────────────────────────────────
# Validate against workflow/schemas/config.schema.yaml (also fills in defaults
# for any omitted parameters). Path is relative to this file (workflow/rules/).
validate(config, "../schemas/config.schema.yaml")

# ── Sample sheet ────────────────────────────────────────────────────────
# Columns: sample_id, condition, replicate, input_control, igg_control, peak_mode, notes.
#   * peak_mode empty            -> control-only sample (Input/IgG): aligned + bigWig,
#                                   usable as a MACS2 -c control, but NOT peak-called.
#   * peak_mode narrow | broad   -> IP sample: peaks called in that mode.
#   * input_control              -> sample_id of the matched Input control.
#   * igg_control                -> sample_id of the matched IgG control.
#   * condition                  -> replicate group (all IP rows sharing a condition
#                                   are replicates; reproducibility follows the count).
# The config `control_type` (input|igg, default input) selects which of the two an
# IP uses as its MACS2 -c control; the pipeline falls back to the other column when
# the selected one is empty for that sample, and to no control if both are empty.
samples_df = pd.read_csv(config["samples_table"], dtype=str).fillna("")
for col in ("sample_id", "condition", "replicate", "input_control", "igg_control", "peak_mode", "notes"):
    if col not in samples_df.columns:
        samples_df[col] = ""
    samples_df[col] = samples_df[col].astype(str).str.strip()
samples_df["peak_mode"] = samples_df["peak_mode"].str.lower()

SAMPLES = samples_df["sample_id"].tolist()

# Per-sample maps
PEAK_MODE  = dict(zip(samples_df["sample_id"], samples_df["peak_mode"]))
CONDITION  = dict(zip(samples_df["sample_id"], samples_df["condition"]))
INPUT_CTRL = {s: (c or None) for s, c in zip(samples_df["sample_id"], samples_df["input_control"])}
IGG_CTRL   = {s: (c or None) for s, c in zip(samples_df["sample_id"], samples_df["igg_control"])}

# Which control drives peak calling: "input" (default) or "igg". The chosen column
# wins; if it is empty for a sample, fall back to the other; else no control.
CONTROL_TYPE = str(config["control_type"]).lower()

def _resolve_control(sample):
    inp, igg = INPUT_CTRL.get(sample), IGG_CTRL.get(sample)
    primary, secondary = (inp, igg) if CONTROL_TYPE == "input" else (igg, inp)
    return primary or secondary or None

# Effective MACS2 control per sample (respects control_type + fallback)
CONTROL = {s: _resolve_control(s) for s in SAMPLES}

# IP (peak-called) vs. control-only samples
IP_SAMPLES      = [s for s in SAMPLES if PEAK_MODE.get(s) in ("narrow", "broad")]
CONTROL_SAMPLES = [s for s in SAMPLES if PEAK_MODE.get(s) not in ("narrow", "broad")]
NARROW_SAMPLES  = [s for s in IP_SAMPLES if PEAK_MODE[s] == "narrow"]
BROAD_SAMPLES   = [s for s in IP_SAMPLES if PEAK_MODE[s] == "broad"]
# IP samples with a resolved control get an IP-over-control (log2) ratio track
RATIO_SAMPLES   = [s for s in IP_SAMPLES if CONTROL.get(s)]

# ── Sample-sheet validation (fail fast) ─────────────────────────────────
_dups = sorted({s for s in SAMPLES if SAMPLES.count(s) > 1})
if _dups:
    raise ValueError(f"samples.csv: duplicate sample_id(s): {_dups}")

_bad_modes = sorted({PEAK_MODE[s] for s in SAMPLES if PEAK_MODE[s] not in ("narrow", "broad", "")})
if _bad_modes:
    raise ValueError(f"samples.csv: invalid peak_mode value(s): {_bad_modes}; "
                     "use 'narrow', 'broad', or leave empty for a control-only sample")

if CONTROL_TYPE not in ("input", "igg"):
    raise ValueError(f"config: control_type must be 'input' or 'igg', got {CONTROL_TYPE!r}")

_ctrl_refs = {c for s in IP_SAMPLES for c in (INPUT_CTRL.get(s), IGG_CTRL.get(s)) if c}
_missing_ctrl = sorted(c for c in _ctrl_refs if c not in SAMPLES)
if _missing_ctrl:
    raise ValueError("samples.csv: input_control/igg_control refers to unknown "
                     f"sample_id(s): {_missing_ctrl}")

if not IP_SAMPLES:
    raise ValueError("samples.csv: no IP samples found (every row has an empty peak_mode). "
                     "At least one row must set peak_mode to 'narrow' or 'broad'.")

# ── Replicate groups (IP samples grouped by condition) ──────────────────
GROUPS = (samples_df[samples_df["sample_id"].isin(IP_SAMPLES)]
          .groupby("condition")["sample_id"].apply(list).to_dict())

# All replicates of a condition must share one peak mode (so per-group
# consensus/IDR is well-defined).
for g, members in GROUPS.items():
    modes = {PEAK_MODE[s] for s in members}
    if len(modes) > 1:
        raise ValueError(f"samples.csv: condition '{g}' mixes peak_mode {sorted(modes)}; "
                         "all replicates of a condition must share one peak_mode")

GROUP_MODE = {g: PEAK_MODE[members[0]] for g, members in GROUPS.items()}   # narrow|broad per group

def _repro_method(members):
    n = len(members)
    if n >= 3:
        return "majority"
    if n == 2:
        return "idr"
    return "single"

GROUP_METHOD  = {g: _repro_method(m) for g, m in GROUPS.items()}
IDR_GROUPS    = [g for g in GROUPS if GROUP_METHOD[g] == "idr"]
NONIDR_GROUPS = [g for g in GROUPS if GROUP_METHOD[g] != "idr"]
IDR_SAMPLES   = [s for g in IDR_GROUPS for s in GROUPS[g]]

# IDR pairs (QC): all within-condition replicate pairs among IP samples
IDR_PAIRS = []
for group, members in GROUPS.items():
    for a, b in combinations(members, 2):
        IDR_PAIRS.append((group, a, b))
# IP samples that appear in at least one QC IDR pair (need relaxed peaks for QC)
PAIR_SAMPLES = sorted({s for _, a, b in IDR_PAIRS for s in (a, b)})

# ── Differential-binding contrasts (DESeq2 over the consensus matrix) ─────
# Each contrast compares two conditions (from `condition`); both must have IP
# samples (i.e. appear in GROUPS). Empty list = no differential binding.
CONTRASTS = config.get("contrasts") or []
_bad_contrasts = [(c.get("name", "?"), c.get(k)) for c in CONTRASTS
                  for k in ("condition_a", "condition_b") if c.get(k) not in GROUPS]
if _bad_contrasts:
    raise ValueError("config: contrasts reference unknown condition(s) "
                     f"{_bad_contrasts}; valid conditions with IP samples: {sorted(GROUPS)}")

# Differential binding requires replicates: DESeq2 has no within-group variance
# to estimate dispersion from in a single-replicate design. Skip (don't run) any
# contrast whose conditions do NOT both have >=2 IP replicates.
def _contrast_replicated(c):
    return len(GROUPS[c["condition_a"]]) >= 2 and len(GROUPS[c["condition_b"]]) >= 2

RUNNABLE_CONTRASTS = [c for c in CONTRASTS if _contrast_replicated(c)]
for c in CONTRASTS:
    if not _contrast_replicated(c):
        na, nb = len(GROUPS[c["condition_a"]]), len(GROUPS[c["condition_b"]])
        print(f"[chipseq] skipping differential-binding contrast '{c['name']}' "
              f"({c['condition_a']}={na} rep, {c['condition_b']}={nb} rep): "
              f"DESeq2 needs >=2 replicates per condition.", file=sys.stderr)

CONTRAST_NAMES   = [c["name"] for c in RUNNABLE_CONTRASTS]
CONTRAST_BY_NAME = {c["name"]: c for c in RUNNABLE_CONTRASTS}

# ── Output directories (all relative to the working dir) ────────────────
RESULT_DIR             = "results"
FASTQC_DIR             = f"{RESULT_DIR}/fastqc"
FASTP_DIR              = f"{RESULT_DIR}/fastp"
ALIGN_DIR              = f"{RESULT_DIR}/aligned"
TMP_DIR                = f"{RESULT_DIR}/tmp"
FILTERED_DIR           = f"{RESULT_DIR}/filtered"
DEDUP_DIR              = f"{RESULT_DIR}/dedup"
BLACKLIST_FILTERED_DIR = f"{RESULT_DIR}/blacklist_filtered"
PEAKS_DIR              = f"{RESULT_DIR}/peaks"
QC_DIR                 = f"{RESULT_DIR}/qc"

# Signal tracks
BIGWIG_DIR       = f"{RESULT_DIR}/bigwig"          # RPGC depth-normalized (all samples)
RATIO_BIGWIG_DIR = f"{RESULT_DIR}/ratio_bigwig"    # log2 IP/input (bamCompare)

# Consensus (Module B) directories
RELAXED_PEAKS_DIR  = f"{RESULT_DIR}/peaks_relaxed"
CONSENSUS_DIR      = f"{RESULT_DIR}/consensus"

# QC-pipeline directories (aliases + QC-only outputs)
RMD_BAM_DIR    = BLACKLIST_FILTERED_DIR   # QC alias: results/blacklist_filtered
PEAK_DIR       = PEAKS_DIR                 # QC alias: results/peaks
BEDGRAPH_DIR   = f"{RESULT_DIR}/bedgraph"
DEEPTOOLS_DIR  = f"{RESULT_DIR}/deeptools"
FRIP_DIR       = f"{RESULT_DIR}/FRiP"
IDR_DIR        = f"{RESULT_DIR}/idr"
RELAXED_DIR    = f"{RESULT_DIR}/qc_relaxed_peaks"
COMPLEXITY_DIR = f"{RESULT_DIR}/library_complexity"
ANNOT_DIR      = f"{RESULT_DIR}/peak_annotation"
CROSSCORR_DIR  = f"{RESULT_DIR}/qc_crosscorr"          # ENCODE NSC/RSC (phantompeakqualtools)
JSD_DIR        = f"{RESULT_DIR}/qc_fingerprint"        # deepTools fingerprint JSD metrics
REPRO_DIR      = f"{RESULT_DIR}/idr_reproducibility"   # ENCODE self/pooled pseudo-replicate IDR

# Downstream-analysis directories (workflow/rules/downstream.smk)
ANNOTATION_DIR = f"{RESULT_DIR}/annotation"            # ChIPseeker peak annotation + GO
MOTIF_DIR      = f"{RESULT_DIR}/motifs"                # HOMER motif enrichment
DIFFBIND_DIR   = f"{RESULT_DIR}/diff_binding"          # DESeq2 differential binding
OVERLAP_DIR    = f"{RESULT_DIR}/peak_overlap"          # peak-set Jaccard / overlap

# ── Reference data / config ─────────────────────────────────────────────
GENOME_2BIT  = os.path.join("ref", "hg38.2bit")   # QC: computeGCBias --genome
GTF_FILE     = config["gtf"]
PROMOTER_BED = config["promoter_bed"]
ENHANCER_BED = config["enhancer_bed"]
EGS          = config["effective_genome_size"]
MACS2_GENOME = config["macs2_genome"]
MACS2_QVALUE = config["macs2_qvalue"]
BROAD_CUTOFF = config["broad_cutoff"]

# ── Helpers ─────────────────────────────────────────────────────────────
def _alt(names):
    """Regex alternation for wildcard_constraints; matches nothing if empty."""
    return "|".join(re.escape(n) for n in names) if names else "a^"

def peak_ext(sample):
    """narrowPeak / broadPeak file extension for an IP sample's mode."""
    return "narrowPeak" if PEAK_MODE[sample] == "narrow" else "broadPeak"

def peak_file(sample):
    """Final MACS2 peak file for an IP sample (extension follows its peak_mode)."""
    return f"{PEAKS_DIR}/{sample}_peaks.{peak_ext(sample)}"

def all_peak_files():
    return [peak_file(s) for s in IP_SAMPLES]

def relaxed_peak_file(sample):
    """Module-B relaxed peak file (IDR input) for an IP sample."""
    return f"{RELAXED_PEAKS_DIR}/{sample}_relaxed.{peak_ext(sample)}"

def group_ext(group):
    return "narrowPeak" if GROUP_MODE[group] == "narrow" else "broadPeak"

def idr_peak_file(group):
    """Module-B IDR consensus peaks for a 2-replicate group."""
    return f"{CONSENSUS_DIR}/idr/{group}.idr_peaks.{group_ext(group)}"

def _group_relaxed_inputs(wildcards):
    return [relaxed_peak_file(s) for s in GROUPS[wildcards.group]]

def control_bam(wildcards):
    """Resolved control BAM (Input or IgG per control_type) for an IP sample (for
    the DAG); [] if the sample has no control."""
    c = CONTROL.get(wildcards.sample)
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else []

def control_arg(wildcards):
    """MACS2 '-c <control>.nobl.bam' string for an IP sample; '' if no control."""
    c = CONTROL.get(wildcards.sample)
    return f"-c {BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else ""

def ratio_input_bam(wildcards):
    """Resolved control BAM for the IP-over-control ratio track (RATIO_SAMPLES only)."""
    c = CONTROL[wildcards.sample]
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam"

# fastp adapter handling: AUTO-DETECT adapters for paired-end reads by default
# (--detect_adapter_for_pe). If adapter sequences are provided in config
# (adapter_r1 / adapter_r2, non-empty), pass them explicitly instead — that
# OVERRIDES auto-detection. Leave them unset/empty to auto-detect.
def _fastp_adapter_args():
    r1 = str(config.get("adapter_r1") or "").strip()
    r2 = str(config.get("adapter_r2") or "").strip()
    if r1:
        args = f"--adapter_sequence {r1}"
        if r2:
            args += f" --adapter_sequence_r2 {r2}"
        return args
    return "--detect_adapter_for_pe"

FASTP_ADAPTER_ARGS = _fastp_adapter_args()

# ── ENCODE pseudo-replicate reproducibility (2-replicate IDR conditions) ──
# For each 2-rep condition we assess reproducibility with self-pseudoreplicates
# (each replicate split in half → N1, N2) and pooled pseudoreplicates (both
# replicates pooled, then split → Np), plus the true-replicate IDR count (Nt),
# following the ENCODE self-consistency / rescue-ratio scheme. A "unit" is one
# thing that gets split into two pseudo-halves: a single replicate (self) or a
# pooled condition.
REPRO_UNITS = [f"self__{s}" for s in IDR_SAMPLES] + [f"pool__{g}" for g in IDR_GROUPS]

def _unit_group(unit):
    return CONDITION[unit[len("self__"):]] if unit.startswith("self__") else unit[len("pool__"):]

def _unit_mode(unit):
    return GROUP_MODE[_unit_group(unit)]

def _unit_rep_sample(unit):
    """Representative sample whose control the pseudo-rep peak calls reuse."""
    return unit[len("self__"):] if unit.startswith("self__") else GROUPS[_unit_group(unit)][0]

REPRO_NARROW_UNITS = [u for u in REPRO_UNITS if _unit_mode(u) == "narrow"]
REPRO_BROAD_UNITS  = [u for u in REPRO_UNITS if _unit_mode(u) == "broad"]

def pseudo_source_bam(wildcards):
    """The BAM a pseudo-rep unit is split from (a replicate's, or the pooled BAM)."""
    u = wildcards.unit
    if u.startswith("self__"):
        return f"{BLACKLIST_FILTERED_DIR}/{u[len('self__'):]}.nobl.bam"
    return f"{REPRO_DIR}/pool/{u[len('pool__'):]}.bam"

def unit_control_bam(wildcards):
    c = CONTROL.get(_unit_rep_sample(wildcards.unit))
    return f"{BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else []

def unit_control_arg(wildcards):
    c = CONTROL.get(_unit_rep_sample(wildcards.unit))
    return f"-c {BLACKLIST_FILTERED_DIR}/{c}.nobl.bam" if c else ""
