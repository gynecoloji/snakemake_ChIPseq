# Downstream ChIP-seq analysis: peak annotation + GO (ChIPseeker/clusterProfiler),
# motif enrichment (HOMER), differential binding (DESeq2 over the consensus count
# matrix), peak-set overlap (bedtools Jaccard), and signal heatmaps / metagene
# profiles (deepTools). Consumes the primary pipeline's peaks, consensus set,
# count matrix and bigWigs. Shared names live in common.smk.
#
# Run alone with:  snakemake --use-conda --cores N downstream_all

def _overlap_targets():
    """Peak-set overlap needs >=2 IP samples to be meaningful."""
    if len(IP_SAMPLES) >= 2:
        return [os.path.join(OVERLAP_DIR, "jaccard_matrix.tsv"),
                os.path.join(OVERLAP_DIR, "jaccard_heatmap.png")]
    return []


rule downstream_all:
    input:
        # Peak annotation + GO (per IP sample)
        expand(os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.annotation.tsv"), sample=IP_SAMPLES),
        expand(os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.GO.tsv"), sample=IP_SAMPLES),
        # Motif enrichment (per IP sample)
        expand(os.path.join(MOTIF_DIR, "{sample}", "knownResults.txt"), sample=IP_SAMPLES),
        # Differential binding (per contrast; empty when no contrasts configured)
        expand(os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.deseq2.tsv"), contrast=CONTRAST_NAMES),
        # Peak-set overlap (Jaccard)
        _overlap_targets(),
        # Signal heatmaps + metagene
        os.path.join(DEEPTOOLS_DIR, "Heatmap_peaks.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_peaks.png"),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_genebody.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_genebody.png")


# ── Peak annotation + GO (ChIPseeker + clusterProfiler) ──────────────────
rule annotate_peaks:
    wildcard_constraints:
        sample = _alt(IP_SAMPLES)
    input:
        peaks = lambda w: peak_file(w.sample)
    output:
        anno = os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.annotation.tsv"),
        featdist = os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.feature_distribution.png"),
        disttss = os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.dist_to_tss.png"),
        go = os.path.join(ANNOTATION_DIR, "{sample}", "{sample}.GO.tsv")
    params:
        name = "{sample}"
    conda:
        "../envs/chipseeker.yaml"
    log:
        "logs/annotate_peaks/{sample}.log"
    script:
        "../scripts/peak_annotation.R"


# ── Motif enrichment (HOMER) ─────────────────────────────────────────────
rule motif_enrichment:
    wildcard_constraints:
        sample = _alt(IP_SAMPLES)
    input:
        peaks = lambda w: peak_file(w.sample),
        genome = config["human_fasta"]
    output:
        known = os.path.join(MOTIF_DIR, "{sample}", "knownResults.txt"),
        html = os.path.join(MOTIF_DIR, "{sample}", "homerResults.html")
    params:
        outdir = os.path.join(MOTIF_DIR, "{sample}"),
        size = lambda w: "200" if PEAK_MODE[w.sample] == "narrow" else "given",
        bed = os.path.join(TMP_DIR, "{sample}.homer.bed")
    threads: 8
    conda:
        "../envs/homer.yaml"
    log:
        "logs/motif_enrichment/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} {TMP_DIR} logs/motif_enrichment
        # HOMER wants a peak/BED file (chrom start end name score strand); MACS2
        # narrowPeak/broadPeak already carry these in columns 1-6.
        cut -f1-6 {input.peaks} > {params.bed}
        findMotifsGenome.pl {params.bed} {input.genome} {params.outdir} \
            -size {params.size} -p {threads} \
            -preparsedDir {params.outdir}/preparsed > {log} 2>&1
        rm -f {params.bed}
        """


# ── Differential binding (DESeq2 over the consensus count matrix) ─────────
rule differential_binding:
    wildcard_constraints:
        contrast = _alt(CONTRAST_NAMES)
    input:
        counts = f"{CONSENSUS_DIR}/consensus_counts.txt",
        peaks = f"{CONSENSUS_DIR}/consensus_peaks.bed"
    output:
        results = os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.deseq2.tsv"),
        sig = os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.significant.tsv"),
        ma = os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.MA.png"),
        volcano = os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.volcano.png"),
        pca = os.path.join(DIFFBIND_DIR, "{contrast}", "{contrast}.PCA.png")
    params:
        name = "{contrast}",
        cond_a = lambda w: CONTRAST_BY_NAME[w.contrast]["condition_a"],
        cond_b = lambda w: CONTRAST_BY_NAME[w.contrast]["condition_b"],
        bams_a = lambda w: [f"{BLACKLIST_FILTERED_DIR}/{s}.nobl.bam"
                            for s in GROUPS[CONTRAST_BY_NAME[w.contrast]["condition_a"]]],
        bams_b = lambda w: [f"{BLACKLIST_FILTERED_DIR}/{s}.nobl.bam"
                            for s in GROUPS[CONTRAST_BY_NAME[w.contrast]["condition_b"]]]
    conda:
        "../envs/chipseeker.yaml"
    log:
        "logs/differential_binding/{contrast}.log"
    script:
        "../scripts/differential_binding.R"


# ── Peak-set overlap (bedtools Jaccard, pairwise over IP samples) ─────────
rule peak_jaccard:
    input:
        peaks = all_peak_files()
    output:
        long = os.path.join(OVERLAP_DIR, "jaccard_long.tsv")
    params:
        samples = IP_SAMPLES,
        peakdir = PEAKS_DIR
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/peak_overlap/jaccard.log"
    shell:
        """
        mkdir -p {OVERLAP_DIR} {TMP_DIR} logs/peak_overlap
        for s in {params.samples}; do
            pk=$(ls {params.peakdir}/${{s}}_peaks.narrowPeak {params.peakdir}/${{s}}_peaks.broadPeak 2>/dev/null | head -1)
            sort -k1,1 -k2,2n "$pk" > {TMP_DIR}/${{s}}.ovl.sorted.bed
        done
        echo -e "sample_a\tsample_b\tjaccard" > {output.long}
        for a in {params.samples}; do
            for b in {params.samples}; do
                j=$(bedtools jaccard -a {TMP_DIR}/$a.ovl.sorted.bed -b {TMP_DIR}/$b.ovl.sorted.bed \
                    | awk 'NR==2{{print $3}}')
                [ -z "$j" ] && j=0
                echo -e "$a\t$b\t$j" >> {output.long}
            done
        done 2> {log}
        for s in {params.samples}; do rm -f {TMP_DIR}/${{s}}.ovl.sorted.bed; done
        """


rule peak_overlap_matrix:
    input:
        long = os.path.join(OVERLAP_DIR, "jaccard_long.tsv")
    output:
        matrix = os.path.join(OVERLAP_DIR, "jaccard_matrix.tsv"),
        plot = os.path.join(OVERLAP_DIR, "jaccard_heatmap.png")
    conda:
        "../envs/deeptools.yaml"   # provides matplotlib
    log:
        "logs/peak_overlap/matrix.log"
    script:
        "../scripts/peak_overlap.py"


# ── Signal heatmap over the consensus peak set (peak-centered) ───────────
rule deeptools_peak_heatmap:
    input:
        bigwigs = expand(os.path.join(BIGWIG_DIR, "{sample}.bw"), sample=SAMPLES),
        peaks = f"{CONSENSUS_DIR}/consensus_peaks.bed"
    output:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix_peaks.mat.gz"),
        heatmap = os.path.join(DEEPTOOLS_DIR, "Heatmap_peaks.png"),
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_peaks.png")
    threads: 16
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_peak_heatmap/peaks.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_peak_heatmap
        computeMatrix reference-point -p {threads} \
            --referencePoint center \
            -S {input.bigwigs} -R {input.peaks} \
            -a 3000 -b 3000 --skipZeros \
            -o {output.matrix} 2> {log}
        plotHeatmap -m {output.matrix} --dpi 300 --heatmapWidth 20 \
            -out {output.heatmap} --plotFileFormat png --sortUsing mean 2>> {log}
        plotProfile -m {output.matrix} --dpi 300 \
            -out {output.profile} --plotFileFormat png 2>> {log}
        """


# ── Metagene profile / heatmap over gene bodies (scale-regions) ──────────
rule deeptools_metagene:
    input:
        bigwigs = expand(os.path.join(BIGWIG_DIR, "{sample}.bw"), sample=SAMPLES),
        gtf = GTF_FILE
    output:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix_genebody.mat.gz"),
        heatmap = os.path.join(DEEPTOOLS_DIR, "Heatmap_genebody.png"),
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_genebody.png")
    threads: 16
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_metagene/genebody.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_metagene
        computeMatrix scale-regions -p {threads} \
            -S {input.bigwigs} -R {input.gtf} \
            --regionBodyLength 5000 -a 2000 -b 2000 \
            --skipZeros --metagene \
            -o {output.matrix} 2> {log}
        plotHeatmap -m {output.matrix} --dpi 300 --heatmapWidth 20 \
            -out {output.heatmap} --plotFileFormat png --sortUsing mean 2>> {log}
        plotProfile -m {output.matrix} --dpi 300 \
            -out {output.profile} --plotFileFormat png 2>> {log}
        """
