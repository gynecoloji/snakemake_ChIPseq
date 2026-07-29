# ChIP-seq QC pipeline: deepTools coverage/correlation/fingerprint/GC/TSS QC,
# FRiP, IDR on relaxed peaks, library complexity, reads-in-annotation, peak
# summaries, a FastQC-only MultiQC report, and a self-contained interactive HTML
# QC report.
#
# Consumes the primary pipeline's results/ outputs (blacklist_filtered/, dedup/,
# filtered/, peaks/, bigwig/). Shared config, samples, directory constants and
# helpers live in common.smk (included first by workflow/Snakefile).

localrules: multiqc_fastqc


# Aggregate target for the QC pipeline. Run it alone (after the primary pipeline
# outputs exist) with:  snakemake --use-conda --cores N qc_all
rule qc_all:
    input:
        # deepTools coverage + QC
        expand(os.path.join(BEDGRAPH_DIR, "{sample}.nobl.RPGC.bedgraph"), sample=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        os.path.join(DEEPTOOLS_DIR, "fragmentsize.txt"),
        os.path.join(DEEPTOOLS_DIR, "chipseq_fingerprint.png"),
        os.path.join(DEEPTOOLS_DIR, "chipseq_fingerprint.tab"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        expand(os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.png"), sample=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "fragment_lengths.txt"),
        os.path.join(DEEPTOOLS_DIR, "correlation_matrix.tab"),
        os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json"),
        # TSS: heatmap/profile + numeric enrichment score
        os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        os.path.join(QC_DIR, "tss_enrichment_mqc.txt"),
        # FRiP (IP samples) + IDR (relaxed) + library complexity
        [os.path.join(FRIP_DIR, f"{s}.frip.txt") for s in IP_SAMPLES],
        [expand(os.path.join(IDR_DIR, "{group}--{rep1}--{rep2}--idr_peaks.{condition}.txt"),
                condition=[group_ext(g)], group=[g], rep1=[r1], rep2=[r2]) for g, r1, r2 in IDR_PAIRS],
        expand(os.path.join(COMPLEXITY_DIR, "{sample}_complexity.txt"), sample=SAMPLES),
        # peak + annotation summary
        os.path.join(QC_DIR, "peak_summary_mqc.txt"),
        os.path.join(ANNOT_DIR, "reads_in_annotations_mqc.txt"),
        # ENCODE QC: strand cross-correlation, fingerprint JSD, IDR reproducibility
        os.path.join(QC_DIR, "cross_correlation_mqc.txt"),
        os.path.join(QC_DIR, "fingerprint_jsd_mqc.txt"),
        os.path.join(QC_DIR, "idr_reproducibility_mqc.txt"),
        # FastQC-only MultiQC report
        os.path.join(QC_DIR, "multiqc_fastqc.html"),
        # Interactive HTML QC report
        os.path.join(QC_DIR, "chipseq_qc_report.html")


# 1. RPGC bedgraph per sample (coverage QC; main pipeline makes bigWigs, not bedgraphs)
rule deeptools_bedgraph:
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        bedgraph = os.path.join(BEDGRAPH_DIR, "{sample}.nobl.RPGC.bedgraph")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_bedgraph/{sample}.log"
    shell:
        """
        mkdir -p {BEDGRAPH_DIR} logs/deeptools_bedgraph
        bamCoverage -p {threads} \
            --outFileFormat bedgraph \
            --effectiveGenomeSize {EGS} \
            --normalizeUsing RPGC \
            --binSize 10 --extendReads \
            --bam {input.bam} -o {output.bedgraph} 2> {log}
        """


# 2. Fragment size distribution
rule deeptools_fragmentsize:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        table = os.path.join(DEEPTOOLS_DIR, "fragmentsize.txt"),
        raw = os.path.join(DEEPTOOLS_DIR, "fragment_lengths.txt")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_fragmentsize/fragmentsize.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_fragmentsize
        bamPEFragmentSize -p {threads} \
            -hist {output.plot} \
            -T "Fragment size of PE ChIP-seq data" \
            --maxFragmentLength 1500 \
            -b {input.bams} \
            --outRawFragmentLengths {output.raw} \
            --table {output.table} 2> {log}
        """


# 3. Fingerprint (signal-to-noise; the canonical ChIP-seq IP-vs-input QC)
rule deeptools_plotfingerprint:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "chipseq_fingerprint.png"),
        table = os.path.join(DEEPTOOLS_DIR, "chipseq_fingerprint.tab")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_plotfingerprint/fingerprint.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_plotfingerprint
        plotFingerprint -p {threads} \
            -b {input.bams} \
            --ignoreDuplicates \
            -T "Fingerprints" \
            --skipZeros \
            --plotFileFormat png \
            -plot {output.plot} \
            --outRawCounts {output.table} 2> {log}
        """


# 4. Correlation / PCA across samples
rule deeptools_cor_multibam:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES)
    output:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz"),
        counts = os.path.join(DEEPTOOLS_DIR, "deeptools_readCounts.tab")
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/multibam.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_correlation
        multiBamSummary bins \
            -bs 5000 \
            --ignoreDuplicates \
            -p {threads} \
            --bamfiles {input.bams} \
            -out {output.npz} \
            --outRawCounts {output.counts} 2> {log}
        """


rule deeptools_cor_scatterplot:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/scatterplot.log"
    shell:
        """
        plotCorrelation --corData {input.npz} \
            --whatToPlot scatterplot \
            --skipZero \
            --plotTitle "Scatterplot" \
            --plotFileFormat png \
            --corMethod spearman \
            --log1p \
            --plotFile {output.plot} 2> {log}
        """


rule deeptools_cor_heatmap:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        cormat = os.path.join(DEEPTOOLS_DIR, "correlation_matrix.tab")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/heatmap.log"
    shell:
        """
        plotCorrelation --corData {input.npz} \
            --whatToPlot heatmap \
            --skipZero \
            --plotTitle "Heatmap" \
            --plotFileFormat png \
            --corMethod spearman \
            --log1p \
            --outFileCorMatrix {output.cormat} \
            --plotFile {output.plot} 2> {log}
        """


rule deeptools_cor_pca:
    input:
        npz = os.path.join(DEEPTOOLS_DIR, "deeptools_multiBAM.out.npz")
    output:
        plot = os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        data = os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.tab")
    threads: 4
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation/pca.log"
    shell:
        """
        plotPCA --corData {input.npz} \
            --plotTitle "PCA" \
            --plotFileFormat png \
            --ntop 1000 \
            --plotFile {output.plot} \
            --outFileNameData {output.data} 2> {log}
        """


# 5. GC bias
rule deeptools_gc_bias:
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        genome = GENOME_2BIT
    output:
        freq = os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.txt"),
        plot = os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.png")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_gc_bias/{sample}.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_gc_bias
        computeGCBias -b {input.bam} \
            --effectiveGenomeSize {EGS} \
            -p {threads} \
            --genome {input.genome} \
            -freq {output.freq} \
            --biasPlot {output.plot} \
            --plotFileFormat png 2> {log}
        """


# 6. TSS signal: heatmap + profile + profile-data table (reuses the main
#    pipeline's RPGC bigWigs instead of regenerating them).
rule deeptools_tss:
    input:
        bigwigs = expand(os.path.join(BIGWIG_DIR, "{sample}.bw"), sample=SAMPLES),
        gtf = GTF_FILE
    output:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix.mat.gz"),
        heatmap = os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        profiledata = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.data.tab")
    threads: 16
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_tss/tss.log"
    shell:
        """
        mkdir -p {DEEPTOOLS_DIR} logs/deeptools_tss
        computeMatrix reference-point \
            -p {threads} \
            --referencePoint TSS \
            -S {input.bigwigs} \
            -R {input.gtf} \
            -a 2000 -b 2000 \
            --skipZeros \
            -o {output.matrix} 2> {log}
        plotHeatmap \
            -m {output.matrix} \
            --dpi 300 \
            --zMin -3 --zMax 3 \
            --heatmapWidth 20 \
            -out {output.heatmap} \
            --plotFileFormat png \
            --sortUsing mean 2>> {log}
        plotProfile \
            -m {output.matrix} \
            --dpi 300 \
            -out {output.profile} \
            --plotFileFormat png \
            --outFileNameData {output.profiledata} 2>> {log}
        """


# 6b. Downsampled TSS heatmap matrix (JSON) for the interactive report's canvas
#     heatmap; reuses the existing matrix.mat.gz (no computeMatrix rerun).
rule deeptools_tss_heatmap_downsample:
    input:
        matrix = os.path.join(DEEPTOOLS_DIR, "matrix.mat.gz")
    output:
        json = os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json")
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/deeptools_tss/downsample.log"
    shell:
        """
        python workflow/scripts/downsample_tss_matrix.py {input.matrix} \
            -o {output.json} --nrows 180 --ncols 80 > {log} 2>&1
        """


# 7. Numeric TSS enrichment score per sample (from the profile-data table)
rule tss_enrichment_score:
    input:
        profile = os.path.join(DEEPTOOLS_DIR, "Profile_TSS.data.tab")
    output:
        tsv = os.path.join(QC_DIR, "tss_enrichment_scores.tsv"),
        mqc = os.path.join(QC_DIR, "tss_enrichment_mqc.txt")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/tss_enrichment_score/tss.log"
    script:
        "../scripts/tss_score.py"


# 8. FRiP (fraction of reads in peaks) per IP sample (peak file follows peak_mode)
rule FRiP:
    wildcard_constraints:
        sample = _alt(IP_SAMPLES)
    input:
        bamfile = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        peakfile = lambda w: peak_file(w.sample)
    output:
        fripfile = os.path.join(FRIP_DIR, "{sample}.frip.txt")
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/FRiP/{sample}.log"
    shell:
        """
        mkdir -p {FRIP_DIR} logs/FRiP
        total=$(samtools view -c {input.bamfile})
        in_peaks=$(bedtools intersect -u -abam {input.bamfile} -b {input.peakfile} | samtools view -c) 2> {log}
        frip=$(echo "scale=4; $in_peaks / $total" | bc)
        echo -e "{wildcards.sample}\t$in_peaks\t$total\t$frip" > {output.fripfile}
        """


# 9. QC-side relaxed MACS2 calls for IDR (IDR expects relaxed peaks, not -q0.05)
rule qc_relaxed_peaks_narrow:
    wildcard_constraints:
        sample = _alt([s for s in PAIR_SAMPLES if PEAK_MODE[s] == "narrow"])
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        control = control_bam
    output:
        peaks = os.path.join(RELAXED_DIR, "{sample}_relaxed.narrowPeak")
    params:
        outdir = RELAXED_DIR,
        name = "{sample}",
        genome = MACS2_GENOME,
        pvalue = config["idr_relaxed_pvalue"],
        top_n = config["idr_top_n_peaks"],
        control_arg = control_arg
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/qc_relaxed_peaks/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/qc_relaxed_peaks
        macs2 callpeak -t {input.bam} {params.control_arg} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_summits.bed \
              {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        """


rule qc_relaxed_peaks_broad:
    wildcard_constraints:
        sample = _alt([s for s in PAIR_SAMPLES if PEAK_MODE[s] == "broad"])
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        control = control_bam
    output:
        peaks = os.path.join(RELAXED_DIR, "{sample}_relaxed.broadPeak")
    params:
        outdir = RELAXED_DIR,
        name = "{sample}",
        genome = MACS2_GENOME,
        pvalue = config["idr_relaxed_pvalue"],
        top_n = config["idr_top_n_peaks"],
        broad_cutoff = BROAD_CUTOFF,
        control_arg = control_arg
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/qc_relaxed_peaks/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/qc_relaxed_peaks
        macs2 callpeak -t {input.bam} {params.control_arg} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel --broad --broad-cutoff {params.broad_cutoff} -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_peaks.gappedPeak \
              {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        """


# 10. IDR on within-condition replicate pairs (relaxed peaks)
rule idr:
    input:
        rep1 = lambda w: os.path.join(RELAXED_DIR, f"{w.rep1}_relaxed.{w.condition}"),
        rep2 = lambda w: os.path.join(RELAXED_DIR, f"{w.rep2}_relaxed.{w.condition}")
    output:
        peaks = os.path.join(IDR_DIR, "{group}--{rep1}--{rep2}--idr_peaks.{condition}.txt")
    conda:
        "../envs/idr.yaml"
    wildcard_constraints:
        condition = "narrowPeak|broadPeak"
    log:
        "logs/idr/{group}--{rep1}--{rep2}--idr_{condition}.log"
    shell:
        """
        mkdir -p {IDR_DIR} logs/idr
        idr --samples {input.rep1} {input.rep2} \
            --input-file-type {wildcards.condition} \
            --rank p.value \
            --output-file {output.peaks} \
            --plot \
            --log-output-file {log}
        """


# 11. Library complexity (NRF / PBC1 / PBC2) on the pre-dedup filtered BAM
rule calculate_library_complexity:
    input:
        bam = os.path.join(FILTERED_DIR, "{sample}.sorted.filtered.bam")
    output:
        txt = os.path.join(COMPLEXITY_DIR, "{sample}_complexity.txt")
    conda:
        "../envs/bedtools.yaml"  # needs samtools, bedtools, bc
    threads: 8
    log:
        "logs/library_complexity/{sample}.log"
    params:
        temp_dir = os.path.join(TMP_DIR, "{sample}_complexity"),
        tmp = os.path.join(TMP_DIR, "{sample}_complexity.bed")
    shell:
        """
        mkdir -p {COMPLEXITY_DIR} logs/library_complexity {params.temp_dir}

        # Extract fragments (BEDPE) from the name-sorted BAM
        echo "Extracting fragments from BAM file..." >> {log}
        samtools sort -n -@ {threads} {input.bam} | \
        bedtools bamtobed -bedpe -i stdin > {params.tmp} 2>> {log}
        # position-only key (chrom, start, end) — NO read name: PCR duplicates share
        # coordinates but have different names, so the name must be dropped for
        # `uniq -c` to collapse them into one location with the right multiplicity.
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1==$4) {{
                start = ($2 < $5) ? $2 : $5;
                end = ($3 > $6) ? $3 : $6;
                print $1, start, end;
            }}
        }}' {params.tmp} | \
        sort -k1,1 -k2,2n -k3,3n > {params.temp_dir}/fragments.bed

        fragment_count=$(wc -l < {params.temp_dir}/fragments.bed)
        echo "Total extracted fragments: $fragment_count" >> {log}

        # Fragment counts by genomic location (for PCR-duplicate stats)
        sort -k1,1 -k2,2n -k3,3n {params.temp_dir}/fragments.bed | \
        uniq -c > {params.temp_dir}/fragment_counts.txt

        unique=$(wc -l < {params.temp_dir}/fragment_counts.txt)
        one_read=$(awk '$1 == 1' {params.temp_dir}/fragment_counts.txt | wc -l)
        two_reads=$(awk '$1 == 2' {params.temp_dir}/fragment_counts.txt | wc -l)
        total_reads=$(samtools view -c {input.bam})

        echo "Unique locations: $unique" >> {log}
        echo "Locations with exactly one fragment: $one_read" >> {log}
        echo "Locations with exactly two fragments: $two_reads" >> {log}
        echo "Total mapped reads: $total_reads" >> {log}

        if [ "$one_read" -eq 0 ]; then
            echo "WARNING: No unique fragments found, setting to 1 to prevent division by zero" >> {log}
            one_read=1
        fi

        nrf=$(echo "scale=6; $unique / $fragment_count" | bc)
        pbc1=$(echo "scale=6; $one_read / $unique" | bc)
        if [ "$two_reads" -eq 0 ]; then
            two_reads=1
            echo "CRITICAL: two_reads is 0 before PBC2 calculation, forcing to 1" >> {log}
        fi
        pbc2=$(echo "scale=6; $one_read / $two_reads" | bc)

        echo "## Library Complexity Metrics for {wildcards.sample} ##" > {output.txt}
        echo -e "Total Reads\t$total_reads" >> {output.txt}
        echo -e "Total Fragments\t$fragment_count" >> {output.txt}
        echo -e "Distinct Fragment Locations (Nd)\t$unique" >> {output.txt}
        echo -e "Locations with 1 Fragment (N1)\t$one_read" >> {output.txt}
        echo -e "Locations with 2 Fragments (N2)\t$two_reads" >> {output.txt}
        echo -e "NRF (Nd/Total)\t$nrf" >> {output.txt}
        echo -e "PBC1 (N1/Nd)\t$pbc1" >> {output.txt}
        echo -e "PBC2 (N1/N2)\t$pbc2" >> {output.txt}

        echo -e "\n## Quality Assessment ##" >> {output.txt}
        if (( $(echo "$nrf > 0.9" | bc -l) )); then echo -e "NRF: $nrf - High complexity (>0.9)" >> {output.txt}
        elif (( $(echo "$nrf > 0.8" | bc -l) )); then echo -e "NRF: $nrf - Good complexity (0.8-0.9)" >> {output.txt}
        elif (( $(echo "$nrf > 0.7" | bc -l) )); then echo -e "NRF: $nrf - Moderate complexity (0.7-0.8)" >> {output.txt}
        else echo -e "NRF: $nrf - Low complexity (<0.7)" >> {output.txt}; fi

        if (( $(echo "$pbc1 > 0.9" | bc -l) )); then echo -e "PBC1: $pbc1 - Near ideal (>0.9)" >> {output.txt}
        elif (( $(echo "$pbc1 > 0.8" | bc -l) )); then echo -e "PBC1: $pbc1 - Good (0.8-0.9)" >> {output.txt}
        elif (( $(echo "$pbc1 > 0.7" | bc -l) )); then echo -e "PBC1: $pbc1 - Moderate (0.7-0.8)" >> {output.txt}
        else echo -e "PBC1: $pbc1 - Severe bottlenecking (<0.7)" >> {output.txt}; fi

        if (( $(echo "$pbc2 > 10" | bc -l) )); then echo -e "PBC2: $pbc2 - Near ideal (>10)" >> {output.txt}
        elif (( $(echo "$pbc2 > 3" | bc -l) )); then echo -e "PBC2: $pbc2 - Good (3-10)" >> {output.txt}
        elif (( $(echo "$pbc2 > 1" | bc -l) )); then echo -e "PBC2: $pbc2 - Moderate (1-3)" >> {output.txt}
        else echo -e "PBC2: $pbc2 - Severe bottlenecking (<1)" >> {output.txt}; fi

        rm -rf {params.temp_dir} {params.tmp}
        echo "Library complexity calculation completed for {wildcards.sample}" >> {log}
        """


# 12. Reads in promoters vs enhancers (signal distribution; all samples)
rule reads_in_annotations:
    input:
        bams = expand(os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"), sample=SAMPLES),
        promoter = PROMOTER_BED,
        enhancer = ENHANCER_BED
    output:
        tsv = os.path.join(ANNOT_DIR, "reads_in_annotations.tsv"),
        mqc = os.path.join(ANNOT_DIR, "reads_in_annotations_mqc.txt")
    params:
        samples = SAMPLES,
        bamdir = RMD_BAM_DIR
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/reads_in_annotations/annot.log"
    shell:
        """
        mkdir -p {ANNOT_DIR} logs/reads_in_annotations
        printf "# id: reads_in_annotations\n# section_name: 'Reads in annotations'\n# description: 'Fraction of reads overlapping promoters / enhancers.'\n# plot_type: 'table'\nSample\tTotal\tIn promoter\tIn enhancer\tPromoter frac\tEnhancer frac\n" > {output.mqc}
        echo -e "sample\ttotal\tin_promoter\tin_enhancer\tpromoter_frac\tenhancer_frac" > {output.tsv}
        for s in {params.samples}; do
            bam={params.bamdir}/$s.nobl.bam
            total=$(samtools view -c $bam)
            prom=$(bedtools intersect -u -abam $bam -b {input.promoter} | samtools view -c)
            enh=$(bedtools intersect -u -abam $bam -b {input.enhancer} | samtools view -c)
            awk -v s=$s -v t=$total -v p=$prom -v e=$enh 'BEGIN{{
                pf=(t>0)?p/t:0; ef=(t>0)?e/t:0;
                printf "%s\\t%d\\t%d\\t%d\\t%.4f\\t%.4f\\n", s, t, p, e, pf, ef
            }}' | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# 13. Peak count / width summary + FRiP, per IP sample (narrow or broad)
rule peak_summary:
    input:
        peaks = all_peak_files(),
        frips = [os.path.join(FRIP_DIR, f"{s}.frip.txt") for s in IP_SAMPLES]
    output:
        tsv = os.path.join(QC_DIR, "peak_summary.tsv"),
        mqc = os.path.join(QC_DIR, "peak_summary_mqc.txt")
    params:
        samples = IP_SAMPLES,
        peakdir = PEAK_DIR,
        fripdir = FRIP_DIR
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/peak_summary/peak_summary.log"
    shell:
        """
        mkdir -p {QC_DIR} logs/peak_summary
        printf "# id: peak_summary\n# section_name: 'Peak summary'\n# description: 'Peak mode, count, width stats, and FRiP per IP sample.'\n# plot_type: 'table'\nSample\tMode\tPeaks\tMean width\tMin width\tMax width\tFRiP\n" > {output.mqc}
        echo -e "sample\tmode\tn_peaks\tmean_width\tmin_width\tmax_width\tFRiP" > {output.tsv}
        for s in {params.samples}; do
            pk=$(ls {params.peakdir}/${{s}}_peaks.narrowPeak {params.peakdir}/${{s}}_peaks.broadPeak 2>/dev/null | head -1)
            mode=narrow; [ "${{pk##*.}}" = "broadPeak" ] && mode=broad
            frip=$(cut -f4 {params.fripdir}/$s.frip.txt)
            awk -v s=$s -v mode=$mode -v frip=$frip 'BEGIN{{mn=""}} {{
                w=$3-$2; n++; sum+=w;
                if(mn==""||w<mn)mn=w; if(w>mx)mx=w
            }} END{{
                printf "%s\\t%s\\t%d\\t%.1f\\t%d\\t%d\\t%s\\n", s, mode, n, (n>0?sum/n:0), (mn==""?0:mn), (mx==""?0:mx), frip
            }}' "$pk" | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── ENCODE: strand cross-correlation (NSC / RSC / fragment length) ───────
# 14. phantompeakqualtools run_spp.R per IP sample. ENCODE cross-correlation is a
#     single-ended tag metric, so we take the first mate of each pair.
rule cross_correlation:
    wildcard_constraints:
        sample = _alt(IP_SAMPLES)
    input:
        bam = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam")
    output:
        metrics = os.path.join(CROSSCORR_DIR, "{sample}.spp.out"),
        plot = os.path.join(CROSSCORR_DIR, "{sample}.spp.pdf")
    params:
        r1 = os.path.join(TMP_DIR, "{sample}.r1.bam")
    threads: 4
    conda:
        "../envs/phantompeakqualtools.yaml"
    log:
        "logs/cross_correlation/{sample}.log"
    shell:
        """
        mkdir -p {CROSSCORR_DIR} {TMP_DIR} logs/cross_correlation
        rm -f {output.metrics}
        samtools view -b -f 64 {input.bam} > {params.r1} 2> {log}
        Rscript "$(which run_spp.R)" -c={params.r1} -savp={output.plot} \
            -out={output.metrics} -p={threads} >> {log} 2>&1
        rm -f {params.r1}
        """


# 15. Cross-correlation summary (NSC>=1.05, RSC>=0.8 are the ENCODE targets)
rule cross_correlation_summary:
    input:
        expand(os.path.join(CROSSCORR_DIR, "{sample}.spp.out"), sample=IP_SAMPLES)
    output:
        tsv = os.path.join(QC_DIR, "cross_correlation.tsv"),
        mqc = os.path.join(QC_DIR, "cross_correlation_mqc.txt")
    params:
        samples = IP_SAMPLES,
        ccdir = CROSSCORR_DIR
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/cross_correlation/summary.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/cross_correlation
        printf "# id: cross_correlation\n# section_name: 'Strand cross-correlation (ENCODE)'\n# description: 'NSC (>=1.05) and RSC (>=0.8) from phantompeakqualtools; est. fragment length and QualityTag.'\n# plot_type: 'table'\nSample\tFragLen\tNSC\tRSC\tQualityTag\n" > {output.mqc}
        echo -e "sample\test_frag_len\tNSC\tRSC\tquality_tag" > {output.tsv}
        for s in {params.samples}; do
            # run_spp.R cols: 1 file, 2 numReads, 3 estFragLen (may be CSV), 4 corr,
            # 5 phantomPeak, 6 corr, 7 argmin, 8 mincorr, 9 NSC, 10 RSC, 11 QualityTag
            awk -v s=$s 'BEGIN{{FS="\t";OFS="\t"}} {{split($3,a,","); print s,a[1],$9,$10,$11}}' \
                {params.ccdir}/$s.spp.out | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── ENCODE: quantitative fingerprint (Jensen-Shannon distance) ───────────
# 16. deepTools plotFingerprint IP-vs-its-control, per IP sample with a control.
rule fingerprint_jsd:
    wildcard_constraints:
        sample = _alt(RATIO_SAMPLES)
    input:
        ip = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam"),
        ip_bai = os.path.join(RMD_BAM_DIR, "{sample}.nobl.bam.bai"),
        control = ratio_input_bam,
        control_bai = lambda w: ratio_input_bam(w) + ".bai"
    output:
        metrics = os.path.join(JSD_DIR, "{sample}.jsd.txt"),
        plot = os.path.join(JSD_DIR, "{sample}.fingerprint.png")
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/fingerprint_jsd/{sample}.log"
    shell:
        """
        mkdir -p {JSD_DIR} logs/fingerprint_jsd
        plotFingerprint -b {input.ip} {input.control} \
            --labels {wildcards.sample} control \
            --JSDsample {input.control} \
            --ignoreDuplicates --skipZeros \
            -p {threads} \
            --outQualityMetrics {output.metrics} \
            --plotFile {output.plot} > {log} 2>&1
        """


# 17. Fingerprint JSD summary (JS distance + % genome enriched, per IP sample)
rule fingerprint_jsd_summary:
    input:
        expand(os.path.join(JSD_DIR, "{sample}.jsd.txt"), sample=RATIO_SAMPLES)
    output:
        tsv = os.path.join(QC_DIR, "fingerprint_jsd.tsv"),
        mqc = os.path.join(QC_DIR, "fingerprint_jsd_mqc.txt")
    params:
        samples = RATIO_SAMPLES,
        jsddir = JSD_DIR
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/fingerprint_jsd/summary.log"
    shell:
        r"""
        mkdir -p {QC_DIR} logs/fingerprint_jsd
        printf "# id: fingerprint_jsd\n# section_name: 'Fingerprint JSD (ENCODE)'\n# description: 'deepTools plotFingerprint IP-vs-control: JS distance and %% genome enriched.'\n# plot_type: 'table'\nSample\tJS distance\tPct genome enriched\tAUC\n" > {output.mqc}
        echo -e "sample\tjs_distance\tpct_genome_enriched\tauc" > {output.tsv}
        for s in {params.samples}; do
            # outQualityMetrics has a header row; map column names -> index, then
            # print the IP-labeled row (JSDsample=control, so JS Distance is IP-vs-control).
            awk -v s="$s" 'BEGIN{{FS="\t";OFS="\t"}}
                NR==1{{for(i=1;i<=NF;i++){{h=$i; gsub(/^ +| +$/,"",h); col[h]=i}} next}}
                $1==s{{print s, $(col["JS Distance"]), $(col["% genome enriched"]), $(col["AUC"])}}' \
                {params.jsddir}/$s.jsd.txt | tee -a {output.tsv} >> {output.mqc}
        done 2> {log}
        """


# ── ENCODE: pseudo-replicate reproducibility (self-consistency + rescue) ─
# 18. Pool the two replicate BAMs of a 2-rep condition.
rule repro_pool:
    wildcard_constraints:
        group = _alt(IDR_GROUPS)
    input:
        lambda w: [os.path.join(RMD_BAM_DIR, f"{s}.nobl.bam") for s in GROUPS[w.group]]
    output:
        os.path.join(REPRO_DIR, "pool", "{group}.bam")
    threads: 4
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/idr_reproducibility/pool_{group}.log"
    shell:
        """
        mkdir -p {REPRO_DIR}/pool logs/idr_reproducibility
        samtools merge -f -@ {threads} {output} {input} 2> {log}
        """


# 19. Split a source BAM into two disjoint pseudo-replicate halves (mates kept
#     together via a deterministic read-name hash).
rule repro_split:
    wildcard_constraints:
        unit = _alt(REPRO_UNITS)
    input:
        bam = pseudo_source_bam
    output:
        pr1 = os.path.join(REPRO_DIR, "split", "{unit}.pr1.bam"),
        pr2 = os.path.join(REPRO_DIR, "split", "{unit}.pr2.bam")
    params:
        hdr = os.path.join(TMP_DIR, "{unit}.hdr.sam"),
        b1 = os.path.join(TMP_DIR, "{unit}.pr1.body.sam"),
        b2 = os.path.join(TMP_DIR, "{unit}.pr2.body.sam")
    threads: 4
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/idr_reproducibility/split_{unit}.log"
    shell:
        r"""
        mkdir -p {REPRO_DIR}/split {TMP_DIR} logs/idr_reproducibility
        samtools view -H {input.bam} > {params.hdr} 2> {log}
        samtools view {input.bam} | awk -v o1={params.b1} -v o2={params.b2} \
            'BEGIN{{CH="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz:._-/"}}
             {{ n=$1; s=0; L=length(n); for(i=1;i<=L;i++){{p=index(CH,substr(n,i,1)); if(p==0)p=1; s+=p}}
                if(s%2==0) print > o1; else print > o2 }}' 2>> {log}
        cat {params.hdr} {params.b1} | samtools sort -@ {threads} -o {output.pr1} - 2>> {log}
        cat {params.hdr} {params.b2} | samtools sort -@ {threads} -o {output.pr2} - 2>> {log}
        rm -f {params.hdr} {params.b1} {params.b2}
        """


# 20. Relaxed peaks on each pseudo-half (mode-aware), reusing the unit's control.
rule repro_relaxed_narrow:
    wildcard_constraints:
        unit = _alt(REPRO_NARROW_UNITS),
        half = "pr1|pr2"
    input:
        bam = os.path.join(REPRO_DIR, "split", "{unit}.{half}.bam"),
        control = unit_control_bam
    output:
        peaks = os.path.join(REPRO_DIR, "peaks", "{unit}.{half}_relaxed.narrowPeak")
    params:
        outdir = os.path.join(REPRO_DIR, "peaks"),
        name = "{unit}.{half}",
        genome = MACS2_GENOME,
        pvalue = config["idr_relaxed_pvalue"],
        top_n = config["idr_top_n_peaks"],
        control_arg = unit_control_arg
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/idr_reproducibility/relaxed_{unit}.{half}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/idr_reproducibility
        macs2 callpeak -t {input.bam} {params.control_arg} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_summits.bed \
              {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        """


rule repro_relaxed_broad:
    wildcard_constraints:
        unit = _alt(REPRO_BROAD_UNITS),
        half = "pr1|pr2"
    input:
        bam = os.path.join(REPRO_DIR, "split", "{unit}.{half}.bam"),
        control = unit_control_bam
    output:
        peaks = os.path.join(REPRO_DIR, "peaks", "{unit}.{half}_relaxed.broadPeak")
    params:
        outdir = os.path.join(REPRO_DIR, "peaks"),
        name = "{unit}.{half}",
        genome = MACS2_GENOME,
        pvalue = config["idr_relaxed_pvalue"],
        top_n = config["idr_top_n_peaks"],
        broad_cutoff = BROAD_CUTOFF,
        control_arg = unit_control_arg
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/idr_reproducibility/relaxed_{unit}.{half}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/idr_reproducibility
        macs2 callpeak -t {input.bam} {params.control_arg} -f BAMPE -g {params.genome} \
            --outdir {params.outdir} -n {params.name}_relaxedtmp \
            --nomodel --broad --broad-cutoff {params.broad_cutoff} -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.broadPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_peaks.gappedPeak \
              {params.outdir}/{params.name}_relaxedtmp_sorted.broadPeak
        """


# 21. IDR between a unit's two pseudo-halves -> count of reproducible peaks.
rule repro_idr_narrow:
    wildcard_constraints:
        unit = _alt(REPRO_NARROW_UNITS)
    input:
        pr1 = os.path.join(REPRO_DIR, "peaks", "{unit}.pr1_relaxed.narrowPeak"),
        pr2 = os.path.join(REPRO_DIR, "peaks", "{unit}.pr2_relaxed.narrowPeak")
    output:
        npeaks = os.path.join(REPRO_DIR, "idr", "{unit}.n.txt")
    params:
        threshold = config["idr_threshold"],
        idr_out = os.path.join(REPRO_DIR, "idr", "{unit}.idr.txt")
    conda:
        "../envs/idr.yaml"
    log:
        "logs/idr_reproducibility/idr_{unit}.log"
    shell:
        """
        mkdir -p {REPRO_DIR}/idr logs/idr_reproducibility
        idr --samples {input.pr1} {input.pr2} \
            --input-file-type narrowPeak --rank p.value \
            --idr-threshold {params.threshold} \
            --output-file {params.idr_out} > {log} 2>&1
        awk -v t={params.threshold} 'BEGIN{{c=-log(t)/log(10)}} $12>=c' {params.idr_out} | wc -l > {output.npeaks}
        """


rule repro_idr_broad:
    wildcard_constraints:
        unit = _alt(REPRO_BROAD_UNITS)
    input:
        pr1 = os.path.join(REPRO_DIR, "peaks", "{unit}.pr1_relaxed.broadPeak"),
        pr2 = os.path.join(REPRO_DIR, "peaks", "{unit}.pr2_relaxed.broadPeak")
    output:
        npeaks = os.path.join(REPRO_DIR, "idr", "{unit}.n.txt")
    params:
        threshold = config["idr_threshold"],
        idr_out = os.path.join(REPRO_DIR, "idr", "{unit}.idr.txt")
    conda:
        "../envs/idr.yaml"
    log:
        "logs/idr_reproducibility/idr_{unit}.log"
    shell:
        """
        mkdir -p {REPRO_DIR}/idr logs/idr_reproducibility
        idr --samples {input.pr1} {input.pr2} \
            --input-file-type broadPeak --rank p.value \
            --idr-threshold {params.threshold} \
            --output-file {params.idr_out} > {log} 2>&1
        awk -v t={params.threshold} 'BEGIN{{c=-log(t)/log(10)}} $11>=c' {params.idr_out} | wc -l > {output.npeaks}
        """


# 22. Reproducibility summary: self-consistency + rescue ratios per condition.
rule idr_reproducibility_summary:
    input:
        true_idr = [idr_peak_file(g) for g in IDR_GROUPS],
        counts = [os.path.join(REPRO_DIR, "idr", f"{u}.n.txt") for u in REPRO_UNITS]
    output:
        tsv = os.path.join(QC_DIR, "idr_reproducibility.tsv"),
        mqc = os.path.join(QC_DIR, "idr_reproducibility_mqc.txt")
    params:
        groups = IDR_GROUPS,
        members = {g: GROUPS[g] for g in IDR_GROUPS},
        true_idr = {g: idr_peak_file(g) for g in IDR_GROUPS},
        repro_idr_dir = os.path.join(REPRO_DIR, "idr")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/idr_reproducibility/summary.log"
    script:
        "../scripts/idr_reproducibility.py"


# 23. FastQC-only MultiQC report (everything else moves to the interactive report)
rule multiqc_fastqc:
    input:
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_001_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_001_fastqc.html"), sample=SAMPLES)
    output:
        html = os.path.join(QC_DIR, "multiqc_fastqc.html")
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/multiqc_fastqc/multiqc.log"
    shell:
        """
        mkdir -p {QC_DIR} logs/multiqc_fastqc
        multiqc -f -m fastqc {FASTQC_DIR}/ \
            --outdir {QC_DIR} \
            --filename multiqc_fastqc.html > {log} 2>&1
        """


# 15. Interactive HTML QC report aggregating both pipelines' numeric QC tables and
#     deepTools QC (drawn client-side as interactive SVG/canvas charts).
rule qc_report:
    input:
        # numeric sources
        expand(os.path.join(ALIGN_DIR, "{sample}.bowtie2.log"), sample=SAMPLES),
        expand(os.path.join(FILTERED_DIR, "{sample}.idxstats.txt"), sample=SAMPLES),
        expand(os.path.join(DEDUP_DIR, "{sample}.dedup.metrics.txt"), sample=SAMPLES),
        expand(os.path.join(COMPLEXITY_DIR, "{sample}_complexity.txt"), sample=SAMPLES),
        os.path.join(QC_DIR, "peak_summary.tsv"),
        os.path.join(QC_DIR, "tss_enrichment_scores.tsv"),
        os.path.join(QC_DIR, "blacklist_filtering_stats.txt"),
        os.path.join(ANNOT_DIR, "reads_in_annotations.tsv"),
        os.path.join(QC_DIR, "cross_correlation.tsv"),
        os.path.join(QC_DIR, "fingerprint_jsd.tsv"),
        os.path.join(QC_DIR, "idr_reproducibility.tsv"),
        os.path.join(CONSENSUS_DIR, "consensus_peaks.bed"),
        # embedded plots (PNG)
        os.path.join(DEEPTOOLS_DIR, "fragmentSize.png"),
        os.path.join(DEEPTOOLS_DIR, "Heatmap_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "Profile_TSS.png"),
        os.path.join(DEEPTOOLS_DIR, "chipseq_fingerprint.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_heatmap.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_PCA.png"),
        os.path.join(DEEPTOOLS_DIR, "deeptools_scatterplot.png"),
        expand(os.path.join(DEEPTOOLS_DIR, "{sample}.gc_content.png"), sample=SAMPLES),
        os.path.join(DEEPTOOLS_DIR, "fragment_lengths.txt"),
        os.path.join(DEEPTOOLS_DIR, "correlation_matrix.tab"),
        os.path.join(DEEPTOOLS_DIR, "tss_heatmap_downsampled.json"),
    output:
        html = os.path.join(QC_DIR, "chipseq_qc_report.html")
    params:
        results = RESULT_DIR,
        samples = ",".join(SAMPLES),
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/qc_report/qc_report.log"
    shell:
        """
        mkdir -p {QC_DIR} logs/qc_report
        python workflow/scripts/build_qc_report.py \
            --results-dir {params.results} \
            --out {output.html} \
            --samples {params.samples} \
            --generated "$(date -u '+%Y-%m-%d %H:%M UTC')" > {log} 2>&1
        """
