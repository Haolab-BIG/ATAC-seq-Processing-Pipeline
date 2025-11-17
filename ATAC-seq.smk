configfile : "config.yaml"

SAMPLES_NEST = config.get("samples", {})

SAMPLE2GROUP = {
    sample: group
    for group, sd in SAMPLES_NEST.items()
    for sample in sd.keys()
}

ALL_SAMPLES = SAMPLE2GROUP.keys()

ALL_GROUPS = set(SAMPLE2GROUP.values())

def fq_path(sample, read):
    grp = SAMPLE2GROUP.get(sample)
    return SAMPLES_NEST[grp][sample][read]

def grp_sample(group):
    SAMPLE=[]
    for sample in SAMPLE2GROUP:
        grp = SAMPLE2GROUP.get(sample)
        if grp==group:
            SAMPLE.append(sample)
    return SAMPLE

def grp_sample_bam(group):
    BAM=[]
    for sample in SAMPLE2GROUP:
        grp = SAMPLE2GROUP.get(sample)
        if grp==group:
            BAM.append(f"sorted/{sample}.sorted.bam")
    return BAM

def make_output():
    output = []
    output += expand("data/{sample}/{sample}_R{Num}_val_{Num}_fastqc.html", sample=ALL_SAMPLES, Num=["1","2"])
    output += expand("bigwig/{sample}.bw", sample=ALL_SAMPLES)
    output += expand("data/multiqc_report/{sample}/multiqc_report.html", sample=ALL_SAMPLES)
    output.append("TFBS/bindetect_results.txt")
    if config["species"]=="hsa":
        output.append("diff/annotated_peak.csv")
    return output

rule all:
    input:
        make_output()

rule qc_before_filtered:
    output:
        html_R1=temp("qc/{sample}/{sample}_R1_fastqc.html"),
        html_R2=temp("qc/{sample}/{sample}_R2_fastqc.html"),
    params:
        FQ1 = lambda wc: fq_path(wc.sample, "R1"),
        FQ2 = lambda wc: fq_path(wc.sample, "R2"),
        threads=8,
    container:
        config["container"]
    shell:
        """
        fastqc -t {params.threads} {params.FQ1} {params.FQ2} -o qc/{wildcards.sample}
        """

rule trim_fastq:
    input:
        html_R1="qc/{sample}/{sample}_R1_fastqc.html",
        html_R2="qc/{sample}/{sample}_R2_fastqc.html",
    output:
        trimmed_FQ1=temp("qc/{sample}/{sample}_R1_val_1.fq.gz"),
        trimmed_FQ2=temp("qc/{sample}/{sample}_R2_val_2.fq.gz"),
    params:
        FQ1 = lambda wc: fq_path(wc.sample, "R1"),
        FQ2 = lambda wc: fq_path(wc.sample, "R2"),
    threads: 
        8
    container:
        config["container"]
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} -o qc/{wildcards.sample} {params.FQ1} {params.FQ2}
        """

rule qc_after_filtered:
    input:
        trimmed_FQ1="qc/{sample}/{sample}_R1_val_1.fq.gz",
        trimmed_FQ2="qc/{sample}/{sample}_R2_val_2.fq.gz",
    output:
        html_R1="data/{sample}/{sample}_R1_val_1_fastqc.html",
        html_R2="data/{sample}/{sample}_R2_val_2_fastqc.html",
    threads: 
        8
    container:
        config["container"]
    shell:
        """
        fastqc -t {threads} qc/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz qc/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz -o data/{wildcards.sample}
        """

rule bwa_mem:
    input:
        trimmed_FQ1="qc/{sample}/{sample}_R1_val_1.fq.gz",
        trimmed_FQ2="qc/{sample}/{sample}_R2_val_2.fq.gz",
    output:
        mapped_bam=temp("mapped/{sample}.bam"),
    params:
        BWAMEM2IDX = config["BWA-MEM2IDX"],
        home=config["home"],
        threads=8
    container:
        config["container"]
    shell:
        """ 
        bwa-mem2 mem -t {params.threads} -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" {params.BWAMEM2IDX}  qc/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz qc/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz | \
        samblaster | sambamba view -S -f bam /dev/stdin -t {params.threads} | sambamba sort /dev/stdin -t {params.threads} -o {output.mapped_bam} -q
        """

rule filter:
    input:
        mapped_bam="mapped/{sample}.bam",
    output:
        filtered_bam="filtered_bam/{sample}.bam",
    params:
        home=config["home"],
        threads=4
    container:
        config["container"]
    shell:
        """ 
        samtools idxstats {input.mapped_bam} | cut -f 1 | grep -v "chrUn" | grep -v "chrM" | \
        xargs samtools view -bh -@ {params.threads} -F 3844 -q 30 {input.mapped_bam} > {output.filtered_bam}
        """

rule index_filtered_bam:
    input:
        filtered_bam="filtered_bam/{sample}.bam",
    output:
        filtered_bam_index=temp("filtered_bam/{sample}.bam.bai"),
    threads:
        4
    container:
        config["container"]
    shell:
        """ 
        samtools index -@ {threads} {input.filtered_bam}
        """

rule shift_reads:
    input:
        filtered_bam_index="filtered_bam/{sample}.bam.bai",
    output:
        shifted_bam_tmp=temp("shifted_bam/{sample}.bam.tmp"),
        short_metrics=temp("shifted_bam/{sample}.short.metrics"),
    threads:
        4
    params:
        filtered_bam="filtered_bam/{sample}.bam",
    container:
        config["container"]
    shell:
        """ 
        alignmentSieve -b {params.filtered_bam} -o shifted_bam/{wildcards.sample}.bam.tmp -p {threads} --filterMetrics shifted_bam/{wildcards.sample}.short.metrics --ATACshift
        """

rule sort_shifted_bam:
    input:
        shifted_bam_tmp="shifted_bam/{sample}.bam.tmp",
    output:
        shifted_bam="shifted_bam/{sample}.bam",
    threads:
        4
    params:
        filtered_bam="filtered_bam/{sample}.bam",
    container:
        config["container"]
    shell:
        """ 
        samtools sort -o {output.shifted_bam} -@ {threads} {input.shifted_bam_tmp} -T shifted_bam/{wildcards.sample}
        """

rule index_shifted_bam:
    input:
        shifted_bam="shifted_bam/{sample}.bam",
    output:
        shifted_bam_index=temp("shifted_bam/{sample}.bam.bai"),
    threads:
        4
    container:
        config["container"]
    shell:
        """ 
        samtools index -@ {threads} {input.shifted_bam}
        """

rule sort_by_name:
    input:
        shifted_bam_index="shifted_bam/{sample}.bam.bai",
    output:
        name_sorted="name_sorted/{sample}.bam",
    threads:
        2
    params:
        shifted_bam="shifted_bam/{sample}.bam",
    container:
        config["container"]
    shell:
        """ 
        if [ ! -d "tmp" ]; then 
        mkdir tmp
        fi
        sambamba sort -n --nthreads {threads} --tmpdir tmp --out {output.name_sorted} {params.shifted_bam}
        """

rule genrich:
    input:
        name_sorted="name_sorted/{sample}.bam",
    output:
        narrowPeak="peaks/genrich/{sample}/{sample}_peaks.narrowPeak",
    container:
        config["container"]
    shell:
        """ 
        chr_exl=$(samtools view -H {input.name_sorted} | grep -E "chrUn|chrM" | cut -d \':\' -f 2| cut -d $\'\\t\' -f 1 | sed -z \'s/\\\\n/,/g\')
        Genrich -t "{input.name_sorted}" -o peaks/genrich/{wildcards.sample}/{wildcards.sample}_peaks.narrowPeak \
        -f peaks/genrich/{wildcards.sample}/{wildcards.sample}.f.bedgraph -k peaks/genrich/{wildcards.sample}/{wildcards.sample}.k.bedgraph \
        -e $chr_exl -j -y -D -d 150 -p 0.01 -g 50
        """

rule sort_shifted_bam2:
    input:
        narrowPeak="peaks/genrich/{sample}/{sample}_peaks.narrowPeak",
    output:
        mapped_sorted="sorted/{sample}.sorted.bam",
    threads:
        10
    params:
        name_sorted="name_sorted/{sample}.bam",
    container:
        config["container"]
    shell:
        """ 
        samtools sort {params.name_sorted} -o {output.mapped_sorted} -@ {threads} -T mapped/{wildcards.sample}
        """

rule index_mapped_bam:
    input:
        mapped_sorted="sorted/{sample}.sorted.bam",
    output:
        mapped_sorted_index="sorted/{sample}.sorted.bam.bai",
    container:
        config["container"]
    shell:
        """ 
        samtools index {input.mapped_sorted}
        """

rule merge_bam:
    input:
        expand(rules.index_mapped_bam.output.mapped_sorted_index, sample=lambda wc: grp_sample(wc.group))
    output:
        mapped_bam="group/{group}.bam",
    threads:
        10
    params:
        mapped_sorted_files_str=lambda wc: " ".join(grp_sample_bam(wc.group))
    container:
        config["container"]
    shell:
        """
        samtools merge -@ {threads} {output.mapped_bam} {params.mapped_sorted_files_str}
        """

rule index_group_bam:
    input:
        mapped_bam="group/{group}.bam",
    output:
        group_bam_index=temp("group/{group}.bam.bai"),
    container:
        config["container"]
    shell:
        """
        samtools index {input.mapped_bam}
        """

rule bamPE_fragment_size:
    input:
        expand(rules.index_group_bam.output.group_bam_index, group=ALL_GROUPS)
    output:
        fragmentSize_metric=temp("deepTools_qc/bamPEFragmentSize/fragmentSize.metric.pdf"),
    params:
        sorted_files_str=lambda w: " ".join(expand("sorted/{sample}.sorted.bam", sample=ALL_SAMPLES)),
        threads=24,
    container:
        config["container"]
    shell:
        """
        bamPEFragmentSize \
        --bamfiles {params.sorted_files_str} \
        --binSize 1000000 \
        -o {output.fragmentSize_metric} \
        --table deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv \
        -p {params.threads}
        """

rule CSAW:
    input:
        fragmentSize_metric="deepTools_qc/bamPEFragmentSize/fragmentSize.metric.pdf",
    output:
        DiffBinding_significant=temp("DiffBinding_significant.bed"),
    container:
        config["container"]
    shell:
        """
        Rscript scripts/CSAW_ATACseq.R samplesheet.txt sorted/ peaks/genrich/ deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv
        """

rule tidy_CSAW:
    input:
        DiffBinding_significant="DiffBinding_significant.bed",
    output:
        tidy_DiffBinding_significant="diff/DiffBinding_significant.bed",
    container:
        config["container"]
    shell:
        """
        mv {input.DiffBinding_significant} {output.tidy_DiffBinding_significant}
        """

rule atacorrect:
    input:
        tidy_DiffBinding_significant="diff/DiffBinding_significant.bed",
    output:
        corrected_bw=temp("bias_correction/{group}_corrected.bw"),
    threads:
        25
    params:
        fa_in=config["fa_in"],
        group_bam="group/{group}.bam",
    container:
        config["container"]
    shell:
        """
        TOBIAS ATACorrect -b {params.group_bam} \
                          -g {params.fa_in} \
                          -p {input.tidy_DiffBinding_significant} \
                          --cores {threads} \
                          --outdir bias_correction \
                          --prefix {wildcards.group}
        """

rule footprinting:
    input:
        corrected_bw="bias_correction/{group}_corrected.bw",
    output:
        footprints_bw=temp("footprinting/{group}_footprints.bw"),
    threads:
        25
    params:
        tidy_DiffBinding_significant="diff/DiffBinding_significant.bed",
    container:
        config["container"] 
    shell:
        """
        TOBIAS FootprintScores --signal {input.corrected_bw} \
                               --regions {params.tidy_DiffBinding_significant} \
                               --output {output.footprints_bw} \
                               --cores {threads}
        """

rule bindetect:
    input:
        expand(rules.footprinting.output.footprints_bw, group=ALL_GROUPS),
    output:
        TFBS="TFBS/bindetect_results.txt",
    threads:
        25
    params:
        tidy_DiffBinding_significant="diff/DiffBinding_significant.bed",
        footprinting=lambda w: " ".join(expand(rules.footprinting.output.footprints_bw, group=ALL_GROUPS)),
        fa_in=config["fa_in"],
        all_motifs=config["all_motifs"],
    container:
        config["container"] 
    shell:
        """
        TOBIAS BINDetect --motifs {params.all_motifs} \
                         --signals {params.footprinting} \
                         --genome {params.fa_in} \
                         --peaks {params.tidy_DiffBinding_significant} \
                         --cores {threads} \
                         --motif-pvalue 0.01 \
                         --outdir TFBS/
            """

rule bamCoverage:
    input:
        mapped_sorted="sorted/{sample}.sorted.bam",
    output:
        bw="bigwig/{sample}.bw",
    params:
        threads=25,
    container:
        config["container"] 
    shell:
        """
        bamCoverage -b {input.mapped_sorted} -o {output.bw} -of bigwig -bs 50 -p {params.threads} --ignoreDuplicates --minMappingQuality 20
        """

if config["species"]=="hsa":
    rule ChIPpeakAnno:
        input:
            tidy_DiffBinding_significant="diff/DiffBinding_significant.bed",
        output:
            annotated_peak="diff/annotated_peak.csv",
        container:
            config["container"] 
        shell:
            """
            Rscript scripts/CHIPPEAKANNO_ATACseq.R {input.tidy_DiffBinding_significant}
            """

rule Final_MultiQC:
    input:
        name_sorted="name_sorted/{sample}.bam",
    output:
        multiqc_report_html="data/multiqc_report/{sample}/multiqc_report.html",
    container:
        config["container"] 
    shell:
        """
        multiqc data/ -o data/multiqc_report/{wildcards.sample} --force
        """