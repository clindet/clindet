rule STAR_1_pass:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz"
        # SJ="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.out.tab"
    output:
        pass1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_pass1.log",
        sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/", # "/"" mustq in the config string
        ref=config['resources'][genome_version]['REFFA'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}"
    threads: 10
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype None --outFileNamePrefix {params.out_dir} \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat
        awk '$1~/chr[1-2XY]/ && $6==0 && $5>0 && $7>0' {params.out_dir}SJ.out.tab > {output.sj_filt}
        touch {output.pass1}
        """



rule STAR_map:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz",
        sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
    output:
        bam=temp("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam"),
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
        um_fq1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R1.fq",
        um_fq2="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R2.fq"
    params:
        out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/", # "/"" must in the config string
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        star_index=config['softwares']['star']['index'][genome_version],
        rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}"
    threads: 10
    conda:
        config['softwares']['star']['conda']
    shell:
        """ 
        STAR --genomeDir {params.star_index} --runThreadN={threads} \
            --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.out_dir} \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline '{params.rg}' --sjdbFileChrStartEnd {input.sj_filt}\
            --sjdbGTFfile {params.gtf} \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c

        mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.out_dir}/Unmapped.out.mate1 {output.um_fq1}
        mv {params.out_dir}/Unmapped.out.mate2 {output.um_fq2}

        touch {output.stamp}
        """


# rule cal_count_matrix:
#     input:
#         gtf=config['resources'][genome_version]['GTF'],
#         index=expand("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",sample = paired_samples)
#     params:
#         bam=expand("{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",sample = paired_samples)
#     conda:
#         config['softwares']['subread']['conda']
#     threads: 60
#     output:
#         "{project}/{genome_version}/results/summary/feature_count_all_{sample}.txt"
#     shell:
#         """ featureCounts -T {threads} \
#            -a {input.gtf} \
#            -t gene -g gene_name \
#            -o {output} {params.bam}"""



## RSEM need re-write, RSEM not work with genome type mapping

# rule STAR_RSEM_map:
#     input:
#         R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
#         R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz",
#         sj_filt='{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.filt'
#     output:
#         bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.rsem.sorted.bam",
#         stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
#         um_fq1="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R1.fq",
#         um_fq2="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_unmapped_R2.fq"
#     params:
#         out_dir="{project}/{genome_version}/results/mapped/STAR/{sample}/", # "/"" must in the config string
#         ref=config['resources'][genome_version]['REFFA'],
#         gtf=config['resources'][genome_version]['GTF'],
#         star_index=config['softwares']['star']['index'][genome_version],
#         rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}"
#     threads: 10
#     conda:
#         config['softwares']['star']['conda']
#     shell:
#         """ 
#         STAR --genomeDir {params.star_index} --runThreadN={threads} \
#             --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.out_dir} \
#             --outReadsUnmapped Fastx \
#             --outSAMattrRGline '{params.rg}' --sjdbFileChrStartEnd {input.sj_filt}\
#             --sjdbGTFfile {params.gtf} \
#             --outFilterMultimapNmax 50 \
#             --peOverlapNbasesMin 10 \
#             --alignSplicedMateMapLminOverLmate 0.5 \
#             --alignSJstitchMismatchNmax 5 -1 5 5 \
#             --chimSegmentMin 10 \
#             --chimOutType WithinBAM HardClip \
#             --chimJunctionOverhangMin 10 \
#             --chimScoreDropMax 30 \
#             --chimScoreJunctionNonGTAG 0 \
#             --chimScoreSeparation 1 \
#             --chimSegmentReadGapMax 3 \
#             --chimMultimapNmax 50 \
#             --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c

#         mv {params.out_dir}/Aligned.sortedByCoord.out.bam {output.bam}
#         mv {params.out_dir}/Unmapped.out.mate1 {output.um_fq1}
#         mv {params.out_dir}/Unmapped.out.mate2 {output.um_fq2}

#         touch {output.stamp}
#         """


rule cal_exp_RSEM:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz"
    params:
        rsem_ref=config['softwares']['rsem']['index'][genome_version],
        result_prefix="{project}/{genome_version}/results/summary/RSEM"
    threads: 10
    conda:
        config['softwares']['star']['conda']
    output:
        genes="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results",
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam"
        # bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results"
    shell:
        """{config[softwares][rsem][cal_exp]}   --paired-end \
            --star -p {threads} --star-output-genome-bam \
            --star-gzipped-read-file  --sort-bam-by-coordinate \
            {input.R1} {input.R2} \
           {params.rsem_ref} {params.result_prefix}/{wildcards.sample}/{wildcards.sample}
           """

rule RSEM_sort_bam:
    input:
        genes="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results",
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.bam",
    threads: 10
    conda:
        'samtools114'
    output:
        sort_bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam"
        # bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.genes.results"
    shell:
        """
        samtools sort -@ {threads} -o {output.sort_bam} {input.bam}
        samtools index {output.sort_bam}
        """

rule RSEM_bam2bigwig:
    input:
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam"
    output:
        bw="{project}/{genome_version}/results/RSEM/bigwig/{sample}/{sample}.bw"
    conda:'deeptools'
    threads:10
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} -p {threads} --normalizeUsing RPKM
        """


### kallisto index 

###### kallisto index 
# rule kallisto_index:
#     input:
#         ref_cdna=config['resources'][genome_version]['CDNA']
#     output:
#         index=config['resources'][genome_version]['CDNA'] + '.idx'
#     threads: 10
#     conda:
#         config['softwares']['kallisto']['conda']
#     shell:
#         """
#         kallisto index {input.ref_cdna}
#         """

           

###### kallisto quanto



### samlom 


# ######  summarys
# rule bed_to_interval_list:
#     input:
#         bed=get_sample_bed,
#         dict=config['resources'][genome_version]['REFFA_DICT'],
#     output:
#         "/public/ClinicalExam/lj_sih/projects/project_QC/data/{sample}/{sample}.interval_list",
#     params:
#         extra="--SORT true",  # sort output interval list before writing
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v1.7.0/bio/picard/bedtointervallist"


# rule picard_collect_hs_metrics:
#     input:
#         bam="{project}/{genome_version}/results/recal/{sample}.bam",
#         reference=config['resources'][genome_version]['REFFA'],
#         bait_intervals="/public/ClinicalExam/lj_sih/projects/project_QC/data/{sample}/{sample}.interval_list",
#         target_intervals="/public/ClinicalExam/lj_sih/projects/project_QC/data/{sample}/{sample}.interval_list",
#     output:
#         "result/stats/hs_metrics/{sample}.txt",
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v1.7.0/bio/picard/collecthsmetrics"


# rule STAR_oneshot:
#     input:
#         R1="{project}/{genome_version}/results/trimmed/{sample}_R1.fastq.gz",
#         R2="{project}/{genome_version}/results/trimmed/{sample}_R2.fastq.gz"
#         # SJ="{project}/{genome_version}/results/mapped/STAR/{sample}/_STARpass1/SJ.out.tab"
#     output:
#         # oneshot="{project}/{genome_version}/results/mapped/STAR_oneshot/{sample}/{sample}_oneshot.log",
#         bam="{project}/{genome_version}/results/mapped/STAR_oneshot/{sample}/{sample}.sorted.bam"
#     params:
#         out_dir="{project}/{genome_version}/results/mapped/STAR_oneshot/{sample}/", # "/"" mustq in the config string
#         ref=config['resources'][genome_version]['REFFA'],
#         star_index=config['softwares']['star']['index'][genome_version],
#         rg=r"ID:{sample} PL:ILLUMINA.NovaSeq LB:RNA-Seq SM:{sample}",
#         gtf=config['resources'][genome_version]['GTF']
#     threads: 10
#     conda:
#         config['softwares']['star']['conda']
#     shell:
#         """ 
#         STAR --genomeDir {params.star_index} --runThreadN={threads} \
#             --outSAMtype None --outFileNamePrefix {params.out_dir} \
#             --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat
#             --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.out_dir} \
#             --genomeLoad NoSharedMemory --outSAMunmapped Within \
#             --outSAMattrRGline '{params.rg}' \
#             --sjdbGTFfile {params.gtf} \
#             --outFilterMultimapNmax 50 \
#             --peOverlapNbasesMin 10 \
#             --alignSplicedMateMapLminOverLmate 0.5 \
#             --alignSJstitchMismatchNmax 5 -1 5 5 \
#             --chimSegmentMin 10 \
#             --chimOutType WithinBAM HardClip \
#             --chimJunctionOverhangMin 10 \
#             --chimScoreDropMax 30 \
#             --chimScoreJunctionNonGTAG 0 \
#             --chimScoreSeparation 1 \
#             --chimSegmentReadGapMax 3 \
#             --chimMultimapNmax 50 \
#             --readFilesIn {input.R1} {input.R2} --readFilesCommand gunzip -c
#         """