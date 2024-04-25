rule map_reads:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}-{group}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}-{group}_R2.fastq.gz"
    output:
        "{project}/{genome_version}/results/mapped/{sample_type}/{sample}-{group}.sorted.bam",
    params:
        ref=config['resources'][genome_version]['REFFA'],
        rg=r"@RG\tID:{sample}_{group}\tPL:Illumia\tSM:{sample}_{group}\tLB:DNA-Seq"
    threads: 30 
    conda:
        config['softwares']['samtools']['conda']
    shell:
        """ {config[softwares][bwa][mem][call]} -t 30 -MR '{params.rg}' \
        {params.ref} \
        {input.R1} {input.R2} | samtools fixmate -O bam - - | \
        samtools sort -@ 30 -O bam -o {output}
        """

rule mark_duplicates:
    input:
        "{project}/{genome_version}/results/mapped/{sample_type}/{sample}-{group}.sorted.bam",
    output:
        bam=temp("{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam"),
        metrics="{project}/{genome_version}/results/qc/dedup/{sample_type}/{sample}-{group}.metrics.txt"
    params:
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        {config[softwares][gatk4][MarkDuplicates][call]} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
        -I {input} \
        -O {output.bam} \
        -M {output.metrics}
        """

rule recalibrate_base_qualities:
    input:
        bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
        #bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
        ref=config['resources'][genome_version]['REFFA'],
        known_1=config['resources']['varanno'][genome_version]['KNOWN_SITES1'],
        known_2=config['resources']['varanno'][genome_version]['KNOWN_SITES2']
    output:
        recal_table="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.grp",
    params:
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """ 
            {config[softwares][gatk4][BaseRecalibrator][call]} --use-original-qualities -R {input.ref} \
            -I {input.bam} \
            -O {output.recal_table} \
            --known-sites {input.known_1} \
            --known-sites {input.known_2}
        """


rule apply_base_quality_recalibration:
    input:
        bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
        #bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam.bai",
        ref=config['resources'][genome_version]['REFFA'],
        recal_table="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.grp",
    output:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
    conda:
        config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """{config[softwares][gatk4][ApplyBQSR][call]} --add-output-sam-program-record -use-original-qualities \
            --bqsr-recal-file {input.recal_table} \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam}  
            samtools index {output.bam}
        """


rule bed_to_interval_list:
    input:
        bed=get_sample_bed,
        dict=config['resources'][genome_version]['REFFA_DICT'],
    output:
        it="{project}/{genome_version}/picard/{sample}.interval_list"
    params:
        extra="--SORT true",  # sort output interval list before writing
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
            export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
            java -jar /public/ClinicalExam/lj_sih/softwares/picard.jar BedToIntervalList --INPUT {input.bed} \
            --SEQUENCE_DICTIONARY {input.dict} {params.extra} \
            --OUTPUT {output.it}
        """


rule picard_collect_wes:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
        ref=config['resources'][genome_version]['REFFA'],
        reference=config['resources'][genome_version]['REFFA'],
        bait_intervals="{project}/{genome_version}/picard/{sample}.interval_list",
        target_intervals="{project}/{genome_version}/picard/{sample}.interval_list",
    output:
        hs="{project}/{genome_version}/results/stats/{sample_type}/wes_metrics/{sample}-{group}.txt"
    resources:
        mem_mb=1024
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
        java -jar /public/ClinicalExam/lj_sih/softwares/picard.jar CollectHsMetrics \
        --INPUT {input.bam} \
        --OUTPUT {output.hs} \
        --REFERENCE_SEQUENCE {input.reference} \
        --BAIT_INTERVALS {input.bait_intervals} \
        --TARGET_INTERVALS {input.target_intervals}
        """

# rule picard_collect_hs_metrics:
#     input:
#         bam="{project}/{genome_version}/results/recal/{sample}-{group}.bam",
#         reference=config['resources'][genome_version]['REFFA'],
#         bait_intervals="/public/ClinicalExam/lj_sih/projects/project_QC/data/{sample}/{sample}.interval_list",
#         target_intervals="/public/ClinicalExam/lj_sih/projects/project_QC/data/{sample}/{sample}.interval_list",
#     output:
#         "result/stats/hs_metrics/{sample}-{group}.txt",
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v1.7.0/bio/picard/collecthsmetrics"