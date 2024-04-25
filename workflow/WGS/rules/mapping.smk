rule map_reads:
    input:
        R1="{project}/{genome_version}/results/trimmed/{sample}-{group}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/trimmed/{sample}-{group}_R2.fastq.gz"
    output:
        temp("{project}/{genome_version}/results/mapped/{sample_type}/{sample}-{group}.sorted.bam"),
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
        bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
        metrics="{project}/{genome_version}/results/qc/dedup/{sample_type}/{sample}-{group}.metrics.txt"
    params:
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        {config[softwares][gatk4][MarkDuplicates][call]} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
        -I {input} \
        -O {output.bam} \
        -M {output.metrics} --TMP_DIR {params.temp_directory}
        """


### for faster run, may consider not run applyBQSR

if recal:
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
            """  export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
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
            """ export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
            {config[softwares][gatk4][ApplyBQSR][call]} --add-output-sam-program-record -use-original-qualities \
                --bqsr-recal-file {input.recal_table} \
                -R {input.ref} \
                -I {input.bam} \
                -O {output.bam} 
                samtools index {output.bam} 
            """
else:
    rule recal_link:
        input:
            bam="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bam",
            bai="{project}/{genome_version}/results/dedup/{sample_type}/{sample}-{group}.sorted.bai"
        output:
            bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
            bai="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam.bai"
        params:
            temp_directory=config['params']['java']['temp_directory']
        shell:
            """
            ln -s $(realpath {input.bam}) {output.bam}
            ln -s $(realpath {input.bai}) {output.bai}
            """

# rule rename_chr_bam:
#     input:
#         bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#     output:
#         bam="{project}/{genome_version}/results/recal_rename/{sample_type}/{sample}-{group}.bam"
#     conda:
#         config['softwares']['samtools']['conda']
#     params:
#         name_map=config['softwares']['jvarkit']['bamrenamechr']['hg19_remove_chr']
#     shell:
#         """ {config[softwares][jvarkit][bamrenamechr][call]} \
#             --bamcompression 9 -i \
#             --samoutputformat BAM \
#             -R {input.ref} -f {params.name_map} \
#             -o {output.bam}  {input.bam}
#             samtools index {output.bam}
#         """

rule telomerecat:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        len="{project}/{genome_version}/results/telomerecat/{sample_type}/{sample}-{group}.csv"
    conda:
        config['softwares']['samtools']['conda']
    threads:8
    params:
        outdir="{project}/{genome_version}/results/telomerecat/{sample_type}/{sample}-{group}"
    shell:
        """ telomerecat bam2length -p{threads} -v 2 \
            --temp_dir /public/ClinicalExam/lj_sih/tmp \
            -r {input.ref} --outbam_dir {params.outdir} \
            --output {output.len}  {input.bam}
        """


rule picard_collect_wgs:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        "{project}/{genome_version}/results/stats/{sample_type}/wgs_metrics/{sample}-{group}.txt"
    resources:
        mem_mb=1024
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && \
        java -jar {config[softwares][PICARD]} CollectWgsMetrics \
        -I {input.bam} \
        -O {output} -R {input.ref}
        """

rule picard_flength_wgs:
    input:
        bam="{project}/{genome_version}/results/recal/{sample_type}/{sample}-{group}.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        txt="{project}/{genome_version}/results/stats/{sample_type}/wgs_metrics/{sample}-{group}-fragment-length.txt",
        pdf="{project}/{genome_version}/results/stats/{sample_type}/wgs_metrics/{sample}-{group}-fragment-length.pdf"
    resources:
        mem_mb=1024
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory']
    shell:
        """
        {config[softwares][gatk4][call]} CollectInsertSizeMetrics \
        -I {input.bam} --TMP_DIR {params.temp_directory} \
        -O {output.txt} -R {input.ref} -H {output.pdf} -M 0.5
        """