rule arriba_fusion:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log"
    output:
        tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"
        # dis_tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion_discarded.tsv"
    # conda:
    #     config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        blacklist=config['softwares']['arriba'][genome_version]['blacklist'],
        known_fus=config['softwares']['arriba'][genome_version]['known_fusions'],
        pro_dom=config['softwares']['arriba'][genome_version]['protein_domains']
    singularity:
        config['singularity']['arriba']['sif']
    shell:
        """
            /arriba_v2.4.0/arriba \
            -x {input.bam}  -o {output.tsv} \
            -a {params.ref} -g {params.gtf} \
            -b {params.blacklist} -k {params.known_fus} \
            -t {params.known_fus} \
            -p {params.pro_dom}    
        """

rule arriba_draw:
    input:
        # bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log",
        tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv"
    output:
        # tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.tsv",
        pdf="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion.pdf"
        # dis_tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion_discarded.tsv"
    # conda:
    #     config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['resources'][genome_version]['REFFA'],
        gtf=config['resources'][genome_version]['GTF'],
        c_band=config['softwares']['arriba'][genome_version]['cytobands'],
        blacklist=config['softwares']['arriba'][genome_version]['blacklist'],
        known_fus=config['softwares']['arriba'][genome_version]['known_fusions'],
        pro_dom=config['softwares']['arriba'][genome_version]['protein_domains']
    singularity:
        config['singularity']['arriba']['sif']
    shell:
        """   
            /arriba_v2.4.0/draw_fusions.R \
            --fusions={input.tsv}\
            --alignments={input.bam} \
            --output={output.pdf} \
            --annotation={params.gtf} \
            --cytobands={params.c_band} \
            --proteinDomains={params.pro_dom}
        """



rule easyfuse_link:
    input:
        unpack(get_rna_fastq)
    output:
        link_dir=directory("{project}/{genome_version}/results/easyfuse/link_dir/{sample}"),
        R1="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_R2.fastq.gz",
        log="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_link.log"
        # dis_tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion_discarded.tsv"
    # conda:
    #     config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory'],
    singularity:
        config['singularity']['easyfuse']['sif']
    shell:
        """
            ln -s $(realpath {input.R1}) {output.R1}
            ln -s $(realpath {input.R2}) {output.R2}
            touch {output.log}
        """

rule easyfuse_call:
    input:
        link_dir="{project}/{genome_version}/results/easyfuse/link_dir/{sample}",
        R1="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_R1.fastq.gz",
        R2="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_R2.fastq.gz",
        log="{project}/{genome_version}/results/easyfuse/link_dir/{sample}/{sample}_link.log"
    output:
        out_dir=directory("{project}/{genome_version}/results/easyfuse/out/{sample}/FusionSummary"),
        # dis_tsv="{project}/{genome_version}/results/fusion/{sample}_arriba_fusion_discarded.tsv"
    # conda:
    #     config['softwares']['samtools']['conda']
    params:
        temp_directory=config['params']['java']['temp_directory'],
        out_dir="{project}/{genome_version}/results/easyfuse/out/{sample}"
    singularity:
        config['singularity']['easyfuse']['sif']
    shell:
        """
            python /code/easyfuse/processing.py \
            -i {input.link_dir} \
            -o {output.out_dir} 
        """

rule TRUST4_TBCR:
    input:
        bam="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}.sorted.bam",
        stamp="{project}/{genome_version}/results/mapped/STAR/{sample}/{sample}_star.log"
    output:
        "{project}/{genome_version}/results/IG/TRUST4/{sample}_report.tsv"
    conda:
        config['softwares']['trust4']['conda']
    threads:10
    params:
        temp_directory=config['params']['java']['temp_directory'],
        ref=config['softwares']['trust4'][genome_version]['ref'],
        f=config['softwares']['trust4'][genome_version]['f'],
        oref="{project}/{genome_version}/results/IG/TRUST4/{sample}"
    shell:
        """
        {config[softwares][trust4][call]} \
            -b {input.bam}  -f {params.f} \
            -o {params.oref} --ref {params.ref}  -t {threads}
        """


# rule arriba_fusion_oneshot:
#     input:
#         bam="{project}/{genome_version}/results/mapped/STAR_oneshot/{sample}/{sample}.sorted.bam"
#     output:
#         tsv="{project}/{genome_version}/results/fusion_oneshot/{sample}_arriba_fusion.tsv"
#         # dis_tsv="{project}/{genome_version}/results/fusion_oneshot/{sample}_arriba_fusion_discarded.tsv"
#     conda:
#         config['softwares']['samtools']['conda']
#     params:
#         temp_directory=config['params']['java']['temp_directory'],
#         ref=config['resources'][genome_version]['REFFA'],
#         gtf=config['resources'][genome_version]['GTF'],
#         blacklist=config['softwares']['arriba'][genome_version]['blacklist'],
#         known_fus=config['softwares']['arriba'][genome_version]['known_fusions'],
#         pro_dom=config['softwares']['arriba'][genome_version]['protein_domains']
#     singularity:
#         config['singularity']['arriba']['sif']
#     shell:
#         """
#             /arriba_v2.4.0/arriba \
#             -x {input.bam}  -o {output.tsv} \
#             -a {params.ref} -g {params.gtf} \
#             -b {params.blacklist} -k {params.known_fus} \
#             -t {params.known_fus} \
#             -p {params.pro_dom}    
#         """

