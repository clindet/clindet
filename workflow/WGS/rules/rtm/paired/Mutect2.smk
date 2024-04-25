# rule M2_ST:
#     input:
#         log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam"
#     output:
#         # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
#         table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table"
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         ref=config['resources'][genome_version]['REFFA'],
#         temp_directory=config['params']['java']['temp_directory'],
#         af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         GetPileupSummaries -R {params.ref}  \
#         -I {input.Tum} \
#         -V {params.af_vcf} \
#         -O {output.table}
#         """

# rule M2_SNC:
#     input:
#         log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
#     output:
#         # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
#         table="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         ref=config['resources'][genome_version]['REFFA'],
#         temp_directory=config['params']['java']['temp_directory'],
#         af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         GetPileupSummaries -R {params.ref}  \
#         -I {input.NC} \
#         -V {params.af_vcf} \
#         -O {output.table}
#         """


# rule M2_contam:
#     input:
#         log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
#         T="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
#         NC="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
#     output:
#         # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
#         seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
#         ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         ref=config['resources'][genome_version]['REFFA'],
#         temp_directory=config['params']['java']['temp_directory']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         CalculateContamination \
#         -I {input.T} \
#         -matched {input.NC} \
#         -tumor-segmentation {output.seg} \
#         -O {output.ctam} 
#         """


# rule call_variants_q2:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#         pon="{project}/{genome_version}/results/vcf/paired/PoN/MM_pon.vcf.gz",
#         pon_log='{project}/{genome_version}/logs/paired/Mutect2_PoNVCF_MM.log'
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         temp_directory=config['params']['java']['temp_directory'],
#         af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         Mutect2 -R {input.ref} \
#         --native-pair-hmm-threads {threads} \
#         -I {input.Tum} \
#         -I {input.NC} \
#         -O {output.vcf} \
#         -normal {wildcards.sample}_NC \
#         -pon {input.pon} \
#         --germline-resource {params.af_vcf}
#         """

# rule M2_filter:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf",
#         table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
#         seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
#         ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_filter.vcf"
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         temp_directory=config['params']['java']['temp_directory'],
#         af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         FilterMutectCalls \
#         -R {input.ref} \
#         -V {input.vcf} \
#         --contamination-table {input.ctam} \
#         --stats {input.vcf}.stats \
#         --tumor-segmentation {input.seg} \
#         -O {output.vcf}
#         """


rule mutect2:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        dbsnp=config['resources'][genome_version]['DBSNP']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        Mutect2 \
        --reference {input.ref} \
        --input {input.NC} \
        --input {input.Tum} \
        --normal-sample {wildcards.sample}-NC \
        --tumor-sample {wildcards.sample}-T \
        --dbsnp {params.dbsnp} \
        --output {output.vcf}
        """


rule mutect2_filter:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        dbsnp=config['resources'][genome_version]['DBSNP']
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        FilterMutectCalls \
        --variant {input.vcf} \
        --output {output.vcf}
        """
