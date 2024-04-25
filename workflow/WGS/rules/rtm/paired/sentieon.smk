# rule pon_GB_pseudo:
#     input:
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA'],
#         dbsnp=config['resources'][genome_version]['DBSNP']
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz"
#     params:
#     threads: 10
#     shell:
#         """
#         /public/astra/bin/_sentieon/bin/sentieon driver -t {threads} -r {input.ref} -i {input.NC} \
#         --algo TNhaplotyper --detect_pon \
#         --cosmic COSMIC_VCF --dbsnp {input.dnsnp} {output.vcf}
#         """

# rule call_variants_pon:
#     input:
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         ref=config['resources'][genome_version]['REFFA']
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_Mutect2.vcf",
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         temp_directory=config['params']['java']['temp_directory'],
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         Mutect2 -R {input.ref} \
#         --native-pair-hmm-threads {threads} \
#         -I {input.NC} --max-mnp-distance 0 \
#         -O {output.vcf} \
#         """

# rule pon_GB:
#     input:
#         expand("{project}/{genome_version}/results/vcf/paired/PoN/{MM_sample}_Mutect2.vcf",MM_sample = paired_samples,project = project,genome_version = genome_version)
#     output:
#         log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_{sample}.log'
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         ref=config['resources'][genome_version]['REFFA'],
#         temp_directory=config['params']['java']['temp_directory'],
#         vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input),
#         bed=config['resources'][genome_version]['GENOME_BED']
#     threads: 10
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         GenomicsDBImport -R {params.ref} \
#         --merge-input-intervals true --sites-only-vcf-output true \
#         --genomicsdb-workspace-path /public/ClinicalExam/lj_sih/resource/mutFilter/pon_db_WGS \
#         {params.vcfs}
#         touch {output.log}
#         """


# # don't use -L params, will stack
# rule M2_CSPN:
#     input:
#         '{project}/{genome_version}/logs/paired/Mutect2_PoNDB_{sample}.log'
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
#         log='{project}/{genome_version}/logs/paired/Mutect2_PoNVCF_{sample}.log'
#     params:
#         gatk4=config['softwares']['gatk4']['call'],
#         ref=config['resources'][genome_version]['REFFA'],
#         af_vcf=config['resources'][genome_version]['MUTECT2_VCF'],
#         temp_directory=config['params']['java']['temp_directory'],
#         vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input)
#     threads: 20
#     shell:
#         """
#         export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
#         CreateSomaticPanelOfNormals -R {params.ref} \
#         --germline-resource {params.af_vcf} \
#         -V gendb://../../../resource/mutFilter/pon_db_WGS \
#         --sites-only-vcf-output true \
#         -O {output.vcf}
#         touch {output.log}
#         """

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


rule call_variants_sentieon:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        pon=config['resources'][genome_version]['WGS_PON']
    output:
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
        ori_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.ori_data"
        # contam_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.contam_data"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF'],
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
    threads: 10
    shell:
        """
        /public/astra/bin/_sentieon/bin/sentieon driver -t {threads} -r {input.ref} \
        -i {input.Tum}  \
        -i {input.NC} \
        --algo TNhaplotyper2 --tumor_sample {wildcards.sample}_T \
        --normal_sample {wildcards.sample}_NC \
        --pon {input.pon} \
        {params.temp_vcf} \
        --algo OrientationBias --tumor_sample {wildcards.sample}_T \
        {output.ori_data}
        """


rule filter_sentieon:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        temp_vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sent_temp.vcf",
        ori_data="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.ori_data"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/sentieon.vcf"
    threads: 10
    shell:
        """
        /public/astra/bin/_sentieon/bin/sentieon driver -t {threads} \
        -r  {input.ref}  \
        --algo TNfilter --tumor_sample {wildcards.sample}_T \
        --normal_sample {wildcards.sample}_NC \
        -v {input.temp_vcf}   \
        --orientation_priors {input.ori_data}  \
        {output.vcf}
        """
# --vcf GERMLINE_RESOURCE \
        # --algo ContaminationModel --tumor_sample {wildcards.sample}_T \
        # --normal_sample {wildcards.sample}_NC \
        # --tumor_segments {output.contam_data}.segments \
        # {output.contam_data} 
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
