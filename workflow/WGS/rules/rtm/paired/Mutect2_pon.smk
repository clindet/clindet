rule call_variants_pon:
    input:
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_Mutect2.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        Mutect2 -R {input.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.NC} --max-mnp-distance 0 \
        -O {output.vcf} \
        """

rule pon_GB:
    input:
        expand("{project}/{genome_version}/results/vcf/paired/PoN/{MM_sample}_Mutect2.vcf",MM_sample = paired_samples,project = project,genome_version = genome_version)
    output:
        log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_{sample}.log'
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input),
        bed=config['resources'][genome_version]['GENOME_BED']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GenomicsDBImport -R {params.ref} \
        --merge-input-intervals true --sites-only-vcf-output true \
        --genomicsdb-workspace-path /public/ClinicalExam/lj_sih/resource/mutFilter/pon_db_WGS \
        {params.vcfs}
        touch {output.log}
        """


# don't use -L params, will stack
rule M2_CSPN:
    input:
        '{project}/{genome_version}/logs/paired/Mutect2_PoNDB_{sample}.log'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        log='{project}/{genome_version}/logs/paired/Mutect2_PoNVCF_{sample}.log'
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF'],
        temp_directory=config['params']['java']['temp_directory'],
        vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input)
    threads: 20
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        CreateSomaticPanelOfNormals -R {params.ref} \
        --germline-resource {params.af_vcf} \
        -V gendb://../../../resource/mutFilter/pon_db_WGS \
        --sites-only-vcf-output true \
        -O {output.vcf}
        touch {output.log}
        """

rule M2_ST:
    input:
        log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {params.ref}  \
        -I {input.Tum} \
        -V {params.af_vcf} \
        -O {output.table}
        """

rule M2_SNC:
    input:
        log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {params.ref}  \
        -I {input.NC} \
        -V {params.af_vcf} \
        -O {output.table}
        """


rule M2_contam:
    input:
        log='{project}/{genome_version}/logs/paired/Mutect2_PoNDB_MM.log',
        T="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        NC="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        CalculateContamination \
        -I {input.T} \
        -matched {input.NC} \
        -tumor-segmentation {output.seg} \
        -O {output.ctam} 
        """


rule call_variants_q2:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        pon="{project}/{genome_version}/results/vcf/paired/PoN/MM_pon.vcf.gz",
        pon_log='{project}/{genome_version}/logs/paired/Mutect2_PoNVCF_MM.log'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        Mutect2 -R {input.ref} \
        --native-pair-hmm-threads {threads} \
        -I {input.Tum} \
        -I {input.NC} \
        -O {output.vcf} \
        -normal {wildcards.sample}_NC \
        -pon {input.pon} \
        --germline-resource {params.af_vcf}
        """

rule M2_filter:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        af_vcf=config['resources'][genome_version]['MUTECT2_VCF']
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        FilterMutectCalls \
        -R {input.ref} \
        -V {input.vcf} \
        --contamination-table {input.ctam} \
        --stats {input.vcf}.stats \
        --tumor-segmentation {input.seg} \
        -O {output.vcf}
        """
