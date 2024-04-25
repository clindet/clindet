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
        log='logs/paired/Mutect2_PoNDB_{sample}.log'
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input),
        bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed'
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GenomicsDBImport -R {params.ref} -L {params.bed} \
        --merge-input-intervals true --sites-only-vcf-output true \
        --genomicsdb-workspace-path /public/ClinicalExam/lj_sih/resource/mutFilter/pon_db \
        {params.vcfs}
        touch {output.log}
        """
# don't use -L params, will stack
rule M2_CSPN:
    input:
        'logs/paired/Mutect2_PoNDB_MM.log'
    output:
        vcf=expand("{project}/{genome_version}/results/vcf/paired/PoN/{batch}_pon.vcf.gz",batch = plat_batch,project = project,genome_version = genome_version),
        log=expand('logs/paired/Mutect2_PoNVCF_{batch}.log',batch = plat_batch,)
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed',
        vcfs=lambda wildcards, input: ' -V ' + ' -V '.join(input)
    threads: 20
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        CreateSomaticPanelOfNormals -R {params.ref} \
        --germline-resource /public/ClinicalExam/lj_sih/resource/af-only-gnomad.raw.sites.hg19.vcf.gz \
        -V gendb://../../../resource/mutFilter/pon_db \
        --sites-only-vcf-output true \
        -O {output.vcf}
        touch {output.log}
        """

rule M2_ST:
    input:
        log='logs/paired/Mutect2_PoNDB_MM.log',
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed'
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {params.ref}  \
        -I {input.Tum} \
        -V /public/ClinicalExam/lj_sih/resource/af-only-gnomad.raw.sites.hg19.vcf.gz \
        -L {params.bed} \
        -O {output.table}
        """

rule M2_SNC:
    input:
        log='logs/paired/Mutect2_PoNDB_MM.log',
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed'
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {params.ref}  \
        -I {input.NC} \
        -V /public/ClinicalExam/lj_sih/resource/af-only-gnomad.raw.sites.hg19.vcf.gz \
        -L {params.bed} \
        -O {output.table}
        """


rule M2_contam:
    input:
        log='logs/paired/Mutect2_PoNDB_MM.log',
        T="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        NC="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        ref=config['resources'][genome_version]['REFFA'],
        temp_directory=config['params']['java']['temp_directory'],
        bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed'
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
        bed=get_sample_bed,
        pon="{project}/{genome_version}/results/vcf/paired/PoN/MM_pon.vcf.gz",
        pon_log='logs/paired/Mutect2_PoNVCF_MM.log'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
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
        --germline-resource /public/ClinicalExam/lj_sih/resource/af-only-gnomad.raw.sites.hg19.vcf.gz \
        --intervals {input.bed}
        """

rule M2_filter:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed,
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2.vcf",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Mutect2_filter.vcf"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
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
