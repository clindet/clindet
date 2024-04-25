rule M2_ST:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        bed=config['resources'][genome_version]['WES_BED'],
        vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],#/public/ClinicalExam/lj_sih/resource/af-only-gnomad.raw.sites.hg19.vcf.gz
        ref=config['resources'][genome_version]['REFFA']
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory']
        # bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed',
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {input.ref}  \
        -I {input.Tum} \
        -V  {input.vcf} \
        -L {input.bed} \
        -O {output.table}
        """

rule M2_SNC:
    input:
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=config['resources'][genome_version]['WES_BED'],
        vcf=config['resources'][genome_version]['MUTECT2_germline_vcf']
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        table="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
        # bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed',
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GetPileupSummaries -R {input.ref}  \
        -I {input.NC} \
        -V {input.vcf} \
        -L {input.bed} \
        -O {output.table}
        """


rule M2_contam:
    input:
        T="{project}/{genome_version}/results/recal/{sample}/{sample}-T_pileupsummaries.table",
        NC="{project}/{genome_version}/results/recal/{sample}/{sample}-NC_pileupsummaries.table",
        ref=config['resources'][genome_version]['REFFA'],
        bed=config['resources'][genome_version]['WES_BED']
    output:
        # vcf="{project}/{genome_version}/results/vcf/paired/PoN/{sample}_pon.vcf.gz",
        seg="{project}/{genome_version}/results/recal/{sample}/{sample}_segments.table",
        ctam="{project}/{genome_version}/results/recal/{sample}/{sample}_calculatecontamination.table"
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory']
        # bed='/public/home/lijf/projects/clinmm/meta/bed/rnaseq/hg19/hg19.exon.bed'
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
        pon=config['resources'][genome_version]['WES_PON'],
        germ_vcf=config['resources'][genome_version]['MUTECT2_germline_vcf'],
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
        --germline-resource {input.germ_vcf} \
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
