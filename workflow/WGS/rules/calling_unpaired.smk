rule unpaired_coverageBed:
    input:
        a=get_sample_bed,
        b="results/recal/unpaired/{sample}-{group}.bam"
    output:
        bam="results/stats/unpaired/{sample}-{group}.cov"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "v1.7.0/bio/bedtools/coveragebed"

rule unpaired_coverageBed_2:
    input:
        a=get_sample_bed,
        b="results/mapped/unpaired/{sample}-{group}.sorted.bam"
    output:
        bam="results/stats_raw/{sample}-{group}.cov"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "v1.7.0/bio/bedtools/coveragebed"

rule unpaired_call_variants_mutect2:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        vcf="results/vcf/unpaired/{sample}/Mutect2.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call']
    threads: 8
    shell:
        """
        {params.gatk4} \
        Mutect2 -R {input.ref} \
        -I {input.Tum} \
        -O {output.vcf} \
        --intervals {input.bed} 
        """


rule unpaired_call_variants_HaplotypeCaller:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        vcf="results/vcf/unpaired/{sample}/HaplotypeCaller.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call']
    shell:
        """
        {params.gatk4} \
        HaplotypeCaller -R {input.ref} \
        -I {input.Tum} \
        -O {output.vcf} \
        --intervals {input.bed} \
        --native-pair-hmm-threads 1 --annotate-with-num-discovered-alleles -A UniqueAltReadCount -A ReferenceBases \
        -A PossibleDeNovo -A Coverage -A DepthPerAlleleBySample -A DepthPerSampleHC -A StrandBiasBySample -A StrandOddsRatio &
        """

rule unpaired_call_variants_UnifiedGenoTyper:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        vcf="results/vcf/unpaired/{sample}/UnifiedGenoTyper.vcf",
    params:
        gatk3=config['softwares']['gatk3']['call']
    shell:
        """
        {params.gatk3} \
        --unsafe -T UnifiedGenotyper -R {input.ref} \
        -I {input.Tum} \
        -o {output.vcf} \
        --intervals {input.bed}  -nt 24 -dcov 5000 --unsafe -glm BOTH
        """

rule unpaired_call_config_strelka:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("results/vcf/unpaired/{sample}/Manta"),
        tamp="results/vcf/unpaired/{sample}-Manta.log"
    params:
    conda:
        "strelka"
    shell:
        """
            configManta.py \
            --tumorBam {input.Tum} \
            --runDir {output.dir} \
            --exome \
            --callRegions {input.bed} \
            --referenceFasta {input.ref}
            {output.dir}/runWorkflow.py
            touch {output.tamp}
        """

rule unpaired_call_strelka:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("results/vcf/unpaired/{sample}/Strelka"),
        tamp="results/vcf/unpaired/{sample}-Strelka.log"
    params:
    conda:
        "strelka"
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.Tum} \
        --referenceFasta {input.ref} \
        --exome \
        --runDir {output.dir} \
        --callRegions {input.bed}
        {output.dir}/runWorkflow.py -m local 
        touch {output.tamp}
        """

rule unpaired_call_strelka_somatic:
    input:
        Tum="results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("results/vcf/unpaired/{sample}/StrelkaSomatic"),
        tamp="results/vcf/unpaired/{sample}-StrelkaSomatice.log"
    params:
    conda:
        "strelka"
    shell:
        """
        configureStrelkaSomaticWorkflow.py \
        --tumorBam {input.Tum} \
        --referenceFasta {input.ref} \
        --runDir {output.dir} \
        --callRegions {input.bed}
        {output.dir}/runWorkflow.py -m local 
        touch {output.tamp}
        """

rule unpaired_merge_strelka:
    input:
        tamp="results/vcf/unpaired/{sample}-Strelka.log"
    output:
        "results/vcf/unpaired/{sample}/strelka.vcf"
    params:
        snp="results/vcf/unpaired/{sample}/Strelka/results/variants/variants.vcf.gz"
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | grep 'PASS' >> {output}
        """

rule unpaired_merge_strelka_somatic:
    input:
        tamp="results/vcf/unpaired/{sample}-StrelkaSomatice.log"
    output:
        "results/vcf/unpaired/{sample}/strelkasomatic.vcf"
    params:
        snp="results/vcf/unpaired/{sample}/StrelkaSomatic/results/variants/somatic.snvs.vcf.gz",
        indel="results/vcf/unpaired/{sample}/StrelkaSomatic/results/variants/somatic.indels.vcf.gz",
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | grep 'PASS' >>  {output}
        zcat {params.indel} | grep -v '#' | grep 'PASS' >>  {output}
        sed -i "s;TUMOR;{wildcards.sample}_T;"  {output}
        sed -i "s;NORMAL;{wildcards.sample}_NC;"  {output}
        """


rule unpaired_vardict_single_mode:
    input:
        reference=config['resources']['hg19']['REFFA'],
        regions=get_sample_bed,
        bam="results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="results/vcf/unpaired/{sample}/vardict.vcf"
    params:
        extra="",
        bed_columns="-c 1 -S 2 -E 3 -g 4",  # Optional, default is -c 1 -S 2 -E 3 -g 4
        allele_frequency_threshold="0.01",  # Optional, default is 0.01
    threads: 1
    log:
        "logs/unpaired/varscan_{sample}_s_.log",
    wrapper:
        "v1.10.0/bio/vardict"

rule unpaired_freebayes:
    input:
        ref=config['resources']['hg19']['REFFA'],
        samples="results/recal/unpaired/{sample}-T.bam",
        indexes="results/recal/unpaired/{sample}-T.bam.bai",
        regions=get_sample_bed
    output:
        "results/vcf/unpaired/{sample}/freebayes.vcf"
    log:
        "logs/freebayes/unpaired/{sample}.log",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    wrapper:
        "v1.7.0/bio/freebayes"