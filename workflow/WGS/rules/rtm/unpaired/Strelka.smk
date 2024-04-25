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

