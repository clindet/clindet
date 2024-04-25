rule unpaired_call_config_strelka:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/unpaired/{sample}/Manta"),
        tamp="{project}/{genome_version}/results/vcf/unpaired/{sample}-Manta.log"
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
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/unpaired/{sample}/Strelka"),
        tamp="{project}/{genome_version}/results/vcf/unpaired/{sample}-Strelka.log"
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
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/unpaired/{sample}/StrelkaSomatic"),
        tamp="{project}/{genome_version}/results/vcf/unpaired/{sample}-StrelkaSomatice.log"
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
        tamp="{project}/{genome_version}/results/vcf/unpaired/{sample}-Strelka.log"
    output:
        "{project}/{genome_version}/results/vcf/unpaired/{sample}/strelka.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/unpaired/{sample}/Strelka/results/variants/variants.vcf.gz"
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | grep 'PASS' >> {output}
        """

rule unpaired_merge_strelka_somatic:
    input:
        tamp="{project}/{genome_version}/results/vcf/unpaired/{sample}-StrelkaSomatice.log"
    output:
        "{project}/{genome_version}/results/vcf/unpaired/{sample}/strelkasomatic.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/unpaired/{sample}/StrelkaSomatic/results/variants/somatic.snvs.vcf.gz",
        indel="{project}/{genome_version}/results/vcf/unpaired/{sample}/StrelkaSomatic/results/variants/somatic.indels.vcf.gz",
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

