rule call_config_strelka:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/paired/{sample}/Manta"),
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-Manta.log",
        indelcd="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/candidateSmallIndels.vcf.gz"
    params:
    conda:
        "strelka"
    shell:
        """
            configManta.py \
            --tumorBam {input.Tum} \
            --bam {input.NC} \
            --runDir {output.dir} \
            --exome \
            --callRegions {input.bed}.gz \
            --referenceFasta {input.ref}
            {output.dir}/runWorkflow.py
            touch {output.tamp}
        """

# for Manta candidate
rule call_strelka_manta:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        indelcd="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/candidateSmallIndels.vcf.gz",
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaManta"),
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaManta.log"
    params:
    threads:8
    conda:
        "strelka"
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.Tum} \
        --bam {input.NC} \
        --referenceFasta {input.ref} \
        --indelCandidates {input.indelcd} \
        --exome \
        --runDir {output.dir} \
        --callRegions {input.bed}.gz
        {output.dir}/runWorkflow.py -m local -j {threads}
        touch {output.tamp}
        """

rule call_strelka_somatic_manta:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        indelcd="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/candidateSmallIndels.vcf.gz",
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomaticManta"),
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaSomaticeManta.log"
    params:
    threads:8
    conda:
        "strelka"
    shell:
        """
        configureStrelkaSomaticWorkflow.py \
        --tumorBam {input.Tum} \
        --normalBam {input.NC} \
        --indelCandidates {input.indelcd} \
        --referenceFasta {input.ref} \
        --runDir {output.dir} \
        --callRegions {input.bed}.gz
        {output.dir}/runWorkflow.py -m local -j {threads}
        touch {output.tamp}
        """

rule merge_strelka_manta:
    input:
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaManta.log"
    output:
        "{project}/{genome_version}/results/vcf_germline/paired/{sample}/strelkamanta.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaManta/results/variants/variants.vcf.gz"
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | {{ grep 'PASS' || true;}} >> {output}
        """
    # or the uppon line can be rewrite as :zcat {params.snp} | grep -v '#' | grep 'PASS'| cat >> {output}

rule merge_strelka_somatic_manta:
    input:
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaSomaticeManta.log"
    output:
        "{project}/{genome_version}/results/vcf/paired/{sample}/strelkasomaticmanta.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomaticManta/results/variants/somatic.snvs.vcf.gz",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomaticManta/results/variants/somatic.indels.vcf.gz",
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | {{ grep 'PASS' || true; }} >>  {output}
        zcat {params.indel} | grep -v '#' | {{ grep 'PASS' || true; }} >>  {output}
        sed -i "s;TUMOR;{wildcards.sample}_T;"  {output}
        sed -i "s;NORMAL;{wildcards.sample}_NC;"  {output}
        """
### pure strelka
rule call_strelka:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/paired/{sample}/Strelka"),
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-Strelka.log"
    params:
    conda:
        "strelka"
    shell:
        """
        configureStrelkaGermlineWorkflow.py \
        --bam {input.Tum} \
        --bam {input.NC} \
        --referenceFasta {input.ref} \
        --exome \
        --runDir {output.dir} \
        --callRegions {input.bed}.gz
        {output.dir}/runWorkflow.py -m local 
        touch {output.tamp}
        """

rule call_strelka_somatic:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        dir=directory("{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomatic"),
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaSomatice.log"
    params:
    conda:
        "strelka"
    shell:
        """
        configureStrelkaSomaticWorkflow.py \
        --tumorBam {input.Tum} \
        --normalBam {input.NC} \
        --referenceFasta {input.ref} \
        --runDir {output.dir} \
        --callRegions {input.bed}.gz
        {output.dir}/runWorkflow.py -m local 
        touch {output.tamp}
        """

rule merge_strelka:
    input:
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-Strelka.log"
    output:
        "{project}/{genome_version}/results/vcf/paired/{sample}/strelka.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/Strelka/results/variants/variants.vcf.gz"
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | {{ grep 'PASS' || true;}} >> {output}
        """
    # or the uppon line can be rewrite as :zcat {params.snp} | grep -v '#' | grep 'PASS'| cat >> {output}

rule merge_strelka_somatic:
    input:
        tamp="{project}/{genome_version}/results/logs/strelka/paired/{sample}-StrelkaSomatice.log"
    output:
        "{project}/{genome_version}/results/vcf/paired/{sample}/strelkasomatic.vcf"
    params:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomatic/results/variants/somatic.snvs.vcf.gz",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/StrelkaSomatic/results/variants/somatic.indels.vcf.gz",
    conda:
        "strelka"
    shell:
        """
        zcat {params.snp} | grep '#' > {output}
        zcat {params.snp} | grep -v '#' | {{ grep 'PASS' || true; }} >>  {output}
        zcat {params.indel} | grep -v '#' | {{ grep 'PASS' || true; }} >>  {output}
        sed -i "s;TUMOR;{wildcards.sample}_T;"  {output}
        sed -i "s;NORMAL;{wildcards.sample}_NC;"  {output}
        """
