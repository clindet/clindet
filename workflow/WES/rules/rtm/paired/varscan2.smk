## varscan2 call #####
rule varscan2_mpileup:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        tumor="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        regions=get_sample_bed,
    output:
        normal=temp("{project}/{genome_version}/results/recal/paired/{sample}-NC.mp"),
        tumor=temp("{project}/{genome_version}/results/recal/paired/{sample}-T.mp")
    threads: 2
    conda:"samtools114"
    shell:
        """
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.normal} > {output.normal}
        samtools mpileup -q 1 -Q 1 -f {input.ref} -l {input.regions} {input.tumor} > {output.tumor}
        """

rule varscan2_call:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.mp",
        tumor="{project}/{genome_version}/results/recal/paired/{sample}-T.mp",
        regions=get_sample_bed,
    output:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.vcf"
    log:
        "{project}/{genome_version}/logs/varscan2/paired/{sample}.log",
    shell:
        """
        {config[softwares][varscan2][call]} somatic {input.normal} {input.tumor} --output-snp {output.snp} --output-indel {output.indel}
        """

rule varscan2_merge:
    input:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.snp.vcf",
        indel="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.indel.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.vcf",
        log="{project}/{genome_version}/logs/varscan2/paired/{sample}.log"
    threads: 1
    shell:
        """
        python scripts/varscan2vcf.py {input.snp} --sample_name {wildcards.sample} > {output.vcf}
        python scripts/varscan2vcf.py {input.indel} --sample_name {wildcards.sample} | grep -v '#' >> {output.vcf}
        touch {output.log}
        """

rule varscan2_filter_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/varscan2.vcf"
    output:
        vcf="{project}/{genome_version}/results/filter_vcf/paired/{sample}/varscan2.vcf"
    threads: 1
    params:
        caller='varscan2'
    script:
        """
        ../../scripts/vcf_filter_somtic.R
        """
