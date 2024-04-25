rule deepvariant_call:
    input:
        # reference=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant/{sample}.deepvariant.vcf",
        gvcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant/{sample}.deepvariant.gvcf",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant",
    threads: 10
    singularity: config['singularity']['deepvariant']['sif']
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref={params.ref} \
        --reads={input.Tum} \
        --output_vcf={output.vcf} \
        --output_gvcf={output.gvcf} \
        --intermediate_results_dir {params.out_prefix} \
        --num_shards={threads}
        """

rule deepvariant_somatic_call:
    input:
        # reference=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic/{sample}.deepvariant.vcf"
        # gvcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic/{sample}.deepvariant.gvcf",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic",
    threads: 10
    singularity: config['singularity']['deepvariant_somatic']['sif']
    shell:
        """
        /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=WGS \
        --ref={params.ref} \
        --reads_normal={input.NC} \
        --reads_tumor={input.Tum} \
        --output_vcf={output.vcf} \
        --sample_name_tumor={wildcards.sample}_T \
        --sample_name_normal={wildcards.sample}_NC \
        --num_shards={threads} \
        --logging_dir={params.out_prefix}/logs \
        --intermediate_results_dir  {params.out_prefix} 
        """

rule deepvariant_filter_germline:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic/{sample}.deepvariant.vcf",
        # gvcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic/{sample}.deepvariant.gvcf",
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/deepvariant.vcf"
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant",
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="GERMLINE"'  {input.vcf} > {output.vcf} 
        """

rule deepvariant_filter_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant_somatic/{sample}.deepvariant.vcf",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant.vcf"
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant",
    threads: 1
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf} 
        """