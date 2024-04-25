rule deepvariant_call:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.vcf",
        gvcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant/{sample}.deepvariant.gvcf",
        # dir="{project}/{genome_version}/results/vcf/paired/{sample}/deepvariant,
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/unpaired/{sample}/deepvariant",
    threads: 10
    singularity: config['singularity']['deepvariant']['sif']
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref={input.reference} \
        --reads={input.Tum} \
        --output_vcf={output.vcf} \
        --output_gvcf={output.gvcf} \
        --intermediate_results_dir {params.out_prefix} \
        --num_shards={threads}
        """
