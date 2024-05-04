rule rnaindel:
    input:
        bam="{project}/{genome_version}/results/summary/RSEM/{sample}/{sample}.STAR.genome.sort.bam",
        ref=config['resources'][genome_version]['REFFA'],
    output:
        vcf="{project}/{genome_version}/results/vcf/{sample}.vcf.gz"
    params:
        data_dir=config['singularity']['rnaindel']['sif']
    singularity: config['singularity']['rnaindel'][genome_version]['data_dir']
    threads:8
    shell:
        """
        rnaindel PredictIndels -i {input.bam} \
        -o {output.vcf} \
        -r {input.ref} -d {params.data_dit} -p {threads}
        """
