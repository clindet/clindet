rule lofreq_call_up:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/lofreq.vcf"
    threads: 10
    singularity: config['singularity']['lofreq']['sif']
    shell:
        """
        lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output.vcf} {input.Tum}
        """
