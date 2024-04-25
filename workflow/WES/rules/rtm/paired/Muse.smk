rule muse_call:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        snp="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.MuSE.txt"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_prefix="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}",
    threads: 10
    singularity: config['singularity']['muse']['sif']
    shell:
        """
        /MuSE/bin/MuSE call -f {params.ref} -O {params.out_prefix} -n {threads} {input.Tum} {input.NC}
        """

rule muse_sump:
    input:
        txt="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}.MuSE.txt"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/muse.vcf"
    params:
        dbsnp=config['resources'][genome_version]['DBSNP_GZ']
    threads: 10
    singularity: config['singularity']['muse']['sif']
    shell:
        """
        /MuSE/bin/MuSE sump -I {input.txt} -O {output.vcf} -E -n {threads} -D {params.dbsnp}
        """
