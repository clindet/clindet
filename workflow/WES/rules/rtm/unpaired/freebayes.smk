rule unpaired_freebayes:
    input:
        ref=config['resources'][genome_version]['REFFA'],
        samples="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        regions=get_sample_bed
    output:
        "{project}/{genome_version}/results/vcf/unpaired/{sample}/freebayes.vcf"
    log:
        "logs/freebayes/unpaired/{sample}.log",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    wrapper:
        "v1.7.0/bio/freebayes"