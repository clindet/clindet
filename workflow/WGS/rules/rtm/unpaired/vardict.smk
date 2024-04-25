rule unpaired_vardict_single_mode:
    input:
        reference=config['resources']['hg19']['REFFA'],
        regions=get_sample_bed,
        bam="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/vardict.vcf"
    params:
        extra="",
        bed_columns="-c 1 -S 2 -E 3 -g 4",  # Optional, default is -c 1 -S 2 -E 3 -g 4
        allele_frequency_threshold="0.01",  # Optional, default is 0.01
    threads: 1
    log:
        "logs/unpaired/varscan_{sample}_s_.log",
    wrapper:
        "v1.10.0/bio/vardict"