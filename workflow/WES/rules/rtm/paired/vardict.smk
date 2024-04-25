# rule vardict_paired_mode:
#     input:
#         reference=config['resources'][genome_version]['REFFA'],
#         regions=get_sample_bed,
#         bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
#     params:
#         extra="",
#     threads: 10
#     log:
#         "{project}/{genome_version}/logs/paired/varscan_{sample}_tn.log",
#     wrapper:
#         "v1.7.0/bio/vardict"


rule vardict_paired_mode:
    input:
        reference=config['resources'][genome_version]['REFFA'],
        regions=get_sample_bed,
        bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    params:
        extra="",
    threads: 10
    # log:
    #     "{project}/{genome_version}/logs/paired/varscan_{sample}_tn.log",
    wrapper:
        "v1.7.0/bio/vardict"

# rule vardict_paired_mode:
#     input:
#         reference=config['resources'][genome_version]['REFFA'],
#         regions=get_sample_bed,
#         bam="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         normal="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
#     params:
#         extra="",
#     threads: 10
#     # log:
#     #     "{project}/{genome_version}/logs/paired/varscan_{sample}_tn.log",
#     wrapper:
#         "v1.7.0/bio/vardict"


rule vardict_filter_somatic:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict_filter.vcf"
    threads: 1
    params:
        caller='vardict'
    script:
        "../../../scripts/vcf_filter_somtic.R"

rule vardict_filter_germline:
    input:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/vardict.vcf"
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/vardict_germline.vcf"
    threads: 1
    params:
        caller='vardict'
    script:
        "../../../scripts/vcf_filter_germline.R"