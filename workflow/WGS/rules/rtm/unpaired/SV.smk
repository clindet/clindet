rule SV_unp_delly:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/unpaired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
        """

rule delly_unp_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}_filter.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        bcftools filter 'FILTER="PASS"'  {input.vcf} > {output.vcf}
        """


rule SV_sansa_snp_annodelly:
    input:
        vcf="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        anno="{project}/{genome_version}/results/sv/unpaired/{sample}/SV_anno_{sample}.bcf",
        query="{project}/{genome_version}/results/sv/unpaired/{sample}/query_{sample}.tsv.gz"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        db=config['softwares']['sansa'][genome_version]['db'],
        g=config['softwares']['sansa'][genome_version]['g'],
        t=10000
    shell:
        """
        {config[softwares][sansa][call]} annotate -i Name  -g {params.g} -t {params.t} \
        -a {output.anno} -o {output.query} {input.vcf} 
        """


# rule msi: 
#     input:
#         bed=config['bed_top100'],
#         list=config['list'],
#         tumor='bam/tumor.sort.mkdup.bam',
#     output:
#         'msi/sample.msi.result'
#     params:
#         '-l 1 -p 1 -q 1 -s 1 -b'
#     threads:
#         16
#     log:
#         'logs/msi.log'
#     singularity:
#         '/data_sas/ypu/git_repository/HEMECDx/test_data/tumorOnly/msisensor_v0.6.sif'
#     shell:
#         "msisensor msi {params} {threads} -e {input.bed} -d {input.list} "
#         "-t {input.tumor} -o {output} 1>{log} 2>&1 "