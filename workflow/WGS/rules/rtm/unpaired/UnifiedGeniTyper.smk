rule unpaired_call_variants_UnifiedGenoTyper:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/UnifiedGenoTyper.vcf",
    params:
        gatk3=config['softwares']['gatk3']['call']
    shell:
        """
        {params.gatk3} \
        --unsafe -T UnifiedGenotyper -R {input.ref} \
        -I {input.Tum} \
        -o {output.vcf} \
        --intervals {input.bed}  -nt 24 -dcov 5000 --unsafe -glm BOTH
        """