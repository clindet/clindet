rule unpaired_call_variants_HaplotypeCaller:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        ref=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/HaplotypeCaller.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call']
    shell:
        """
        {params.gatk4} \
        HaplotypeCaller -R {input.ref} \
        -I {input.Tum} \
        -O {output.vcf} \
        --native-pair-hmm-threads 1 --annotate-with-num-discovered-alleles -A UniqueAltReadCount -A ReferenceBases \
        -A PossibleDeNovo -A Coverage -A DepthPerAlleleBySample -A DepthPerSampleHC -A StrandBiasBySample -A StrandOddsRatio &
        """