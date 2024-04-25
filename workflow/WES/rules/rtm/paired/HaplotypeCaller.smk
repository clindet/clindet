rule call_variants_HaplotypeCaller:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/HaplotypeCaller.vcf",
    params:
        gatk4=config['softwares']['gatk4']['call'],
        temp_directory=config['params']['java']['temp_directory'],
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        HaplotypeCaller -R {input.ref} \
        -I {input.Tum} \
        -I {input.NC} \
        -O {output.vcf} \
        --intervals {input.bed} \
        --native-pair-hmm-threads 1 --annotate-with-num-discovered-alleles -A UniqueAltReadCount -A ReferenceBases \
        -A PossibleDeNovo -A Coverage -A DepthPerAlleleBySample -A DepthPerSampleHC -A StrandBiasBySample -A StrandOddsRatio
        """

