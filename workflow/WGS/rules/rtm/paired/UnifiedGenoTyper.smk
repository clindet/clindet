rule call_variants_UnifiedGenoTyper:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources']['hg19']['REFFA'],
        bed=get_sample_bed
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/UnifiedGenoTyper.vcf",
    params:
        gatk3=config['softwares']['gatk3']['call'],
        #temp_directory=config['params']['java']['temp_directory'],
        temp_directory="{project}/{genome_version}/results/tmp/paired/{sample}/"
    resources:
        tmpdir="{project}/{genome_version}/results/tmp/paired/{sample}/"
    shell:
        """
        mkdir -p {params.temp_directory}
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk3} \
        --unsafe -T UnifiedGenotyper -R {input.ref} \
        -I {input.Tum} \
        -I {input.NC} \
        -o {output.vcf} \
        --intervals {input.bed}  -nt 24 -dcov 5000 --unsafe -glm BOTH
        """