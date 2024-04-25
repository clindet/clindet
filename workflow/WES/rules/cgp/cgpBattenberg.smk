### caveman_normal_pane
rule battenberg_call: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam"
    output:
        # out_dir=directory('results/vcf/paired/{sample}/caveman'),
        log='logs/paired/cgpbattenberg_{sample}.log'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        igbed=config['singularity']['cgpbattenberg']['ignorebed'],
        gc=config['singularity']['cgpbattenberg']['gc'],
        impute=config['singularity']['cgpbattenberg']['impute'],
        prob_loci=config['singularity']['cgpbattenberg']['prob_loci'],
        loci_1000=config['singularity']['cgpbattenberg']['loci_1000'],
        out_dir='results/vcf/paired/{sample}/battenberg',
        log='results/vcf/paired/{sample}/battenberg/log'
        # gender=get_gender
    singularity:
        config['singularity']['cgpbattenberg']['sif']
    shell:
        """
        battenberg.pl -outdir {params.out_dir} \
        -r {params.ref}.fai \
        -tb {input.Tum} -prob-loci {params.prob_loci} \
        -nb {input.NC} \
        -ge XY \
        -impute-info {params.impute} \
        -thousand-genomes-loc  {params.loci_1000}\
        -ignore-contigs-file {params.igbed} \
        -gc-correction-loc {params.gc} \
        -species Human \
        -assembly 37 \
        -t {threads}
        touch {output.log}
        """