### caveman_normal_pane
rule CM_cnv:
    input:
        rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    output:
        Tcnv="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_Tcnv.bed",
        NCcnv="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_NCcnv.bed"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 1
    script:
        "../../scripts/caveman/cnv_bed.R"

### caveman_normal_pane
rule CM_call: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        Tcnv="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_Tcnv.bed",
        NCcnv="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_NCcnv.bed"
    output:
        # out_dir=directory('results/vcf/paired/{sample}/caveman'),
        log='{project}/{genome_version}/logs/paired/caveman_{sample}.log'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        igbed=config['singularity']['caveman']['ignorebed'],
        out_dir='{project}/{genome_version}/results/vcf/paired/{sample}/caveman'
    singularity:
        config['singularity']['cgpwgs']['sif']
    shell:
        """
        caveman.pl \
        -outdir {params.out_dir} \
        -reference {params.ref}.fai \
        -tc {input.Tcnv} -nc {input.NCcnv} \
        -tumour-bam {input.Tum} -normal-bam {input.NC} \
        -ignore-file {params.igbed} \
        -s HUMAN -sa GRCh37 -seqType WXS -t {threads} -no-flagging
        touch {output.log}
        """


rule CM_flag: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        log='{project}/{genome_version}/logs/paired/caveman_{sample}.log'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/caveman.vcf"
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        c=config['singularity']['caveman'][genome_version]['flag']['c'],
        v=config['singularity']['caveman'][genome_version]['flag']['v'],
        u=config['singularity']['caveman'][genome_version]['flag']['u'],
        b=config['singularity']['caveman'][genome_version]['flag']['b'],
        ab=config['singularity']['caveman'][genome_version]['flag']['ab'],
        g=config['singularity']['caveman'][genome_version]['flag']['g'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/caveman/{sample}_T_vs_{sample}_NC.muts.ids.vcf.gz",
    singularity:
        config['singularity']['caveman']['sif']
    shell:
        """
        cgpFlagCaVEMan.pl \
        -i {params.vcf} \
        -o {output} \
        -m {input.Tum} -n {input.NC} \
        -s Human  -t genome \
        -ref {params.ref}.fai \
        -c  {params.c} \
        -v  {params.v} \
        -u  {params.u} \
        -g  {params.g} \
        -b  {params.b} \
        -ab {params.ab}
        """



rule CM_germ_flag: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        log='{project}/{genome_version}/logs/paired/caveman_{sample}.log'
    output:
        vcf="{project}/{genome_version}/results/vcf_germline/paired/{sample}/caveman.vcf"
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        c=config['singularity']['caveman'][genome_version]['flag']['c'],
        v=config['singularity']['caveman'][genome_version]['flag']['v'],
        u=config['singularity']['caveman'][genome_version]['flag']['u'],
        b=config['singularity']['caveman'][genome_version]['flag']['b'],
        ab=config['singularity']['caveman'][genome_version]['flag']['ab'],
        g=config['singularity']['caveman'][genome_version]['flag']['g'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/caveman/{sample}_T_vs_{sample}_NC.snps.ids.vcf.gz",
    singularity:
        config['singularity']['caveman']['sif']
    shell:
        """
        cgpFlagCaVEMan.pl \
        -i {params.vcf} \
        -o {output} \
        -m {input.Tum} -n {input.NC} \
        -s Human  -t genome \
        -ref {params.ref}.fai \
        -c  {params.c} \
        -v  {params.v} \
        -u  {params.u} \
        -g  {params.g} \
        -b  {params.b} \
        -ab {params.ab}
        """


# # rule CM_UM: 
# #     input:
# #         vcfs=expand("analysis/pindel_normal/log/{sample}_pindel_NC.log",sample = paired_samples)
# #     output:
# #         'analysis/normalPanel/pindel_{sample}.gff3.gz'
# #     threads: 20
# #     params:
# #         ref=config['resources'][genome_version]['REFFA'],
# #         gff3='analysis/normalPanel/pindel_{sample}',
# #     singularity:
# #         config['singularity']['caveman']['sif']
# #         # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
# #     shell:
# #         """  
# #         """