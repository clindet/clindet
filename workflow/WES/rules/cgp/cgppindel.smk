# # pindel.pl -noflag \
# # -reference /public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta \
# # -simrep /public/ClinicalExam/lj_sih/resource/mutFilter/simpleRepeats.bed.gz \
# # -genes /public/ClinicalExam/lj_sih/resource/mutFilter/hg19.exon.bed.gz \
# # -filter /public/ClinicalExam/lj_sih/resource/mutFilter/WES/human/hg19/cgpPindel/perl/rules/pulldownRules.lst \
# # -assembly GRCh37 \
# # -species Human \
# # -seqtype WXS \
# # -tumour /public/ClinicalExam/lj_sih/resource/mutFilter/hg19_FAKE.bam \
# # -normal /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/recal/paired/MM-015-NC.bam \
# # -outdir /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/pindel \
# # -cpus 20


# ### pindel_normal_pane
# rule PI_NP: 
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         NC="/public/ClinicalExam/lj_sih/resource/mutFilter/hg19_FAKE.bam",
#     output:
#         log="analysis/pindel_normal/log/{sample}_pindel_NC.log"
#     threads: 20
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#         out_dir='analysis/pindel_normal/{sample}',
#     singularity:
#         config['singularity']['cgppindel']['sif']
#         # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
#     shell:
#         """
#         pindel.pl -noflag \
#         -reference /public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta \
#         -simrep /public/ClinicalExam/lj_sih/resource/mutFilter/simpleRepeats.bed.gz \
#         -genes /public/ClinicalExam/lj_sih/resource/mutFilter/hg19.exon.bed.gz \
#         -assembly GRCh37 \
#         -species Human \
#         -seqtype WXS \
#         -tumour {input.NC} \
#         -normal {input.Tum} \
#         -outdir {params.out_dir} \
#         -cpus {threads} > {output.log}
#         """


# ### pindel_normal_pane
# rule PI_UM: 
#     input:
#         vcfs=expand("analysis/pindel_normal/log/{sample}_pindel_NC.log",sample = paired_samples)
#     output:
#         'analysis/normalPanel/pindel_{sample}.gff3.gz'
#     threads: 20
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#         gff3='analysis/normalPanel/pindel_{sample}'
#     singularity:
#         config['singularity']['cgppindel']['sif']
#         # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
#     shell:
#         """
#         pindel_np_from_vcf.pl -o {params.gff3} -samp_id NORMAL analysis/pindel_normal/*/*.vcf.gz        
#         """

# ### pindel_normal_pane
# rule PI_call: 
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         out_dir=directory('results/vcf/paired/{sample}/cgppindel'),
#         log='logs/paired/cgppindel_{sample}.log'
#     threads: 20
#     params:
#         ref=config['resources'][genome_version]['REFFA'],
#     singularity:
#         config['singularity']['cgppindel']['sif']
#         # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
#     shell:
#         """
#         pindel.pl \
#         -reference /public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta \
#         -simrep /public/ClinicalExam/lj_sih/resource/mutFilter/simpleRepeats.bed.gz \
#         -genes /public/ClinicalExam/lj_sih/resource/mutFilter/hg19.exon.bed.gz \
#         -exclude chrUn% \
#         -unmatched analysis/normalPanel/cgppindel.gff3.gz \
#         -filter /public/ClinicalExam/lj_sih/softwares/cgp/cgpPindel/perl/rules/pulldownRules.lst \
#         -softfil /public/ClinicalExam/lj_sih/softwares/cgp/cgpPindel/perl/rules/softRules.lst \
#         -assembly GRCh37 \
#         -species Human \
#         -seqtype WXS \
#         -tumour {input.Tum} \
#         -normal {input.NC} \
#         -outdir {output.out_dir} \
#         -cpus {threads} > {output.log}
#         """


# ### pindel_normal_pane
# rule PI_ggz: 
#     input:
#         log='logs/paired/cgppindel_{sample}.log',
#         germ_bed='results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.germline.bed'
#     output:
#         log='logs/paired/germline_bed_{sample}.log'
#     threads: 20
#     singularity:
#         config['singularity']['cgppindel']['sif']
#         # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
#     shell:
#         """
#         bgzip {input.germ_bed}
#         bgzip {input.germ_bed}.gz
#         tabix -p {input.germ_bed}.gz
#         touch {output.log}
#         """
# ## filter an format DP and AD tag
# rule cgppindel_filter_somatic:
#     input:
#         vcf='logs/paired/cgppindel_{sample}.log'
#     output:
#         vcf="{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel_filter.vcf"
#     threads: 1
#     params:
#         caller='cgppindel',
#         vcf='results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
#     script:
#         "../../scripts/vcf_filter_somtic.R"

# pindel.pl -noflag \
# -reference /public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta \
# -simrep /public/ClinicalExam/lj_sih/resource/mutFilter/simpleRepeats.bed.gz \
# -genes /public/ClinicalExam/lj_sih/resource/mutFilter/hg19.exon.bed.gz \
# -filter /public/ClinicalExam/lj_sih/resource/mutFilter/WES/human/hg19/cgpPindel/perl/rules/pulldownRules.lst \
# -assembly GRCh37 \
# -species Human \
# -seqtype WXS \
# -tumour /public/ClinicalExam/lj_sih/resource/mutFilter/hg19_FAKE.bam \
# -normal /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/recal/paired/MM-015-NC.bam \
# -outdir /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/pindel \
# -cpus 20


### pindel_normal_pane
rule PI_NP: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        NC="/public/ClinicalExam/lj_sih/resource/mutFilter/hg19_FAKE.bam",
    output:
        log="{project}/{genome_version}/analysis/pindel_normal/log/{sample}_pindel_NC.log"
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        out_dir='{project}/{genome_version}/analysis/pindel_normal/{sample}',
        simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
        genes=config['singularity']['cgppindel'][genome_version]['genes']
    singularity:
        config['singularity']['cgppindel']['sif']
        # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
    shell:
        """
        pindel.pl -noflag \
        -reference {params.ref} \
        -simrep {params.simrep} \
        -genes {params.genes} \
        -assembly {genome_version} \
        -species Human \
        -seqtype WXS \
        -tumour {input.NC} \
        -normal {input.Tum} \
        -outdir {params.out_dir} \
        -cpus {threads} > {output.log}
        """


### pindel_normal_pane
rule PI_UM: 
    input:
        vcfs=expand("{project}/{genome_version}/analysis/pindel_normal/log/{sample}_pindel_NC.log",sample = paired_samples,project = project,genome_version = genome_version)
    output:
        '{project}/{genome_version}/analysis/normalPanel/pindel_{sample}.gff3.gz'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gff3='analysis/normalPanel/pindel_{sample}'
    singularity:
        config['singularity']['cgppindel']['sif']
        # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
    shell:
        """
        pindel_np_from_vcf.pl -o {params.gff3} -samp_id NORMAL analysis/pindel_normal/*/*.vcf.gz        
        """

### pindel_normal_pane
rule PI_call: 
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        NP_gff3=config['singularity']['cgppindel'][genome_version]['WES']['normal_panel']
    output:
        out_dir=directory('{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel'),
        log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log'
    threads: 20
    params:
        ref=config['resources'][genome_version]['REFFA'],
        simrep=config['singularity']['cgppindel'][genome_version]['simrep'],
        genes=config['singularity']['cgppindel'][genome_version]['genes']
    singularity:
        config['singularity']['cgppindel']['sif']
        # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
    shell:
        """
        pindel.pl \
        -reference {params.ref} \
        -simrep /public/ClinicalExam/lj_sih/resource/mutFilter/simpleRepeats.bed.gz \
        -genes /public/ClinicalExam/lj_sih/resource/mutFilter/hg19.exon.bed.gz \
        -exclude chrUn% \
        -unmatched {input.NP_gff3} \
        -filter /public/ClinicalExam/lj_sih/softwares/cgp/cgpPindel/perl/rules/pulldownRules.lst \
        -softfil /public/ClinicalExam/lj_sih/softwares/cgp/cgpPindel/perl/rules/softRules.lst \
        -assembly GRCh37 \
        -species Human \
        -seqtype WXS \
        -tumour {input.Tum} \
        -normal {input.NC} \
        -outdir {output.out_dir} \
        -cpus {threads} > {output.log}
        """


### pindel_normal_pane
rule PI_ggz: 
    input:
        log='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log',
        germ_bed='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.germline.bed'
    output:
        log='{project}/{genome_version}/logs/paired/germline_bed_{sample}.log'
    threads: 20
    singularity:
        config['singularity']['cgppindel']['sif']
        # '/public/ClinicalExam/lj_sih/softwares/pindel.sif'
    shell:
        """
        bgzip {input.germ_bed}
        bgzip {input.germ_bed}.gz
        tabix -p {input.germ_bed}.gz
        touch {output.log}
        """
## filter an format DP and AD tag
rule cgppindel_filter_somatic:
    input:
        vcf='{project}/{genome_version}/results/logs/paired/cgppindel_{sample}.log'
    output:
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel_filter.vcf"
    threads: 1
    params:
        caller='cgppindel',
        vcf='{project}/{genome_version}/results/vcf/paired/{sample}/cgppindel/{sample}_T_vs_{sample}_NC.flagged.vcf.gz'
    script:
        "../../scripts/vcf_filter_somtic.R"
