rule CNA_ASCAT:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        #ascat_pre_dir= directory(config['resources'][genome_version]['ASCAT_WES_PreDIR'])
    output:
        rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        alleles_prefix="/public/ClinicalExam/lj_sih/projects/project_clindet/resources/human/hg19/WES/ASCAT_pre",
        loci_prefix="/public/ClinicalExam/lj_sih/projects/project_clindet/resources/human/hg19/WES/ASCAT_pre",
        # gender=,
        # GCcontentfile=,
        # replictimingfile=,
        sample_index = lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ASCAT.R"

rule CNA_ASCAT_sc:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
    output:
        rdata="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}/{sample}_ASCATsc.rdata"
    params:
        wd="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ASCATsc.R"
        
### PoN of factesCH
rule CNA_snp_pileup_nor:
    input:
        # Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        pileup="{project}/{genome_version}/results/cnv/db/pon_pileups/{sample}_pileup.tsv.gz"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        dbsnp=config['resources'][genome_version]['DBSNP'],
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 2
    shell:
        "snp-pileup -g -A -q15 -r20 -Q20 -P100 -v {params.dbsnp} {output.pileup} {input.NC}"

rule pon_facetsCH:
    input:
        expand("{project}/{genome_version}/results/cnv/db/pon_pileups/{MM_sample}_pileup.tsv.gz",MM_sample = paired_samples,project = project,genome_version = genome_version)
    output:
        log='logs/paired/facets_PoNDB_{sample}.log'
    conda: 'sankemake'
    params:
        wd="{project}/{genome_version}/results/cnv/db/pon_pileups"
    threads: 10
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && {params.gatk4} \
        GenomicsDBImport -R {params.ref} -L {params.bed} \
        --merge-input-intervals true --sites-only-vcf-output true \
        --genomicsdb-workspace-path /public/ClinicalExam/lj_sih/resource/mutFilter/pon_db \
        {params.vcfs}
        touch {output.log}
        """

### call facetsCH
rule CNA_snp_pileup_tum:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        # NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        pileup="{project}/{genome_version}/results/cnv/db/paired_facetsCH/{sample}/{sample}_pileup.tsv.gz",
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        dbsnp=config['resources'][genome_version]['DBSNP'],
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 2
    shell:
        "snp-pileup -g -A -q15 -r20 -Q20 -P100 -v {params.dbsnp} {output.pileup} {input.Tum}"

# rule CNA_facetsCH:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
#     params:
#         wd="{project}/{genome_version}/results/cnv/paired/{sample}",
#         # gender=,
#         sample_index= lambda wildcards: wildcards.sample
#     threads: 8
#     script:
#         "../scripts/ASCAT.R"

rule CNA_exomedepth:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        bed=get_sample_bed,
        ref=config['resources'][genome_version]['REFFA']
    output:
        rds="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_exomedepth.rds",
        tsv="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_exomedepth.tsv"
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ExomeDepth.R"

rule SM_check:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        tsv="{project}/{genome_version}/results/stats/paired/sample_check/{sample}/{sample}_check.tsv",
    params:
        wd="{project}/{genome_version}/results/stats/paired/sample_check/{sample}",
        rdata="{project}/{genome_version}/results/stats/paired/sample_check/{sample}/{sample}_check.rdata",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/sample_check.R"

rule CNA_ABSOLUTE_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata",
        maf="{project}/{genome_version}/results/maf/paired/{sample}/merge/{sample}_filter.maf"
    output:
        # rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
        seg="{project}/{genome_version}/results/cnv/GISTIC2/{sample}/{sample}.seg",
        ndt_seg="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.seg.txt",
        ndt_maf="{project}/{genome_version}/results/clone/PhylogicNDT/{sample}/{sample}.maf",
        absolute_dir=directory("{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}"),
        # absolute_pdf="{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}/DoAbsolute.called.ABSOLUTE.plots.pdf"
    # conda:'snakemake'
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ABSOLUTE.R"



rule CNA_Battenberg:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        # wd=directory("{project}/{genome_version}/results/cnv/Battenberg/{sample}"),
        sclo="{project}/{genome_version}/results/cnv/Battenberg/{sample}/{sample}_copynumber.txt"
    params:
        # wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        wd="{project}/{genome_version}/results/cnv/Battenberg/{sample}",
        script=Path(str(workflow.current_basedir) + '/../../../scripts/subclone.R'),
        # mkdir -p {params.wd}
        # gender=,
        sample_index = lambda wildcards: wildcards.sample
    threads: 8
    shell:
        """
        Rscript {params.script} -t {wildcards.sample}  -n {wildcards.sample}_NC \
        --tb {input.Tum} \
        --nb {input.NC} --sex Male \
        -o {params.wd}
        """

# rule freec_config:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#         #ascat_pre_dir= directory(config['resources'][genome_version]['ASCAT_WES_PreDIR'])
#     output:
#         rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
#     params:
#         wd="{project}/{genome_version}/results/cnv/paired/{sample}",
#         alleles_prefix="/public/ClinicalExam/lj_sih/projects/project_clindet/resources/human/hg19/WES/ASCAT_pre",
#         loci_prefix="/public/ClinicalExam/lj_sih/projects/project_clindet/resources/human/hg19/WES/ASCAT_pre",
#         # gender=,
#         sample_index = lambda wildcards: wildcards.sample
#     threads: 8
#     shell:
#         """
#         """


rule sequenza_bam2seqz:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        ref=config['resources'][genome_version]['REFFA'],
        gc=config['softwares']['sequenza'][genome_version]['gc']
        # gc="{project}/{genome_version}/genome/{genome_version}_gc50.wig.gz"
    output:
        seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.seqz.gz",
    threads: 8
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils bam2seqz \
        --normal {input.NC} \
        --tumor {input.Tum} \
        --fasta {input.ref} -gc {input.gc} \
        --output {output.seqz}
        """

rule sequenza_seqz_binning:
    input:
        seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}.seqz.gz"
    output:
        bin_seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}_bin50.seqz.gz",
    conda:
        config['softwares']['sequenza']['conda']
    shell:
        """
        sequenza-utils seqz_binning -w 50 --seqz {input.seqz} -o {output.bin_seqz}
        """


rule sequenza_call:
    input:
        bin_seqz="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}_bin50.seqz.gz",
    output:
        segment="{project}/{genome_version}/results/cnv/sequenza/{sample}/{sample}_segments.txt",
    params:
        wd="{project}/{genome_version}/results/cnv/sequenza/{sample}",
    script:
        "../../../scripts/sequenza.R"

rule ASCAT_GISTIC:
    input:
        cnv_rdata="{project}/{genome_version}/results/cnv/paired/ascat/{sample}/{sample}_ASCAT.rdata"
    output:
        # rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
        seg="{project}/{genome_version}/results/cnv/GISTIC2_seg/{sample}/{sample}.seg"
        # absolute_pdf="{project}/{genome_version}/results/cnv/ABSOLUTE/{sample}/DoAbsolute.called.ABSOLUTE.plots.pdf"
    # conda:'snakemake'
    params:
        wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        # gender=,
        sample_index= lambda wildcards: wildcards.sample
    threads: 8
    script:
        "../../../scripts/ASCAT_to_GISTIC2.R"