#### delly workflow
rule SV_delly:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
        """

rule delly_filter:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        bcftools filter -i 'FILTER="PASS"'  {input.vcf} > {output.vcf}
        """

rule delly_pre_merge:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_delly.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        bcftools annotate -x INFO/END,INFO/SVMETHOD,INFO/CONSENSUS -i 'SVTYPE="BND"'  {input.vcf} > {output.vcf}
        """

rule SV_sansa_annodelly:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}_filter.vcf"
    output:
        anno="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_anno_{sample}.bcf",
        query="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/query_{sample}.tsv.gz"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        db=config['softwares']['sansa'][genome_version]['db'],
        g=config['softwares']['sansa'][genome_version]['g'],
        t=1000000
    shell:
        """
        {config[softwares][sansa][call]} annotate -i Name  -g {params.g} -t {params.t} \
        -a {output.anno} -o {output.query} {input.vcf} 
        """

#### BRASS work flow
rule SV_brass_bamstat:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas"
    singularity: config['singularity']['brass']['sif']
    shell:
        """
            bam_stats -i {input.Tum} -o {output.Tum}
            bam_stats -i {input.NC} -o {output.NC}
        """

rule brass_ascat:
    input:
        rdata="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}/{sample}_ASCATsc.rdata"
    output:
        ascat="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}.ascat"
    script:
        "../scripts/stats_ascat.R"



rule SV_brass:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        Tumbas="{project}/{genome_version}/results/recal/paired/{sample}-T.bam.bas",
        NCbas="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam.bas",
        ascat="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}.ascat"
    output:
        log="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/{sample}_brass.log"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        gc=config['singularity']['brass']['gc'],
        b=config['singularity']['brass']['b'],
        d=config['singularity']['brass']['d'],
        cb=config['singularity']['brass']['cb'],
        ct=config['singularity']['brass']['ct'],
        vi=config['singularity']['brass']['vi'],
        mi=config['singularity']['brass']['mi'],
        out_dir="{project}/{genome_version}/results/sv/paired/BRASS/{sample}"
    singularity: config['singularity']['brass']['sif']
    threads:20
    shell:
        """
        mkdir -p {params.out_dir}
        brass.pl -s human -as hg19 -pr WGS \
        -b {params.b} \
        -c {threads} -o {params.out_dir} \
        -d {params.d} -g {params.ref} \
        -gc {params.gc} -cb {params.cb} \
        -vi {params.vi} -ct {params.ct} -mi {params.mi}\
        -t {input.Tum} \
        -n {input.NC} -ss {input.ascat}
        touch {output.log}
        """

rule brass_pre_merge:
    input:
        log="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/{sample}_brass.log"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_brass.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/sv/paired/BRASS/{sample}/{sample}_T_vs_{sample}_NC.annot.vcf.gz"
    shell:
        """
        zcat {params.vcf} | bcftools annotate  -i 'SVTYPE="BND"' > {output.vcf}
        """

rule manta_pre_merge:
    input:
        tamp="{project}/{genome_version}/results/vcf/paired/{sample}/{sample}-Manta.log",
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_manta.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        vcf="{project}/{genome_version}/results/vcf/paired/{sample}/Manta/results/variants/somaticSV.vcf.gz"
    shell:
        """
        zcat {params.vcf}  | bcftools annotate -i 'SVTYPE="BND" && FILTER="PASS" ' > {output.vcf}
        """

#### gridss workflow
rule SV_gridss:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
    output:
        sv="{project}/{genome_version}/results/sv/paired/gridss/{sample}/SV_gridss_{sample}.vcf",
    params:
        ref=config['resources'][genome_version]['REFFA'],
    shell:
        """
        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
        """


rule SV_svaba:
    input:
        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
        dbsnp_indel=config['resources'][genome_version]['DBSNP_INDEL'],
    output:
        sv="{project}/{genome_version}/results/sv/paired/svaba/{sample}/{sample}.svaba.somatic.sv.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        wd="{project}/{genome_version}/results/sv/paired/svaba/{sample}"
    threads:20
    singularity: config['singularity']['svaba']['sif']
    shell:
        """
        svaba run -t {input.Tum} -n {input.NC} -p {threads} \
        -D {input.dbsnp_indel} -a {params.wd}/{wildcards.sample} -G {params.ref}
        """

rule jasmine_merge:
    input:
        delly="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_delly.vcf",
        brass="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_brass.vcf",
        manta="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_manta.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_merge.vcf"
    params:
        ref=config['resources'][genome_version]['REFFA']
    shell:
        """
        /public/ClinicalExam/lj_sih/softwares/Jasmine/jasmine --output_genotypes --allow_intrasample  --comma_filelist  file_list={input.brass},{input.manta},{input.delly} out_file={output.vcf}
        """
# echo 'export PATH=/public/ClinicalExam/lj_sih/softwares/go/bin:$PATH' >> ~/.bashrc && \
#   source ~/.bashrc
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

# rule merge_SV:
#     input:
#         rdata="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}/{sample}_ASCATsc.rdata"
#     output:
#         ascat="{project}/{genome_version}/results/cnv/ASCATsc/paired/{sample}.ascat"
#     script:
#         "../scripts/stats_ascat.R"


#rule SV_lumpy:
#    input:
#        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#    output:
        # sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
#    params:
#        ref=config['resources'][genome_version]['REFFA'],
#    shell:
#        """
#        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
#        
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.NC_spl}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.NC_dis}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.T_spl}
#        samtools fixmate -O bam - - | samtools sort -@ 30 -O bam -o {output.T_dis}
#       """


#rule SV_lumpy:
#    input:
#        Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#        NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#    output:
#        sv="{project}/{genome_version}/results/sv/paired/DELLY/{sample}/SV_delly_{sample}.vcf",
#    params:
#        ref=config['resources'][genome_version]['REFFA'],
#    shell:
#        """
#        {config[softwares][delly][call]} call -g {params.ref} {input.Tum} {input.NC} > {output}
#        """



rule SV_merge_final:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_merge.vcf"
    output:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_final.vcf"
    shell:
        """
        bcftools filter -i 'SUPP>=2' {input.vcf} > {output.vcf}
        """

rule SV_sansa_merge:
    input:
        vcf="{project}/{genome_version}/results/sv/paired/merge/{sample}/{sample}_final.vcf"
    output:
        anno="{project}/{genome_version}/results/sv/paired/merge/{sample}/SV_anno_{sample}.bcf",
        query="{project}/{genome_version}/results/sv/paired/merge/{sample}/final_{sample}.tsv.gz"
    params:
        ref=config['resources'][genome_version]['REFFA'],
        db=config['softwares']['sansa'][genome_version]['db'],
        g=config['softwares']['sansa'][genome_version]['g'],
        t=1000000
    shell:
        """
        {config[softwares][sansa][call]} annotate -i Name  -g {params.g} -t {params.t} \
        -a {output.anno} -o {output.query} {input.vcf}
        """