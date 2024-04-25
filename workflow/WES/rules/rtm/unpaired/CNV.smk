### purple
### amber first

rule unpaired_amber:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        pcf="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber/{sample}.amber.baf.pcf",
        output_dir=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber")
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        amber  -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -target_regions_bed {input.bed} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -loci {params.loci}
        """

rule unpaired_cobalt:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        bed=get_sample_bed
    output:
        pcf="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt/{sample}.cobalt.ratio.pcf",
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt"
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt",
        tumor_only_diploid_bed=config['singularity']['hmftools'][genome_version]['cobalt']['tumor_only_diploid_bed'],
        gc_profile=config['singularity']['hmftools'][genome_version]['cobalt']['gc_profile'],
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        cobalt -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -output_dir {params.output_dir} \
        -threads {threads} \
        -pcf_gamma 50 \
        -tumor_only_diploid_bed {params.tumor_only_diploid_bed} \
        -gc_profile {params.gc_profile}
        """

rule unpaired_purple:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        amber=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber"),
        cobalt=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/cobalt"),
        sage_vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/{sample}.sage.vcf.gz",
        ref_genome=config['resources'][genome_version]['REFFA'],
        bed=get_sample_bed
    output:
        qc="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple/{sample}.purple.qc",
        output_dir=directory("{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/purple")
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        tumor_only_diploid_bed=config['singularity']['hmftools'][genome_version]['purple']['tumor_only_diploid_bed'],
        gc_profile=config['singularity']['hmftools'][genome_version]['purple']['gc_profile'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['purple']['ensembl_data_dir'],
        somatic_hotspots=config['singularity']['hmftools'][genome_version]['purple']['somatic_hotspots'],
        driver_gene_panel=config['singularity']['hmftools'][genome_version]['purple']['driver_gene_panel'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['purple']['ref_genome_version'],
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        purple \
        -tumor {wildcards.sample} \
        -amber {input.amber} \
        -cobalt {input.cobalt} \
        -target_regions_bed {input.bed} \
        -gc_profile {params.gc_profile} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -somatic_vcf {input.sage_vcf} \
        -somatic_hotspots {params.somatic_hotspots} \
        -driver_gene_panel {params.driver_gene_panel} \
        -circos /opt/circos-0.69-2/bin/circos \
        -output_dir {output.output_dir}
        """
# amber  -tumor HE13_YXY_D -tumor_bam /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/recal/unpaired/HE13_YXY_D-T.bam \
# -output_dir /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/cnv/unpaired/purple/HE13_YXY_D/amber \
# -threads 16  \
# -target_regions_bed /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/target_test.bed \
# -loci /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/GermlineHetPon.37.vcf.gz 


# cobalt \
# -tumor HE13_YXY_D \
# -tumor_bam /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/recal/unpaired/HE13_YXY_D-T.bam \
# -output_dir /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/cnv/unpaired/purple/HE13_YXY_D/cobalt \
# -threads 16 \
# -target_regions_bed /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/target_test.bed \
# -tumor_only_diploid_bed /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/copy_number/DiploidRegions.37.bed.gz \
# -gc_profile /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/copy_number/GC_profile.1000bp.37.cnp


# purple \
#    -tumor HE13_YXY_D \
#    -amber /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/cnv/unpaired/purple/HE13_YXY_D/amber \
#    -cobalt /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/cnv/unpaired/purple/HE13_YXY_D/cobalt \
#    -gc_profile /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/copy_number/GC_profile.1000bp.37.cnp \
#    -ref_genome_version 37 \
#    -ref_genome /public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/Homo_sapiens_assembly19.fasta \
#    -ensembl_data_dir /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/common/ensembl_data \
#    -threads 16 \
#    -target_regions_bed /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/target_test.bed \
#    -somatic_vcf /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/vcf/unpaired/HE13_YXY_D/HE13_YXY_D.sage.vcf.gz \
#    -somatic_hotspots /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/variants/KnownHotspots.somatic.37.vcf.gz \
#    -driver_gene_panel /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/common/DriverGenePanel.37.tsv \
#    -circos /opt/circos-0.69-2/bin/circos \
#    -output_dir /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/cnv/unpaired/purple/HE13_YXY_D/purple

### control freec



#### cnvkit


### exon
