rule unpaired_sage:
    input:
        Tum="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam",
        indexes="{project}/{genome_version}/results/recal/unpaired/{sample}-T.bam.bai",
        ref_genome=config['resources'][genome_version]['REFFA']
    output:
        vcf="{project}/{genome_version}/results/vcf/unpaired/{sample}/{sample}.sage.vcf.gz",
    params:
        output_dir="{project}/{genome_version}/results/cnv/unpaired/purple/{sample}/amber",
        loci=config['singularity']['hmftools'][genome_version]['amber']['loci'],
        high_confidence_bed=config['singularity']['hmftools'][genome_version]['sage']['high_confidence_bed'],
        ensembl_data_dir=config['singularity']['hmftools'][genome_version]['sage']['ensembl_data_dir'],
        coverage_bed=config['singularity']['hmftools'][genome_version]['sage']['coverage_bed'],
        hotspots=config['singularity']['hmftools'][genome_version]['sage']['hotspots'],
        panel_bed=config['singularity']['hmftools'][genome_version]['sage']['panel_bed'],
        ref_genome_version=config['singularity']['hmftools'][genome_version]['sage']['ref_genome_version'],
    threads: 10
    singularity:config['singularity']['hmftools']['sif']
    shell:
        """
        sage \
        -tumor {wildcards.sample} -tumor_bam {input.Tum} \
        -ref_genome_version {params.ref_genome_version} \
        -ref_genome {input.ref_genome} \
        -ensembl_data_dir {params.ensembl_data_dir} \
        -threads {threads} \
        -coverage_bed {params.coverage_bed} \
        -hotspots  {params.hotspots} \
        -panel_bed {params.panel_bed} \
        -high_confidence_bed {params.high_confidence_bed} \
        -output_vcf {output.vcf}
        """


# sage \
#     -tumor HE13_YXY_D -tumor_bam /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/recal/unpaired/HE13_YXY_D-T.bam \
#     -ref_genome_version 37 \
#     -ref_genome /public/ClinicalExam/lj_sih/projects/project_clindet/reference/b37/Homo_sapiens_assembly19.fasta \
#     -ensembl_data_dir /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/common/ensembl_data \
#     -threads 16 -write_bqr_plot \
#     -coverage_bed /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/variants/CoverageCodingPanel.37.bed.gz \
#     -hotspots /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/variants/KnownHotspots.somatic.37.vcf.gz \
#     -panel_bed /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/variants/ActionableCodingPanel.37.bed.gz \
#     -high_confidence_bed /public/ClinicalExam/lj_sih/projects/project_clindet/data/pipeline/v5_34/ref/37/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz \
#     -output_vcf /public/ClinicalExam/lj_sih/projects/project_pipeline/WES/Breast/b37/results/vcf/unpaired/HE13_YXY_D/HE13_YXY_D.sage.vcf.gz \
#     -vis_output_dir vis_mut -vis_pass_only

