# rule CNA_Batt:
#     input:
#         Tum="{project}/{genome_version}/results/recal_rename/paired/{sample}-T.bam",
#         NC="{project}/{genome_version}/results/recal_rename/paired/{sample}-NC.bam",
#     output:
#         rdata="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
#     params:
#         ref=config['resources']['hg19']['REFFA'],
#         igbed=config['singularity']['caveman']['ignorebed'],
#         out_dir='{project}/{genome_version}/results/subclone/paired/{sample}/battenberg'
#     threads: 8
#     singularity:
#         config['singularity']['cgpbattenberg']['sif']
#     script:
#         """
#         battenberg.pl -outdir {params.out_dir} \
#         -r {params.ref}.fai \
#         -tb {input.Tum} \
#         -nb {input.NC} \
#         -ge XY \
#         -impute-info impute_info.txt \
#         -thousand-genomes-loc 1000genomesloci \
#         -ignore-contigs-file GRCh38.exclude.contigs.txt \
#         -gc-correction-loc battenberg_wgs_gc_correction_1000g_v3 \
#         -species Human \
#         -assembly 37 \
#         -t {threads}

#         Rscript subclone.R -t MM-100 -n NC-100 --tb {input.Tum} \
#         --nb {input.NC} --sex Female -o MM-100
#         """

# rule CNA_snp_pileup_nor:
#     input:
#         Tum="{project}/{genome_version}/results/recal/paired/{sample}-T.bam",
#         # NC="{project}/{genome_version}/results/recal/paired/{sample}-NC.bam",
#     output:
#         n_pileup="{project}/{genome_version}/results/cnv/paired/{sample}/{sample}_ASCAT.rdata",
#     params:
#         wd="{project}/{genome_version}/results/cnv/paired/{sample}",
        
#         # gender=,
#         sample_index= lambda wildcards: wildcards.sample
#     threads: 8
#     script:
#         "snp-pileup -g -A -q15 -r20 -Q20 -P100 -v $DBSNP $OUTPUT_PILEUP $INPUT_BAM"

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