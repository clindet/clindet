library(R.utils)
library(glue)
library(maftools)
library(tidyverse)

genome_version <- snakemake@wildcards[['genome_version']]

input_tumor_bam <- snakemake@input[['Tum']]
input_tumor_bam <- getAbsolutePath(input_tumor_bam)

input_normal_bam <- snakemake@input[['NC']]
input_normal_bam <- getAbsolutePath(input_normal_bam)


output_tsv <- snakemake@output[['tsv']]
output_tsv <- getAbsolutePath(output_tsv)


output_rdata <- snakemake@params[['rdata']]
output_rdata <- getAbsolutePath(output_rdata)

threads = as.numeric(snakemake@threads)
print(threads)

# gender <- snakemake@params[['gender']]
gender = 'XX'
work_dir <- snakemake@params[['wd']]
sample_index <- snakemake@params[['sample_name']]


if(dir.exists(work_dir)){
} else {
  dir.creat(work_dir)
}
setwd(work_dir)


bams = c(
    input_tumor_bam,
    input_normal_bam
)


tryCatch({
    res = maftools::sampleSwaps(bams = bams, build = genome_version)
}, error = function(e){
    report <- data.frame(
        X_bam = basename(input_tumor_bam)  %>% str_remove('.bam'),
        Y_bam = basename(input_normal_bam)  %>% str_remove('.bam'),
        concordant_snps = 0,
        concordant_snps = 0,
        fract_concordant_snps = 0,
        cor_coef = 0,
        XY_possibly_paired = 'No'
    )
    write_tsv(report,output_tsv)
},finally = {
    report <- res$pairwise_comparison
    saveRDS(res, file = output_rdata)
    write_tsv(report,output_tsv)
})

