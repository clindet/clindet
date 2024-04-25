library(R.utils)

input_tumor_bam <- snakemake@input[['Tum']]
input_tumor_bam <- getAbsolutePath(input_tumor_bam)

input_normal_bam <- snakemake@input[['NC']]
input_normal_bam <- getAbsolutePath(input_normal_bam)

output_rdata <- snakemake@output[['rdata']]
output_rdata <- getAbsolutePath(output_rdata)

threads = as.numeric(snakemake@threads)
print(threads)

# gender <- snakemake@params[['gender']]
gender = 'XX'
work_dir <- snakemake@params[['wd']]
sample_index <- snakemake@params[['sample_name']]
####################################################################
library(ASCAT.sc)
####################################################################

####################################################################
# ref_fasta <- "/public/home/lijf/env/genome/broad/hg19/ucsc.hg19.fasta"

####################################################################
setwd(work_dir)
####################################################################

####################################################################
res <- run_sc_sequencing(tumour_bams=input_tumor_bam,
                         allchr=c(1:22), ## might need a "chr" instead of "" if hg38
                         sex="male",
                         chrstring_bam="chr",  ## might need a "chr" instead of "" if hg38
                         purs = seq(0.1, 1, 0.01) , ## purity grid values
                         ploidies = seq(1.7,3, 0.01), ## average ploidy grid values
                         maxtumourpsi=5, ##maximum tumour ploidy
                         binsize=500000, ## bin size - reduce if enough sequencing reads (look at dpb (depth per bin) value in plots can go down to 100bp or even lower)
                         build="hg19", ## either hg19 or hg38 so far
                         projectname=sample_index,
                         MC.CORES=threads, ##number of cores available
                         multipcf=FALSE) ##use multipcf for multi-track segmentation if multi-sample sequencing
####################################################################
save(res,file=output_rdata) ## save results for reuse later on
####################################################################