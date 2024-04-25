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

if(dir.exists(work_dir)){
} else {
  dir.creat(work_dir)
}
setwd(work_dir)

library(ASCAT)
# input_tumor_bam = '/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/recal/paired/MM-014-T.bam'
# input_normal_bam = '/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/recal/paired/MM-014-NC.bam'
ascat.prepareHTS(
  tumourseqfile = input_tumor_bam,
  normalseqfile = input_normal_bam,
  tumourname = "Tumor",
  normalname = "Germline",
  allelecounter_exe = "/public/ClinicalExam/lj_sih/softwares/alleleCount/bin/alleleCounter",
  alleles.prefix = "/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/MMWES/hg19/results/cnv/ASCAT_pre/alleleData/Cleaned/alleleData_chr",
  loci.prefix = "/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/MMWES/hg19/results/cnv/ASCAT_pre/alleleData/Cleaned/loci_chr",
  gender = gender,
  genomeVersion = "hg19",
  nthreads = threads,
  chrom_names = c(1:22),
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")
  

ascat.bc = ascat.loadData(
  Tumor_LogR_file = "Tumor_LogR.txt", 
  Tumor_BAF_file = "Tumor_BAF.txt", 
  Germline_LogR_file = "Germline_LogR.txt", 
  Germline_BAF_file = "Germline_BAF.txt", 
  gender = gender, genomeVersion = "hg19", isTargetedSeq=T)


ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, 
GCcontentfile = "/public/ClinicalExam/lj_sih/resource/ASCAT/hg19/GC_correction.txt", 
replictimingfile = "/public/ClinicalExam/lj_sih/resource/ASCAT/hg19/Rep_correction.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc, penalty=25)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = output_rdata)
