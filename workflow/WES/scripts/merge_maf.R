#!/public/home/lijf/env/miniconda3/envs/clindet/bin/Rscript

library(stringr)
library(data.table)
library(tools)
library(tidyverse)
input_maf_dir <- snakemake@params[['dir']]
output_merge_maf <- snakemake@output[['maf']]
output_filter_merge_maf <- snakemake@output[['filter_maf']]

sample_name <- paste0(file_path_sans_ext(basename(output_merge_maf)),'_T')

# input_maf_dir <- '2022-BGI'
# output_merge_maf <- '2022-BGI/merge.maf'


# callers_vcf <- list.files(input_maf_dir, ".maf", full.names = TRUE)
callers_vcf <- snakemake@input[['maf1']]
callers_vcf <- callers_vcf[!str_detect(callers_vcf,'strelkasomatic.vcf.maf')]
# print(callers_vcf)
# saveRDS(callers_vcf,'test.rds')


final <- NULL
#callers <- c("Mutect2", "HaplotypeCaller", "Strelka", "VarDict", "Pindel", "UnifiedGenoTyper", "Freebayes", "Varscan2", "Lofreq")
#callers <- c("Mutect2", "HaplotypeCaller", "Strelka", "StrelkaSomatic", "UnifiedGenoTyper")

index_name <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")

final.tmp <- NULL
for(caller_file in callers_vcf) {
  caller <- file_path_sans_ext(file_path_sans_ext(basename(caller_file)))
  dat <- data.table::fread(caller_file, data.table = FALSE)
  if ("Fusion" %in% colnames(dat)) {
    dat <- dat[,-which(colnames(dat) %in% c("Fusion", "Method", "Frame", "CONSENSUS"))]
  }
  if (caller == "Lofreq") {
    i <- ncol(dat) - 1
    i2 <- ncol(dat)
    dat$t_depth <- dat[,i]
    dat$t_ref_count <- round(dat[,i] * (1 - dat[,i2]))
    dat$t_alt_count <- round(dat[,i] * dat[,i2], 0)
    dat <- dat[,-c(i, i2)]
  }
  
  if (!"Chromosome" %in% colnames(dat)) next
  if (nrow(dat) == 0) next
  dat$caller <- caller
  dat$t_vaf <- as.numeric(dat$t_alt_count) / as.numeric(dat$t_depth)
  dat$n_vaf <- as.numeric(dat$n_alt_count) / as.numeric(dat$n_depth)
  
  # fil1 <- !dat$Variant_Classification %in% c("3'UTR", "3'Flank", "5'Flank", "Intron",
  #                                            "Silent", "5'UTR", "RNA", "IGR",
  #                                            "Splice_Region", "Translation_Start_Site", "Targeted_Region")
  # fil2 <- str_detect(dat$all_effects, "missense_variant|stop|splice|inframe|frameshift")
  fil1 <- !dat$Variant_Classification %in% c("3'UTR", "3'Flank", "5'Flank", "Intron",
                                             "5'UTR", "RNA", "IGR", "Targeted_Region")
  fil2 <- str_detect(dat$all_effects, "missense_variant|stop|splice|inframe|frameshift|synonymous")
  dat <- dat[fil1 | fil2,]
  dat <- dat[is.na(dat$gnomAD_AF) | dat$gnomAD_AF <= 0.1,]
  dat <- dat[is.na(dat$AF) | dat$AF <= 0.1,]
  dat$FILTER <- paste0(caller,":",dat$FILTER,"|")
  # print(caller)
  # print(colnames(dat))R
  final.tmp <- rbind(final.tmp, dat)
  final.tmp <- final.tmp[order(final.tmp[, index_name[1]], final.tmp[, index_name[2]], final.tmp[, index_name[3]],
                               final.tmp[, index_name[4]], final.tmp[, index_name[5]]),]
  idx <- which(duplicated(final.tmp[,index_name]))
  for (i in idx) {
    if (!str_detect(final.tmp[i, "caller"], final.tmp[i - 1, "caller"]))
      final.tmp[i - 1, "caller"] <- paste0(c(final.tmp[i - 1, "caller"], final.tmp[i, "caller"]), collapse = ",")
      final.tmp[i - 1, "FILTER"] <- paste0(c(final.tmp[i - 1, "FILTER"], final.tmp[i, "FILTER"]), collapse = "")
    if (is.na(final.tmp[i, "t_vaf"]) || final.tmp[i, "t_vaf"] == "") {
      next
    }
    if (is.na(final.tmp[i-1,"t_vaf"]) || final.tmp[i-1,"t_vaf"] == "" || final.tmp[i,"t_vaf"] > final.tmp[i - 1, "t_vaf"])  {
      final.tmp[i, "caller"] <- final.tmp[i - 1, "caller"]
      final.tmp[i, "FILTER"] <- final.tmp[i - 1, "FILTER"]
      final.tmp[i -1,] <- final.tmp[i, ] 
    }
  }
  final.tmp <- final.tmp[!duplicated(final.tmp[,index_name]),]
  final.tmp$Tumor_Sample_Barcode <- sample_name
}
vaf <- final.tmp$t_vaf
vaf2 <- final.tmp$n_vaf
fil1 <- vaf >= 0.01 | is.na(vaf)
final.tmp <- final.tmp[fil1,]
final.tmp <- final.tmp[!is.na(final.tmp$Chromosome),]
chr <- final.tmp$Chromosome
start <- final.tmp$Start_Position
end <- final.tmp$End_Position
ref <- final.tmp$Reference_Allele
alt <- final.tmp$Tumor_Seq_Allele2

data.table::fwrite(final.tmp,output_merge_maf, sep = "\t", row.names = FALSE, quote = FALSE)

