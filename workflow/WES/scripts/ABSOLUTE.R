library(tidyverse)
library(GenomicRanges)

input_cnv <- snakemake@input[['cnv_rdata']]
input_cnv <- R.utils::getAbsolutePath(input_cnv)

input_maf <- snakemake@input[['maf']]
input_maf <- R.utils::getAbsolutePath(input_maf)

threads = snakemake@threads


output_ndt_seg <- snakemake@output[['ndt_seg']]
output_ndt_maf <- snakemake@output[['ndt_maf']]
output_seg <- snakemake@output[['seg']]

out_absolute_dir <- snakemake@output[['absolute_dir']]

sample_name <- basename(input_cnv) %>% str_remove('_ASCAT.rdata')

load(input_cnv)

snp_df <- ascat.bc$SNPpos
snp_df$end <- snp_df$Position + 1

snp_gr <- makeGRangesFromDataFrame(
  snp_df,seqnames.field = 'Chromosome',start.field ='Position',
  end.field = 'end'
)

segment_gr <- makeGRangesFromDataFrame(
  ascat.output$segments,
  seqnames.field = 'chr',start.field = 'startpos',end.field = 'endpos'
)

seg <- ascat.output$segments
ov_reg <- findOverlaps(snp_gr,segment_gr)
snp_sum <- as.data.frame(ov_reg) %>% count(subjectHits)

seg$Segment_Mean <- log2(seg$nMajor + seg$nMinor) -1
seg$Num_Probes <- 0
seg$Num_Probes[snp_sum$subjectHits] <- snp_sum$n
seg$Sample = sample_name

Seg <- seg %>% dplyr::select(Sample,chr,startpos,endpos,Num_Probes,Segment_Mean)
colnames(Seg) <- c("Sample","Chromosome","Start","End", "Num_Probes","Segment_Mean")


Maf <- data.table::fread(input_maf)
Maf$Chromosome <- Maf$Chromosome %>% str_remove('chr')
Maf$Tumor_Sample_Barcode <- Maf$Tumor_Sample_Barcode %>% str_remove('_T')

Maf_gr <- makeGRangesFromDataFrame(Maf,
  seqnames.field  = 'Chromosome',start.field = 'Start_Position',end.field = 'End_Position'
)

ov_maf <- findOverlaps(Maf_gr,segment_gr)
maf_sum <- as.data.frame(ov_maf)
Maf$local_cn_a1 <- NA
Maf$local_cn_a1[maf_sum$queryHits] <- seg$nMinor[maf_sum$subjectHits]
Maf$local_cn_a2 <- NA
Maf$local_cn_a2[maf_sum$queryHits] <- seg$nMajor[maf_sum$subjectHits]


ndt_seg <- seg %>% dplyr::select(chr,startpos,endpos,nMinor,nMajor)
colnames(ndt_seg) <- c('Chromosome','Start','End','A1.Seg.CN','A2.Seg.CN')

write_tsv(Seg,output_seg)
write_tsv(ndt_seg,output_ndt_seg)
write_tsv(Maf,output_ndt_maf)

#运行
library(DoAbsolute)
DoAbsolute(Seg = Seg, Maf = Maf,
           platform = "SNP_6.0",
           copy.num.type = "total",
           results.dir = out_absolute_dir,
           nThread = threads,
           #nThread用于计算的核数cores,默认值为1
           keepAllResult = TRUE,
           #keepAllResult为TRUE（默认）是清除所有结果，否则清除结果目录并保留最重要的结果
           verbose = TRUE)

