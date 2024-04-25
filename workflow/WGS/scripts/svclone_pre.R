library(tidyverse)
# sv_tsv <- '/public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/MMWGSPE300_RJ/hg19/results/sv/paired/merge/MM-213/final_MM-213.tsv.gz'

genome_version <- snakemake@wildcards[['genome_version']]

sv_tsv  <- snakemake@input[['txt']]
# input_tumor_bam <- getAbsolutePath(input_tumor_bam)


output_txt <- snakemake@output[['txt']]
# output_rdata <- getAbsolutePath(output_rdata)

sv_info <- read_tsv(sv_tsv)
sv_info_sub <- sv_info %>% select(query.chr,query.start,query.chr2,query.end)#,query.svtype,query.ct)
colnames(sv_info_sub) <- c('chr1','pos1','chr2','pos2')#,'classification')
# sv_info_sub <- sv_info_sub  %>% mutate(
#     dir1 = '+',
#     dir2 = '+'
# )
write_tsv(sv_info_sub,output_txt)
# chr1	pos1	dir1	chr2	pos2	dir2	classification
# 22	18240676		-	22	18232335	-	INV
# 22	19940482		-	22	19937820    +	DEL
# 22	21383572		-	22	21382745	+	DUP
# 22	21395573		+	22	21395746	+	INV

# svclone annotate -i MM-213-SV.tsv \
#  -b /public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/MMWGSPE300_RJ/hg19/results/recal/paired/MM-213-T.bam \
#  -s MM-213 --sv_format simple -o MM-213-anno.tsv \
#  -cfg /public/ClinicalExam/lj_sih/softwares/SVclone/svclone_config.ini


# svclone count -i /public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/scripts/MM-213-anno.tsv/MM-213_svin.txt -b /public/ClinicalExam/lj_sih/projects/project_pipeline/WGS/MMWGSPE300_RJ/hg19/results/recal/paired/MM-213-T.bam \
#  -s MM-213 -o MM-213-count \
#  -cfg /public/ClinicalExam/lj_sih/softwares/SVclone/svclone_config.ini