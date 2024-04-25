library(VariantAnnotation)
library(tidyverse)
# vcf<-readVcf("~/Desktop/MM-014_T_vs_MM-014_NC.flagged.vcf.gz")


output_vcf <- snakemake@output[['vcf']]
caller <- snakemake@params[['caller']]

# cgppindel vardict varscan2 need filter
if(caller == 'vardict') {
    input_vcf <- snakemake@input[['vcf']]
    # input_vcf <- '/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/vcf/paired/MM-014/vardict.vcf'
    vcf <- readVcf(input_vcf)

    status_filter <- c('LikelySomatic','StrongSomatic')
    tag_filter <- c('PASS')
    filter_index <- (vcf@fixed$FILTER %in% tag_filter) & (vcf@info$STATUS %in% status_filter)
    filter_vcf <- vcf[filter_index]
    writeVcf(filter_vcf,output_vcf,index = F)

} else if (caller == 'cgppindel') {
    input_vcf <- snakemake@params[['vcf']]
    vcf <- readVcf(input_vcf)
    geno(header(vcf)) <- rbind(geno(header(vcf)),
        data.frame(row.names = 'DP',
                    Number = '1',
                    Type = 'Integer',
                    Description = 'Read depth'))

    geno(vcf)$DP <- geno(vcf)$FD
    info(vcf)
    AD_tag <- unlist(str_c(geno(vcf)$FD - geno(vcf)$FC,',',geno(vcf)$FC))
    geno(vcf)$AD <- matrix(AD_tag,ncol = 2, byrow = F,dimnames = dimnames(geno(vcf)$FD))
    filter_vcf <- vcf[VariantAnnotation::fixed(vcf)$FILTER == 'PASS']
    writeVcf(filter_vcf,output_vcf,index = F)
} else if (caller == 'varscan2') {
    input_vcf <- snakemake@input[['vcf']]
    # input_vcf <- '/public/ClinicalExam/lj_sih/projects/project_pipeline/WES/results/vcf/paired/MM-014/varscan2.snp.vcf'
    vcf <- readVcf(input_vcf)

} else {
    stop('Not support caller')
}


