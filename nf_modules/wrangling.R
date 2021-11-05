## Uploading libraries
library("dplyr")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "CA518-1-full_variant_table.vcf"

## Place args into named object
file_name <- args[1]

## Uploading data
vcf_1 <- read.table(file = file_name, header = T, 
                  sep = "\t", stringsAsFactors = F)

## Changing name just for strings
colnames(vcf_1) <- c("CHROM", "POS", "ID", "REF", "ALT", "GT")

## Script
vcf_2 <- vcf_1 %>% 
  mutate(GT = gsub(x = GT,
                               pattern = ":.*",
                               replacement = ""))

## Adding genotype info
vcf_3 <- vcf_2 %>% 
  mutate(HETHOM = case_when(vcf_2$GT == "1/0" | vcf_2$GT == "0/1" ~ "Heterozygous",
                              vcf_2$GT == "1/1" | vcf_2$GT == "0/0" ~ "Homozygous"))

## Getting allele 1 and 2
vcf_4 <- vcf_3 %>% 
  mutate(ALLELE_1 = gsub(x = GT,
                         pattern = "/[0-1]",
                         replacement = ""))

vcf_5 <- vcf_4 %>% 
  mutate(ALLELE_2 = gsub(x = GT,
                         pattern = "[0-1]/",
                         replacement = ""))

## Getting genotiping of allele
vcf_6 <- vcf_5 %>% 
  mutate (ALLELE_1 = ifelse(vcf_4$ALLELE_1 == 0, vcf_4$ALLELE_1 <- vcf_4$REF, vcf_4$ALLELE_1 <- vcf_4$ALT ))
vcf_7 <- vcf_6 %>%
  mutate (ALLELE_2 = ifelse(vcf_5$ALLELE_2 == 0, vcf_5$ALLELE_1 <- vcf_5$REF, vcf_5$ALLELE_1 <- vcf_5$ALT ))

## Selecting columns
vcf_final <- vcf_7 %>% 
  select("CHROM","POS", "REF", "ALT", "GT", "ALLELE_1", "ALLELE_2", "HETHOM")

## Saving as tsv
## Getting name
final_name <- gsub(pattern = "*.vcf.test", replacement = ".tsv", x = file_name)

## Saving
write.table(vcf_final, file = final_name, 
            sep = "\t", col.names = T, row.names = F)
