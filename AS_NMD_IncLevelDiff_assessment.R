# to investigate whether ILDs of NMD-target SE events are different
# in (putatively) NMD-inhibited conditions vs. splicing alteration conditions
# based on output from rMATS/SpliceTools

library("dplyr")
library("stringr") # for str_split_fixed
library("ggplot2")

#SETUP ####
ann <- read.table("../3_NMD_and_expression/NMD_and_expression/reference_files/hg38_plus_Akata_inverted.bed.converted.bed")
ebv_ann <- ann %>%
  filter(V1 == "chrEBV_Akata_inverted")
ebv_ann[13:16] <- str_split_fixed(ebv_ann$V4, "_", 4)
ebv_genes <- unique(ebv_ann$V15)

set_tpm <- 10
