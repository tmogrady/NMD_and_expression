# to investigate whether ILDs of NMD-target SE events are different
# in (putatively) NMD-inhibited conditions vs. splicing alteration conditions
# based on output from rMATS/SpliceTools

library("dplyr")
library("stringr") # for str_split_fixed
library("ggplot2")

#SETUP ####
ann <- read.table("reference_files/hg38_plus_Akata_inverted.bed.converted.bed")
ebv_ann <- ann %>%
  filter(V1 == "chrEBV_Akata_inverted")
ebv_ann[13:16] <- str_split_fixed(ebv_ann$V4, "_", 4)
ebv_genes <- unique(ebv_ann$V15)

set_tpm <- 10


#INPUTS ####
#get SE
ds1_se <- read.table("../../1_NMD_inhibitor_RNA-seq/CC115/Mutu_1uM/1_JCEC/Mutu_CC115_1uM_vs_cntl_test_cntl_SE.MATS.JCEC.txt",
                     header = TRUE)
ds2_se <- read.table("../../../Splicing_EBV/1_RNASeq_data_and_analyses/3_reactivation/9b_Mutu_I_polyA_Zta/2_rMATS/1_JCEC/Mutu_polyA_Zta_test_cntl_SE.MATS.JCEC.txt",
                     header = TRUE)
#identify NMD-targeting SE
ds1_nmd_bed_neg <- read.table("../../1_NMD_inhibitor_RNA-seq/CC115/Mutu_1uM/3_SpliceTools/SETranslateNMD_Mutu_CC115_1uM_vs_cntl_test_cntl_SE.MATS.JCEC_FDR_0.0005/6_NMD/SE_NMD_lists.bed/SE_NMD.neg_IncDiff.bed",
                              header = TRUE)
ds2_nmd_bed_neg <- read.table("../../../Splicing_EBV/1_RNASeq_data_and_analyses/2022-12-21_SpliceTools/3_output/9b_Mutu_Zta/SETranslateNMD_Mutu_polyA_Zta_test_cntl_SE.MATS.JCEC_FDR_0.0005/6_NMD/SE_NMD_lists.bed/SE_NMD.neg_IncDiff.bed",
                              header = TRUE)

#ANALYSIS #####
#get ILDs of NMD-targeting SEs
ds1_nmd_neg_ild <- left_join(ds1_nmd_bed_neg, ds1_se, 
                         by = c("chr", "strand", "SE_donor" = "upstreamEE", "SE_acceptor" = "downstreamES"),
                         multiple = "all")
ds2_nmd_neg_ild <- left_join(ds2_nmd_bed_neg, ds2_se, 
                             by = c("chr", "strand", "SE_donor" = "upstreamEE", "SE_acceptor" = "downstreamES"),
                             multiple = "all")
#some are quasi-duplicates. Pick the biggest one in the right direction
ds1_nmd_neg_ild_slim <- ds1_nmd_neg_ild %>%
  filter(IncLevelDifference < 0) %>%
  group_by(chr, strand, SE_donor, SE_acceptor) %>%
  summarise(ILD = min(IncLevelDifference))
ds2_nmd_neg_ild_slim <- ds2_nmd_neg_ild %>%
  filter(IncLevelDifference < 0) %>%
  group_by(chr, strand, SE_donor, SE_acceptor) %>%
  summarise(ILD = min(IncLevelDifference))



