# to investigate expression levels of genes putatively targeted to NMD
# based on output from rMATS/SpliceTools

library("dplyr")
library("stringr") # for str_split_fixed
library("ggplot2")

#INPUTS ####
exp <- read.table("../1_NMD_inhibitor_RNA-seq/CC115/Akata_1uM/2_DESeq/Test.deseq2.results.txt", header=TRUE, sep="\t")
nmd_genes <- read.table("../1_NMD_inhibitor_RNA-seq/CC115/Akata_1uM/3_SpliceTools/SETranslateNMD_Akata_CC115_1uM_vs_cntl_test_cntl_SE.MATS.JCEC_FDR_0.0005/6_NMD/gene_lists/SE_NMD.neg_IncDiff.txt",
                        header = FALSE)
set_tpm <- 10

#ANALYSIS ####
#remove EBV genes from DE file:
ann <- read.table("/reference_files/hg38_plus_Akata_inverted.bed.converted.bed")
ebv_ann <- ann %>%
  filter(V1 == "chrEBV_Akata_inverted")
ebv_ann[13:16] <- str_split_fixed(ebv_ann$V4, "_", 4)
ebv_genes <- unique(ebv_ann$V15)

exp <- exp %>%
  filter(!gene %in% ebv_genes)
#57152

#identify NMD-targeted genes in DE file
exp_nmd_genes <- exp %>%
  mutate(NMD = ifelse(gene %in% nmd_genes$V1, "yes", "no"))

#first plot ####
ggplot(exp_nmd_genes, aes(x = log2FC, colour = NMD)) +
  stat_ecdf() +
  xlim(-7.5, 7.5)

#expression filter: UPDATE WITH COLUMN NAMES ####
#to improve: make this smarter
exp_tpm <- exp_nmd_genes %>%
  filter(A0A..TPMs. > set_tpm) %>%
  filter(A0B..TPMs. > set_tpm) %>%
  filter(A0C..TPMs. > set_tpm)

#pdf("2022-12-09_NMD_log2FC_Akata_CC115_1uM_ks_tpm10_reference.pdf", width = 2.5, height = 3.5)
ggplot(exp_tpm, aes(x = log2FC, colour = NMD)) +
  stat_ecdf() +
  scale_colour_manual(values = c("darkgrey", "#fbb4ae"),
                      labels = c("no NMD-AS", "NMD-AS")) +
  #geom_vline(xintercept = 0, colour = "grey") +
  xlim(-2,2) +
  theme_classic() + theme(legend.position = "bottom") +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    size = 0.5))
#dev.off()

#pdf("2022-12-19_NMD_log2FC_Akata_CC115_1uM_ks_tpm10.pdf", width = 2, height = 2.5)
ggplot(exp_tpm, aes(x = log2FC, colour = NMD)) +
  stat_ecdf(size=1.5) +
  scale_colour_manual(values = c("darkgrey", "#b2df8a"),
                      labels = c("no NMD-AS", "NMD-AS")) +
  xlim(-2,2) +
  theme_classic() + 
  guides(colour = "none") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        #legend.title = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    size = 0.5))
#dev.off()

#STATS ####
nmd_fc <- exp_tpm %>%
  filter(NMD == "yes") %>%
  pull(log2FC)
no_nmd_fc <- exp_tpm %>%
  filter(NMD == "no") %>%
  pull(log2FC)

ks.test(nmd_fc, no_nmd_fc)

