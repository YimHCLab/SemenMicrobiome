#Load library----

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library('qiime2R')
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(dunn.test)

#Add CFU data into Metadata----

#omit the S13 data because only 4 byte in raw sequences

metadata3 <- as.data.frame(metadata3)

metadata4 <- merge(metadata3, CFU_data, by = "sample-ID", all = TRUE)

write.table(metadata4, file = "metadata4.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

#plot the CFU on Control vs infertile

metadata5 <- subset(metadata4, metadata4$`sample-ID` != "S1" )
metadata5 <- subset(metadata5, metadata5$`sample-ID` != "S13" )
metadata5 <- subset(metadata5, metadata5$`sample-ID` != "SNTC" )
metadata5 <- subset(metadata5, metadata5$`sample-ID` != "81S150" )

qqnorm(metadata5$CFU)
qqline(metadata5$CFU)

#data are not normal

# Mann-Whitney U test

metadata6 <- metadata5

metadata6$Fertility <- factor(metadata6$Fertility, levels = c("Control", "Infertile"))
metadata6$logCFU <- log10(metadata6$CFU)

ggboxplot(metadata6, x = "Fertility", y = "logCFU", color = "Fertility", 
          add = "jitter", ylab = "Log10(CFU in 300ul seminal fluid)",
          title = "Estimated CFU of seminal microbiome") +
  stat_compare_means(comparisons = list(c("Control", "Infertile")), label = "p.signif",
                     method = "wilcox.test", 
                     size = 4, vjust = -0.5)

# Group G1 G2 G3 comparison of CFU

metadata6$Group <- factor(metadata6$Group, levels = c("G1", "G2", "G3"))

library(rstatix)  
stat.test <- dunn_test(logCFU ~ Group, data = metadata6)

ggboxplot(metadata6, x = "Group", y = "logCFU", color = "Group", 
          add = "jitter", ylab = "Log10(CFU in 300ul seminal fluid)",
          title = "Estimated CFU of seminal microbiome") + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif",
     y.position = c(7.5, 7, 8)) 


#Import data ----
##import data from QIIM2 output in the hf_diver folder, that is after host decontamination 
#and phylogenetic tree construction
#table is not normalised with sequencing depth

Oct23_16S <-qza_to_phyloseq(
  features="table-normalized.qza",
  tree="rooted-tree.qza", "taxonomy.qza",
  metadata = "metadata4.txt"
)
Oct23_16S

#omitting sample 81S150 because of the missing CFU data
#https://github.com/joey711/phyloseq/issues/618

Oct23_16S.f <- subset_samples(Oct23_16S, sample_names(Oct23_16S) != "81S150")

#making is negative or not

sample_data(Oct23_16S.f)$is.neg <- sample_data(Oct23_16S.f)$Sample_or_Control == "Control Sample"

#threshold of 0.56
Oct23f.contamdf.prev056 <- isContaminant(Oct23_16S.f, method="prevalence", neg="is.neg", threshold=0.56)
table(Oct23f.contamdf.prev056$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
Oct23_16S.f.pa <- transform_sample_counts(Oct23_16S.f, function(abund) 1*(abund>0))
Oct23_16S.f.pa.neg <- prune_samples(sample_data(Oct23_16S.f.pa)$Sample_or_Control == "Control Sample", Oct23_16S.f.pa)
Oct23_16S.f.pa.pos <- prune_samples(sample_data(Oct23_16S.f.pa)$Sample_or_Control == "True Sample", Oct23_16S.f.pa)

# Make data.frame of prevalence in positive and negative samples (threshold = 0.56)
Oct23_16S.f.pa.056 <- data.frame(Oct23_16S.f.pos=taxa_sums(Oct23_16S.f.pa.pos), Oct23_16S.f.neg=taxa_sums(Oct23_16S.f.pa.neg),
                        contaminant = Oct23f.contamdf.prev056$contaminant)

plot <- ggplot(data=Oct23_16S.f.pa.056, aes(x=Oct23_16S.f.neg, y=Oct23_16S.f.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") 

plot + ggtitle("The number of taxa were observed in negative controls and positive samples (threshold = 0.56)")

Oct23_16S.f.noncontam056 <- prune_taxa(!Oct23f.contamdf.prev056$contaminant, Oct23_16S.f)

Oct23_16S.f.noncontam056

table(Oct23f.contamdf.prev056$contaminant)

