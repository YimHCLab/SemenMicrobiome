#Load Library ----
library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(agricolae)
library(FSA)
library(dplyr)
library(decontam); packageVersion("decontam")
library('qiime2R')
packages <- c("MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "picante", "pheatmap", "igraph", "rgexf", 
              "ggalluvial", "ggh4x", "FSA", "gridExtra", "aplot", "NST", "GGally")
# Now check or install
for(x in packages){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

#Import Ellis Fok threshold = 0.56 & Remove Neg----

#change to Microtable format
Oct23_16S.f.con.micro <- phyloseq2meco(Oct23_16S.f.noncontam056)

#*-Tran_Venn from Decontam to see Negative
#merge samples as one community for each group
Oct23_16S.f.con.micro1 <- Oct23_16S.f.con.micro$merge_samples(use_group = "Fertility")

##create trans_venn object in ratio====
t3 <- trans_venn$new(Oct23_16S.f.con.micro1, ratio ="seqratio")
t3$plot_venn()

#this decided whether no feature from negative control present in samples

##Perform absolute quantification and rarefaction (QMP)----
#https://github.com/raeslab/QMP/blob/master/QMP.R

#removing 2 negative controls from phyloseq object
Oct23_16S.f.noncontam056

Oct23_16S.056 <- subset_samples(Oct23_16S.f.noncontam056, sample_names(Oct23_16S.f.noncontam056) != "SNTC")
Oct23_16S.056 <- subset_samples(Oct23_16S.056, sample_names(Oct23_16S.056) != "S1")

##create trans_venn object in ratio====
t3 <- trans_venn$new(Oct23_16S.f.con.micro1, ratio ="seqratio")
t3$plot_venn()

#this decided whether no feature from negative control present in samples

##Perform absolute quantification and rarefaction (QMP)----
#https://github.com/raeslab/QMP/blob/master/QMP.R

#removing 2 negative controls from phyloseq object
Oct23_16S.f.noncontam056

Oct23_16S.056 <- subset_samples(Oct23_16S.f.noncontam056, sample_names(Oct23_16S.f.noncontam056) != "SNTC")
Oct23_16S.056 <- subset_samples(Oct23_16S.056, sample_names(Oct23_16S.056) != "S1")

###Function for QMP----
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[row.names(cnv_corrected_abundance_table),]) # make sure the order of the samples is the same in both files    
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}

###Performing QMP using Nature paper ----
Oct23_16S.056.otu <- as.data.frame(otu_table(Oct23_16S.056 ))
Oct23_16S.056.otu <- t(Oct23_16S.056.otu)

Oct23_16S.056.sam <- sample_data(Oct23_16S.056)

Oct23_16S.056.CFU <- Oct23_16S.056.sam[,16]
Oct23_16S.056.CFU <- as.matrix(Oct23_16S.056.CFU)
colnames(Oct23_16S.056.CFU) <- "#"

rarefy_even_sampling_depth(Oct23_16S.056.otu,Oct23_16S.056.CFU)
Oct23_16S.056.rar <- rarefy_even_sampling_depth(Oct23_16S.056.otu,Oct23_16S.056.CFU)

###Importing this absolute abundance back to Phyloseq object----

Oct23_16S.056.rar.t <- as.data.frame(t(Oct23_16S.056.rar))
Oct23_16S.056.rar.t2 <- Oct23_16S.056.rar.t

#Convert number into integer

Oct23_16S.056.rar.t2[,1:87] <- lapply(Oct23_16S.056.rar.t2[,1:87], function(x) as.integer(round(x)))

Oct23_16S.056.a <- phyloseq(otu_table(Oct23_16S.056.rar.t2, taxa_are_rows = TRUE), 
                            sample_data(Oct23_16S.056), 
                            tax_table(Oct23_16S.056), phy_tree(Oct23_16S.056))

#Convert Phyloseq to Microtable of Microeco---- 

Oct23_16S.056.a.micro <- phyloseq2meco(Oct23_16S.056.a)

Oct23_16S.056.a.micro.c <- clone(Oct23_16S.056.a.micro)

Oct23_16S.056.a.micro.c$sample_table <- subset(Oct23_16S.056.a.micro.c$sample_table, is.neg == "FALSE")

Oct23_16S.056.a.micro.c$tidy_dataset()

Oct23_16S.056.a.micro.c

#Calculate Abundance, Alpa and Beta diversity----
# use default parameters, calculate the taxa abundance at each taxonomic rank using cal_abund()
Oct23_16S.056.a.micro.c$cal_abund(rel = FALSE)

# return dataset$taxa_abund
class(Oct23_16S.056.a.micro.c$taxa_abund)

# show part of the relative abundance at Phylum level
Oct23_16S.056.a.micro.c$taxa_abund$Phylum[1:5, 1:5]

#save the taxa abundance file to a local place easily
Oct23_16S.056.a.micro.c$save_abund(dirpath = "taxa_abund")

#calculate the alpha diversity; with phylogenetic

Oct23_16S.056.a.micro.c$cal_alphadiv(measures = c("Shannon", "Simpson", 
                                                  "Chao1","InvSimpson", "PD", "Observed",
                                                  "Coverage"))

Oct23_16S.056.a.micro.c$cal_alphadiv(measures = c("Shannon", "ACE"), PD = TRUE)

class(Oct23_16S.056.a.micro.c$alpha_diversity)

# save dataset$alpha_diversity to a directory
Oct23_16S.056.a.micro.c$save_alphadiv(dirpath = "alpha_diversity")

# calculate beta diversity
Oct23_16S.056.a.micro.c$cal_betadiv(unifrac = TRUE)

class(Oct23_16S.056.a.micro.c$beta_diversity)

# save dataset$beta_diversity to a directory
Oct23_16S.056.a.micro.c$save_betadiv(dirpath = "beta_diversity")

##use 10 Phyla with the highest abundance====
t1 <- trans_abund$new(dataset = Oct23_16S.056.a.micro.c, taxrank = "Phylum",
                      ntaxa = 10)

t1$plot_heatmap(facet = "Fertility", 
           xtext_keep = FALSE, withmargin = FALSE,
           color_values = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))+
  labs(fill = "CFU in 300ul Seminal Fluid")

##Diversity class====
##Alpha Diversity====

###Fertility ====
#there is a bug, ACE is NA, and so need to do this.

Oct23_16S.056.a.micro.c$cal_alphadiv(measures = c("Shannon","Simpson",
                                                  "Observed","Chao1"), PD = TRUE)

t5.f <- trans_alpha$new(dataset = Oct23_16S.056.a.micro.c, group = "Fertility")
t5.f$data_stat <- na.omit(t5.f$data_stat)
t5.f$data_alpha <- na.omit(t5.f$data_alpha)

t5.f$cal_diff(method = "t.test")

t5.f$res_diff

t5.f$plot_alpha(measure = "Shannon", shape = "Fertility")
t5.f$plot_alpha(measure = "PD", shape = "Fertility")
t5.f$plot_alpha(measure = "Chao1", shape = "Fertility")
t5.f$plot_alpha(measure = "Observed", shape = "Fertility")
t5.f$plot_alpha(measure = "Simpson", shape = "Fertility")


##Beta Diversity====
###Fertility====
###Bray ====
#Control vs infertile

t6.f <- trans_beta$new(dataset = Oct23_16S.056.a.micro.c, group = "Fertility", measure = "bray")
t6.f$cal_ordination(ordination = "PCoA")
class(t6.f$res_ordination)
t6.f$plot_ordination(plot_color = "Fertility", plot_shape = "Fertility", plot_type = c("point", "ellipse"))+
  labs(title = "Bray Curtis")

#Permanova Overall
t6.f$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.f$res_manova

#Permanova pairwise
t6.f$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.f$res_manova
                                      
#Permutation Test
t6.f$cal_betadisper()
t6.f$res_betadisper

###jaccard ====
#Control vs infertile
t6.jaccard.f <- trans_beta$new(dataset = Oct23_16S.056.a.micro.c, group = "Fertility", measure = "jaccard")

t6.jaccard.f$cal_ordination(ordination = "PCoA")
class(t6.jaccard.f$res_ordination)
t6.jaccard.f$plot_ordination(plot_color = "Fertility", plot_shape = "Fertility", plot_type = c("point", "ellipse"))+
  labs(title = "Jaccard")

#Permanova Overall
t6.jaccard.f$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.jaccard.f$res_manova

#Permanova pairwise
t6.jaccard.f$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.jaccard.f$res_manova
                                      
#Permutation Test
t6.jaccard.f$cal_betadisper()
t6.jaccard.f$res_betadisper

###unwei_unifrac ====

#Control vs Infertility
t6.unwei_unifrac.f <- trans_beta$new(dataset = Oct23_16S.056.a.micro.c, group = "Fertility", measure = "unwei_unifrac")

t6.unwei_unifrac.f$cal_ordination(ordination = "PCoA")
class(t6.unwei_unifrac.f$res_ordination)
t6.unwei_unifrac.f$plot_ordination(plot_color = "Fertility", plot_shape = "Fertility", plot_type = c("point", "ellipse"))+
                                       labs(title = "Unweighted Unifrac")

#Permanova Overall
t6.unwei_unifrac.f$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.unwei_unifrac.f$res_manova

#Permanova Pairwise
t6.unwei_unifrac.f$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.unwei_unifrac.f$res_manova
                                      
#Permutation Test
t6.unwei_unifrac.f$cal_betadisper()
t6.unwei_unifrac.f$res_betadisper

###wei_unifrac ====
#Control vs Infertility

t6.wei_unifrac.f <- trans_beta$new(dataset = Oct23_16S.056.a.micro.c, group = "Fertility", measure = "wei_unifrac")

t6.wei_unifrac.f$cal_ordination(ordination = "PCoA")
class(t6.wei_unifrac.f$res_ordination)
t6.wei_unifrac.f$plot_ordination(plot_color = "Fertility", plot_shape = "Fertility", plot_type = c("point", "ellipse"))+
                                      labs(title = "Weighted Unifrac")

#Permanova Overall
t6.wei_unifrac.f$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.wei_unifrac.f$res_manova

#Permanova pairwise
t6.wei_unifrac.f$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
t6.wei_unifrac.f$res_manova
                                      
#Permutation Test
t6.wei_unifrac.f$cal_betadisper()
t6.wei_unifrac.f$res_betadisper

###Network====
##Co-abundance----
###Overall----
#Spearman correlation using WGCNA package

library(WGCNA)

####Control----
Control <- clone(Oct23_16S.056.a.micro.c)

Control$sample_table <- subset(Control$sample_table, Fertility == "Control")
Control$tidy_dataset()

#####Network----
t9.con <- trans_network$new(dataset = Control, cor_method = "spearman", 
                            use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)
# construct network; require igraph package
t9.con$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t9.con$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t9$res_network

# # modularity for this undirected network with igraph cluster_fast_greedy function
t9.con$cal_module(method = "cluster_fast_greedy")

# save network
# open network.gexf file using Gephi(https://gephi.org/)
# require rgexf package
t9.con$save_network(filepath = "network_Con.gexf")

# calculate network attributes
t9.con$cal_network_attr()
t9.con$res_network_attr

# get node properties
t9.con$get_node_table(node_roles = TRUE)
# return t9$res_node_table

# get edge properties
t9.con$get_edge_table()
# return t9$res_edge_table 

t9.con$get_adjacency_matrix()
# return t9$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
t9.con$plot_taxa_roles(use_type = 1, roles_color_background = TRUE)

# plot node roles with phylum information
t9.con$plot_taxa_roles(use_type = 2) + ylim(-3,3)

t9.con$cal_eigen()
# return t9$res_eigen



####Infertile----
Infertile <- clone(Oct23_16S.056.a.micro.c)

Infertile$sample_table <- subset(Infertile$sample_table, Fertility == "Infertile")
Infertile$tidy_dataset()

#####Network----
t9.In <- trans_network$new(dataset = Infertile, cor_method = "spearman", 
                           use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)
# construct network; require igraph package
t9.In$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
t9.In$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return t9$res_network

# # modularity for this undirected network with igraph cluster_fast_greedy function
t9.In$cal_module(method = "cluster_fast_greedy")

# save network
# open network.gexf file using Gephi(https://gephi.org/)
# require rgexf package
t9.In$save_network(filepath = "network_In.gexf")

# calculate network attributes
t9.In$cal_network_attr()
t9.In$res_network_attr

# get node properties
t9.In$get_node_table(node_roles = TRUE)
# return t9$res_node_table

# get edge properties
t9.In$get_edge_table()
# return t9$res_edge_table 

t9.In$get_adjacency_matrix()
# return t9$res_adjacency_matrix

# add_label = TRUE can be used to directly add text labels for points
t9.In$plot_taxa_roles(use_type = 1, roles_color_background = TRUE)

# plot node roles with phylum information
t9.In$plot_taxa_roles(use_type = 2) + ylim(-3,3)

t9.In$cal_eigen()
# return t9.In$res_eigen

##Explainable Classes====

###Sample Data table cloning and selecting the column----

new_test <- clone(Oct23_16S.056.a.micro.c)

new_test$sample_table <- data.frame(new_test$sample_table, env_data_16S[rownames(new_test$sample_table), ])
# now new_test$sample_table has the whole data

#fill blank of the sample table with NA

new_test$sample_table[new_test$sample_table == ""] <- NA 

new_test

library(mice)
##Correlation between metadata----
#there is missing data in some patients, so use mice package to interpolate the data

t10 <- trans_env$new(dataset = new_test, env_cols = 8:13, character2numeric = TRUE, complete_na = TRUE)

t10$cal_autocor()
t10$cal_autocor(group = "Fertility")

##RDA analysis----

t10$cal_ordination(method = "RDA", taxa_level = "Genus")

# As the main results of RDA are related with the projection and angles between different arrows,
# we adjust the length of the arrow to show them clearly using several parameters.

t10$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

t10$plot_ordination(plot_color = "Fertility", plot_shape = "Fertility")

t10.table <- t10$res_ordination_axis


#Use 0.56 threshold Grouping by defb119#########
#Import Ellis Fok threshold = 0.56 & Remove Neg----

library(file2meco)
library(microeco)
library(phyloseq)
library(tidyverse)

##Copy Microtable ====
Oct23_16S.056.a.micro2.c <- clone(Oct23_16S.056.a.micro.c)

Oct23_16S.056.a.micro2.c$sample_table <- subset(Oct23_16S.056.a.micro2.c$sample_table, is.neg == "FALSE")

Oct23_16S.056.a.micro2.c$tidy_dataset()

Oct23_16S.056.a.micro2.c

# use default parameters, calculate the taxa abundance at each taxonomic rank using cal_abund()
Oct23_16S.056.a.micro2.c$cal_abund(rel = FALSE)

# return dataset$taxa_abund
class(Oct23_16S.056.a.micro2.c$taxa_abund)

# show part of the relative abundance at Phylum level
Oct23_16S.056.a.micro2.c$taxa_abund$Phylum[1:5, 1:5]

#save the taxa abundance file to a local place easily
Oct23_16S.056.a.micro2.c$save_abund(dirpath = "taxa_abund")

#calculate the alpha diversity; with phylogenetic

Oct23_16S.056.a.micro2.c$cal_alphadiv(measures = c("Shannon", 
                                                   "Simpson", "Observed"),
                                      PD = TRUE)

class(Oct23_16S.056.a.micro2.c$alpha_diversity)

# save dataset$alpha_diversity to a directory
Oct23_16S.056.a.micro2.c$save_alphadiv(dirpath = "alpha_diversity")

# calculate beta diversity
Oct23_16S.056.a.micro2.c$cal_betadiv(unifrac = TRUE)
class(Oct23_16S.056.a.micro2.c$beta_diversity)
# save dataset$beta_diversity to a directory
Oct23_16S.056.a.micro2.c$save_betadiv(dirpath = "beta_diversity")

##use 10 Phyla with the highest abundance====
g1.phyla <- trans_abund$new(dataset = Oct23_16S.056.a.micro2.c, taxrank = "Phylum",
                      ntaxa = 10)

g1.phyla$plot_heatmap(facet = "Group", 
                xtext_keep = FALSE, withmargin = FALSE,
                color_values = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))+
  labs(fill = "CFU")
                                      
##show 40 taxa at Genus level====
g2 <- trans_abund$new(dataset = Oct23_16S.056.a.micro2.c, 
                      taxrank = "Genus", ntaxa = 40)
g2$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)+
  labs(fill = "CFU")

##Diversity class====
##Alpha Diversity====

Oct23_16S.056.a.micro2.c$cal_alphadiv(measures = c("Shannon","Simpson",
                                                  "Observed","Chao1"), PD = TRUE)

g5 <- trans_alpha$new(dataset = Oct23_16S.056.a.micro.c, group = "Group")
g5$data_stat <- na.omit(g5$data_stat)
g5$data_alpha <- na.omit(g5$data_alpha)

g5$cal_diff(method = "anova")

g5$res_diff

library(writexl)

write_xlsx(g5$res_diff,"alpha_diversity_anova.xlsx")

#https://myaseen208.com/agricolae/articles/MultipleComparisons.html
#Means with the same letter are not significantly different.

g5$plot_alpha(measure = "Shannon", shape = "Group")+ labs(title = "Shannon diversity among G1-G3 group",
                       caption = "ANOVA with duncan test
            Means with the same letter are not significantly different")

PD.plot.g <- g5$plot_alpha(measure = "PD", shape = "Group")

PD.plot.g + labs(title = "PD diversity among G1-G3 group",
                 caption = "ANOVA with duncan test
            Means with the same letter are not significantly different")


Chao1.plot.g <- g5$plot_alpha(measure = "Chao1", shape = "Group")

Chao1.plot.g + labs(title = "Chao1 diversity among G1-G3 group",
                    caption = "ANOVA with duncan test
            Means with the same letter are not significantly different")



Simpson.plot.g <- g5$plot_alpha(measure = "Simpson", shape = "Group")
Simpson.plot.g  + labs(title = "Simpson diversity among G1-G3 group",
                       caption = "ANOVA with duncan test
            Means with the same letter are not significantly different")

Observed.plot.g <- g5$plot_alpha(measure = "Observed", shape = "Group")
Observed.plot.g + labs(title = "Observed taxa among G1-G3 group",
                       caption = "ANOVA with duncan test
            Means with the same letter are not significantly different")

##Beta Diversity====
###Bray ====

g6 <- trans_beta$new(dataset = Oct23_16S.056.a.micro2.c, 
                     group = "Group", measure = "bray")

g6$cal_ordination(ordination = "PCoA")
class(g6$res_ordination)

g6$plot_ordination(plot_color = "Group", plot_shape = "Group", 
                   plot_type = c("point", "ellipse"))+ 
  labs(title = "Bray Curtis")

#Permanova Overall
g6$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6$res_manova  


#Permanova pairwise
g6$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6$res_manova     
                                      

#Permutation Test
g6$cal_betadisper()
g6$res_betadisper

###jaccard ====
g6.jaccard <- trans_beta$new(dataset = Oct23_16S.056.a.micro2.c, 
                             group = "Group", measure = "jaccard")

g6.jaccard$cal_ordination(ordination = "PCoA")
class(g6.jaccard$res_ordination)
g6.jaccard$plot_ordination(plot_color = "Group", plot_shape = "Group", 
                           plot_type = c("point", "ellipse"))+
  labs(title = "Jaccard")

#Permanova Overall
g6.jaccard$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.jaccard$res_manova  

#Permanova Pairwise
g6.jaccard$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.jaccard$res_manova 
                                      
#Permutation Test
g6.jaccard$cal_betadisper()
g6.jaccard$res_betadisper

###unwei_unifrac ====
g6.unwei_unifrac <- trans_beta$new(dataset = Oct23_16S.056.a.micro2.c, 
                                   group = "Group", measure = "unwei_unifrac")

g6.unwei_unifrac$cal_ordination(ordination = "PCoA")
class(g6.unwei_unifrac$res_ordination)
g6.unwei_unifrac$plot_ordination(plot_color = "Group", plot_shape = "Group", 
                                 plot_type = c("point", "ellipse"))+
  labs(title = "Unweighted Unifrac")

#Permanova Overall
g6.unwei_unifrac$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.unwei_unifrac$res_manova  

#Permanova Pairwise
g6.unwei_unifrac$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.unwei_unifrac$res_manova  

#Permutation Test
g6.unwei_unifrac$cal_betadisper()
g6.unwei_unifrac$res_betadisper

###wei_unifrac ====
g6.wei_unifrac <- trans_beta$new(dataset = Oct23_16S.056.a.micro2.c, 
                                 group = "Group", measure = "wei_unifrac")

g6.wei_unifrac$cal_ordination(ordination = "PCoA")
class(g6.wei_unifrac$res_ordination)
g6.wei_unifrac$plot_ordination(plot_color = "Group", plot_shape = "Group", 
                               plot_type = c("point", "ellipse"))+
  labs(title = "Weighted Unifrac")


#Permanova Overall
g6.wei_unifrac$cal_manova(
  manova_all = TRUE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.wei_unifrac$res_manova  


#Permanova pairwise
g6.wei_unifrac$cal_manova(
  manova_all = FALSE,
  manova_set = NULL,
  group = NULL,
  p_adjust_method = "fdr"
)
g6.wei_unifrac$res_manova  

#Permutation Test
g6.wei_unifrac$cal_betadisper()
g6.wei_unifrac$res_betadisper


##Subgrouping the dataset----
###G1----
G1 <- clone(Oct23_16S.056.a.micro2.c)

G1$sample_table <- subset(G1$sample_table, Group == "G1")
G1$tidy_dataset()
###G2----
G2 <- clone(Oct23_16S.056.a.micro2.c)
G2$sample_table <- subset(G2$sample_table, Group == "G2")
G2$tidy_dataset()
###G3----
G3 <- clone(Oct23_16S.056.a.micro2.c)
G3$sample_table <- subset(G3$sample_table, Group == "G3")
G3$tidy_dataset()

                                      
##Network====
###Co-abundance----

# Spearman correlation using WGCNA package

library(WGCNA)
###G1====
g9.G1 <- trans_network$new(dataset = G1, 
                           cor_method = "spearman",
                           use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001,
                           add_taxa_name = c("Phylum", "Genus"))


# construct network; require igraph package
g9.G1$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, complete_na = TRUE)
# use arbitrary coefficient threshold to contruct network
g9.G1$cal_network(COR_p_thres = 0.01, COR_cut = 0.7, complete_na = TRUE)
# return t9$res_network

# # modularity for this undirected network with igraph cluster_fast_greedy function
g9.G1$cal_module(method = "cluster_fast_greedy")

g9.G1$cal_module()
nodes <- V(g9.G1$res_network)$name
V(g9.G1$res_network)$RelativeAbundance <- apply(g9.G1$data_abund[, nodes], 2, sum) * 100/sum(dataset$otu_table)
g9.G1$res_network
g9.G1$save_network("network_G1.gexf")

# save network
# open network.gexf file using Gephi(https://gephi.org/)
# require rgexf package
g9.G1$save_network(filepath = "network_G1.gexf")

# calculate network attributes
g9.G1$cal_network_attr()
g9.G1$res_network_attr

# get node properties
g9.G1$get_node_table(node_roles = TRUE)
# return t9$res_node_table

# get edge properties
g9.G1$get_edge_table()
# return t9$res_edge_table 

g9.G1$get_adjacency_matrix()
# return t9$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
g9.G1$plot_taxa_roles(use_type = 1, roles_color_background = TRUE)

# plot node roles with phylum information
g9.G1$plot_taxa_roles(use_type = 2)

g9.G1$cal_eigen()
# return t9$res_eigen
                                      

###G2====
g9.G2 <- trans_network$new(dataset = G2, 
                           cor_method = "spearman",
                           use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001,
                           add_taxa_name = c("Phylum", "Genus"))


# construct network; require igraph package
g9.G2$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, complete_na = TRUE)
# use arbitrary coefficient threshold to contruct network
g9.G2$cal_network(COR_p_thres = 0.01, COR_cut = 0.7, complete_na = TRUE)
# return t9$res_network

# # modularity for this undirected network with igraph cluster_fast_greedy function
g9.G2$cal_module(method = "cluster_fast_greedy")

g9.G2$cal_module()
nodes <- V(g9.G2$res_network)$name
V(g9.G2$res_network)$RelativeAbundance <- apply(g9.G2$data_abund[, nodes], 2, sum) * 100/sum(g9.G2$otu_table)
g9.G2$res_network
g9.G2$save_network("network_G2.gexf")

# open network.gexf file using Gephi(https://gephi.org/)

# calculate network attributes
g9.G2$cal_network_attr()
g9.G2$res_network_attr


# get node properties
g9.G2$get_node_table(node_roles = TRUE)
# return t9$res_node_table

# get edge properties
g9.G2$get_edge_table()
# return t9$res_edge_table 

g9.G2$get_adjacency_matrix()
# return t9$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
g9.G2$plot_taxa_roles(use_type = 1, roles_color_background = TRUE)

# plot node roles with phylum information
g9.G2$plot_taxa_roles(use_type = 2)

g9.G2$cal_eigen()
g9.G2$res_eigen
# return t9$res_eigen

###G3====
g9.G3 <- trans_network$new(dataset = G3, 
                           cor_method = "spearman",
                           use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001,
                           add_taxa_name = c("Phylum", "Genus"))


# construct network; require igraph package
g9.G3$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE, complete_na = TRUE)
# use arbitrary coefficient threshold to contruct network
g9.G3$cal_network(COR_p_thres = 0.01, COR_cut = 0.7, complete_na = TRUE)
# return t9$res_network

# # modularity for this undirected network with igraph cluster_fast_greedy function
g9.G3$cal_module(method = "cluster_fast_greedy")

g9.G3$cal_module()
nodes <- V(g9.G3$res_network)$name
V(g9.G3$res_network)$RelativeAbundance <- apply(g9.G3$data_abund[, nodes], 2, sum) * 100/sum(g9.G3$otu_table)
g9.G3$res_network
g9.G3$save_network("network_G3.gexf")

# open network.gexf file using Gephi(https://gephi.org/)

# calculate network attributes
g9.G3$cal_network_attr()
g9.G3$res_network_attr


# get node properties
g9.G3$get_node_table(node_roles = TRUE)
# return t9$res_node_table

# get edge properties
g9.G3$get_edge_table()
# return t9$res_edge_table 

g9.G3$get_adjacency_matrix()
# return t9$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
g9.G3$plot_taxa_roles(use_type = 1, roles_color_background = TRUE)

# plot node roles with phylum information
g9.G3$plot_taxa_roles(use_type = 2)

g9.G3$cal_eigen()
g9.G3$res_eigen
# return t9$res_eigen


                                      
