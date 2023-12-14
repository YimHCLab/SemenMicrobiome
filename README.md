Quantitative Microbiome Profiling (QMP) of Seminal Fluid Microbiome

This analysis is divided into 2 parts:

Part 1 - Raw sequences were analyzed using https://github.com/biofuture/microbiome_singulairty_version

Part 2 - The output from Part 1 were further analyzed using Rstudio
---------------------------------------------------------------------------------
Part 1:




Part 2:

R scripts for Decontamination, Quantitative Microbiome Profiling (QMP) and Mediation Analyses as described in "DEFB119 stratifies dysbiosis with distorted networks in the seminal microbiome associated with male infertility" manuscript.

Output files from QIIME2 including table.qza, rooted-tree.qza, and taxonomy.qza were run through the following R Scripts in the following order:

1. Decontam.R
2. QMP_Microeco_analysis.R
3. Mediation_Analysis.R
