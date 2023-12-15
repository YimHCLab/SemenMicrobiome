Quantitative Microbiome Profiling (QMP) of Seminal Fluid Microbiome for the "DEFB119 stratifies dysbiosis with distorted networks in the seminal microbiome associated with male infertility" manuscript.
------------------------------------------------------------------------------------------------------

This analysis is divided into 2 parts:

Part 1 - Raw sequences Analysis using QIIME2
------------------------------------------------------------------------------------------------------
The pipeline for analysing the raw sequences can be found at https://github.com/biofuture/microbiome_singulairty_version and https://figshare.com/articles/software/MRC_microbiome_16S_analysis_pipeline_singularity_version/24784080

Linux Scripts for performing the singularity analysis are in the file Singularity.pipeline.script.txt


Part 2 - The output from Part 1 was further analyzed using Rstudio
-------------------------------------------------------------------------------------------------------

R scripts for Decontamination, Quantitative Microbiome Profiling (QMP) and Mediation Analyses

Output files from QIIME2 including table.qza, rooted-tree.qza, and taxonomy.qza were run through the following R Scripts in the following order:

1. Decontam.R
2. QMP_Microeco_analysis.R
3. Mediation_Analysis.R
