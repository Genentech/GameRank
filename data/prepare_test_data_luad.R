#
# Script to generate test dataset pulled from TCGA 
#
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )
library( tibble )
library( labelled )

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("curatedTCGAData")
library( curatedTCGAData )
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAutils")
library( TCGAutils )

#   Study Abbreviation Study Name
# ---------------------------------------------------------------------------
# 1  ACC  Adrenocortical Carcinoma
# 2  BLCA Bladder Urothelial Carcinoma
# 3  BRCA Breast Invasive Carcinoma
# 4  CESC Cervical Squamous Cell Carcinoma And Endocervical Adenocarcinoma
# 5  CHOL Cholangiocarcinoma
# 6  CNTL Controls
# 7  COAD Colon Adenocarcinoma
# 8  DLBC Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
# 9  ESCA Esophageal Carcinoma
# 10 GBM  Glioblastoma Multiforme
# 11 HNSC Head And Neck Squamous Cell Carcinoma
# 12 KICH Kidney Chromophobe
# 13 KIRC Kidney Renal Clear Cell Carcinoma
# 14 KIRP Kidney Renal Papillary Cell Carcinoma
# 15 LAML Acute Myeloid Leukemia
# 16 LGG  Brain Lower Grade Glioma
# 17 LIHC Liver Hepatocellular Carcinoma
# 18 LUAD Lung Adenocarcinoma
# 19 LUSC Lung Squamous Cell Carcinoma
# 20 MESO Mesothelioma
# 21 OV   Ovarian Serous Cystadenocarcinoma
# 22 PAAD Pancreatic Adenocarcinoma
# 23 PCPG Pheochromocytoma And Paraganglioma
# 24 PRAD Prostate Adenocarcinoma
# 25 READ Rectum Adenocarcinoma
# 26 SARC Sarcoma
# 27 SKCM Skin Cutaneous Melanoma
# 28 STAD Stomach Adenocarcinoma
# 29 TGCT Testicular Germ Cell Tumors
# 30 THCA Thyroid Carcinoma
# 31 THYM Thymoma
# 32 UCEC Uterine Corpus Endometrial Carcinoma
# 33 UCS  Uterine Carcinosarcoma
# 34 UVM  Uveal Melanoma

# ExperimentList data types   Description
# ----------------------------------------------------------------------------
# SummarizedExperiment*
# RNASeqGene                RSEM TPM gene expression values
# RNASeq2GeneNorm           Upper quartile normalized RSEM TPM gene expression values
# miRNAArray                Probe-level  miRNA expression values
# miRNASeqGene              Gene-level log2 RPM miRNA expression values
# mRNAArray                 Unified gene-level mRNA expression values
# mRNAArray_huex            Gene-level mRNA expression values from Affymetrix Human Exon Array
# mRNAArray_TX_g4502a       Gene-level mRNA expression values from Agilent 244K Array
# mRNAArray_TX_ht_hg_u133a  Gene-level mRNA expression values from Affymetrix Human Genome U133 Array
# GISTIC_AllByGene          Gene-level GISTIC2 copy number values
# GISTIC_ThresholdedByGene  Gene-level GISTIC2 thresholded discrete copy number values
# RPPAArray                 Reverse Phase Protein Array normalized protein expression values
#
# RangedSummarizedExperiment
# GISTIC_Peaks              GISTIC2 thresholded discrete copy number values in recurrent peak regions
#
# SummarizedExperiment with HDF5Array DelayedMatrix
# Methylation_methyl27      Probe-level methylation beta values from Illumina HumanMethylation 27K      BeadChip
# Methylation_methyl450     Probe-level methylation beta values from Infinium HumanMethylation 450K     BeadChip
#
# RaggedExperiment
# CNASNP                    Segmented somatic Copy Number Alteration calls from SNP array
# CNVSNP                    Segmented germline Copy Number Variant calls from SNP Array
# CNASeq                    Segmented somatic Copy Number Alteration calls from low pass DNA Sequencing
# Mutation*                 Somatic mutations calls
# CNACGH_CGH_hg_244a        Segmented somatic Copy Number Alteration calls from CGH Agilent Microarray 244A
# CNACGH_CGH_hg_415k_g4124a Segmented somatic Copy Number Alteration calls from CGH Agilent Microarray 415K
# * All can be converted to RangedSummarizedExperiment (except RPPAArray) with TCGAutils

#
# Use Copy Number Aberrations (CNA) and Singular Nucleotide Variations (SNV) for features
#

dir_out <- "~/GameRank/data/"

# Lung ----
file_rdata <- file.path( dir_out, "tcga_luad_cna_cnv.Rdata" )
dc <- c("LUAD")
ay <- c("CNASNP","CNVSNP")
dorun <- TRUE
mae <- curatedTCGAData( diseaseCode = dc, assays = ay, dry.run = !dorun )

var_clin <- getClinicalNames( diseaseCode = dc )
cod <- colData( mae ) %>% as_tibble 
cod <- cod %>% 
  dplyr::select( all_of( intersect( cod %>% colnames, c("patientID", var_clin ) ) ) )
cod 
mae # N = 518

# Simplyfy the RaggedExperiments to RangedSummarizedExperiments
# RaggedExperiment mutation objects become a genes by patients RangedSummarizedExperiment object containing '1' if there is a 
# non-silent mutation somewhere in the gene, and '0' otherwise. "CNA" and "CNV" segmented copy number are reduced using a weighted
# mean in the rare cases of overlapping (non-disjoint) copy number regions.
# 
# These functions rely on 'TxDb.Hsapiens.UCSC.hg19.knownGene' and 'org.Hs.eg.db' to map to the 'hg19' NCBI build.
mae <- simplifyTCGA( mae )
df_cna <- longFormat( mae[["LUAD_CNASNP-20160128_simplified"]] ) %>%  as_tibble %>% 
  transmute( patientID = colname, rowname = sprintf( "%s_cna", rowname ), value ) %>%
  pivot_wider( names_from = "rowname", values_from = "value", values_fill = list( value = NA_real_ ) )
df_cna 

df_cnv <- longFormat( mae[["LUAD_CNVSNP-20160128_simplified"]] ) %>% as_tibble %>%
  transmute( patientID = colname, rowname = sprintf( "%s_cnv", rowname ), value ) %>%
  pivot_wider( names_from = "rowname", values_from = "value", values_fill = list( value = NA_real_ ) )
df_cnv

dat <- cod %>%
  dplyr::left_join( df_cna, "patientID" ) %>%
  dplyr::left_join( df_cnv, "patientID" )
dat 

lst_keys <- c("patientID")
lst_outcomes <- c("days_to_death", "vital_status" )
lst_meta <- c("days_to_last_followup","date_of_initial_pathologic_diagnosis","tumor_tissue_site","histological_type")
lst_vars <- var_clin %>%
  union( union( df_cna %>% colnames, df_cnv %>% colnames ) ) %>% 
  setdiff( lst_keys ) %>%
  setdiff( lst_outcomes ) %>%
  setdiff( lst_meta )
lst_vars

dat <- dat %>% dplyr::select( all_of( intersect(c(lst_keys, lst_meta, lst_outcomes, lst_vars), dat %>% colnames ) ) )
save( dat, lst_keys, lst_outcomes, lst_meta, lst_vars, file = file_rdata )

