# Create Britroc sample plots and rds files from original paper
# June 15th, 2021
# 
# Input:
# All data is grabbed straight from the paper.
# This code is all copy and pasted from the manuscript repo. 
# Just placed in the needed order.

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(NMF)
  library(flexmix)
  library(YAPSA)
  library(pheatmap)
  library(viridis)
  library(hrbrthemes)
})

source(file = 'scripts/main_functions.R')
source(file = 'scripts/helper_functions.R')

# Read in ploidy and purity estimates obtained from Battenberg plots.
deep_metrics <- read.csv("manuscript_Rmarkdown/data/britroc_deepWGS_ploidy_purity_49.csv", header=T, as.is = T)
deep_metrics$purity <- deep_metrics$purity/100
deep_metrics <- deep_metrics %>%
  mutate(name= sub("[0-9]+_tumor_","",x=wgs_ID)) %>%
  dplyr::select(name,ploidy, purity) %>%
  tidyr::gather("metric", "deep", -name)

# Get ploidy for sWGS samples
CN <- readRDS("data/input_copy_numbers/britroc_absolute_copynumber.rds")
Biobase::pData(CN)$ploidy <- getPloidy(CN) %>%
  .$out

shallow_metrics <- Biobase::pData(CN) %>%
  dplyr::select(name, ploidy, purity) %>%
  filter(name %in% deep_metrics$name) %>%
  tidyr::gather("metric", "shallow", -name)

# Annotate sWGS samples with star_rating
samp_annotation_all <- read.csv("manuscript_Rmarkdown/data/britroc_sample_data.csv", as.is=T)
samp_annotation <- samp_annotation_all %>% 
  filter(IM.JBLAB_ID %in% shallow_metrics$name & Failed !="Y") %>%
  dplyr::select(IM.JBLAB_ID , star_rating)

all_CN<-readRDS("data/input_copy_numbers/britroc_absolute_copynumber.rds")
# all_CN<-all_CN[,colnames(all_CN)%in%samp_annotation[!samp_annotation$Failed=="Y","IM.JBLAB_ID"]]

#extract 3 star CN from each case
result<- samp_annotation_all %>% 
  filter(!samp_annotation_all$Failed=="Y") %>%
  dplyr::select(Britroc_No,IM.JBLAB_ID,star_rating) %>%
  dplyr::group_by(Britroc_No) %>%
  dplyr::slice(which.max(star_rating))
ids<-result %>% filter(star_rating==3)
ids<-ids$IM.JBLAB_ID
hq_CN<-all_CN[,colnames(all_CN)%in%ids]
ids<-result %>% filter(star_rating==2)
ids<-ids$IM.JBLAB_ID
lq_CN<-all_CN[,colnames(all_CN)%in%ids]

#extract copy-number features
CN_features<-extractCopynumberFeatures(hq_CN)
britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features)
component_by_signature <- generateSignatures(britroc_sample_component_matrix, 7)
saveRDS(component_by_signature, file = paste0("data/pilot_output/QDNA_60kb_britroc_component_by_signature.rds"))

pdf(file = paste0("pilot_plots/compareSignatures/CS_heatmap_QDNA_britroc.pdf"), 7, 7, onefile = FALSE)
basismap(component_by_signature, Rowv = NA, main = "Signature x Component matrix")
dev.off()
