# Compare Copy Number Signatures Across Cohorts
# May 19th, 2021
# 
# Input:
# Change from a list of Component-by-Signature Matrices from NMFfit objects to the NMFfit objects themselves.
# This change will permit for easier plotting.
# This means that NMFfit objects must be saved in the cn_signatures_pilot.R script.

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

source(file = '~/repositories/cnsignatures/scripts/main_functions.R')
source(file = '~/repositories/cnsignatures/scripts/helper_functions.R')


compare_signature_by_component <- function (cs_matrix_list, mat_names, saveFile) {
  
  nsig <- 7
  # britroc feat_sig matrix
  feat_sig_mat<-basis(cs_matrix_list[[1]])
  colnames(cs_matrix_list[[1]]@fit@W) <- c('s7', 's1', 's6', 's4', 's3', 's2', 's5')
  reord_1 <- as.integer(c(2,6,5,4,7,3,1))
  names(reord_1)<-paste0("s",1:7)
  feat_sig_mat<-feat_sig_mat[,reord_1]
  colnames(feat_sig_mat)<-paste0("s",1:nsig)
  sig_feat_mat<-t(feat_sig_mat)
  
  # Second feat_sig matrix
  feat_sig_mat_2<-basis(cs_matrix_list[[2]])
  colnames(cs_matrix_list[[2]]@fit@W) <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')
  sig_feat_mat_2<-t(feat_sig_mat_2)
  colnames(feat_sig_mat_2)<-paste0("s",1:nsig)
  
  # Third feat_sig matrix
  feat_sig_mat_3<-basis(cs_matrix_list[[3]])
  colnames(cs_matrix_list[[3]]@fit@W) <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')
  sig_feat_mat_3<-t(feat_sig_mat_3)
  colnames(feat_sig_mat_3)<-paste0("s",1:nsig)
  
  # browser()  
  # determine matching signatures and their correlation
  # reord_2<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_2,2,cor,x,method="s"))})
  # sigcor_2<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_2,2,cor,x,method="s"))})
  reord_2 <- as.integer(c(2,3,4,5,6,7,1))
  
  # reord_3<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_3,2,cor,x,method="s"))})
  # sigcor_3<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_3,2,cor,x,method="s"))})
  reord_3 <- as.integer(c(1,5,2,7,3,6,4))
  
  # plot the feat_sig matrices side by side
  pdf(file = saveFile, 10, 7, onefile = FALSE)
  par(mfrow=c(1,3))
  basismap(cs_matrix_list[[1]],Rowv=NA,Colv=reord_1,main=mat_names[1],tracks=NA)
  basismap(cs_matrix_list[[2]],Rowv=NA,Colv=reord_2,main=mat_names[2],tracks=NA)
  basismap(cs_matrix_list[[3]],Rowv=NA,Colv=reord_3,main=mat_names[3],tracks=NA)
  dev.off()
}

# Save File
output_image = "plots/britroc_plots/compareSignatures/rel_to_abs_CN_signatures_comparison_heatmap.pdf"

# Input
mat1 <- readRDS(file = "data/britroc_output/QDNA_30kb_component_by_signature_absbritroc.rds")
mat2 <- readRDS(file = "data/britroc_output/QDNA_30kb_custSigs_7sigs_component_by_signature_rbritroc91.rds")
mat3 <- readRDS(file = "data/britroc_output/QDNA_30kb_custSigs_7sigs_component_by_signature_rbritroc199.rds")

# Names of Input
mat_names <- c('abs_Britroc_3star', 'rel_Britroc_3star', 'rel_Britroc_2&3star')

compare_signature_by_component(list(mat1, mat2, mat3), mat_names, saveFile = output_image)
