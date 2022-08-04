# get_acns_from_vafs_rascal.R
#
# Calculate Absolute Copy Numbers (ACNs) given Variant Allele Frequencies (Vafs) for several genes.
# Use the rascal package in R to do this transformation.
# Based on instructions in the vignette:
# https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(QDNAseq)
  library(CGHcall)
  library(EnsDb.Hsapiens.v75)
  library(readr)
})


###########################
### Functions
###########################

calculate_vaf_acns <- function (relative_cns, vafs) {
  
  qdnaseq_segs <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_comCNVfilt.tsv'
  qdnaseq_segs <- read.table(file = qdnaseq_segs, header = TRUE)
  qdnaseq_segments <- gather(qdnaseq_segs, sample, segmented, `CC.CHM.1341`:`YW.EC052`, factor_key=TRUE)    # Convert to long
  segments <- copy_number_segments(qdnaseq_segments)                            # Collapse to continuous segments
  
  rascal_batch_solutions <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_comCNVfilt_solutions_mad_optimal.csv", sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")
  
  vafs_cels <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/vafs_and_cellularities.tsv'
  vafs_cels <- data.table::fread(file = vafs_cels, header = TRUE, sep = "\t", fill = TRUE)
  vafs_cels$sample_id <- str_replace_all(vafs_cels$sample_id, "-", ".")
  # To simplify matters lets sum all unique vafs for a given gene/sample and treat them as a single vaf
  vafs_cels[,8:13] <- as.data.frame(apply(vafs_cels[,8:13], 
                                          MARGIN = c(1,2), 
                                          function(x) sum_delimited_elements(x, ';')))
  
  chosenSegmentTablesList <- list()
  j <- 1
  for (i in unique(rascal_batch_solutions$sample)) {
    
    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i)
    vafs <- vafs_cels %>% dplyr::filter(sample_id == i)
    if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN
    
    meta <- vafs %>% dplyr::select(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity)
    vafs <- vafs %>% dplyr::select(-c(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity))
    vaf_idx <- which.max(as.double(vafs[1,]))                                   # Choose max vaf
    if (length(vaf_idx) == 0) next                                              # If there isn't a vaf for the sample skip finding an ACN
    vaf <- as.double(vafs[1, ..vaf_idx])                                        # vaf
    vaf_name <- names(vafs[1, ..vaf_idx])                                         
    vaf_name <- toupper(str_split(vaf_name, pattern = '\\.')[[1]][1])           # vaf name
    
    rcn_obj <- sample_segments
    vaf_gene_rcn <- get_segmented_rcn_for_gene(rcn_obj, vaf_name)
    if (vaf_gene_rcn == FALSE) { message(paste('No segmented CN call found for', rcn_obj$sample[1], 'VAF gene.', sep = ' ')); next}                                          # If there isn't a CN for the VAF skip finding an ACN
    solution_set <- solutions %>%                                               # Get from running find best fit solution
      dplyr::select(ploidy, cellularity) %>%
      dplyr::mutate(tp53_absolute_copy_number = relative_to_absolute_copy_number(vaf_gene_rcn, ploidy, cellularity)) %>%
      dplyr::mutate(tp53_tumour_fraction = tumour_fraction(tp53_absolute_copy_number, cellularity))
    
    sol_idx <- which.min(abs(as.double(vaf) - (solution_set$tp53_tumour_fraction*100)))
    solution_set <- solution_set[,]
    absolute_segments <- mutate(sample_segments,
                                   copy_number = relative_to_absolute_copy_number(copy_number, 
                                          solution_set[sol_idx,]$ploidy, 
                                          solution_set[sol_idx,]$cellularity))
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome, start=start, end=end, segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j = j+1
  }
  
  # names(chosenSegmentTablesList) <- rascal_batch_solutions$sample
  saveRDS(chosenSegmentTablesList, file = "/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalVAF.rds")
}


# 1. Look up vaf gene start and end in genome reference
# 2. Filter rcn obj (ex. QDNAseq obj) for seg. CN at the vaf gene location
get_segmented_rcn_for_gene <- function (rcn_obj, gene) {

  granges_gene <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_name == gene)
  grRCN = makeGRangesFromDataFrame(rcn_obj)
  hits = findOverlaps(grRCN, granges_gene)
  
  if (length(hits) == 0) {
    return(FALSE)
  } else {
    overlaps <- pintersect(grRCN[queryHits(hits)], granges_gene[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(granges_gene[subjectHits(hits)])
    idx <- queryHits(hits)[which.max(percentOverlap)]
    return(rcn_obj$copy_number[idx])
  }
}

sum_delimited_elements <- function (element, delimiter) {
  elements <- str_split(element, pattern = delimiter)[[1]]
  element <- sum(as.numeric(elements))
  return(element)
}



# Add highest MAD column to rascal solutions table
# input_file: input path 
getOptimalMADSolutions <- function (input_file, rcn_file, segs_file, save_file = TRUE) {
  rascal_batch_solutions <- read.table(file = input_file, sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(input_file$sample, "-", ".")
  
  temp <- rascal_batch_solutions %>% 
    group_by(sample) %>% 
    mutate(mad_optimal = (min(distance) == distance)) 
  temp <- temp[temp$mad_optimal,] %>% 
    mutate(mad_optimal = (max(cellularity) == cellularity)) %>% 
    ungroup()
  rascal_batch_solutions <- left_join(rascal_batch_solutions, temp, by = c('sample', 'ploidy', 'cellularity', 'distance')) %>% 
    replace_na(list(mad_optimal = FALSE))
  
  if (save_file) {
    write.csv(rascal_batch_solutions, file = str_replace(input_file, 'solutions.csv', 'rascal_solutions.csv'))
  }
  
  qdnaseq_segs <- rcn_file
  qdnaseq_segs <- read.table(file = qdnaseq_segs, header = TRUE)
  qdnaseq_segments <- gather(qdnaseq_segs, sample, segmented, `CC.CHM.1341`:`YW.EC052`, factor_key=TRUE)
  segments <- copy_number_segments(qdnaseq_segments)                            # Collapse to continuous segments
  
  rascal_batch_solutions <- rascal_batch_solutions %>% dplyr::filter(mad_optimal == TRUE)
  chosenSegmentTablesList <- list()
  j <- 1
  for (i in rascal_batch_solutions$sample) {
    sample_segments <- dplyr::filter(segments, sample == i)
    absolute_segments <- dplyr::mutate(sample_segments, 
                                       copy_number = relative_to_absolute_copy_number(copy_number, 
                                                                                      rascal_batch_solutions$ploidy[j], 
                                                                                      rascal_batch_solutions$cellularity[j]))
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome, start=start, end=end, segVal=copy_number)
    chosenSegmentTablesList[[j]] <- as.data.frame(absolute_segments)
    j = j+1
  }
  names(chosenSegmentTablesList) <- rascal_batch_solutions$sample
  saveRDS(chosenSegmentTablesList, file = segs_file)
  
   
}


