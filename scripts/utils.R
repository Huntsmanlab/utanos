# sWGS Utils
# Functions used in the analysis of shallow WGS data
# combined functions from:
# get_acns_from_vafs_rascal.R
# filter_segmented_CN_calls.R
# 

###########################
### Libraries to load
###########################
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
# calculate_vaf_acns
# getOptimalMADSolutions
# get_segmented_rcn_for_gene
# sum_delimited_elements
# minmax_column
#
# gene_cr
# annotation_cr
# segments_to_copy_number
# StringToGRanges
# GetGRangesFromEnsDb
#
# compare_bin_CNs
# removeBlacklist

#####

# Calculate Absolute Copy Numbers (ACNs) given Variant Allele Frequencies (Vafs) for several genes.
# Use the rascal package in R to do this transformation.
# Based on instructions in the vignette:
# https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd
# Parameters:
# relative_cns - 
# vafs - 
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

minmax_column <- function (df_col) {
  preproc2 <- preProcess(as.data.frame(df_col), method=c("range"))
  minmaxed <- predict(preproc2, as.data.frame(df_col))
  return(minmaxed[,1])
}

gene_cr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    if (sum(toupper(queryset$Gene) %in% toupper(line)) > 0) {
      queryset_matches <- c(queryset_matches, i)
    }
  }
  return(targetset[unique(queryset_matches),])
}

annotation_cr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    queryset_matches <- c(queryset_matches, which(toupper(queryset$gene_name) %in% toupper(line)))
  }
  return(queryset[queryset_matches,])
}

segments_to_copy_number <- function(segs, bin_size, genome = 'hg19', Xincluded = FALSE) {
  
  # Stop execution if we don't have the required input
  stopifnot(is.list(segs))
  stopifnot(length(segs) != 0)
  stopifnot("chromosome" %in% names(segs[[1]]))
  stopifnot("start" %in% names(segs[[1]]), is.numeric(segs[[1]]$start))
  stopifnot("end" %in% names(segs[[1]]), is.numeric(segs[[1]]$end))
  stopifnot("segVal" %in% names(segs[[1]]), is.numeric(segs[[1]]$segVal))
  # Create template binned genome 
  genome_chrs <- 22
  if (Xincluded) { genome_chrs <- 23}
  if (genome == 'hg38') {
    chroms <- getChromInfoFromUCSC("hg38", map.NCBI=TRUE) %>% 
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  } else {
    chroms <- getChromInfoFromUCSC("hg19") %>% 
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size)) 
  }
  chroms$chrom <- sub('chr', '', chroms$chrom)
  genome_template <- data.frame(chromosome = rep(chroms$chrom, chroms$nbins), 
                                start = rep(rep(1,dim(chroms)[1]), chroms$nbins),
                                end = rep(chroms$size, chroms$nbins)) %>% 
    dplyr::group_by(chromosome,start,end) %>% 
    dplyr::mutate(id = row_number()) %>% 
    dplyr::mutate(start = start + ((id-1) * bin_size), 
                  end = start + bin_size - 1) %>% 
    dplyr::select(-id) %>% 
    dplyr::ungroup()
  # Build dataframe of per-bin copy numbers
  out <- c()
  for (name in names(segs)) {
    sample <- segs[[name]] %>% dplyr::mutate(size = end - start + 1) %>% 
      dplyr::mutate(nbins = size/bin_size)
    sample <- as.data.frame(lapply(sample, rep, sample$nbins)) %>% 
      group_by(start,end,segVal) %>% 
      dplyr::mutate(id = row_number()) %>% 
      dplyr::mutate(chromosome = as.character(chromosome),
                    start = (floor(start/bin_size) * bin_size + 1) + ((id-1) * bin_size), 
                    end = start + bin_size - 1) %>% 
      dplyr::select(-c(size,nbins,id)) %>% 
      dplyr::ungroup()
    per_bin_copy_numbers <- left_join(genome_template, sample, by = c('chromosome', 'start', 'end'))
    per_bin_copy_numbers$sample_id <- name
    out <- rbind(out, per_bin_copy_numbers)
  }
  names(out) <- c("chromosome", "start", "end", "state", "sample_id")
  return(out)
}

StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

GetGRangesFromEnsDb <- function(
    ensdb,
    standard.chromosomes = TRUE,
    biotypes = c(listGenebiotypes(ensdb)),
    verbose = TRUE
) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }
  
  # extract genes from each chromosome
  if (verbose) {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
    })
  }
  
  # combine
  tx <- do.call(what = c, args = tx)
  tx <- tx[tx$gene_biotype %in% biotypes]
  return(tx)
}

# function for extracting the copy number for a given sample
# can handle QDNAseqCopyNumbers object or a copy number data frame
compare_bin_CNs <- function(objs, sample, bin_area) {
  outlist <- vector(mode = "list", length = length(objs))
  for (i in 1:length(objs)) {
    obj <- objs[[i]][,sample]
    copy_number_values <- Biobase::assayDataElement(obj, "copynumber")[,1]
    segmented_values <- Biobase::assayDataElement(obj, "segmented")[,1]
    df <- Biobase::fData(obj) %>%
      rename_with(tolower) %>%
      tibble::rownames_to_column(var = "id") %>%
      as_tibble() %>%
      dplyr::select(id, chromosome, start, end) %>%
      dplyr::mutate_at(vars(start, end), as.integer) %>%
      dplyr::mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
      dplyr::mutate(sample = sample) %>%
      dplyr::mutate(copy_number = copy_number_values) %>%
      dplyr::mutate(segmented = segmented_values) %>%
      dplyr::select(sample, chromosome, start, end, copy_number, segmented)
    outlist[[i]] <- df %>% dplyr::filter((start < (bin_area+75000)) & (start > (bin_area-75000)))
  }
  return(outlist)
}

### Remove the blacklist regions from hg19 per-bin dataframe
# DESCRIPTION
# Parameters:
# data - 
removeBlacklist <- function(data) {
  # Read in blacklist file
  blacklist = as.data.frame(data.table::fread(file = "~/repos/cnsignatures/data/external_datasets/binBlacklist.txt", sep = ' ', header = TRUE))
  
  # Convert blacklist and data to GRanges objects and find indices of overlaps
  grBL = makeGRangesFromDataFrame(blacklist)
  grData = makeGRangesFromDataFrame(data)
  overlaps = findOverlaps(grData, grBL)
  
  # Replace state of data rows with chr and start overlapping with blacklist
  data$state = replace(data$state, overlaps@from, NA)
  
  # Return original dataframe with state of blacklisted segments of genome identified as low mapability
  return(data)
}
