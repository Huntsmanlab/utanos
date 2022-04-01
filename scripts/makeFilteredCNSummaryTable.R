# This script is used to combine CN profiles from a QDNAseq or CGHcall object.
# Then convert this table into a useable/helpful output format for wetlab researchers.

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(CGHcall)
})

# Collapse down to segments
collapse_table_to_segments <- function(df) {
  
  stopifnot(is.data.frame(df))
  stopifnot("chromosome" %in% names(df))
  stopifnot("start" %in% names(df), is.numeric(df$start))
  stopifnot("end" %in% names(df), is.numeric(df$end))
  stopifnot("mean.cn" %in% names(df), is.numeric(df$mean.cn))
  
  df %>%
    dplyr::filter(!is.na(mean.cn)) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(chromosome, start) %>%
    dplyr::mutate(new_segment = row_number() == 1 | !(chromosome == lag(chromosome) & mean.cn == lag(mean.cn))) %>%
    dplyr::mutate(segment = cumsum(new_segment)) %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      mean_copy_number = dplyr::first(mean.cn),
      gain_probability = dplyr::first(gain.freq),
      loss_probability = dplyr::first(loss.freq),
      gain_count = dplyr::first(gain.count),
      loss_count = dplyr::first(loss.count),
      bin_count = n(),
      sum_of_bin_lengths = sum(length)
    )
}

# Function that does the work
make_summary_table <- function(CNobj, snames, lowT, highT, pL, pG, prop, save_path) {
  
  # Retrieve just the samples of interest
  CNobj <- CNobj[,sampleNames(CNobj) %in% snames]
  
  df <- data.frame(chromosome = chromosomes(CNobj),
                   start = bpstart(CNobj),
                   end = bpend(CNobj))
  
  cns <- segmented(CNobj)
  com_cns <- as.data.frame(cns) %>% mutate(mean = rowMeans(across(where(is.numeric))))
  df$mean.cn <- com_cns$mean
  
  uni.chrom <- unique(df$chrom)
  
  nclass <-3
  if (!is.null(probamp(CNobj))) nclass <- nclass+1 
  if (!is.null(probdloss(CNobj))) nclass <- nclass+1 
  if(nclass==3) {df$loss.freq <- rowMeans(probloss(CNobj)); df$gain.freq <- rowMeans(probgain(CNobj))}
  if(nclass==4) {df$loss.freq <- rowMeans(probloss(CNobj)); df$gain.freq <- rowMeans(probgain(CNobj))+rowMeans(probamp(CNobj))}
  if(nclass==5) {df$loss.freq <- rowMeans(probloss(CNobj))+rowMeans(probdloss(CNobj)); df$gain.freq <- rowMeans(probgain(CNobj))+rowMeans(probamp(CNobj))}
  
  bin_cns_loss <- as.data.frame(cns < lowT) %>% dplyr::mutate(total = rowSums(across(where(is.logical))))
  bin_cns_gain <- as.data.frame(cns > highT) %>% dplyr::mutate(total = rowSums(across(where(is.logical))))
  df$loss.count <- bin_cns_loss$total
  df$gain.count <- bin_cns_gain$total
  browser()
  df <- df %>% dplyr::filter((loss.freq > pL) | (gain.freq > pG),
                      (loss.count/ncol(cns) > prop) | (gain.count/ncol(cns) > prop),
                      (mean.cn < lowT) | (mean.cn > highT)
                      )
  df <- collapse_table_to_segments(df)
  
  df <- df %>% transform(chromosome = as.character(chromosome))
  df$chromosome[df$chromosome == 23] <- 'X'
  df$coordinates <- paste0('chr', df$chromosome, ':', as.character(df$start), '-', as.character(df$end))
  df$segment <- NULL
  
  write.table(df, file = save_path, sep = '\t', col.names = TRUE, row.names = FALSE)
}


# Input
ccnvFILT_obj <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/15kb_rCN_comCNVfilt.rds')
# Output
save_path <- '~/Documents/projects/cn_signatures_shallowWGS/data/summaryTables/NSMP_samples_batches1-13.tsv'
# List of samples for which data is desired
# Pay attention to the formatting! Its important!
samples_of_interest <- c('CC-VGH-1176', 'CC-VGH-1187', 
                         'CC-VGH-1190', 'CC-VGH-1196', 
                         'CC-VGH-1239', 'CC-VGH-1264',
                         'CC-VGH-1293', 'CC-VGH-1302',
                         'CC-VGH-0033a', 'CC-VGH-1181')

# Probability of gain or loss threshold - i.e. declare the minimum mean probability of loss or gain across samples
# Expl. For region Chr1:2850001-2865000, and probabilities of loss for 10 samples of:
# 0.814 0 0 0 0 0.125 0.975 0 0 0
# The mean probability of loss would be 0.1914
prob_loss <- 0.2
prob_gain <- 0.2
# Log space relative copy number thresholds
low_threshold <- -0.2
high_threshold <- 0.2
# Presence threshold - i.e. In what minimum proportion of samples does the gain or loss need to be present?
proportion_threshold <- 0.2

# Run command
make_summary_table(ccnvFILT_obj, samples_of_interest, low_threshold, high_threshold, prob_loss, prob_gain, proportion_threshold, save_path)
