# This script is used to combine CN profiles from a QDNAseq or CGHcall object.
# Then convert this table into a useable/helpful output format for wetlab researchers.

#' Collapses down table to segmented copy numbers
#'
#' @description
#'
#' CollapseTableToSegments transforms relative copy-number calls to segment tables.
#'
#' @param df A dataframe. A dataframe with copy number calls (with columns 'chromosome', 'start', 'end', 'mean.cn')
#' @returns A dataframe. A dataframe of segment summaries of various characteristics (derived from copy-number calls)
#'
#' @export
CollapseTableToSegments <- function(df) {

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
      gain_samples = dplyr::first(gain.samples),
      loss_samples = dplyr::first(loss.samples),
      bin_count = n(),
      sum_of_bin_lengths = sum(length)
    )
}

#' Makes summary tables
#'
#' @description
#'
#' MakeSummaryTable takes in a Copy Number object and writes a TSV table containing all the summaries.
#' It is the function that does the work.
#'
#' @param CNobj An object. A a ballgown object containing Copy Number data.
#' @param snames A vector. A vector of sample IDs that are of interest.
#' @param lowT A whole number. A threshold below which there is copy number loss.
#' @param highT A whole number. A threshold above which there is copy number gain.
#' @param pL A float. A probability threshold (between 0 and 1), to which loss frequency is compared.
#' @param pG A float. A probability threshold (between 0 and 1), to which gain frequency is compared.
#' @param prop A float. A proportion threshold (between 0 and 1), to which loss and gain occurences are compared (separately).
#' @param save_path A string. A path (directory) to where segment tables should be saved. ex. '~/Documents/test_project'.
#' @returns Nothing.
#'
#' @export
MakeSummaryTable <- function(CNobj, snames, lowT, highT, pL, pG, prop, save_path) {

  # Retrieve just the samples of interest
  CNobj <- CNobj[,sampleNames(CNobj) %in% snames]

  df <- data.frame(chromosome = chromosomes(CNobj),
                   start = bpstart(CNobj),
                   end = bpend(CNobj))

  cns <- segmented(CNobj)
  com_cns <- as.data.frame(cns) %>% mutate(mean = rowMeans(across(where(is.numeric))))
  df$mean.cn <- com_cns$mean

  uni.chrom <- unique(df$chrom)

  # Add gain / loss metrics per bin to output table
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

  # Add per-gain or loss sample names to output table
  mat <- as.data.frame(t(colnames(bin_cns_loss)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:n(), each = dim(bin_cns_loss)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_cns_loss[,c(1:dim(cns)[2])])
  mat[!mask] <- NA
  mat <- as.data.frame(mat) %>% tidyr::unite("names", everything(), sep = ',', na.rm = TRUE)
  df$loss.samples <- mat$names
  mat <- as.data.frame(t(colnames(bin_cns_gain)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:n(), each = dim(bin_cns_gain)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_cns_gain[,c(1:dim(cns)[2])])
  mat[!mask] <- NA
  mat <- as.data.frame(mat) %>% tidyr::unite("names", everything(), sep = ',', na.rm = TRUE)
  df$gain.samples <- mat$names
  df <- df %>% dplyr::filter((loss.freq > pL) | (gain.freq > pG),
                      (loss.count/ncol(cns) > prop) | (gain.count/ncol(cns) > prop),
                      (mean.cn < lowT) | (mean.cn > highT)
                      )

  # Collapse the data frame to genomic segments
  df <- CollapseTableToSegments(df)
  df <- df %>% transform(chromosome = as.character(chromosome))
  df$chromosome[df$chromosome == 23] <- 'X'
  df$coordinates <- paste0('chr', df$chromosome, ':', as.character(df$start), '-', as.character(df$end))
  df$segment <- NULL

  # Add genes present in the region to the output table
  mycoords.gr = df %>% dplyr::select(chromosome, start, end) %>%
                    mutate(chromosome=paste0('chr', chromosome)) %>%
                    makeGRangesFromDataFrame
  overlaps <- findOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                 single.strand.genes.only=FALSE), mycoords.gr)
  entrez_ids <- names(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                            single.strand.genes.only=FALSE)[overlaps@from])
  symbols <- getSYMBOL(entrez_ids, data='org.Hs.eg')
  genes_table <- data.frame(symbols = symbols, entrez_id = entrez_ids, entry = overlaps@to) %>%
                      dplyr::filter(!is.na(symbols))
  genes_table <- genes_table %>%
                    group_by(entry) %>%
                    summarise(features_in_region = paste0(symbols, collapse = ","))
  df$entry <- 1:dim(df)[1]
  df <- df %>% dplyr::left_join(genes_table, by = c('entry'))
  df$entry <- NULL

  # Save the output table
  write.table(df, file = save_path, sep = '\t', col.names = TRUE, row.names = FALSE)
}

