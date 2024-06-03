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
    dplyr::mutate(new_segment = dplyr::row_number() == 1 | !(chromosome == dplyr::lag(chromosome) &
                                                               mean.cn == dplyr::lag(mean.cn))) %>%
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
      gain_proportion = dplyr::first(gain.prop),
      loss_proportion = dplyr::first(loss.prop),
      gain_samples = dplyr::first(gain.samples),
      loss_samples = dplyr::first(loss.samples),
      bin_count = dplyr::n(),
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
#' @param lowT A whole number. A threshold below which there is copy number loss.
#' @param highT A whole number. A threshold above which there is copy number gain.
#' @param pL A float. A probability threshold (between 0 and 1), to which loss frequency is compared.
#' @param pG A float. A probability threshold (between 0 and 1), to which gain frequency is compared.
#' @param prop A float. A proportion threshold (between 0 and 1), to which loss and gain occurrences are compared (separately).
#' @param find_peaks Boolean. If TRUE, returns a table
#' @param snames A vector. A vector of sample IDs that are of interest.
#' @param ref_genome A string. The reference genome used for alignment.
#' @param save_path A string. A path (directory) to where segment tables should be saved. ex. '~/Documents/test_project'.
#' @returns Nothing.
#'
#' @export
MakeSummaryTable <- function(CNobj,
                             lowT, highT,
                             pL = 0.1, pG = 0.1,
                             prop = 0.5,
                             find_peaks = FALSE,
                             snames = FALSE,
                             ref_genome = "hg19",
                             save_path = FALSE) {

  # Initial checks
  if (!("EnsDb.Hsapiens.v75" %in% installed.packages()) &&
      !("EnsDb.Hsapiens.v86" %in% installed.packages()) &&
      !("EnsDb.Mmusculus.v79" %in% installed.packages())) {
    stop(
      "An annotation package must be installed to use this function. Ex. 'EnsDb.Hsapiens.v75'",
      call. = FALSE)
  }

  output <- list()

  # Retrieve just the samples of interest
  if (snames != FALSE) {
    CNobj <- CNobj[,sampleNames(CNobj) %in% snames]
  }

  df <- data.frame(chromosome = factor(QDNAseq::chromosomes(CNobj),
                                       levels = unique(QDNAseq::chromosomes(CNobj))),
                   start = QDNAseq::bpstart(CNobj),
                   end = QDNAseq::bpend(CNobj))

  cns <- log2(CGHbase::segmented(CNobj))
  com_cns <- as.data.frame(cns) %>% dplyr::mutate(mean = rowMeans(dplyr::across(where(is.numeric))))
  df$mean.cn <- com_cns$mean

  # Add gain / loss metrics per bin to output table
  nclass <- 3
  if (!is.null(CGHbase::probamp(CNobj))) { nclass <- nclass + 1 }
  if (!is.null(CGHbase::probdloss(CNobj))) { nclass <- nclass + 1 }
  if(nclass==3) {
    df$loss.freq <- rowMeans(CGHbase::probloss(CNobj))
    df$gain.freq <- rowMeans(CGHbase::probgain(CNobj)) }
  if(nclass==4) {
    df$loss.freq <- rowMeans(CGHbase::probloss(CNobj))
    df$gain.freq <- rowMeans(CGHbase::probgain(CNobj)) + rowMeans(CGHbase::probamp(CNobj)) }
  if(nclass==5) {
    df$loss.freq <- rowMeans(CGHbase::probloss(CNobj)) + rowMeans(CGHbase::probdloss(CNobj))
    df$gain.freq <- rowMeans(CGHbase::probgain(CNobj)) + rowMeans(CGHbase::probamp(CNobj)) }

  bin_cns_loss <- as.data.frame(cns < lowT) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
  bin_cns_gain <- as.data.frame(cns > highT) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
  df$loss.count <- bin_cns_loss$total
  df$gain.count <- bin_cns_gain$total
  df$loss.prop <- bin_cns_loss$total/dim(cns)[2]
  df$gain.prop <- bin_cns_gain$total/dim(cns)[2]

  # Add per-gain or loss sample names to output table
  mat <- as.data.frame(t(colnames(bin_cns_loss)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:dplyr::n(), each = dim(bin_cns_loss)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_cns_loss[,c(1:dim(cns)[2])])
  mat[!mask] <- NA
  mat <- as.data.frame(mat) %>% tidyr::unite("names", everything(), sep = ',', na.rm = TRUE)
  df$loss.samples <- mat$names
  mat <- as.data.frame(t(colnames(bin_cns_gain)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:dplyr::n(), each = dim(bin_cns_gain)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_cns_gain[,c(1:dim(cns)[2])])
  mat[!mask] <- NA
  mat <- as.data.frame(mat) %>% tidyr::unite("names", everything(), sep = ',', na.rm = TRUE)
  df$gain.samples <- mat$names
  df <- df %>% dplyr::filter((loss.freq > pL) | (gain.freq > pG),
                      (loss.count/ncol(cns) > prop) | (gain.count/ncol(cns) > prop),
                      (mean.cn < lowT) | (mean.cn > highT)
                      )

  # Find the 'peaks' in these proportions.
  # They can be interpreted as shortest overlap regions (SORs) for multiple samples.
  if (find_peaks != FALSE) {
    if (!("gsignal" %in% installed.packages())) {
      stop("The gsignal package must be installed to use the peaks option.",
           call. = FALSE)
    }
    gain_peaks <- gsignal::findpeaks(data = df$gain.freq,
                                     MinPeakDistance = 1, MinPeakHeight = pG)
    loss_peaks <- gsignal::findpeaks(data = df$loss.freq,
                                     MinPeakDistance = 1, MinPeakHeight = pL)
    peaks <- df[c(gain_peaks$loc, loss_peaks$loc),]
    output[['peaks']] <- peaks
  }

  # Collapse the data frame to genomic segments
  df <- CollapseTableToSegments(df)
  df$segment <- NULL

  # Final data re-formatting and add annotations
  if (ref_genome == "hg19") {
    genome_annotation <- annotables::grch37
    gr_genes <- GenomicFeatures::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  } else if (ref_genome == "hg38") {
    genome_annotation <- annotables::grch38
    gr_genes <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  } else if (ref_genome == "mm10") {
    genome_annotation <- annotables::grcm38
    gr_genes <- GenomicFeatures::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  } else {
    stop("Please provide an accepted reference genome.")
  }
  df <- df %>% dplyr::mutate(chromosome = as.character(chromosome))
  df$coordinates <- paste0('chr', df$chromosome, ':',
                           as.character(df$start), '-', as.character(df$end))
  mycoords.gr <- df %>% dplyr::select(chromosome, start, end)
  mycoords.gr <- GenomicRanges::makeGRangesFromDataFrame(mycoords.gr)
  overlaps <- IRanges::findOverlaps(gr_genes, mycoords.gr)
  genes_table <- data.frame(ENSEMBL = gr_genes@ranges@NAMES[overlaps@from], entry = overlaps@to)
  genes_table <- genes_table %>% dplyr::left_join(genome_annotation,
                                                  by = c("ENSEMBL" = "ensgene")) %>%
                    dplyr::select(ENSEMBL, entry, symbol) %>%
                    dplyr::filter(!is.na(symbol)) %>%
                    dplyr::group_by(entry) %>%
                    dplyr::summarise(features_in_region = paste0(symbol, collapse = ","))
  df$entry <- 1:dim(df)[1]
  df <- df %>% dplyr::left_join(genes_table, by = c('entry'))
  df$entry <- NULL
  output[['summary_table']] <- df

  # Optionally, save the output table
  if (save_path != FALSE) {
    write.table(df, file = save_path, sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  return(output)
}

