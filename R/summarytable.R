# Summary Table Functions
# Functions used to :
# 1. Build output tables that summarize gain/loss events
# 2. Build output tables in useable/helpful output formats for non-programmatic purposes.

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

#' Makes Relative Copy-Number Gain/Loss Summary Tables
#'
#' @description
#'
#' MakeSummaryTable takes in a QDNAseq CN-object and returns tables containing agglomerated RCN gain/loss information across the provided samples. \cr
#' Prior to summarizing the gains/losses, the function applies a log2 transformation to center the data at 0. \cr
#' Optionally it also writes out a TSV table containing the summary.
#'
#' @param CNobj An S4 object of type QDNAseqCopyNumbers. This object must contain slots for: \cr
#' copynumber, segmented (standard CN slots) \cr
#' probgain, probdloss, probamp, probnorm, calls, probloss (from running CGHcall)
#' @param lowT Numeric (negative). Ignore gains/losses between this value and 0. ex. -0.1 \cr
#' This parameter helps filter out noise.
#' @param highT Numeric (positive). Ignore gains/losses between this value and 0. ex. 0.1 \cr
#' This parameter helps filter out noise.
#' @param pL (optional) Numeric. A probability threshold (between 0 and 1). Mask loss events where the probability does not rise to this threshold.
#' @param pG (optional) Numeric. A probability threshold (between 0 and 1). Mask gain events where the probability does not rise to this threshold.
#' @param prop (optional) Numeric. A proportion threshold (between 0 and 1). Exclude gain/loss events that occur in less than this proportion of samples.
#' @param find_peaks (optional) Logical. If TRUE, returns a table of peaks found in the summarized data.
#' @param snames (optional) Character vector. A vector of sample IDs that are of interest.
#' @param ref_genome (optional) String. The reference genome used for alignment. \cr
#' Options: 'hg19', 'hg38', and 'mm10'
#' @param save_path (optional) String. A path (directory) to where segment tables should be saved. ex. '~/Documents/test_project'.
#' @param use_calls (optional) Logical. If TRUE, gain/loss calls (e.g. CGHcall) are used to calculate gain/loss counts and proportions.
#' @returns A list of data frames. Always the summary_table, optionally a peaks data frame.
#'
#' @export
MakeSummaryTable <- function(CNobj,
                             lowT, highT,
                             pL = 0.1, pG = 0.1,
                             prop = 0.5,
                             find_peaks = FALSE,
                             snames = FALSE,
                             ref_genome = "hg19",
                             save_path = FALSE,
                             use_calls = FALSE) {

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

  # Determine mean copy number
  cns <- log2(CGHbase::segmented(CNobj))
  com_cns <- as.data.frame(cns) %>% dplyr::mutate(mean = rowMeans(dplyr::across(where(is.numeric))))
  df$mean.cn <- com_cns$mean
  
  # Calculate gains and losses by calls or segmented copy number
  if (use_calls) {
    bin_calls <- CGHbase::calls(CNobj)
    bin_loss <- as.data.frame(bin_calls < 0) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
    bin_gain <- as.data.frame(bin_calls > 0) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
  } else {
    bin_loss <- as.data.frame(cns < lowT) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
    bin_gain <- as.data.frame(cns > highT) %>% dplyr::mutate(total = rowSums(dplyr::across(where(is.logical))))
  }
  
  # Add gains/losses to output table
  df$loss.count <- bin_loss$total
  df$gain.count <- bin_gain$total
  df$loss.prop <- bin_loss$total/dim(cns)[2]
  df$gain.prop <- bin_gain$total/dim(cns)[2]

  # Add per-gain or loss sample names to output table
  mat <- as.data.frame(t(colnames(bin_loss)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:dplyr::n(), each = dim(bin_loss)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_loss[,c(1:dim(cns)[2])])
  mat[!mask] <- NA
  mat <- as.data.frame(mat) %>% tidyr::unite("names", everything(), sep = ',', na.rm = TRUE)
  df$loss.samples <- mat$names
  mat <- as.data.frame(t(colnames(bin_gain)[1:dim(cns)[2]])) %>% dplyr::slice(rep(1:dplyr::n(), each = dim(bin_gain)[1]))
  mat <- as.matrix(mat)
  mask <- as.matrix(bin_gain[,c(1:dim(cns)[2])])
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


#' Generate Human Readable Relative Copy-Number Profiles
#'
#' @description
#'
#' This function takes as input a QDNAseq or CGHcall copy-number object and gives back a long-format table with several useful columns.
#' These columns include; 'sample_id', 'chromosome', 'start', 'end', 'gain_probability',
#' 'loss_probability', 'relative_copy_number', 'bin_count', 'sum_of_bin_lengths',
#' 'cytobands', 'coordinates', and 'size'
#'
#' @param object S4 copy-number object - QDNAseq or CGHcall object
#' @param binsize The binsize used in the copy number object. ex. '30kb', '100kb'
#' @param ref_genome One of the common reference genomes: ex. 'hg19', 'mm10', or 'hg38'
#' @param save_dir (optional) The directory where the tables should be saved. ex. '~/Documents/test_project'
#' @returns Segment tables in long format (by sample id) ready to be written out to a table file (ex. tsv, csv).
#'
#' @export
GenHumanReadableRcnProfile <- function(object, binsize,
                                       ref_genome, save_dir = FALSE) {

  # Create collapsed segments table
  # Add Chromosome cytoband, coordinates (in bp), length of region, and gain or loss tag to each entry
  connection <- RMySQL::dbConnect(RMySQL::MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", dbname=ref_genome)
  cytobands <- DBI::dbGetQuery(conn=connection, statement="SELECT chrom, chromStart, chromEnd, name FROM cytoBand")
  cyto_ranges <- GenomicRanges::makeGRangesFromDataFrame(cytobands)

  segmented <- tidyr::gather(as.data.frame(CGHbase::segmented(object)), sample_id, segmented)
  if (class(object)[1] == 'cghCall') {
    segmented$gainP <- tidyr::gather(as.data.frame(CGHbase::probgain(object)), sample_id, probgain)$probgain
    segmented$ampP <- tidyr::gather(as.data.frame(CGHbase::probamp(object)), sample_id, probamp)$probamp
    segmented$lossP <- tidyr::gather(as.data.frame(CGHbase::probloss(object)), sample_id, probloss)$probloss
    segmented$dlossP <- tidyr::gather(as.data.frame(CGHbase::probdloss(object)), sample_id, probdloss)$probdloss
    segmented$chromosome <- rep(object@featureData@data[["Chromosome"]], dim(object)[2])
    segmented$start <- rep(object@featureData@data[["Start"]], dim(object)[2])
    segmented$end <- rep(object@featureData@data[["End"]], dim(object)[2])
  } else {
    segmented$chromosome <- rep(object@featureData@data[["chromosome"]], dim(object)[2])
    segmented$start <- rep(object@featureData@data[["start"]], dim(object)[2])
    segmented$end <- rep(object@featureData@data[["end"]], dim(object)[2])
    if ('probgain' %in% names(object@assayData)) {
      segmented$gainP <- tidyr::gather(as.data.frame(CGHbase::probgain(object)), sample_id, probgain)$probgain
      segmented$ampP <- tidyr::gather(as.data.frame(CGHbase::probamp(object)), sample_id, probamp)$probamp
      segmented$lossP <- tidyr::gather(as.data.frame(CGHbase::probloss(object)), sample_id, probloss)$probloss
      segmented$dlossP <- tidyr::gather(as.data.frame(CGHbase::probdloss(object)), sample_id, probdloss)$probdloss
    }
  }

  # collapse copy numbers down to segments
  collapsed_segs <- CopyNumberSegments(segmented)
  collapsed_segs$weight <- NULL

  collapsed_segs <- collapsed_segs %>% transform(chromosome = as.character(chromosome))
  collapsed_segs$chromosome[collapsed_segs$chromosome == 23] <- 'X'
  ranges <- GenomicRanges::makeGRangesFromDataFrame(collapsed_segs)
  GenomeInfoDb::seqlevelsStyle(ranges) <-'UCSC'
  hits <- GenomicRanges::findOverlaps(ranges, cyto_ranges)
  temp <- data.frame(ranges = hits@from, cyto_ranges = hits@to, cytobands = cytobands$name[hits@to])
  temp <- temp[,c('ranges','cytobands')]

  # collapse cytobands to a single cell
  temp <- temp %>%
    dplyr::group_by(ranges) %>%
    dplyr::summarise(cytobands = paste(cytobands, collapse = ",")) %>%
    dplyr::ungroup()
  for (i in c(1:dim(temp)[1])) {
    entry <- stringr::str_split(temp$cytobands[i], pattern = ',')[[1]]
    if ( length(entry) > 1 ) {
      temp$cytobands[i] <- paste0(entry[1], '-', entry[length(entry)])
    }
  }
  collapsed_segs$cytobands <- temp$cytobands
  collapsed_segs$coordinates <- paste0(GenomicRanges::seqnames(ranges), ':', GenomicRanges::ranges(ranges))
  collapsed_segs <- collapsed_segs %>% dplyr::mutate(size = end - start + 1)
  collapsed_segs$segment <- NULL
  colnames(collapsed_segs) <- c('sample_id', 'chromosome', 'start', 'end', 'gain_probability',
                                'loss_probability', 'relative_copy_number', 'bin_count',
                                'sum_of_bin_lengths', 'cytobands', 'coordinates', 'size')
  # save segment tables
  if ( save_dir != FALSE ) {
    dir.create(file.path(save_dir, binsize))
    for (i in unique(collapsed_segs$sample_id)) {
      temp <- collapsed_segs[collapsed_segs$sample_id == i, ]
      file_name <- paste0(save_dir, '/', binsize, '/', i, '_', binsize, '_RsegsCytobandsTable.tsv')
      write.table(temp, file = file_name, sep = '\t', col.names = TRUE, row.names = FALSE)
    }
  }

  DBI::dbDisconnect(connection)
  return(collapsed_segs)
}


#' Generate Human Readable Absolute Copy-Number Profiles
#'
#' @description
#'
#' This function takes as input a QDNAseq or CGHcall copy-number object and gives back Cytoband .tsv tables from absolute copy-number calls (in a human readable format).
#'
#' @param object S4 copy-number object - QDNAseq or CGHcall object
#' @param save_path A string. A path (directory) to where segment tables should be saved. ex. '~/Documents/test_project'
#' @returns A table. Collapsed Segment tables in long format (by sample id) ready to be written out to a table file (ex. tsv, csv).
#'
#' @export
GenHumanReadableAcnProfile <- function(object, save_path) {

  # Get Chromosome cytobands, coordinates (in bp), and lengths of regions
  connection <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", dbname="hg19")
  cytobands <- dbGetQuery(conn=connection, statement="SELECT chrom, chromStart, chromEnd, name FROM cytoBand")
  cyto_ranges <- GenomicRanges::makeGRangesFromDataFrame(cytobands)

  # Add cytobands
  collapsed_segs <- plyr::ldply(object, data.frame, .id = 'sample_id')
  collapsed_segs <- collapsed_segs %>% transform(chromosome = as.character(chromosome))
  collapsed_segs$chromosome[collapsed_segs$chromosome == 23] <- 'X'
  ranges <- GenomicRanges::makeGRangesFromDataFrame(collapsed_segs)
  seqlevelsStyle(ranges) <-'UCSC'
  hits <- findOverlaps(ranges, cyto_ranges)
  temp <- data.frame(ranges = hits@from, cyto_ranges = hits@to, cytobands = cytobands$name[hits@to])
  temp <- temp[,c('ranges','cytobands')]

  # collapse cytobands
  temp <- temp %>%
    dplyr::group_by(ranges) %>%
    dplyr::summarise(cytobands = paste(cytobands, collapse = ",")) %>%
    dplyr::ungroup()
  for (i in c(1:dim(temp)[1])) {
    entry <- str_split(temp$cytobands[i], pattern = ',')[[1]]
    if ( length(entry) > 1 ) {
      temp$cytobands[i] <- paste0(entry[1], '-', entry[length(entry)])
    }
  }

  # Re-format for output
  collapsed_segs$cytobands <- temp$cytobands
  collapsed_segs$coordinates <- paste0(seqnames(ranges), ':', ranges(ranges))
  collapsed_segs <- collapsed_segs %>% mutate(size = end - start + 1)
  collapsed_segs$segment <- NULL
  colnames(collapsed_segs) <- c('sample_id', 'chromosome', 'start', 'end',
                                'absolute_copy_number', 'cytobands',
                                'coordinates', 'size')
  collapsed_segs <- collapsed_segs %>% mutate(chromosome = replace(chromosome,
                                                                   chromosome=="X", '23')) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::arrange(as.integer(chromosome), .by_group = TRUE) %>%
    mutate(chromosome = replace(as.character(chromosome),
                                chromosome=='23', 'X'))
  # save segment tables
  dir.create(save_path)
  for (i in unique(collapsed_segs$sample_id)) {
    temp <- collapsed_segs[collapsed_segs$sample_id == i, ]
    write.table(temp, file = paste0(save_path, i, '_segsCytobandsTable.tsv'),
                sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  dbDisconnect(connection)
  return(collapsed_segs)
}

