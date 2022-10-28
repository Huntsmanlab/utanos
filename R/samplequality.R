########################### Sample quality functions ################################
#'
#' @description Function to classify the quality of relative copy number profile calls from QDNAseq or WiseCondorX  
#' @details The sample quality function supports an input of an  QDNAseq object from relative copy number callers for shallow whole genome sequencing data (such as QDNASeq or WiseCondorX) and classifies relative copy number profiles as "High" or "Low" quality
#' @param x *QDNASeq object* containing the relative copy number calls as well as the segmented relative copy number calls
#' 

#' @description Function to format chromosome and position information 
ChromosomeSplitPos <- function(x) {
  x$chr <- str_split_fixed(x$chr, ":", n = 2)
  x$chromosome <- x$chr[, 1]
  x$pos <- str_split_fixed(x$chr[, 2], "-", n = 2)
  x$start <- x$pos[, 1]
  x$end <- x$pos[, 2]
  x$position <- x$chr
  x$chr <- NULL
  x$pos <- NULL
  x[, c("chromosome", "start", "end")] <- lapply(x[, c("chromosome", "start", "end")], as.numeric)
  return(x)
}

#' @description Function to create a dataframe containing both the relative copy numbers and the segmented calls from a QDNASeq object
CopySegFlat <- function(x) {
  copy_number <- x@assayData[["copynumber"]]
  segmented <- x@assayData[["segmented"]]
  copy_number <- as.data.frame(copy_number); segmented <- as.data.frame(segmented)
  copy_number <- copy_number %>%
    dplyr::rownames_to_column("chr") %>%
    pivot_longer(cols = 2:ncol(copy_number), names_to = "sample")
  segmented <- segmented %>%
    dplyr::rownames_to_column("chr") %>%
    pivot_longer(cols = 2:ncol(segmented), names_to = "sample")
  copy_number <- ChromosomeSplitPos(copy_number)
  segmented <- ChromosomeSplitPos(segmented)
  colnames(copy_number)[2] <- "copy_number"
  colnames(segmented)[2] <- "segmented"
  copy_number$pos <- with(copy_number, interaction(sample, position))
  segmented$pos <- with(segmented, interaction(sample, position))
  copy_number <- data.table::setDT(copy_number)
  segmented <- data.table::setDT(segmented)
  comb_table <- copy_number[segmented, on = .(pos)]
  comb_table <- comb_table %>%
    select(chromosome, sample, copy_number, segmented, start = i.start, end = i.end) 
  return(comb_table)
}

#' @description Collapsed segmented dataframe
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number

CollapsedSegs <- function(x) {
  stopifnot(is.data.frame(x))
  stopifnot("sample" %in% names(x))
  stopifnot("chromosome" %in% names(x))
  stopifnot("start" %in% names(x), is.numeric(x))
  stopifnot("end" %in% names(x), is.numeric(x))
  stopifnot("segmented" %in% names(x), is.numeric(x))
  x %>%
    dplyr::filter(!is.na(segmented)) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(sample, chromosome, start) %>%
    dplyr::mutate(new_segment = row_number() == 1 | !(sample == lag(sample) & chromosome == lag(chromosome) & segmented == lag(segmented))) %>%
    dplyr::mutate(segment = cumsum(new_segment))
      return(x)
}

#' @description Calculate segment sizes 
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, copy_number and new_segment (indicator variable for the segment the copy_number call belongs to)

GetSegSizes <- function(x) {
  segmented_sizes <- x %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(
      sample = dplyr::first(sample),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      copy_number = dplyr::first(segmented))
  segmented_sizes <- segmented_sizes %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(seg_sizes = n())
  return(segment_sizes)
}

#' @description Calculate median segment-level variances per sample 
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number
#' 
MedSegVar <- function(x) {
  median_vars <- x %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(
      sample = dplyr::first(sample),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      copy_number = dplyr::first(segmented))
  median_vars <- median_vars %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median_var = var(copy_number, na.rm = TRUE))
  return(median_vars)
}

#' @description Output sample quality decision ("Low" or "High" quality sample)
#' @param x *QDNASeq object* containing the relative copy number calls as well as the segmented relative copy number calls
GetSampleQualityDecision <- function(x, metric = "IQR") {
  comb_dat <- CopySegFlat(x)
  comb_collapsed <- CollapsedSegs(comb_dat)
  seg_sizes <- GetSegSizes(comb_collapsed)
  median_vars <- MedSegVar(comb_collapsed)
  param_dat <- merge(seg_sizes, median_vars, by = "sample")
  if (metric == "SegmedIQR") {
    segment_cutoff <- quantile(param_dat$seg_sizes, 0.75) + 1.5*IQR(param_dat$seg_sizes)
    median_cutoff <- quantile(param_dat$median_var, 0.75) + 1.5*IQR(param_dat$median_var)
    param_dat <- param_dat %>%
      mutate(ifelse(seg_sizes > segment_cutoff & median_var > median_cutoff, "Low", "High"))
  }
  else if (metric == "SegIQR") {
    segment_cutoff <- quantile(param_dat$seg_sizes, 0.75) + 1.5*IQR(param_dat$seg_sizes)
    param_dat <- param_dat %>%
      mutate(ifelse(seg_sizes > segment_cutoff, "Low", "High"))
  }
  else if(metric == "DecisionTree") {
    iso = solitude::isolationForest$new()
    scores_train = param_dat %>%
      solitude::iso$predict() %>%
      arrange(desc(anomaly_score))
    param_dat <- param_dat %>%
      mutate(row_num = row_number())
    param_dat <- merge(param_dat, scores_train, by = c("row_num" = "id"))
    param_dat <- param_dat %>%
      mutate(decision = ifelse(anomaly-score > 0.60, "Low", "High"))
  }
  return(param_dat)
}



