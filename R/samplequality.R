########################### Sample quality functions ################################

#' Make a DF of CNs and Segments from S4 QDNAseq object
#'
#' @description
#' Function to create a dataframe containing both the relative copy numbers and the segmented calls from a QDNAseq object.
#'
CopySegFlat <- function(x) {
  copy_number <- x@assayData[["copynumber"]]
  segmented <- x@assayData[["segmented"]]
  copy_number <- as.data.frame(copy_number); segmented <- as.data.frame(segmented)
  copy_number <- copy_number %>%
    tibble::rownames_to_column("chr") %>%
    pivot_longer(cols = 2:(ncol(copy_number)+1), names_to = "sample")
  segmented <- segmented %>%
    tibble::rownames_to_column("chr") %>%
    pivot_longer(cols = 2:(ncol(segmented)+1), names_to = "sample")
  copy_number <- ChromosomeSplitPos(copy_number)
  segmented <- ChromosomeSplitPos(segmented)
  colnames(copy_number)[2] <- "copy_number"
  colnames(segmented)[2] <- "segmented"
  copy_number <- data.table::setDT(copy_number)
  segmented <- data.table::setDT(segmented)
  comb_table <- merge(copy_number, segmented, by = c("position", "sample"))
  comb_table <- comb_table %>%
    dplyr::select(chromosome = chromosome.x, sample, copy_number, segmented, start = start.y, end = end.y)
  return(comb_table)
}


#' Collapse segments
#'
#' @description Collapsed segmented dataframe
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number
#'
CollapsedSegs <- function(x) {
  x <- as.data.frame(x)
  stopifnot(is.data.frame(x))
  stopifnot("sample" %in% names(x))
  stopifnot("chromosome" %in% names(x))
  x[, c("start", "end", "segmented")] <- lapply(x[, c("start", "end", "segmented")], as.numeric)
  x <- x %>%
    dplyr::filter(!is.na(segmented)) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(sample, chromosome, start) %>%
    dplyr::mutate(new_segment = row_number() == 1 | !(sample == lag(sample) & chromosome == lag(chromosome) & segmented == lag(segmented))) %>%
    dplyr::mutate(segment = cumsum(new_segment))
  return(x)
}


#' Calculate segment sizes
#'
#' @description Calculate segment sizes
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, copy_number and new_segment (indicator variable for the segment the copy_number call belongs to)
#'
GetSegCounts <- function(x) {
  x$segment <- as.character(x$segment)
  x <- x %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(
      sample = dplyr::first(sample),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      copy_number = dplyr::first(segmented))
  x <- x %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(seg_counts = n())
  return(x)
}


#' Calculate median segment-level variance per sample
#'
#' @description Calculate median segment-level variances per sample (median absolute deviance)
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number and new_segment (indicator variable for the segment the copy_number call belongs to)
#'
MedSegVar <- function(x) {
  x <- x %>%
    dplyr::group_by(sample, segment) %>%
    dplyr::summarize(med_dev = median(abs(copy_number - segmented), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(median_sd = median(med_dev, na.rm = TRUE))
  return(x)
}


#' Extract sample grouping
#'
#' @description Extract sample grouping
#' @param x param_dat containing the following columns: sample and sample quality parameters
#'
SampleGrouping <- function(x) {
  x <- x %>%
    mutate(sample_group = case_when(
      str_detect(sample, "CC-CHM") ~ "CC-CHM",
      str_detect(sample, "CC-HAM") ~ "CC-HAM",
      str_detect(sample, "CC-JGH") ~ "CC-JGH",
      str_detect(sample, "CC-LAV") ~ "CC-LAV",
      str_detect(sample, "CC-NSH") ~ "CC-NSH",
      str_detect(sample, "CC-RJH") ~ "CC-RJH",
      str_detect(sample, "CC-SUN") ~ "CC-SUN",
      str_detect(sample, "CC-SSK") ~ "CC-SSK",
      str_detect(sample, "CC-VGH") ~ "CC-VGH",
      str_detect(sample, "CC-WPG") ~ "CC-WPG",
      str_detect(sample, "EC") ~ "EC",
      str_detect(sample, "VOA") ~ "VOA",
      str_detect(sample, "VS") ~ "VS",
      str_detect(sample, "YW") ~ "YW"
    ))
}


#' Calculate and make a quality call for relative copy number profiles
#'
#' @description
#' Function to classify the quality of relative copy number profile calls from QDNAseq or WiseCondorX
#' @param x *QDNASeq object* containing the relative copy number calls as well as the segmented relative copy number calls
#' @details
#' Expects the input of a QDNAseq object containing relative copy-number calls (such as from QDNASeq or WiseCondorX). \cr
#' Classifies relative copy number profiles as "High" or "Low" quality.
#' @return
#' A dataframe of sample quality decisions ("Low" or "High" quality sample)
#'
#' @export
GetSampleQualityDecision <- function(x, metric = "quantile", cutoff = 0.95) {
  comb_dat <- CopySegFlat(x)
  comb_collapsed <- CollapsedSegs(comb_dat)
  seg_counts <- GetSegCounts(comb_collapsed)
  median_vars <- MedSegVar(comb_collapsed)
  param_dat <- merge(seg_counts, median_vars, by = "sample")

  stopifnot(is.character(metric))

  if (is.character(metric)) {
    if (metric %in% c("quantile")) {
      seg_cut <- quantile(param_dat$seg_counts, cutoff)
      med_cut <- quantile(param_dat$median_sd, cutoff)
      param_iqr_dat <- param_dat %>%
        dplyr::mutate(decision = ifelse((seg_counts > seg_cut) | (median_sd > med_cut), "Low", "High"),
                      seg_cut = seg_cut,
                      med_cut = med_cut)
      # param_iqr_dat <- param_dat %>%
      #   dplyr::mutate(decision = ifelse(seg_counts > seg_cut & median_sd > med_cut, "Low", "High"))
      return(param_iqr_dat)
    } else {
      if (metric %in% c("DecisionTree")) {
        iso = solitude::isolationForest$new(sample_size = nrow(param_dat) - 1)
        iso$fit(param_dat[, c(2, 3)])
        scores_train = param_dat %>%
          iso$predict() %>%
          arrange(desc(anomaly_score))
        param_dat <- param_dat %>%
          dplyr::mutate(row_num = row_number())
        param_dat <- full_join(param_dat, scores_train, by = c("row_num" = "id"))
        param_iso_dat <- param_dat %>%
          dplyr::select(sample, seg_counts, median_sd, anomaly_score) %>%
          dplyr::mutate(decision = ifelse(anomaly_score > 0.59, "Low", "High"))
        return(param_iso_dat)
      }
    }
  }
}
