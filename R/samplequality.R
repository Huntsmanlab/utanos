########################### Sample quality functions ################################

#' Make a DF of CNs and Segments from S4 QDNAseq object
#'
#' @description
#' Function to create a dataframe containing both the relative copy numbers and the segmented calls from a QDNAseq object.
CopySegFlat <- function(x) {
  copy_number <- x@assayData[["copynumber"]]
  segmented <- x@assayData[["segmented"]]
  copy_number <- as.data.frame(copy_number); segmented <- as.data.frame(segmented)
  copy_number <- copy_number %>%
    tibble::rownames_to_column("chr") %>%
    tidyr::pivot_longer(cols = 2:(ncol(copy_number)+1), names_to = "sample_id")
  segmented <- segmented %>%
    tibble::rownames_to_column("chr") %>%
    tidyr::pivot_longer(cols = 2:(ncol(segmented)+1), names_to = "sample_id")
  copy_number <- ChromosomeSplitPos(copy_number)
  segmented <- ChromosomeSplitPos(segmented)
  colnames(copy_number)[2] <- "copy_number"
  colnames(segmented)[2] <- "segmented"
  copy_number <- data.table::setDT(copy_number)
  segmented <- data.table::setDT(segmented)
  comb_table <- merge(copy_number, segmented, by = c("position", "sample_id"))
  comb_table <- comb_table %>%
    dplyr::select(chromosome = chromosome.x, sample_id,
                  copy_number, segmented, start = start.y, end = end.y)
  return(comb_table)
}

#' Collapse segments
#'
#' @description Collapsed segmented dataframe
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number
CollapsedSegs <- function(x) {
  x <- as.data.frame(x)
  stopifnot(is.data.frame(x))
  stopifnot("sample_id" %in% names(x))
  stopifnot("chromosome" %in% names(x))
  x[, c("start", "end", "segmented")] <- lapply(x[, c("start", "end", "segmented")], as.numeric)
  x <- x %>%
    dplyr::filter(!is.na(segmented)) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(sample_id, chromosome, start) %>%
    dplyr::mutate(new_segment = dplyr::row_number() == 1 |
                    !(sample_id == dplyr::lag(sample_id) &
                        chromosome == dplyr::lag(chromosome) &
                        segmented == dplyr::lag(segmented))) %>%
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
      sample_id = dplyr::first(sample_id),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      copy_number = dplyr::first(segmented))
  x <- x %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(seg_counts = dplyr::n())
  return(x)
}

#' Calculate median segment-level variance per sample
#'
#' @description Calculate median segment-level variances per sample (median absolute deviance)
#' @param x Dataframe containing the following columns: chromosome, start, stop, segmented, and copy_number and new_segment (indicator variable for the segment the copy_number call belongs to)
#'
MedSegVar <- function(x) {
  x <- x %>%
    dplyr::group_by(sample_id, segment) %>%
    dplyr::summarize(med_dev = median(abs(copy_number - segmented), na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(median_sd = median(med_dev, na.rm = TRUE))
  return(x)
}

#' Extract sample grouping
#'
#' @description Extract sample grouping
#' @param x param_dat containing the following columns: sample_id and sample quality parameters
#'
SampleGrouping <- function(x) {
  x <- x %>%
    dplyr::mutate(sample_group = case_when(
                  str_detect(sample_id, "CC-CHM") ~ "CC-CHM",
                  str_detect(sample_id, "CC-HAM") ~ "CC-HAM",
                  str_detect(sample_id, "CC-JGH") ~ "CC-JGH",
                  str_detect(sample_id, "CC-LAV") ~ "CC-LAV",
                  str_detect(sample_id, "CC-NSH") ~ "CC-NSH",
                  str_detect(sample_id, "CC-RJH") ~ "CC-RJH",
                  str_detect(sample_id, "CC-SUN") ~ "CC-SUN",
                  str_detect(sample_id, "CC-SSK") ~ "CC-SSK",
                  str_detect(sample_id, "CC-VGH") ~ "CC-VGH",
                  str_detect(sample_id, "CC-WPG") ~ "CC-WPG",
                  str_detect(sample_id, "EC") ~ "EC",
                  str_detect(sample_id, "VOA") ~ "VOA",
                  str_detect(sample_id, "VS") ~ "VS",
                  str_detect(sample_id, "YW") ~ "YW"
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
  param_dat <- merge(seg_counts, median_vars, by = "sample_id")

  stopifnot(is.character(metric))
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
      if(metric %in% c("DecisionTree")) {
        iso = solitude::isolationForest$new(sample_size = nrow(param_dat) - 1)
        iso$fit(param_dat[, c(2, 3)])
        scores_train = param_dat %>%
          iso$predict() %>%
          arrange(desc(anomaly_score))
        param_dat <- param_dat %>%
          dplyr::mutate(row_num = dplyr::row_number())
        param_dat <- dplyr::full_join(param_dat, scores_train, by = c("row_num" = "id"))
        param_iso_dat <- param_dat %>%
          dplyr::select(sample_id, seg_counts, median_sd, anomaly_score) %>%
          dplyr::mutate(decision = ifelse(anomaly_score > 0.59, "Low", "High"))
        return(param_iso_dat)
      }
    }
}


#' Remove or set a mask for relative CN bins (QDNAseq object)
#'
#' @description
#' Remove or set a mask for relative CN bins
#'
#' @param cnobj *QDNASeq object* containing the relative copy number calls as well as the segmented relative copy number calls
#' @param filter_by Optional *Dataframe* containing genomic regions to filter/mask. One region expected per row. \cr
#' Expected columns: \cr
#' chr start end etc. etc.
#' @param genome A character *string* of the reference genome used. ex "hg38"
#' @param minimum_overlap An *integer* designating the number minimum number of overlapping bases for a bin to be removed.
#' @param removebins A *Logical*. If TRUE then remove bins overlapping the indicated regions.
#' @param maskgaps A *Logical*. If TRUE then remove bins overlapping the centromeric and telomeric regions. Be sure to indicate the correct reference genome.
#' @param maskcomCNVs A *Logical*. If TRUE then remove bins overlapping commonly copy-number variant regions in this organism. Be sure to indicate the correct reference genome.
#' @param maskname A character *string* for the name to be used in case of setting a mask. Set this field if providing a DF to the filter_by parameter.
#'
#' @details
#' Expects the input of a QDNAseq object containing relative copy-number calls (such as from QDNASeq or WiseCondorX). \cr
#'
#' @return
#' A QDNAseq object filtered or with masks set.
#'
#' @export
FilterCNs <- function (cnobj,
                       filter_by = NULL,
                       genome = "hg19",
                       minimum_overlap = 5000,
                       removebins = FALSE,
                       maskgaps = FALSE,
                       maskcomCNVs = FALSE,
                       maskname = "mask1") {

  stopifnot(dim(cnobj) > 0)
  stopifnot(!missing(filter_by) | maskgaps | maskcomCNVs)
  stopifnot(is.logical(removebins))

  if (!missing(filter_by)) {
    stopifnot("start" %in% colnames(filter_by))
    stopifnot("end" %in% colnames(filter_by))

    cnobj_ranges <- StringToGRanges(featureNames(cnobj), sep = c(':', '-'))
    GenomeInfoDb::seqlevelsStyle(cnobj_ranges) <- 'UCSC'
    GenomeInfoDb::genome(cnobj_ranges) <- genome
    filter_by <- filter_by %>% dplyr::filter(chr != "")
    gfilter_by <- GenomicRanges::makeGRangesFromDataFrame(filter_by)
    GenomeInfoDb::seqlevelsStyle(gfilter_by) <- 'UCSC'
    GenomeInfoDb::genome(gfilter_by) <- genome
    hits <- GenomicRanges::findOverlaps(cnobj_ranges, gfilter_by,
                                        minoverlap = minimum_overlap,
                                        ignore.strand = TRUE)
    if (removebins) {
      cnobj <- cnobj[setdiff(1:length(cnobj_ranges), unique(hits@from)),]
    } else {
      cnobj@featureData@data[[maskname]] <- !(1:length(cnobj_ranges) %in% unique(hits@from))
      new_row <- data.frame(labelDescription = "A user-set mask.",
                            stringsAsFactors = FALSE,
                            row.names = maskname)
      cnobj@featureData@varMetadata <- rbind(cnobj@featureData@varMetadata, new_row)
    }
  }

  if (maskgaps) {
    if (genome == "hg19") { regions_to_mask <- gaps.hg19 }
    if (genome == "hg38") { regions_to_mask <- gaps.hg38 }
    cnobj_ranges <- StringToGRanges(featureNames(cnobj), sep = c(':', '-'))
    GenomeInfoDb::seqlevelsStyle(cnobj_ranges) <- 'UCSC'
    GenomeInfoDb::genome(cnobj_ranges) <- genome
    gr_mask <- GenomicRanges::makeGRangesFromDataFrame(regions_to_mask)
    GenomeInfoDb::genome(gr_mask) <- genome
    hits <- GenomicRanges::findOverlaps(cnobj_ranges, gr_mask,
                                        minoverlap = minimum_overlap,
                                        ignore.strand = TRUE)
    if (removebins) {
      cnobj <- cnobj[setdiff(1:length(cnobj_ranges), unique(hits@from)),]
    } else {
      cnobj@featureData@data[["centro.telo.mask"]] <- !(1:length(cnobj_ranges) %in% unique(hits@from))
      new_row <- data.frame(labelDescription = "A mask for centromeres and telomeres.",
                            stringsAsFactors = FALSE,
                            row.names = "centro.telo.mask")
      cnobj@featureData@varMetadata <- rbind(cnobj@featureData@varMetadata, new_row)
    }
  }

  if (maskcomCNVs) {
    if (genome == "hg19") { regions_to_mask <- hg19.comCNVs.P01S10KB }
    if (genome == "hg38") { regions_to_mask <- hg38.comCNVs.P01S10KB }
    cnobj_ranges <- StringToGRanges(featureNames(cnobj), sep = c(':', '-'))
    GenomeInfoDb::seqlevelsStyle(cnobj_ranges) <- 'NCBI'
    GenomeInfoDb::genome(cnobj_ranges) <- genome
    regions_to_mask <- regions_to_mask %>% dplyr::filter(chr != "")
    gr_mask <- GenomicRanges::makeGRangesFromDataFrame(regions_to_mask)
    GenomeInfoDb::genome(gr_mask) <- genome
    hits <- GenomicRanges::findOverlaps(cnobj_ranges, gr_mask,
                                        minoverlap = minimum_overlap,
                                        ignore.strand = TRUE)
    if (removebins) {
      cnobj <- cnobj[setdiff(1:length(cnobj_ranges), unique(hits@from)),]
    } else {
      cnobj@featureData@data[["comCNV.mask"]] <- !(1:length(cnobj_ranges) %in% unique(hits@from))
      new_row <- data.frame(labelDescription = "A mask for regions commonly CN-aberrant in humans.",
                            stringsAsFactors = FALSE,
                            row.names = "comCNV.mask")
      cnobj@featureData@varMetadata <- rbind(cnobj@featureData@varMetadata, new_row)
    }
  }

  return(cnobj)
}


