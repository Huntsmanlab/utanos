# Functions for scaling relative copy number to absolute and finding the
# best fitting solutions based on distance metrics.
#
# Many functions in this file are copied from another respository (rascal) with
# modifications.
# Original code can be found here: https://github.com/crukci-bioinformatics/rascal

###########################
### Functions
###########################
# CalculateACNs
# GenVafAcns
# GenKploidyAcns
# GenMadAcns
# GetSegmentedRcnForGene
# ReplaceQDNAseqAssaySlots
# AbsoluteCopyNumberDistance
# AbsoluteCopyNumberDistanceGrid
# IsAcceptablePloidyAndCellularity
# FindBestFitSolutions
# AbsoluteCopyNumberScalingFactor
# RelativeToAbsoluteCopyNumber
# TumourFraction
# SumDelimitedElements
# MinmaxColumn



#'  Find best fitting solutions for each sample
#'
#' @description
#' FindRascalSolutions() finds the best ploidy and cellularity pairs from the relative copy number profiles of each sample.
#' Makes use of the find_best_fit_solutions() function from rascal.
#' 
#' @param cnobj An S4 object of type QDNAseqCopyNumbers.
#' @param min_ploidy,max_ploidy the range of ploidies.
#' @param ploidy_step the stepwise increment of ploidies along the grid.
#' @param min_cellularity,max_cellularity the range of cellularities.
#' @param cellularity_step the stepwise increment of cellularities along the
#' grid.
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @param distance_filter_scale_factor the distance threshold above which
#' solutions will be discarded as a multiple of the solution with the smallest
#' distance.
#' @param max_proportion_zero the maximum proportion of fitted absolute copy
#' number values in the zero copy number state.
#' @param min_proportion_close_to_whole_number the minimum proportion of fitted
#' absolute copy number values sufficiently close to a whole number.
#' @param max_distance_from_whole_number the maximum distance from a whole
#' number that a fitted absolute copy number can be to be considered
#' sufficiently close.
#' @param solution_proximity_threshold how close two solutions can be before one
#' will be filtered; reduces the number of best fit solutions where there are
#' many minima in close proximity.
#' @param keep_all set to \code{TRUE} to return all solutions but with
#' additional \code{best_fit} column to indicate which are the local minima that
#' are acceptable solutions (may be useful to avoid computing the distance grid
#' twice)
#' 
#' @return a dataframe containing the calculated solutions. Each solution includes:
#' sample, ploidy, celluarity, and distance.
#'
#' @export
FindRascalSolutions <- function(cnobj,
                                min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
                                min_cellularity = 0.2, max_cellularity = 1.0, cellularity_step = 0.01,
                                distance_function = c("MAD", "RMSD"),
                                distance_filter_scale_factor = 1.25,
                                max_proportion_zero = 0.05,
                                min_proportion_close_to_whole_number = 0.5,
                                max_distance_from_whole_number = 0.15,
                                solution_proximity_threshold = 5,
                                keep_all = FALSE) {
  # convert obj to dataframe
  rcn = ExportBinsQDNAObj(cnobj, type = "segments") %>%
    subset(select = -feature)
  row.names(rcn) <- NULL
  
  pivoted <- tidyr::pivot_longer(rcn, 
                                 cols = !(chromosome:end),
                                 names_to = "sample",
                                 values_to = "segmented", 
                                 values_drop_na = TRUE) %>%
    dplyr::select(sample, tidyr::everything())
  
  # centering
  relative_bin_and_segments_centered <- pivoted %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(dplyr::across(c(segmented), ~ . / stats::median(segmented, na.rm = TRUE))) %>%
    dplyr::ungroup()
  
  # collapse bins
  segments <- CopyNumberSegments(relative_bin_and_segments_centered)
  
  # find best fit soln
  solutions <- segments %>% dplyr::group_by(sample) %>% 
    dplyr::group_modify(~ FindBestFitSolutions(
      .x$copy_number, .x$weight, 
      min_ploidy = min_ploidy, max_ploidy = max_ploidy, ploidy_step = ploidy_step, 
      min_cellularity = min_cellularity, max_cellularity = max_cellularity, cellularity_step = cellularity_step, 
      max_proportion_zero = max_proportion_zero, distance_function = distance_function, 
      solution_proximity_threshold = solution_proximity_threshold,
      min_proportion_close_to_whole_number = min_proportion_close_to_whole_number, 
      keep_all = keep_all
    )) %>% 
    dplyr::ungroup() %>%
    as.data.frame()

  return(solutions)
}


#' Calculate Absolute Copy Numbers
#'
#' @description
#'
#' CalculateACNs() calculates absolute copy numbers (ACNs) from the relative copy number profiles of one or more samples.
#' There are several included options by which to do this.
#' Note: If not providing a table of variant allele frequencies (VAFs) then 'mad' is the only method available. \cr
#'
#' This function makes use of the rascal package in R and instructions can be found in the vignette: \cr
#' https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd  \cr
#'
#' @param cnobj An S4 object of type QDNAseqCopyNumbers. This object must contain a copynumber slot and a segmented slot.
#' @param rascal_sols A tsv or dataframe of the calculated rascal solutions.
#' @param acnmethod The method by which to calculate ACNs. Can be one of: \cr
#' "maxvaf" - Calculate ACNs assuming the maximum discovered VAF for the sample is an appropriate representation for the tumour fraction. \cr
#' char vector - Same as above but rather than using the max vaf provide a character vector of the genes from which to pull VAFs.
#' The genes are assumed to be in order of decreasing precedence. ex. c('TP53', 'KRAS', 'PTEN') \cr
#' "mad" - Calculate ACNs using the mean absolute difference (MAD) column from the solutions table. \cr
#' If using variant allele frequencies from targeted panel sequencing or some other technology: \cr
#' - The variants must must be in a datatable/dataframe. \cr
#' - Required columns: sample_id, chromosome, start, end, gene_name, ref, alt, vaf.
#' - Each row of said table must correspond to a unique variant. \cr
#' - Each variant must have an associated variant allele frequency. \cr
#' - Each row must also be associated with a specific sample. \cr
#' A dataframe of known ploidies per sample - Calculate ACNs at each bin given we already know the ploidy of the samples. \cr
#' @param variants (optional) Dataframe or string. A table of variants containing variant allele frequencies per gene and per sample.
#' @param acn_save_path (optional) String. The output path (absolute path recommended) where to save the result.
#' @param return_sols (optional) Logical. Return the selected rascal solution.
#' @param return_S4 (optional) Logical.
#' @param distance_decimal_places number of decimal places to compare to when finding
#' the minimum distance.
#' 
#' @returns A list containing: \cr
#' 1. A list of dataframes (one for each sample) OR a QDNAseq S4 object. \cr
#' 2. Optionally, a dataframe of the rascal solutions. \cr
#' This function returns ACNs when a rascal solution is found, if one isn't, the sample is skipped.
#' @details
#' ```
#' solutions <- "~/Documents/.../rascal_solutions.csv"
#' rcn.obj <- "~/repos/utanosmodellingdata/sample_copy_number_data/sample_filtered_cn_data.rds"
#' variants <- "~/Documents/.../allvariants.clinvar.cosmic.exons.csv"
#' save_path <- "~/Documents/.../rascal_ACN_segments.rds"
#' variants <- data.table::fread(file = variants, sep = ',')
#' variants <- variants %>% dplyr::rename(chromosome = chr,
#'                                        gene_name = genecode_gene_name)
#' outputs <- CalculateACNs(rcn.obj,
#'                          acn.method = 'maxvaf',
#'                          rascal_sols = solutions,
#'                          variants = variants,
#'                          return_sols = TRUE,
#'                          return_S4 = TRUE)
#' output <- CalculateACNs(rcn.obj,
#'                         rascal_sols = solutions,
#'                         variants = variants,
#'                         acnmethod = c('TP53', 'KRAS', 'BRCA1',
#'                                       'BRCA2', 'PIK3CA', 'PTEN'),
#'                         acn_save_path = save_path)
#' output <- CalculateACNs(rcn.obj,
#'                         rascal_sols = solutions,
#'                         acnmethod = 'mad')
#'
#' ```
#'
#' @export

CalculateACNs <- function (cnobj, acnmethod,
                           rascal_sols = NULL,
                           variants = NULL,
                           acn_save_path = FALSE,
                           return_sols = FALSE,
                           return_S4 = FALSE,
                           distance_decimal_places = 7) {

  output <- list()
  stopifnot(dim(cnobj) > 0)
  stopifnot(inherits(cnobj, c("QDNAseqCopyNumbers", "QDNAseqReadCounts")))

  # If needed, read-in rascal solutions
  if (!is.null(rascal_sols)) {
    if (inherits(rascal_sols, 'character')) {
      rascal_sols <- read.table(file = rascal_sols, sep = ',', header = TRUE)
    } else if (inherits(rascal_sols, 'data.frame')) {
      NULL
    } else {
      stop("Invalid value passed to the 'rascal_sols' parameter. \n
           Please supply either a path or dataframe.")
    }
  }

  # Pull out cns to dataframes, both bin-wise and segmented
  wide_segs <- ExportBinsQDNAObj(cnobj, type = "segments", filter = FALSE)
  long_segs <- wide_segs %>% tidyr::gather(sample, segmented,
                                           5:dim(.)[2], factor_key=TRUE)
  segments <- CopyNumberSegments(long_segs)
  wide_cns <- ExportBinsQDNAObj(cnobj, type = "copynumber", filter = FALSE)
  long_cns <- wide_cns %>% tidyr::gather(sample, copy_number,
                                         5:dim(.)[2], factor_key=TRUE)
  long_cns <- long_cns %>%
    dplyr::mutate(segmented = long_segs$segmented)          # Create combined data.frame

  # If needed, read-in variants
  if (!is.null(variants)) {
    if ('data.frame' %in% class(variants)) {
      NULL
    } else if (class(variants)[1] == 'character') {
      variants <- data.table::fread(file = variants, header = TRUE,
                                    sep = ",", fill = TRUE)
    } else {
      stop("Invalid value passed to the 'variants' parameter. \n
           Please supply either a path or dataframe.")
    }
    variants <- variants %>%
      dplyr::select(chromosome, start, end, sample_id, gene_name,
                    starts_with('ref.'), starts_with('alt.'),
                    starts_with('vaf')) %>%
      dplyr::mutate(vaf = dplyr::select(., starts_with('vaf')) %>%
                      rowSums(na.rm = TRUE)) %>%
      dplyr::group_by(sample_id, gene_name) %>%
      dplyr::summarise(sample_id = dplyr::first(sample_id),
                       gene_name = dplyr::first(gene_name),
                       vaf = max(vaf))
  }

  # Control-flow for acn scaling method
  if ((length(acnmethod) == 1) && (acnmethod == 'maxvaf')) {
    variants <- variants %>%
      dplyr::group_by(sample_id) %>% dplyr::filter(vaf == max(vaf))
    variants <- variants[!duplicated(variants$sample_id),]
    output <- GenVafAcns(segments, rascal_sols, variants,
                         return_sols, long_cns, return_S4)

  } else if ('data.frame' %in% class(acnmethod)) {
    if (any(!(acnmethod$sample_id %in% unique(relative_segs$sample))) == TRUE) {
      if (sum(!(acnmethod$sample_id %in% unique(relative_segs$sample))) == dim(acnmethod)[1]) {
        stop("Sample IDs do not match between ploidy and segment tables.")
      } else {
        warning("There appear to be ploidy sample IDs for which there are no segment tables.")
      }
    }
    output <- GenKploidyAcns(segments, acnmethod, return_sols, long_cns, return_S4)

  } else if (length(acnmethod) > 1) {
    variants <- variants %>% dplyr::group_by(sample_id) %>%
      dplyr::filter(gene_name %in% acnmethod | vaf == max(vaf, na.rm = TRUE)) %>%
      dplyr::filter(gene_name == c(intersect(acnmethod, gene_name),
                                   setdiff(gene_name, acnmethod))[1])
    output <- GenVafAcns(segments, rascal_sols, variants,
                         return_sols, long_cns, return_S4)

  } else if (acnmethod == 'mad') {
    output <- GenMadAcns(segments, rascal_sols, return_sols,
                         long_cns, return_S4, distance_decimal_places)

  } else {
    stop("Invalid value passed to the 'acnmethod' parameter. \n
         Please supply an accepted option.")
  }

  if (return_S4) {
    cnobj <- ReplaceQDNAseqAssaySlots(cnobj, output[["acns"]], output[["acns_segs"]])
    output[['acn.obj']] <- cnobj
    output[["acns"]] <- NULL
    output[["acns_segs"]] <- NULL
  }

  if (acn_save_path != FALSE) {
    saveRDS(output[[1]], file = acn_save_path)
  }
  return(output)
}


# Function that calculates absolute copy numbers using the rascal package and vafs.
GenVafAcns <- function (segments,
                        rascal_batch_solutions,
                        variants,
                        return_sols = FALSE,
                        cns = FALSE,
                        return_S4 = FALSE) {

  # Initialize variables
  output <- list()
  chosenSegmentTablesList <- list()
  j <- 1
  acns <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  acns_segs <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  sols <- data.frame(sample_id = character(),
                     ploidy = double(),
                     cellularity = double(),
                     tumour_fraction = double(),
                     vaf = double(),
                     stringsAsFactors = FALSE)

  # Loop through samples
  for (i in unique(rascal_batch_solutions$sample)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i)
    vafs <- variants %>% dplyr::filter(sample_id == i)
    if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN
    rcn_obj <- sample_segments
    suppressWarnings({
      vaf_gene_rcn <- GetSegmentedRcnForGene(rcn_obj, vafs$gene_name)
    })
    if (vaf_gene_rcn == FALSE) { message(paste('No segmented CN call found for',
                                               rcn_obj$sample[1], 'VAF gene.',
                                               sep = ' ')); next}               # If there isn't a CN for the VAF skip finding an ACN
    solution_set <- solutions %>%                                               # Get from running find best fit solution
      dplyr::select(ploidy, cellularity) %>%
      dplyr::mutate(absolute_copy_number = RelativeToAbsoluteCopyNumber(vaf_gene_rcn, ploidy, cellularity)) %>%
      dplyr::mutate(tumour_fraction = TumourFraction(absolute_copy_number, cellularity))

    sol_idx <- which.min(abs(as.double(vafs$vaf) - (solution_set$tumour_fraction*100)))
    solution_set <- solution_set[,]
    absolute_segments <- dplyr::mutate(sample_segments,
                                copy_number = RelativeToAbsoluteCopyNumber(copy_number,
                                                                               solution_set[sol_idx,]$ploidy,
                                                                               solution_set[sol_idx,]$cellularity))
    sols <- rbind(sols, c(i,
                  as.numeric(solution_set[sol_idx,]$ploidy),
                  as.numeric(solution_set[sol_idx,]$cellularity),
                  round(as.numeric(solution_set[sol_idx,]$tumour_fraction),
                        digits = 3),
                  as.numeric(vafs$vaf)))

    # Optionally, grab needed values to create an S4 QDNAseq object
    if (return_S4) {
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- RelativeToAbsoluteCopyNumber(absolute_cns$copy_number,
                                                                           solution_set[sol_idx,]$ploidy,
                                                                           solution_set[sol_idx,]$cellularity)
      absolute_cns$segmented <- RelativeToAbsoluteCopyNumber(absolute_cns$segmented,
                                                                           solution_set[sol_idx,]$ploidy,
                                                                           solution_set[sol_idx,]$cellularity)
      acns <- cbind(acns, absolute_cns$copy_number)
      acns_segs <- cbind(acns_segs, absolute_cns$segmented)
    }
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome,
                                                             start=start,
                                                             end=end,
                                                             segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j <- j+1
  }

  # Return collapsed acn segments list or an S4 object if wanted
  if (return_S4) {
    colnames(acns) <- unique(rascal_batch_solutions$sample)
    colnames(acns_segs) <- unique(unique(rascal_batch_solutions$sample))
    rownames(acns) <- cns$feature[1:dim(acns)[1]]
    rownames(acns_segs) <- cns$feature[1:dim(acns)[1]]
    output[["acns"]] <- acns
    output[["acns_segs"]] <- acns_segs
  } else {
    output[["acn_segment_tables"]] <- chosenSegmentTablesList
  }

  # Return selected ploidy/cellularity/vaf combos if wanted
  if (return_sols) {
    colnames(sols) <- c('sample_id', 'ploidy', 'cellularity',
                        'tumour_fraction', 'vaf')
    output[['rascal_solutions']] <- as.data.frame(sols)
  }
  return(output)
}


# Function that calculates absolute copy numbers using the rascal package and a DF of known ploidies.
GenKploidyAcns <- function (segments,
                            ploidies,
                            return_sols = FALSE,
                            cns = FALSE,
                            return_S4 = FALSE) {

  chosenSegmentTablesList <- list()
  j <- 1
  plots <- list()
  output <- list()
  cellularity <- 1.0
  acns <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  acns_segs <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  ploidies$ploidy <- as.numeric(ploidies$ploidy)

  for (i in unique(ploidies$sample_id)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    ploidy <- dplyr::filter(ploidies, sample_id == i)
    rcn_obj <- sample_segments
    absolute_segments <- dplyr::mutate(sample_segments,
                                       copy_number = RelativeToAbsoluteCopyNumber(copy_number,
                                                                                  ploidy$ploidy,
                                                                                  cellularity))
    # Optionally, grab needed values to create an S4 QDNAseq object
    if (return_S4) {
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- RelativeToAbsoluteCopyNumber(absolute_cns$copy_number,
                                                                           ploidy$ploidy,
                                                                           cellularity)
      absolute_cns$segmented <- RelativeToAbsoluteCopyNumber(absolute_cns$segmented,
                                                                         ploidy$ploidy,
                                                                         cellularity)
      acns <- cbind(acns, absolute_cns$copy_number)
      acns_segs <- cbind(acns_segs, absolute_cns$segmented)
    }
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome,
                                                             start=start,
                                                             end=end,
                                                             segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j <- j+1
  }

  # Return collapsed acn segments list or an S4 object if wanted
  if (return_S4) {
    colnames(acns) <- unique(rascal_batch_solutions$sample)
    colnames(acns_segs) <- unique(unique(rascal_batch_solutions$sample))
    rownames(acns) <- cns$feature[1:dim(acns)[1]]
    rownames(acns_segs) <- cns$feature[1:dim(acns)[1]]
    output[["acns"]] <- acns
    output[["acns_segs"]] <- acns_segs
  } else{
    output[["acn_segment_tables"]] <- chosenSegmentTablesList
  }

  # return the known ploidies
  if (return_sols) {
    output[["ploidies"]] <- ploidies
  }

  return(output)
}


# Add highest MAD column to rascal solutions table
# Then create and save list of segment tables
# Corollary to calculate_vaf_acns function.
# input_file: input path
#' @export
GenMadAcns <- function (segments,
                        rascal_batch_solutions,
                        return_sols = FALSE,
                        cns = FALSE,
                        return_S4 = FALSE,
                        distance_decimal_places) {

  chosenSegmentTablesList <- list()
  j <- 1
  plots <- list()
  output <- list()
  acns <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  acns_segs <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)

  sols <- data.frame()

  for (i in unique(rascal_batch_solutions$sample)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    # the distance is rounded just before filtering for min distance. This is to ensure distances that only differ 
    # due to floating-point error are treated as equal, and subsequently decided based on max cellularity.
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i) %>%
                      dplyr::mutate(distance = round(distance, distance_decimal_places)) %>%
                      dplyr::filter(distance == min(distance)) %>%
                      dplyr::filter(cellularity == max(cellularity))
    sols <- rbind(sols, solutions)
    absolute_segments <- dplyr::mutate(sample_segments,
                                       copy_number = RelativeToAbsoluteCopyNumber(copy_number,
                                                                                  solutions$ploidy,
                                                                                  solutions$cellularity))
    # Optionally, grab needed values to create an S4 QDNAseq object
    if (return_S4) {
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- RelativeToAbsoluteCopyNumber(absolute_cns$copy_number, 
                                                                           solutions$ploidy,
                                                                           solutions$cellularity)
      absolute_cns$segmented <- RelativeToAbsoluteCopyNumber(absolute_cns$segmented,
                                                                         solutions$ploidy,
                                                                         solutions$cellularity)
      acns <- cbind(acns, absolute_cns$copy_number)
      acns_segs <- cbind(acns_segs, absolute_cns$segmented)
    }
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome,
                                                             start=start,
                                                             end=end,
                                                             segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j <- j+1
  }

  # Return collapsed acn segments list or an S4 object if wanted
  if (return_S4) {
    colnames(acns) <- unique(rascal_batch_solutions$sample)
    colnames(acns_segs) <- unique(unique(rascal_batch_solutions$sample))
    rownames(acns) <- cns$feature[1:dim(acns)[1]]
    rownames(acns_segs) <- cns$feature[1:dim(acns)[1]]
    output[["acns"]] <- acns
    output[["acns_segs"]] <- acns_segs
  } else {
    output[["acn_segment_tables"]] <- chosenSegmentTablesList
  }

  # Return the selected ploidies and cellularities
  if (return_sols) {
    output[["rascal_solutions"]] <- sols
  }

  return(output)
}


# 1. Look up vaf gene start and end in genome reference
# 2. Filter rcn obj (ex. QDNAseq obj) for seg. CN at the vaf gene location
GetSegmentedRcnForGene <- function (rcn_obj, gene) {
  granges_gene <- GenomicFeatures::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, filter = ~ gene_name == gene)
  grRCN <- GenomicRanges::makeGRangesFromDataFrame(rcn_obj)
  hits <- IRanges::findOverlaps(grRCN, granges_gene)

  if (length(hits) == 0) {
    return(FALSE)
  } else {
    overlaps <- GenomicRanges::pintersect(grRCN[S4Vectors::queryHits(hits)],
                                          granges_gene[S4Vectors::subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(granges_gene[S4Vectors::subjectHits(hits)])
    idx <- S4Vectors::queryHits(hits)[which.max(percentOverlap)]
    return(rcn_obj$copy_number[idx])
  }
}

ReplaceQDNAseqAssaySlots <- function(cnobj, new_cns, new_segs) {

  # In case some samples need to be dropped, set a mask
  missing_samples <- cnobj@phenoData@data$name %in% colnames(new_cns)

  if (any(!missing_samples)) {
    warning("Not all samples had ACN solutions.
    These samples are excluded from the returned acn object.
    The probability of loss/gain values in the returned acn object were calculated using the relative CNs object.")
  }

  # Assemble a new QDNAseq object
  new.cnobj <- new('QDNAseqCopyNumbers',
                    bins = cnobj@featureData@data,
                    copynumber = new_cns,
                    phenodata = cnobj@phenoData@data[missing_samples,])
  Biobase::assayDataElement(new.cnobj, "segmented") <- new_segs

  if (length(assayData(cnobj)) > 2) {
    slotstoadd <- names(assayData(cnobj))
    slotstoadd <- slotstoadd[!slotstoadd %in% c("segmented", "copynumber")]
    for (i in slotstoadd)
    Biobase::assayDataElement(new.cnobj, i) <- cnobj@assayData[[i]][,missing_samples]
  }
  return(new.cnobj)
}


#' Compute distance function for fitting relative copy numbers to absolute
#'
#' Computes a distance function on fitting relative copy numbers to whole number
#' absolute values with the given ploidy and cellularity.
#' Copied from rascal.
#'
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' AbsoluteCopyNumberDistance(3, 0.67, copy_number$segmented)
#' @export
AbsoluteCopyNumberDistance <- function(ploidy, cellularity,
                                          relative_copy_numbers, weights = NULL,
                                          distance_function = c("MAD", "RMSD")) {
  
  distance <- 1e10
  if (cellularity < 0.0 || cellularity > 1.0) return(distance)
  
  absolute_copy_numbers <- RelativeToAbsoluteCopyNumber(relative_copy_numbers, ploidy, cellularity)
  
  absolute_copy_number_steps <- round(absolute_copy_numbers)
  
  differences <- abs(absolute_copy_numbers - absolute_copy_number_steps)
  # differences[which(absolute_copy_numbers > 0 & differences > 0.5)] <- 0.5
  
  if (is.null(weights)) weights <- rep(1, length(relative_copy_numbers))
  
  if (is.character(distance_function)) {
    if (distance_function[1] == "MAD") {
      distance <- stats::weighted.mean(differences, weights, na.rm = TRUE)
    } else if (distance_function[1] == "RMSD") {
      distance <- sqrt(stats::weighted.mean(differences ^ 2, weights, na.rm = TRUE))
    }
  }
  
  distance
}

#' Compute distance function for fitting relative copy numbers for a grid of
#' ploidies and cellularities
#'
#' Computes a distance function on fitting relative copy numbers to whole number
#' absolute values for a grid of ploidies and cellularities.
#' Copied from rascal.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param min_ploidy,max_ploidy the range of ploidies.
#' @param ploidy_step the stepwise increment of ploidies along the grid.
#' @param min_cellularity,max_cellularity the range of cellularities.
#' @param cellularity_step the stepwise increment of cellularities along the
#' grid.
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' segments <- CopyNumberSegments(copy_number)
#' distances <- AbsoluteCopyNumberDistanceGrid(segments$copy_number, segments$weight)
#' @export
AbsoluteCopyNumberDistanceGrid <- function(relative_copy_numbers, weights = NULL,
                                               min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
                                               min_cellularity = 0.01, max_cellularity = 1.0, cellularity_step = 0.01,
                                               distance_function = c("MAD", "RMSD")) {
  
  grid <- tidyr::expand_grid(
    ploidy = seq(min_ploidy, max_ploidy, ploidy_step),
    cellularity = seq(min_cellularity, max_cellularity, cellularity_step)
  )

  distances <- purrr::pmap_dbl(grid,
                        AbsoluteCopyNumberDistance,
                        relative_copy_numbers = relative_copy_numbers,
                        weights = weights,
                        distance_function = distance_function)
  
  dplyr::mutate(grid, distance = distances)
}

#' Determine if a given ploidy and cellularity for fitting copy numbers is
#' acceptable
#'
#' Determines if the given ploidy and cellularity are compatible with the
#' provided relative copy numbers based on the proportion of copy number values
#' that correspond to or are below the zero copy number state and on the
#' proportion that are sufficiently close to a whole number copy number state.
#' Copied from rascal.
#'
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param max_proportion_zero the maximum proportion of fitted absolute copy
#' number values in the zero copy number state.
#' @param min_proportion_close_to_whole_number the minimum proportion of fitted
#' absolute copy number values sufficiently close to a whole number.
#' @param max_distance_from_whole_number the maximum distance from a whole
#' number that a fitted absolute copy number can be to be considered
#' sufficiently close.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' IsAcceptablePloidyAndCellularity(3, 0.67, copy_number$segmented)
#' @export
IsAcceptablePloidyAndCellularity <- function(ploidy, cellularity,
                                                 relative_copy_numbers, weights = NULL,
                                                 max_proportion_zero = 0.05,
                                                 min_proportion_close_to_whole_number = 0.5,
                                                 max_distance_from_whole_number = 0.15) {
  
  absolute_copy_numbers <- RelativeToAbsoluteCopyNumber(relative_copy_numbers, ploidy, cellularity)
  absolute_copy_number_steps <- round(absolute_copy_numbers)
  
  if (is.null(weights)) weights <- rep(1, length(absolute_copy_numbers))
  sum_of_weights <- sum(weights)
  
  # filter based on proportion at or below copy number state 0
  prop_zero <- (sum(weights[which(absolute_copy_number_steps <= 0)]) / sum_of_weights)
  if (prop_zero > max_proportion_zero) return(FALSE)
  
  # filter based on proportion not sufficiently close to whole number copy number state
  prop_state <- (sum(weights[which(abs(absolute_copy_numbers - absolute_copy_number_steps) < max_distance_from_whole_number)]) / sum_of_weights)
  if (prop_state < min_proportion_close_to_whole_number) return(FALSE)

  return(TRUE)
}

#' Find best fit solutions for fitting copy numbers through grid-based search
#'
#' Find best fit solutions for a grid-based search through the given ranges of
#' ploidies and cellularities.
#' Copied from rascal.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param min_ploidy,max_ploidy the range of ploidies.
#' @param ploidy_step the stepwise increment of ploidies along the grid.
#' @param min_cellularity,max_cellularity the range of cellularities.
#' @param cellularity_step the stepwise increment of cellularities along the
#' grid.
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @param distance_filter_scale_factor the distance threshold above which
#' solutions will be discarded as a multiple of the solution with the smallest
#' distance.
#' @param max_proportion_zero the maximum proportion of fitted absolute copy
#' number values in the zero copy number state.
#' @param min_proportion_close_to_whole_number the minimum proportion of fitted
#' absolute copy number values sufficiently close to a whole number.
#' @param max_distance_from_whole_number the maximum distance from a whole
#' number that a fitted absolute copy number can be to be considered
#' sufficiently close.
#' @param solution_proximity_threshold how close two solutions can be before one
#' will be filtered; reduces the number of best fit solutions where there are
#' many minima in close proximity.
#' @param keep_all set to \code{TRUE} to return all solutions but with
#' additional \code{best_fit} column to indicate which are the local minima that
#' are acceptable solutions (may be useful to avoid computing the distance grid
#' twice)
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#'
#' segments <- CopyNumberSegments(copy_number)
#'
#' solutions <- FindBestFitSolutions(
#'   segments$copy_number, segments$weight,
#'   min_ploidy = 1.5, max_ploidy = 5.5,
#'   distance_function = "MAD")
#' @export
FindBestFitSolutions <- function(relative_copy_numbers, weights = NULL,
                                    min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
                                    min_cellularity = 0.2, max_cellularity = 1.0, cellularity_step = 0.01,
                                    distance_function = c("MAD", "RMSD"),
                                    distance_filter_scale_factor = 1.25,
                                    max_proportion_zero = 0.05,
                                    min_proportion_close_to_whole_number = 0.5,
                                    max_distance_from_whole_number = 0.15,
                                    solution_proximity_threshold = 5,
                                    keep_all = FALSE) {
  
  # compute distance function for grid-based range of ploidies and cellularities
  distances <- AbsoluteCopyNumberDistanceGrid(
    relative_copy_numbers, weights,
    min_ploidy, max_ploidy, ploidy_step,
    min_cellularity, max_cellularity, cellularity_step,
    distance_function)
  
  # add identifier for each solution to help with marking up the best fit
  # solutions if returning all solutions
  distances <- distances %>%
    dplyr::mutate(id = dplyr::row_number())
  
  # find minima, i.e. those grid points which have a smaller distance than all
  # adjacent points
  distances <- distances %>%
    dplyr::mutate(x = dplyr::dense_rank(cellularity)) %>%
    dplyr::mutate(y = dplyr::dense_rank(ploidy))
  
  solutions <- dplyr::select(distances, id, x, y, distance)
  
  for (xdelta in -1:1) {
    for (ydelta in -1:1) {
      if (xdelta != 0 || ydelta != 0) {
        solutions <- solutions %>%
          dplyr::mutate(xc = x + xdelta, yc = y + ydelta) %>%
          dplyr::left_join(dplyr::select(distances, xc = x, yc = y, dc = distance), by = c("xc", "yc")) %>%
          dplyr::filter(is.na(dc) | distance <= dc) %>%
          dplyr::select(id, x, y, distance)
      }
    }
  }
  
  distances <- dplyr::select(distances, -x, -y)
  
  solutions <- solutions %>%
    dplyr::select(id) %>%
    dplyr::left_join(distances, by = "id")
  
  # only retain solutions with distances no more than
  # distance_filter_scale_factor times the the minimum value
  if (is.numeric(distance_filter_scale_factor) && nrow(solutions) > 1) {
    solutions <- solutions %>%
      dplyr::filter(distance < distance_filter_scale_factor * min(distance))
  }
  
  # retain solutions with acceptable ploidies and cellularities
  solutions <- solutions %>%
    dplyr::rowwise() %>%
    dplyr::filter(
      IsAcceptablePloidyAndCellularity(
        ploidy, cellularity,
        relative_copy_numbers, weights,
        max_proportion_zero = max_proportion_zero,
        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
        max_distance_from_whole_number = max_distance_from_whole_number
      )
    ) %>%
    dplyr::ungroup()
  
  # remove redundant solutions that have very similar ploidy and cellularity to
  # another solution with a smaller distance
  if (is.numeric(solution_proximity_threshold) && nrow(solutions) > 1) {
    
    solution_proximity_threshold2 <- solution_proximity_threshold ^ 2
    
    solutions <- solutions %>%
      dplyr::arrange(distance) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(x = cellularity * 100) %>%
      dplyr::mutate(y = (ploidy - min(ploidy)) * 100 / (max(ploidy) - min(ploidy)))
    
    pairs <- tidyr::expand_grid(rank1 = 1:nrow(solutions), rank2 = 1:nrow(solutions)) %>%
      dplyr::filter(rank2 > rank1) %>%
      dplyr::left_join(dplyr::select(solutions, rank1 = rank, x1 = x, y1 = y, distance1 = distance), by = "rank1") %>%
      dplyr::left_join(dplyr::select(solutions, rank2 = rank, x2 = x, y2 = y, distance2 = distance), by = "rank2") %>%
      dplyr::mutate(xy2 = (x2 - x1) ^ 2 + (y2 - y1) ^ 2)
    
    to_remove <- dplyr::tibble(rank = integer(0))
    
    while (nrow(pairs) > 0) {
      
      rank <- min(pairs$rank1)
      
      proximal <- pairs %>%
        dplyr::filter(rank1 == rank, xy2 <= solution_proximity_threshold2) %>%
        dplyr::select(rank = rank2)
      
      to_remove <- dplyr::bind_rows(to_remove, proximal)
      
      pairs <- pairs %>%
        dplyr::filter(rank1 != rank) %>%
        dplyr::anti_join(proximal, by = c("rank1" = "rank"))
    }
    
    solutions <- solutions %>%
      dplyr::anti_join(to_remove, by = "rank") %>%
      dplyr::select(id, ploidy, cellularity, distance)
  }
  
  # return either the entire distance grid or just the selected solutions
  if (keep_all) {
    distances %>%
      dplyr::mutate(best_fit = id %in% solutions$id) %>%
      dplyr::select(ploidy, cellularity, distance, best_fit)
  } else {
    dplyr::select(solutions, ploidy, cellularity, distance)
  }
}

# Compute the scaling factor for converting relative copy numbers to absolute
# copy numbers for the given ploidy and cellularity.
# Copied from rascal.
AbsoluteCopyNumberScalingFactor <- function(ploidy, cellularity) {
  ploidy + (2 / cellularity) - 2
}

#' Convert relative copy numbers to absolute copy numbers
#'
#' Convert relative copy numbers to absolute copy numbers based on the given
#' ploidy and cellularity.
#' Copied from rascal.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are
#' from the tumour.
#' @return a numeric vector of absolute copy numbers.
#' @examples
#' RelativeToAbsoluteCopyNumber(c(0.98, 1.6, 1.23), 4.01, 0.77)
#' @export
RelativeToAbsoluteCopyNumber <- function(relative_copy_numbers, ploidy, cellularity) {
  
  stopifnot(is.numeric(relative_copy_numbers))
  stopifnot(is.numeric(ploidy))
  stopifnot(is.numeric(cellularity))
  
  ploidy + AbsoluteCopyNumberScalingFactor(ploidy, cellularity) * (relative_copy_numbers - 1)
}

#' Compute tumour DNA fraction for the given absolute copy number and
#' cellularity
#'
#' Compute the tumour DNA fraction for the given absolute copy number(s) in the
#' tumour and the cellularity, i.e. the fraction of all cells that are tumour
#' cells.
#' Copied from rascal.
#'
#' @param absolute_copy_number the absolute copy number(s) (a numeric value).
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @return the fraction of DNA that originates from tumour cells.
#' @examples
#' tumour_fraction(3, 0.82)
#' @export
TumourFraction <- function(absolute_copy_number, cellularity) {

  stopifnot(is.numeric(absolute_copy_number))
  stopifnot(is.numeric(cellularity))

  tumour <- absolute_copy_number * cellularity
  normal <- 2 * (1 - cellularity)
  tumour / (tumour + normal)
}

SumDelimitedElements <- function (element, delimiter) {
  elements <- str_split(element, pattern = delimiter)[[1]]
  element <- sum(as.numeric(elements))
  return(element)
}

MinmaxColumn <- function (df_col) {
  preproc2 <- preProcess(as.data.frame(df_col), method=c("range"))
  minmaxed <- predict(preproc2, as.data.frame(df_col))
  return(minmaxed[,1])
}