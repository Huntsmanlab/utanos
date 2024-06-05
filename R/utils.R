# sWGS Utils
# Functions used in the analysis of shallow WGS data

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
#
# GenHumanReadableAcnProfile
# GenHumanReadableRcnProfile
#

###

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
                           return_S4 = FALSE) {

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
    output <- GenKploidyAcns(segments, acnmethod, )

  } else if (length(acnmethod) > 1) {
    variants <- variants %>% dplyr::group_by(sample_id) %>%
      dplyr::filter(gene_name %in% acnmethod | vaf == max(vaf, na.rm = TRUE)) %>%
      dplyr::filter(gene_name == c(intersect(acnmethod, gene_name),
                                   setdiff(gene_name, acnmethod))[1])
    output <- GenVafAcns(segments, rascal_sols, variants,
                         return_sols, long_cns, return_S4)

  } else if (acnmethod == 'mad') {
      output <- GenMadAcns(segments, rascal_sols, return_sols,
                           long_cns, return_S4)

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
      dplyr::mutate(absolute_copy_number = rascal::relative_to_absolute_copy_number(vaf_gene_rcn, ploidy, cellularity)) %>%
      dplyr::mutate(tumour_fraction = rascal::tumour_fraction(absolute_copy_number, cellularity))

    sol_idx <- which.min(abs(as.double(vafs$vaf) - (solution_set$tumour_fraction*100)))
    solution_set <- solution_set[,]
    absolute_segments <- dplyr::mutate(sample_segments,
                                copy_number = rascal::relative_to_absolute_copy_number(copy_number,
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
      absolute_cns$copy_number <- rascal::relative_to_absolute_copy_number(absolute_cns$copy_number,
                                                                           solution_set[sol_idx,]$ploidy,
                                                                           solution_set[sol_idx,]$cellularity)
      absolute_cns$segmented <- rascal::relative_to_absolute_copy_number(absolute_cns$segmented,
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
                                       copy_number = rascal::relative_to_absolute_copy_number(copy_number,
                                                                                              ploidy$ploidy,
                                                                                              cellularity))
    # Optionally, grab needed values to create an S4 QDNAseq object
    if (return_S4) {
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- rascal::relative_to_absolute_copy_number(absolute_cns$copy_number,
                                                                           ploidy$ploidy,
                                                                           cellularity)
      absolute_cns$segmented <- rascal::relative_to_absolute_copy_number(absolute_cns$segmented,
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
                        return_S4 = FALSE) {

  chosenSegmentTablesList <- list()
  j <- 1
  plots <- list()
  output <- list()
  acns <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)
  acns_segs <- matrix(nrow = sum(cns$sample == cns$sample[1]), ncol = 0)

  for (i in unique(rascal_batch_solutions$sample)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i) %>%
                      dplyr::filter(distance == min(distance)) %>%
                      dplyr::filter(cellularity == max(cellularity))
    absolute_segments <- dplyr::mutate(sample_segments,
                                       copy_number = rascal::relative_to_absolute_copy_number(copy_number,
                                                                                              solutions$ploidy,
                                                                                              solutions$cellularity))
    # Optionally, grab needed values to create an S4 QDNAseq object
    if (return_S4) {
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- rascal::relative_to_absolute_copy_number(absolute_cns$copy_number,
                                                                           solutions$ploidy,
                                                                           solutions$cellularity)
      absolute_cns$segmented <- rascal::relative_to_absolute_copy_number(absolute_cns$segmented,
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
  } else{
    output[["acn_segment_tables"]] <- chosenSegmentTablesList
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

GeneCr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- stringr::str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    if (sum(toupper(queryset$Gene) %in% toupper(line)) > 0) {
      queryset_matches <- c(queryset_matches, i)
    }
  }
  return(targetset[unique(queryset_matches),])
}

AnnotationCr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- stringr::str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    queryset_matches <- c(queryset_matches, which(toupper(queryset$gene_name) %in% toupper(line)))
  }
  return(queryset[queryset_matches,])
}

#' Transforms segment tables into per-bin copy-number tables.
#'
#' @description
#'
#' SegmentsToCopyNumber() transforms segment tables into per-bin copy-number tables.
#' It is an expansion of the calls into per-bin style, where the bin size is user defined.
#' It is the inverse of the CopyNumberSegments function.
#'
#' @param segs A list of dataframes, where each dataframe is a dataframe of the segmented copy-numbers.
#' @param bin_size A natural number. The binsize used in the copy number object, ex. 30 -> 30kb, 100 -> 100kb.
#' @param genome A string. Refers to the reference genome, common reference genomes are: 'hg19', 'mm10', or 'hg38'.
#' @param Xincluded A boolean. Sex chromosomes as a Boolean (e.g. FALSE if not included).
#' @returns A dataframe of per-bin copy numbers (bounded for each sample).
#'
#' @export
SegmentsToCopyNumber <- function(segs, bin_size,
                                 genome = 'hg19', Xincluded = FALSE) {

  # Stop execution if we don't have the required input
  if (is.list(segs) && !is.data.frame(segs)) {
    stopifnot(length(segs) != 0)
    stopifnot("chromosome" %in% names(segs[[1]]))
    stopifnot("start" %in% names(segs[[1]]), is.numeric(segs[[1]]$start))
    stopifnot("end" %in% names(segs[[1]]), is.numeric(segs[[1]]$end))
    stopifnot("segVal" %in% names(segs[[1]]), is.numeric(segs[[1]]$segVal))
  } else if (is.data.frame(segs)) {
    stopifnot(dim(segs)[1] != 0)
    stopifnot("sample" %in% colnames(segs))
    stopifnot("chromosome" %in% colnames(segs))
    stopifnot("start" %in% colnames(segs), is.numeric(segs$start))
    stopifnot("end" %in% colnames(segs), is.numeric(segs$end))
    stopifnot("segVal" %in% colnames(segs), is.numeric(segs$segVal))
    segslist <- list()
    for (samplename in unique(segs$sample)) {
      tempvar <- segs %>% dplyr::filter(sample == samplename) %>%
                    dplyr::select(-sample)
      segslist[[samplename]] <- as.data.frame(tempvar)
    }
    segs <- segslist
  } else {
    stop("Please provide either a list of dataframes or a single dataframe of
         copy-numbers on which to perform this transformation.")
  }

  # Create template binned genome
  chroms <- GetBinnedChromosomes(genome, Xincluded, bin_size)
  genome_template <- data.frame(chromosome = rep(chroms$chrom, chroms$nbins),
                                start = rep(rep(1,dim(chroms)[1]), chroms$nbins),
                                end = rep(chroms$size, chroms$nbins)) %>%
    dplyr::group_by(chromosome,start,end) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
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
      dplyr::group_by(start,end,segVal) %>%
      dplyr::mutate(id = 1:dplyr::n()) %>%
      dplyr::mutate(chromosome = as.character(chromosome),
                    start = (floor(start/bin_size) * bin_size + 1) + ((id-1) * bin_size),
                    end = start + bin_size - 1) %>%
      dplyr::select(-c(size,nbins,id)) %>%
      dplyr::ungroup()
    per_bin_copy_numbers <- dplyr::left_join(genome_template, sample,
                                             by = c('chromosome', 'start', 'end'))
    per_bin_copy_numbers$sample_id <- name
    out <- rbind(out, per_bin_copy_numbers)
  }

  names(out) <- c("chromosome", "start", "end", "segmented", "sample_id")
  return(out)
}

GetBinnedChromosomes <- function(genome, Xincluded = FALSE, bin_size) {

  if (genome == 'hg38') {
    if (Xincluded) genome_chrs <- 23 else genome_chrs <- 22
    chroms <- GenomeInfoDb::getChromInfoFromUCSC("hg38", map.NCBI=TRUE) %>%
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  } else if (genome == 'mm10') {
    if (Xincluded) genome_chrs <- 20 else genome_chrs <- 19
    chroms <- GenomeInfoDb::getChromInfoFromUCSC("mm10") %>%
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  } else {
    if (Xincluded) genome_chrs <- 23 else genome_chrs <- 22
    chroms <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  }
  chroms$chrom <- sub('chr', '', chroms$chrom)
  return(chroms)
}

# Convert a string of genome ranges into GRanges object
#' @export
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- tidyr::separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- GenomicRanges::makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

#' @export
GetGRangesFromEnsDb <- function(
    ensdb,
    standard.chromosomes = TRUE,
    biotypes = c(ensembldb::listGenebiotypes(ensdb)),
    verbose = TRUE
) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = GenomeInfoDb::seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- GenomeInfoDb::keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
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
#' @export
CompareBinCNs <- function(objs, sample, bin_area) {
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
RemoveBlacklist <- function(data) {
  # Read in blacklist file
  blacklist = hg19.blacklistBins

  # Convert blacklist and data to GRanges objects and find indices of overlaps
  grBL = GenomicRanges::makeGRangesFromDataFrame(blacklist)
  grData = GenomicRanges::makeGRangesFromDataFrame(data)
  overlaps = suppressWarnings(IRanges::findOverlaps(grData, grBL))

  # Replace state of data rows with chr and start overlapping with blacklist
  data$state = replace(data$state, overlaps@from, NA)

  # Return original dataframe with state of blacklisted segments of genome identified as low mapability
  return(data)
}


#' Collapses relative copy-number calls to segment tables
#'
#' @description
#'
#' CopyNumberSegments() transforms relative copy-number calls to segment tables.
#' Inverse of the SegmentsToCopyNumber function.
#'
#' @param copy_number A dataframe. A dataframe with copy number calls (with columns 'sample', 'chromosome', 'start', 'end', 'segmented')
#' @returns A dataframe. A dataframe of summaries of various characteristics (derived from copy-number calls)
#'
#' @export
CopyNumberSegments <- function(copy_number) {

  stopifnot(is.data.frame(copy_number))
  stopifnot("sample" %in% names(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("segmented" %in% names(copy_number), is.numeric(copy_number$segmented))

  copy_number %>%
    dplyr::filter(!is.na(segmented)) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(sample, chromosome, start) %>%
    dplyr::mutate(new_segment = dplyr::row_number() == 1 |
                    !(sample == dplyr::lag(sample) &
                        chromosome == dplyr::lag(chromosome) &
                        segmented == dplyr::lag(segmented))) %>%
    dplyr::mutate(segment = cumsum(new_segment)) %>%
    dplyr::group_by(segment) %>%
    dplyr::summarise(
      sample = dplyr::first(sample),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      gain_probability = if ("gainP" %in% colnames(.)) dplyr::first(gainP) + dplyr::first(ampP) else NA,
      loss_probability = if ("lossP" %in% colnames(.)) dplyr::first(lossP) + dplyr::first(dlossP) else NA,
      copy_number = dplyr::first(segmented),
      bin_count = dplyr::n(),
      sum_of_bin_lengths = sum(length),
      weight = sum(length) / median(length)
    )
}


BestRelativeProfile <- function (selection, robjects) {

  if (length(robjects) < 2) {
    stop("Please pass a minimum of two QDNAseq objects to this function.")
  }
  feature_lengths <- sapply(robjects, nrow)
  idx <- which.max(feature_lengths)
  n <- dim(selection)[1]
  dim_names <- rownames(robjects[[idx]]@assayData[["copynumber"]])
  tbins <- as.data.frame(stringr::str_split(dim_names, pattern = ':|-', simplify = TRUE))
  colnames(tbins) <- c('chromosome', 'start', 'end')
  tbins <- tbins %>% dplyr::filter(chromosome != 'Y')
  dim_names <- paste0(tbins$chromosome, ':', tbins$start, '-', tbins$end)
  cns <- tbins
  segs <- tbins
  objstats <- c()

  for (i in names(robjects)) {
    # Pull what needed out of each object in turn to fill the new slots for:
    # copy-number
    sample_ids_mask <- selection$sample_id[i == selection$selection]
    tempmat <- robjects[[i]]@assayData[["copynumber"]][,sample_ids_mask]
    tempmat <- cbind(as.data.frame(tempmat),
                     as.data.frame(stringr::str_split(rownames(tempmat),
                                                      pattern = ':|-',
                                                      simplify = TRUE)))
    cns <- cns %>% dplyr::left_join(tempmat, by = c('chromosome' = 'V1',
                                                    'start' = 'V2',
                                                    'end' = 'V3'))
    # segments
    tempmat <- robjects[[i]]@assayData[["segmented"]][,sample_ids_mask]
    tempmat <- cbind(as.data.frame(tempmat),
                     as.data.frame(stringr::str_split(rownames(tempmat),
                                                      pattern = ':|-',
                                                      simplify = TRUE)))
    segs <- segs %>% dplyr::left_join(tempmat, by = c('chromosome' = 'V1',
                                                      'start' = 'V2',
                                                      'end' = 'V3'))
    # stats
    objstats <- rbind(objstats,
                      robjects[[i]]@phenoData@data[sample_ids_mask,
                                                   c('name',
                                                     'total.reads',
                                                     'expected.variance')])
  }

  # Create copy-numbers array
  cns <- cns[,selection$sample_id]
  copynumbers <- matrix(NA_real_, nrow=nrow(cns), ncol=n,
                        dimnames=list(dim_names, selection$sample_id))
  copynumbers[,1:n] <- as.matrix(cns[,1:(dim(cns)[2])])
  # Create segments array
  segs <- segs[,selection$sample_id]
  segments <- matrix(NA_real_, nrow=nrow(segs), ncol=n,
                     dimnames=list(dim_names, selection$sample_id))
  segments[,1:n] <- as.matrix(segs[,1:dim(segs)[2]])
  # Create bins annotated dataframe
  bins <- matrix(NA_integer_, nrow=nrow(cns), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(tbins[,1:3])
  bins[,4] <- (complete.cases(segments) & complete.cases(copynumbers))
  bins <- Biobase::AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  bins@data$start <- as.integer(bins@data$start)
  bins@data$end <- as.integer(bins@data$end)
  bins@data$use <- as.logical(bins@data$use)
  # Create stats dataframe
  stats <- matrix(NA_integer_, nrow=nrow(selection), ncol=3,
                  dimnames=list(selection$sample_id, c('name', 'total.reads', 'expected.variance')))
  temp_stats <- selection[,1] %>% dplyr::left_join(objstats,
                                                   by = c('sample_id' = 'name')) %>%
    dplyr::rename(name = sample_id)
  stats[,1:3] <- as.matrix(temp_stats[,1:3])
  # Assemble QDNAseq object
  newobj <- new('QDNAseqCopyNumbers',
                bins = bins,
                copynumber = copynumbers,
                phenodata = as.data.frame(stats))
  Biobase::assayDataElement(newobj, "segmented") <- segments
  newobj@featureData@data[["chromosome"]] <- as.character(newobj@featureData@data[["chromosome"]])
  # Add last bit of metadata
  metadata <-  matrix(NA_integer_, nrow = 3, ncol = 1,
                      dimnames=list(c('name', 'total.reads', 'expected.variance')))
  colnames(metadata) <- "labelDescription"
  newobj@phenoData@varMetadata <- as.data.frame(metadata)
  return(newobj)
}

# Input:
# Change from a list of Component-by-Signature Matrices from NMFfit objects to the NMFfit objects themselves.
# This change will permit for easier plotting.
# This means that NMFfit objects must be saved in the cn_signatures_pilot.R script.

#' @export
CompareSignatureByComponent <- function (cs_matrix_list, mat_names, saveFile) {

  nsig <- 7
  # britroc feat_sig matrix
  feat_sig_mat<-basis(cs_matrix_list[[1]])
  colnames(cs_matrix_list[[1]]@fit@W) <- c('s7', 's1', 's6', 's4', 's3', 's2', 's5')
  reord_1 <- as.integer(c(2,6,5,4,7,3,1))
  names(reord_1)<-paste0("s",1:7)
  feat_sig_mat<-feat_sig_mat[,reord_1]
  colnames(feat_sig_mat)<-paste0("s",1:nsig)
  sig_feat_mat<-t(feat_sig_mat)

  # Second feat_sig matrix
  feat_sig_mat_2<-basis(cs_matrix_list[[2]])
  colnames(cs_matrix_list[[2]]@fit@W) <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')
  sig_feat_mat_2<-t(feat_sig_mat_2)
  colnames(feat_sig_mat_2)<-paste0("s",1:nsig)

  # Third feat_sig matrix
  feat_sig_mat_3<-basis(cs_matrix_list[[3]])
  colnames(cs_matrix_list[[3]]@fit@W) <- c('s1', 's2', 's3', 's4', 's5', 's6', 's7')
  sig_feat_mat_3<-t(feat_sig_mat_3)
  colnames(feat_sig_mat_3)<-paste0("s",1:nsig)

  # browser()
  # determine matching signatures and their correlation
  # reord_2<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_2,2,cor,x,method="s"))})
  # sigcor_2<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_2,2,cor,x,method="s"))})
  reord_2 <- as.integer(c(2,3,4,5,6,7,1))

  # reord_3<-apply(feat_sig_mat,2,function(x){which.max(apply(feat_sig_mat_3,2,cor,x,method="s"))})
  # sigcor_3<-apply(feat_sig_mat,2,function(x){max(apply(feat_sig_mat_3,2,cor,x,method="s"))})
  reord_3 <- as.integer(c(1,5,2,7,3,6,4))

  # plot the feat_sig matrices side by side
  pdf(file = saveFile, 10, 7, onefile = FALSE)
  par(mfrow=c(1,3))
  basismap(cs_matrix_list[[1]],Rowv=NA,Colv=reord_1,main=mat_names[1],tracks=NA)
  basismap(cs_matrix_list[[2]],Rowv=NA,Colv=reord_2,main=mat_names[2],tracks=NA)
  basismap(cs_matrix_list[[3]],Rowv=NA,Colv=reord_3,main=mat_names[3],tracks=NA)
  dev.off()
}

#' @export
CreateShallowHRDInput <- function (input_obj, temp_path, output_path) {

  dir.create(file.path(output_path), showWarnings = FALSE, recursive = TRUE)
  exportBins(input_obj, file = paste0(temp_path, "/logtransCNs.txt"), type = 'copynumber', format="tsv")
  exportBins(input_obj, file = paste0(temp_path, "/logtransSegs.txt"), type = 'segments', format="tsv")
  CNs <- data.table::fread(file = paste0(temp_path, "/logtransCNs.txt"), sep = "\t", header = TRUE)
  segs <- data.table::fread(file = paste0(temp_path, "/logtransSegs.txt"), sep = "\t", header = TRUE)

  CNs <- gather(CNs, sample_id, ratio, `CC-CHM-1341`:`YW-EC052`, factor_key=TRUE)
  segs <- gather(segs, sample_id, ratio_median, `CC-CHM-1341`:`YW-EC052`, factor_key=TRUE)

  CNs = cbind(CNs, segs$ratio_median)
  CNs = CNs[,-4]
  CNs = CNs[,-1]

  samples <- unique(CNs$sample_id)
  for (i in samples) {
    X <- CNs %>% dplyr::filter(sample_id == i, !is.na(ratio), !is.na(`segs$ratio_median`))
    X <- X[,-3]
    colnames(X) <- c("chromosome", "start", "ratio", "ratio_median")
    X[,1] <- gsub('X','23',X[,1])
    X[,1] = as.numeric(as.character(X[,1]))

    write.table(X, file = paste0(output_path,"/", i,".bam_ratio.txt"), sep = "\t", row.names = FALSE)
  }

  file.remove(paste0(temp_path, "/logtransCNs.txt"))
  file.remove(paste0(temp_path, "/logtransSegs.txt"))
}

GetChromosomeLengths <- function(build) {
  build <- as.integer(gsub('[^0-9]', '', build))
  if (build == 34 || build == 16) {
    chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
  } else if (build == 35 || build == 17) {
    chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
  } else if (build == 36 || build == 18) {
    chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if (build == 37 || build == 19)  {
    chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  } else { # hg38
    chromosome.lengths <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  }
  names(chromosome.lengths) <- as.character(c(1:22, 'X', 'Y'))
  chromosome.lengths
}


#'  Split Chromosome Name
#'
#' @description
#' Function to split compounded chromosome-position names into separate columns.
#' Distinct columns for chromosome number, start, and ending positions are the result.
#'
#' @param x *dataframe* containing the copy number segments in long format with columns: \cr
#' chr sample segVal
#'
#' @details
#' A QDNAseq object generally contains the chromosome and position information as such: \cr
#' chromosome:start-end.
#' This is inconvenient for plotting and other data wrangling.
#'
#' @return dataframe containing three separate columns for the chromosome number, sample and segVal
#'
ChromosomeSplitPos <- function(x) {
  x$position <- x$chr
  x$chr <- stringr::str_split_fixed(x$chr, ":", n = 2)
  x$chromosome <- x$chr[, 1]
  x$pos <- stringr::str_split_fixed(x$chr[, 2], "-", n = 2)
  x$start <- x$pos[, 1]
  x$end <- x$pos[, 2]
  x$chr <- NULL
  x$pos <- NULL
  x[, c("value", "start", "end")] <- lapply(x[, c("value", "start", "end")], as.numeric)
  return(x)
}


#'  Export the CN bins of a QDNAseq Object to a wide dataframe
#'
#' @description
#' Function to export the bins of a QDNAseq object to a wide dataframe.
#' Can be any of the "copynumber", "segments", or "calls" slots.
#' Unlike the similar function 'exportBins' in the QDNAseq package it does not log-normalize.
#'
#' @param object *QDNAseq Object* containing the copy number segments in long format with columns: \cr
#' chr sample segVal
#' @param type *character* Type of data to export, options are "copynumber" (corrected or uncorrected read counts), "segments", or "calls".
#' @param filter *logical* If @TRUE, bins are filtered, otherwise not.
#' @param digits *numeric* The number of digits to round to.
#' @param chromosomeReplacements *named char vector* A named character vector of chromosome name replacements to be done.
#' Only used when `object` is of class: \cr
#' `cghRaw`, `cghSeg`, `cghCall`, or `cghRegions`, \cr
#' since these classes store chromosome names as integers, whereas all QDNAseq object types use character vectors.
#' Defaults to `c("23"="X", "24"="Y", "25"="MT")` for human.
#'
#' @return A wide dataframe containing a column for each sample as well as 4 additional columns: \cr
#' `feature, chromosome, start, end`
#'
#' @export
ExportBinsQDNAObj <- function(object,
                              type=c("copynumber", "segments", "calls"),
                              filter=TRUE, digits=3,
                              chromosomeReplacements=c("23"="X", "24"="Y", "25"="MT"), ...) {

  type <- match.arg(type)

  if (inherits(object, "QDNAseqSignals")) {
    if (filter) {
      object <- object[binsToUse(object), ]
    }
    feature <- featureNames(object)
    chromosome <- fData(object)$chromosome
    start <- fData(object)$start
    end <- fData(object)$end
    if (inherits(object, "QDNAseqReadCounts")) {
      if (type != "copynumber")
        warning("Ignoring argument 'type' and returning read counts.")
      dat <- assayDataElement(object, "counts")
      type <- "read counts"
    } else {
      if (type == "copynumber") {
        dat <- assayDataElement(object, "copynumber")
      } else if (type == "segments") {
        if (!"segmented" %in% assayDataElementNames(object))
          stop("Segments not found, please run segmentBins() first.")
        dat <- assayDataElement(object, "segmented")
      } else if (type == "calls") {
        if (!"calls" %in% assayDataElementNames(object))
          stop("Calls not found, please run callBins() first.")
        dat <- assayDataElement(object, "calls")
      }
    }
  } else if (inherits(object, c("cghRaw", "cghSeg", "cghCall",
                                "cghRegions"))) {

    feature <- featureNames(object)
    chromosome <- as.character(chromosomes(object))
    for (chromosomeReplacement in names(chromosomeReplacements)) {
      chromosome[chromosome == chromosomeReplacement] <-
        chromosomeReplacements[chromosomeReplacement]
    }
    start <- bpstart(object)
    end <- bpend(object)
    if (inherits(object, c("cghRaw", "cghSeg", "cghCall"))) {
      if (type == "copynumber") {
        dat <- copynumber(object)
      } else if (type == "segments") {
        if (!"segmented" %in% assayDataElementNames(object))
          stop("Segments not found, please run segmentData() first.")
        dat <- segmented(object)
      } else if (type == "calls") {
        if (!"calls" %in% assayDataElementNames(object))
          stop("Calls not found, please run CGHcall() first.")
        dat <- calls(object)
      }
    } else if (inherits(object, "cghRegions")) {
      if (type != "calls")
        warning("Ignoring argument 'type' and returning calls.")
      dat <- regions(object)
    }
  }

  if (is.numeric(digits)) {
    dat <- round(dat, digits=digits)
  }
  oopts2 <- options(scipen=15)
  on.exit(options(oopts2))

  out <- data.frame(feature=feature, chromosome=chromosome, start=start,
                    end=end, dat, check.names=FALSE, stringsAsFactors=FALSE)

  return(out)
}

