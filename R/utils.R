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
#' CalculateACNs() calculates the absolute copy numbers (ACNs) from the relative copy number profiles of one or more samples.
#' There are several included options by which to do this.
#' Note: If not providing a table of variant allele frequencies (VAFs) then 'mad' is the only method available. \cr
#'
#' This function makes use of the rascal package in R and instructions can be found in the vignette: \cr
#' https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd  \cr
#'
#' @param relative_segs A tsv of the relative segmented copy-numbers.
#' @param rascal_sols A tsv of the calculated rascal solutions.
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
#' @param variants A dataframe of the variants including variant allele frequencies per gene and per sample.
#' @param relative_cns A tsv of the relative copy-numbers. Required if 'addplots' param is set to true.
#' @param addplots (optional) Logical. Indicates whether or not plots should be returned alongside the ACNs.
#' @param acn_save_path (optional) String. The output path (absolute path recommended) where to save the result.
#' @param return_sols (optional) Logical. Return the selected rascal solution.
#' @returns A list of dataframes. One DF for each sample where an Absolute Copy Number profile was successfully found. \cr
#' Optionally a second list of plot objects (ggplot); one for each sample where a profile was found.
#' @details
#' ```
#' solutions <- "~/Documents/.../rascal_solutions.csv"
#' rcn_segs <- "~/Documents/.../rCN_segs.tsv"
#' variants <- "~/Documents/.../allvariants.clinvar.cosmic.exons.csv"
#' save_path <- "~/Documents/.../rascal_ACN_segments.rds"
#' variants <- data.table::fread(file = variants, sep = ',')
#' variants <- variants %>% dplyr::rename(chromosome = chr,
#'                                        gene_name = genecode_gene_name)
#' variants$sample_id <- stringr::str_replace_all(variants$sample_id, "-", ".")
#' output <- CalculateACNs(relative_segs = rcn_segs,
#'                         rascal_sols = solutions,
#'                         variants = variants,
#'                         acnmethod = 'maxvaf')
#' output <- CalculateACNs(relative_segs = rcn_segs,
#'                         rascal_sols = solutions,
#'                         variants = variants,
#'                         acnmethod = c('TP53', 'KRAS', 'BRCA1',
#'                                       'BRCA2', 'PIK3CA', 'PTEN'),
#'                         acn_save_path = save_path)
#' ```
#'
#' @export
CalculateACNs <- function (relative_segs, rascal_sols, acnmethod,
                           variants = NULL, relative_cns = FALSE,
                           addplots = FALSE, acn_save_path = FALSE,
                           return_sols = FALSE) {

  relative_segs <- read.table(file = relative_segs, header = TRUE, sep = '\t')
  relative_segs <- relative_segs %>% tidyr::gather(sample, segmented,
                                                   5:dim(.)[2], factor_key=TRUE)    # Convert segs to long
  segments <- CopyNumberSegments(relative_segs)                                     # Collapse segs to continuous segments
  rascal_sols <- read.table(file = rascal_sols,
                            sep = ',', header = TRUE)
  rascal_sols$sample <- stringr::str_replace_all(rascal_sols$sample, "-", ".")

  if (addplots == TRUE) {
    relative_cns <- read.table(file = relative_cns, header = TRUE, sep = '\t')
    relative_cns <- relative_cns %>% tidyr::gather(sample,
                                                   copy_number,
                                                   5:dim(.)[2], factor_key=TRUE)    # Convert cns to long
    relative_cns <- relative_cns %>%
                        dplyr::mutate(segmented = relative_segs$segmented)          # Create combined dataframe
  }

  if (suppressWarnings({!is.null(variants)})) {
    if('data.frame' %in% class(variants)) {                                     # do nothing
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
    variants$sample_id <- stringr::str_replace_all(variants$sample_id, "-", ".")
  }

  if ((length(acnmethod) == 1) && (acnmethod == 'maxvaf')) {
    variants <- variants %>% dplyr::group_by(sample_id) %>%
                                  dplyr::filter(vaf == max(vaf))
    variants <- variants[!duplicated(variants$sample_id),]
    if (addplots == TRUE) {
      output <- GenVafAcns(segments, rascal_sols, variants, return_sols, cns = relative_cns)
    } else {
      output <- GenVafAcns(segments, rascal_sols, variants, return_sols)
    }
  } else if (length(acnmethod) > 1) {
    variants <- variants %>% dplyr::group_by(sample_id) %>%
      dplyr::filter(gene_name %in% acnmethod | vaf == max(vaf, na.rm = TRUE)) %>%
      dplyr::filter(gene_name == c(intersect(acnmethod, gene_name),
                                   setdiff(gene_name, acnmethod))[1])
    if (addplots == TRUE) {
      output <- GenVafAcns(segments, rascal_sols, variants, return_sols, cns = relative_cns)
    } else {
      output <- GenVafAcns(segments, rascal_sols, variants, return_sols)
    }
  } else if (acnmethod == 'mad') {

  } else {
    stop("Invalid value passed to the 'acnmethod' parameter. \n
         Please supply an accepted option.")
  }

  if (acn_save_path != FALSE) {
    saveRDS(output[[1]], file = acn_save_path)
  }
  return(output)
}


# Function that calculates absolute copy numbers using the rascal package and vafs.
GenVafAcns <- function (segments, rascal_batch_solutions, variants,
                        return_sols = FALSE, cns = FALSE) {

  chosenSegmentTablesList <- list()
  j <- 1
  plots <- list()
  sols <- data.frame(sample_id = character(),
                     ploidy = double(),
                     cellularity = double(),
                     vaf_optimal = logical(),
                     stringsAsFactors = FALSE)

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
                  solution_set[sol_idx,]$ploidy,
                  solution_set[sol_idx,]$cellularity,
                  TRUE))
    if (!is.logical(cns)) {                                                     # If cns AND segments passed in, make plots!
      chr_ord <- c(as.character(1:22), 'X')
      absolute_cns <- dplyr::filter(cns, sample == i)
      absolute_cns$copy_number <- rascal::relative_to_absolute_copy_number(absolute_cns$copy_number,
                                                                           solution_set[sol_idx,]$ploidy,
                                                                           solution_set[sol_idx,]$cellularity)
      absolute_cns$segmented <- rascal::relative_to_absolute_copy_number(absolute_cns$segmented,
                                                                           solution_set[sol_idx,]$ploidy,
                                                                           solution_set[sol_idx,]$cellularity)
      absolute_cns$chromosome <- factor(absolute_cns$chromosome, levels = chr_ord)
      absolute_segments$chromosome <- factor(absolute_segments$chromosome, levels = chr_ord)
      plots[[i]] <- rascal::genome_copy_number_plot(absolute_cns, absolute_segments,
                              min_copy_number = 0, max_copy_number = 12,
                              copy_number_breaks = 0:12,
                              ylabel = "Absolute Copy Number",
                              xlabel = "Chromosome") +
                ggplot2::theme(panel.grid.major.y =
                                 ggplot2::element_line(colour = "grey60"),
                               plot.title = ggplot2::element_text(size = 20)) +
                ggplot2::ggtitle(paste0(i, ' - est. ploidy:',
                                        solution_set[sol_idx,]$ploidy,
                                        ' - est. cellularity:',
                                        solution_set[sol_idx,]$cellularity))
    }
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome,
                                                             start=start,
                                                             end=end,
                                                             segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j <- j+1
  }
  if (!is.logical(cns)) {
    output <- list(segment_tables = chosenSegmentTablesList, plots = plots)
  } else{
    output <- list(segment_tables = chosenSegmentTablesList)
  }
  if (return_sols) {
    # output <- c(output, rascal_solutions = sols)
    colnames(sols) <- c('sample_id', 'ploidy', 'cellularity', 'vaf_optimal')
    output[['rascal_solutions']] <- as.data.frame(sols)
  }
  return(output)
}


# Add highest MAD column to rascal solutions table
# Then create and save list of segment tables
# Corollary to calculate_vaf_acns function.
# input_file: input path
#' @export
#' GenMadAcns <- function (segments, rascal_batch_solutions, cns = FALSE) {
GetOptimalMadSolutions <- function (input_file, rcn_file, segs_file, save_file = TRUE) {

  chosenSegmentTablesList <- list()
  j <- 1
  plots <- list()

  for (i in unique(rascal_batch_solutions$sample)) {
    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i)

    solution_set <- solutions %>%
      dplyr::select(ploidy, cellularity) %>%
      dplyr::mutate(absolute_copy_number = rascal::relative_to_absolute_copy_number(ploidy, cellularity)) %>%
      dplyr::mutate(tumour_fraction = rascal::tumour_fraction(absolute_copy_number, cellularity))
  }






  rascal_batch_solutions <- read.table(file = input_file, sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")
  temp <- rascal_batch_solutions %>%
    dplyr::group_by(sample) %>%
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
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
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
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
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
  stopifnot(is.list(segs))
  stopifnot(length(segs) != 0)
  stopifnot("chromosome" %in% names(segs[[1]]))
  stopifnot("start" %in% names(segs[[1]]), is.numeric(segs[[1]]$start))
  stopifnot("end" %in% names(segs[[1]]), is.numeric(segs[[1]]$end))
  stopifnot("segVal" %in% names(segs[[1]]), is.numeric(segs[[1]]$segVal))

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
  names(out) <- c("chromosome", "start", "end", "state", "sample_id")
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
  granges <- GenomicRanges::makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

#' @export
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
  overlaps = IRanges::findOverlaps(grData, grBL)

  # Replace state of data rows with chr and start overlapping with blacklist
  data$state = replace(data$state, overlaps@from, NA)

  # Return original dataframe with state of blacklisted segments of genome identified as low mapability
  return(data)
}


#' Generate Human Readable Relative Copy-Number Profiles
#'
#' @description
#'
#' This function takes as input a QDNAseq or CGHcall copy-number object and gives back a long-format table with several useful columns.
#' These columns include; 'sample', 'chromosome', 'start', 'end', 'gain_probability',
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

  segmented <- tidyr::gather(as.data.frame(CGHbase::segmented(object)), sample, segmented)
  if (class(object)[1] == 'cghCall') {
    segmented$gainP <- tidyr::gather(as.data.frame(CGHbase::probgain(object)), sample, probgain)$probgain
    segmented$ampP <- tidyr::gather(as.data.frame(CGHbase::probamp(object)), sample, probamp)$probamp
    segmented$lossP <- tidyr::gather(as.data.frame(CGHbase::probloss(object)), sample, probloss)$probloss
    segmented$dlossP <- tidyr::gather(as.data.frame(CGHbase::probdloss(object)), sample, probdloss)$probdloss
    segmented$chromosome <- rep(object@featureData@data[["Chromosome"]], dim(object)[2])
    segmented$start <- rep(object@featureData@data[["Start"]], dim(object)[2])
    segmented$end <- rep(object@featureData@data[["End"]], dim(object)[2])
  } else {
    segmented$chromosome <- rep(object@featureData@data[["chromosome"]], dim(object)[2])
    segmented$start <- rep(object@featureData@data[["start"]], dim(object)[2])
    segmented$end <- rep(object@featureData@data[["end"]], dim(object)[2])
    if ('probgain' %in% names(object@assayData)) {
      segmented$gainP <- tidyr::gather(as.data.frame(CGHbase::probgain(object)), sample, probgain)$probgain
      segmented$ampP <- tidyr::gather(as.data.frame(CGHbase::probamp(object)), sample, probamp)$probamp
      segmented$lossP <- tidyr::gather(as.data.frame(CGHbase::probloss(object)), sample, probloss)$probloss
      segmented$dlossP <- tidyr::gather(as.data.frame(CGHbase::probdloss(object)), sample, probdloss)$probdloss
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
  colnames(collapsed_segs) <- c('sample', 'chromosome', 'start', 'end', 'gain_probability',
                                'loss_probability', 'relative_copy_number', 'bin_count',
                                'sum_of_bin_lengths', 'cytobands', 'coordinates', 'size')
  # save segment tables
  if ( save_dir != FALSE ) {
    dir.create(file.path(save_dir, binsize))
    for (i in unique(collapsed_segs$sample)) {
      temp <- collapsed_segs[collapsed_segs$sample == i, ]
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
  collapsed_segs <- plyr::ldply(object, data.frame, .id = 'sample')
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
  colnames(collapsed_segs) <- c('sample', 'chromosome', 'start', 'end',
                                'absolute_copy_number', 'cytobands',
                                'coordinates', 'size')
  collapsed_segs <- collapsed_segs %>% mutate(chromosome = replace(chromosome,
                                                                   chromosome=="X", '23')) %>%
    dplyr::group_by(sample) %>%
    dplyr::arrange(as.integer(chromosome), .by_group = TRUE) %>%
    mutate(chromosome = replace(as.character(chromosome),
                                chromosome=='23', 'X'))
  # save segment tables
  dir.create(save_path)
  for (i in unique(collapsed_segs$sample)) {
    temp <- collapsed_segs[collapsed_segs$sample == i, ]
    write.table(temp, file = paste0(save_path, i, '_segsCytobandsTable.tsv'),
                sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  dbDisconnect(connection)
  return(collapsed_segs)
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

.getChromosomeLengths <- function(build) {
  build <- as.integer(gsub('[^0-9]', '', build))
  if (build == 34 || build == 16) {
    chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
  } else if (build == 35 || build == 17) {
    chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
  } else if (build == 36 || build == 18) {
    chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else {
    chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  }
  names(chromosome.lengths) <- 1:24
  chromosome.lengths
}

