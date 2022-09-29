# sWGS Utils
# Functions used in the analysis of shallow WGS data
# combined functions from:
# get_acns_from_vafs_rascal.R
# filter_segmented_CN_calls.R
#

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
# genHumanReadableACNprofile
# genHumanReadableRCNprofile
#

#####

# Calculate Absolute Copy Numbers (ACNs) given Variant Allele Frequencies (Vafs) for several genes.
# Use the rascal package in R to do this transformation.
# Based on instructions in the vignette:
# https://github.com/crukci-bioinformatics/rascal/blob/master/vignettes/rascal.Rmd
# Parameters:
# relative_cns -
# vafs -

#' Calculate max. VAF-based Absolute Copy Numbers
#'
#' @param rascal_batch_solutions the calculated rascal solutions
#' @param relative_cns the relative copy-numbers
#' @param acn_seg_path the output path where to save the result
#' @param vafs a tsv file of the the variant allele frequencies per gene and per sample
#'
#' @return None
#'
#' @examples
#' calculate_vaf_acns()
#'
#' @export
CalculateVafAcns <- function (rascal_batch_solutions, relative_cns, acn_seg_path, vafs) {

  qdnaseq_segs <- relative_cns
  qdnaseq_segs <- read.table(file = qdnaseq_segs, header = TRUE)
  qdnaseq_segments <- gather(qdnaseq_segs, sample, segmented, `CC.CHM.1341`:`YW.EC052`, factor_key=TRUE)    # Convert to long
  segments <- copy_number_segments(qdnaseq_segments)                            # Collapse to continuous segments

  rascal_batch_solutions <- read.table(file = rascal_batch_solutions, sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")

  vafs_cels <- vafs
  vafs_cels <- data.table::fread(file = vafs_cels, header = TRUE, sep = "\t", fill = TRUE)
  vafs_cels$sample_id <- str_replace_all(vafs_cels$sample_id, "-", ".")
  # To simplify matters lets sum all unique vafs for a given gene/sample and treat them as a single vaf
  vafs_cels[,8:13] <- as.data.frame(apply(vafs_cels[,8:13],
                                          MARGIN = c(1,2),
                                          function(x) sum_delimited_elements(x, ';')))

  chosenSegmentTablesList <- list()
  j <- 1
  for (i in unique(rascal_batch_solutions$sample)) {

    sample_segments <- dplyr::filter(segments, sample == i)
    solutions <- rascal_batch_solutions %>% dplyr::filter(sample == i)
    vafs <- vafs_cels %>% dplyr::filter(sample_id == i)
    if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN

    meta <- vafs %>% dplyr::select(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity)
    vafs <- vafs %>% dplyr::select(-c(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity))
    vaf_idx <- which.max(as.double(vafs[1,]))                                   # Choose max vaf
    if (length(vaf_idx) == 0) next                                              # If there isn't a vaf for the sample skip finding an ACN
    vaf <- as.double(vafs[1, ..vaf_idx])                                        # vaf
    vaf_name <- names(vafs[1, ..vaf_idx])
    vaf_name <- toupper(str_split(vaf_name, pattern = '\\.')[[1]][1])           # vaf name

    rcn_obj <- sample_segments
    vaf_gene_rcn <- get_segmented_rcn_for_gene(rcn_obj, vaf_name)
    if (vaf_gene_rcn == FALSE) { message(paste('No segmented CN call found for', rcn_obj$sample[1], 'VAF gene.', sep = ' ')); next}           # If there isn't a CN for the VAF skip finding an ACN
    solution_set <- solutions %>%                                               # Get from running find best fit solution
      dplyr::select(ploidy, cellularity) %>%
      dplyr::mutate(tp53_absolute_copy_number = relative_to_absolute_copy_number(vaf_gene_rcn, ploidy, cellularity)) %>%
      dplyr::mutate(tp53_tumour_fraction = tumour_fraction(tp53_absolute_copy_number, cellularity))

    sol_idx <- which.min(abs(as.double(vaf) - (solution_set$tp53_tumour_fraction*100)))
    solution_set <- solution_set[,]
    absolute_segments <- mutate(sample_segments,
                                copy_number = relative_to_absolute_copy_number(copy_number,
                                                                               solution_set[sol_idx,]$ploidy,
                                                                               solution_set[sol_idx,]$cellularity))
    absolute_segments <- absolute_segments %>% dplyr::select(chromosome=chromosome, start=start, end=end, segVal=copy_number)
    chosenSegmentTablesList[[i]] <- as.data.frame(absolute_segments)
    j = j+1
  }

  # names(chosenSegmentTablesList) <- rascal_batch_solutions$sample
  saveRDS(chosenSegmentTablesList, file = acn_seg_path)
}

# Add highest MAD column to rascal solutions table
# Then create and save list of segment tables
# Corollary to calculate_vaf_acns function.
# input_file: input path
#' @export
GetOptimalMadSolutions <- function (input_file, rcn_file, segs_file, save_file = TRUE) {

  rascal_batch_solutions <- read.table(file = input_file, sep = ',', header = TRUE)
  rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")
  temp <- rascal_batch_solutions %>%
    group_by(sample) %>%
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

  granges_gene <- genes(EnsDb.Hsapiens.v75, filter = ~ gene_name == gene)
  grRCN = makeGRangesFromDataFrame(rcn_obj)
  hits = findOverlaps(grRCN, granges_gene)

  if (length(hits) == 0) {
    return(FALSE)
  } else {
    overlaps <- pintersect(grRCN[queryHits(hits)], granges_gene[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(granges_gene[subjectHits(hits)])
    idx <- queryHits(hits)[which.max(percentOverlap)]
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

# A true 'utils' function
# Transforms segment tables into per-bin copy-number tables
# It is an expansion of the calls into per-bin style, where the bin size is user defined.
#' @export
SegmentsToCopyNumber <- function(segs, bin_size, genome = 'hg19', Xincluded = FALSE) {

  # Stop execution if we don't have the required input
  stopifnot(is.list(segs))
  stopifnot(length(segs) != 0)
  stopifnot("chromosome" %in% names(segs[[1]]))
  stopifnot("start" %in% names(segs[[1]]), is.numeric(segs[[1]]$start))
  stopifnot("end" %in% names(segs[[1]]), is.numeric(segs[[1]]$end))
  stopifnot("segVal" %in% names(segs[[1]]), is.numeric(segs[[1]]$segVal))
  # Create template binned genome
  genome_chrs <- 22
  if (Xincluded) { genome_chrs <- 23}
  if (genome == 'hg38') {
    chroms <- getChromInfoFromUCSC("hg38", map.NCBI=TRUE) %>%
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  } else {
    chroms <- getChromInfoFromUCSC("hg19") %>%
      head(genome_chrs) %>%
      dplyr::mutate(nbins = ceiling(size/bin_size))
  }
  chroms$chrom <- sub('chr', '', chroms$chrom)
  genome_template <- data.frame(chromosome = rep(chroms$chrom, chroms$nbins),
                                start = rep(rep(1,dim(chroms)[1]), chroms$nbins),
                                end = rep(chroms$size, chroms$nbins)) %>%
    dplyr::group_by(chromosome,start,end) %>%
    dplyr::mutate(id = row_number()) %>%
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
      group_by(start,end,segVal) %>%
      dplyr::mutate(id = row_number()) %>%
      dplyr::mutate(chromosome = as.character(chromosome),
                    start = (floor(start/bin_size) * bin_size + 1) + ((id-1) * bin_size),
                    end = start + bin_size - 1) %>%
      dplyr::select(-c(size,nbins,id)) %>%
      dplyr::ungroup()
    per_bin_copy_numbers <- left_join(genome_template, sample, by = c('chromosome', 'start', 'end'))
    per_bin_copy_numbers$sample_id <- name
    out <- rbind(out, per_bin_copy_numbers)
  }
  names(out) <- c("chromosome", "start", "end", "state", "sample_id")
  return(out)
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
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
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
  blacklist = as.data.frame(data.table::fread(file = "~/repos/cnsignatures/data/external_datasets/binBlacklist.txt", sep = ' ', header = TRUE))

  # Convert blacklist and data to GRanges objects and find indices of overlaps
  grBL = makeGRangesFromDataFrame(blacklist)
  grData = makeGRangesFromDataFrame(data)
  overlaps = findOverlaps(grData, grBL)

  # Replace state of data rows with chr and start overlapping with blacklist
  data$state = replace(data$state, overlaps@from, NA)

  # Return original dataframe with state of blacklisted segments of genome identified as low mapability
  return(data)
}

### Create Cytoband .tsv tables from relative copy-number calls
# DESCRIPTION
# Parameters:
# data -
#' @export
GenHumanReadableRcnProfile <- function(object, binsize) {
  # Expects gl cghcall object

  # Create collapsed segments table
  # Add Chromosome cytoband, coordinates (in bp), length of region, and gain or loss tag to each entry
  connection <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", dbname="hg19")
  cytobands <- dbGetQuery(conn=connection, statement="SELECT chrom, chromStart, chromEnd, name FROM cytoBand")
  cyto_ranges <- makeGRangesFromDataFrame(cytobands)

  segmented <- gather(as.data.frame(segmented(object)), sample, segmented)
  segmented$gainP <- gather(as.data.frame(probgain(object)), sample, probgain)$probgain
  segmented$ampP <- gather(as.data.frame(probamp(object)), sample, probamp)$probamp
  segmented$lossP <- gather(as.data.frame(probloss(object)), sample, probloss)$probloss
  segmented$dlossP <- gather(as.data.frame(probdloss(object)), sample, probdloss)$probdloss
  segmented$chromosome <- rep(object@featureData@data[["Chromosome"]], dim(object)[2])
  segmented$start <- rep(object@featureData@data[["Start"]], dim(object)[2])
  segmented$end <- rep(object@featureData@data[["End"]], dim(object)[2])

  # collapse copy numbers down to segments
  collapsed_segs <- copy_number_segments_new(segmented)
  collapsed_segs$weight <- NULL

  collapsed_segs <- collapsed_segs %>% transform(chromosome = as.character(chromosome))
  collapsed_segs$chromosome[collapsed_segs$chromosome == 23] <- 'X'
  ranges <- makeGRangesFromDataFrame(collapsed_segs)
  seqlevelsStyle(ranges) <-'UCSC'
  hits <- findOverlaps(ranges, cyto_ranges)
  temp <- data.frame(ranges = hits@from, cyto_ranges = hits@to, cytobands = cytobands$name[hits@to])
  temp <- temp[,c('ranges','cytobands')]

  # collapse cytobands to a single cell
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
  collapsed_segs$cytobands <- temp$cytobands
  collapsed_segs$coordinates <- paste0(seqnames(ranges), ':', ranges(ranges))
  collapsed_segs <- collapsed_segs %>% mutate(size = end - start)
  collapsed_segs$segment <- NULL
  colnames(collapsed_segs) <- c('sample', 'chromosome', 'start', 'end', 'gain_probability',
                                'loss_probability', 'relative_copy_number', 'bin_count',
                                'sum_of_bin_lengths', 'cytobands', 'coordinates', 'size')
  # save segment tables
  # collapsed_segs <- collapsed_segs %>% dplyr::group_by(samples)
  save_dir <- '~/Documents/projects/cn_signatures_shallowWGS/data/relativeCN_segs_and_cytoband_tables/'
  dir.create(file.path(save_dir, binsize))
  for (i in unique(collapsed_segs$sample)) {
    temp <- collapsed_segs[collapsed_segs$sample == i, ]
    file_name <- paste0(save_dir, '/', binsize, '/', i, '_', binsize, '_RsegsCytobandsTable.tsv')
    write.table(temp, file = file_name, sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  dbDisconnect(connection)
  return(collapsed_segs)
}

### Create Cytoband .tsv tables from absolute copy-number calls
# DESCRIPTION
# Parameters:
# data -
#' @export
GenHumanReadableAcnProfile <- function(object, save_path) {

  # Get Chromosome cytobands, coordinates (in bp), and lengths of regions
  connection <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", dbname="hg19")
  cytobands <- dbGetQuery(conn=connection, statement="SELECT chrom, chromStart, chromEnd, name FROM cytoBand")
  cyto_ranges <- makeGRangesFromDataFrame(cytobands)

  # Add cytobands
  collapsed_segs <- plyr::ldply(object, data.frame, .id = 'sample')
  collapsed_segs <- collapsed_segs %>% transform(chromosome = as.character(chromosome))
  collapsed_segs$chromosome[collapsed_segs$chromosome == 23] <- 'X'
  ranges <- makeGRangesFromDataFrame(collapsed_segs)
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
  collapsed_segs <- collapsed_segs %>% mutate(size = end - start)
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

### Collapse relative copy-number calls to segment tables
# DESCRIPTION
# Parameters:
# data -
#' @export
CopyNumberSegmentsNew <- function(copy_number) {

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
    dplyr::mutate(new_segment = row_number() == 1 | !(sample == lag(sample) & chromosome == lag(chromosome) & segmented == lag(segmented))) %>%
    dplyr::mutate(segment = cumsum(new_segment)) %>%
    dplyr::group_by(segment) %>%
    dplyr::summarize(
      sample = dplyr::first(sample),
      chromosome = dplyr::first(chromosome),
      start = dplyr::first(start),
      end = dplyr::last(end),
      gain_probability = dplyr::first(gainP) + dplyr::first(ampP),
      loss_probability = dplyr::first(lossP) + dplyr::first(dlossP),
      copy_number = dplyr::first(segmented),
      bin_count = n(),
      sum_of_bin_lengths = sum(length),
      weight = sum(length) / median(length)
    )
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

#'  Relative copy number plotting functions 
#'
#' @description Function to split chromosome names into the separate columns with the chromosome number, start and ending position
#' @details A QDNAseq object generally contains the chromosome and position information in the following format chromosome:start-end. This often needs to be separated into three columns for plotting relative copy numbers and other data wrangling.
#' @param x *dataframe* containing the copy number segments in long format: chr column with format chromosome:start-stop, sample and segVal
#' @return dataframe containing three separate columns for the chromosome number, sample and segVal
#' 
ChromosomeSplit <- function(x) {
  x$chr <- str_split_fixed(x$chr, ":", n = 2)
  x$chromosome <- x$chr[,1] 
  x$pos <- str_split_fixed(x$chr[,2], "-", n = 2)
  x$start <- x$pos[,1]
  x$end <- x$pos[,2]
  x$pos <- NULL
  x$chr <- NULL
  x[, c("segVal", "start", "end")] <- lapply(x[, c("segVal", "start", "end")], as.numeric)
  return(x)
}

plotWisecondorProfiles <- function (input_path) {
  
  # Read in files
  bin_files <- list.files(input_path, pattern = '*_bins.bed', full.names = TRUE)                            # grab just 1
  samples <- list.files(input_path, pattern = '*_bins.bed')
  samples <- sub('_bins.bed', '', samples)
  seg_files <- list.files(input_path, pattern = '*_segments.bed', full.names = TRUE)
  stats_files <- list.files(input_path, pattern = '*_statistics.txt', full.names = TRUE)
  cns <- plyr::ldply(seq_along(bin_files), function(x) {
    df <- data.table::fread(bin_files[[x]], sep = '\t', col.names = c('chromosome', 'start', 'end', 'id', 'ratio', 'zscore'), nThread = 4)
    df <- df[,c(1,2,3,5)]
    df$sample_id <- samples[x]
    colnames(df) = c('chromosome', 'start', 'end', 'state', 'sample_id')
    df})
  # NaNs to NAs
  cns$state[is.nan(cns$state)] <- NA
  # Exponeniate 
  cns$state <- exp(cns$state)
  segs <- lapply(seg_files, function(x) {
    df <- data.table::fread(x, sep = '\t', col.names = c('chromosome', 'start', 'end', 'segVal', 'zscore'))
    df <- df[,c(1,2,3,4)]
    df})
  names(segs) <- samples
  segs <- segs %>% 
    purrr::map(~mutate_at(.x, vars("segVal"), function(x) {exp(x)}))
  stats <- plyr::ldply(seq_along(stats_files), function(x) {
    read_count <- as.integer(sub("Number of reads: ", "", readLines(stats_files[[x]])[26]))
    df <- data.frame(name = samples[x], 
                     total_reads = read_count, 
                     stringsAsFactors=FALSE)
    df})
  rownames(stats) <- samples
  # Create the copynumber object
  cns_wide <- tidyr::spread(cns, sample_id, state)
  dim_names <- paste0(cns_wide$chromosome, ':', cns_wide$start, '-', cns_wide$end)       # create bin dimnames
  copynumbers <- matrix(NA_real_, nrow=nrow(cns_wide), ncol=length(bin_files),
                        dimnames=list(dim_names, samples))
  copynumbers[,1:length(bin_files)] <- as.matrix(cns_wide[,4:(length(bin_files)+3)])
  # Create the segment object
  segs_long <- SegmentsToCopyNumber(segs, bin_size = 30000, genome = 'hg19', Xincluded = TRUE)
  segs_wide <- spread(segs_long, sample_id, state)
  segments <- matrix(NA_real_, nrow=nrow(segs_wide), ncol=length(seg_files),
                     dimnames=list(dim_names, samples))
  segments[,1:length(seg_files)] <- as.matrix(segs_wide[,4:(length(seg_files)+3)])
  # Create the bins
  bins <- matrix(NA_integer_, nrow=nrow(cns_wide), ncol=4,
                 dimnames=list(dim_names, c('chromosome', 'start', 'end', 'use')))
  bins[,1:3] <- as.matrix(cns_wide[,1:3])
  bins[,4] <- complete.cases(segments)
  bins <- Biobase::AnnotatedDataFrame(as.data.frame(bins))
  bins@data$chromosome <- factor(bins@data$chromosome, levels = c(as.character(c(1:22)),'X'))
  bins@data$start <- as.integer(bins@data$start)
  bins@data$end <- as.integer(bins@data$end)
  bins@data$use <- as.logical(bins@data$use)
  
  # Assemble QDNAseq object
  wx_qdnaobj <- new('QDNAseqCopyNumbers', bins=bins, copynumber=copynumbers, phenodata=stats)
  assayDataElement(wx_qdnaobj, "segmented") <- segments
  # Change to total.reads
  colnames(wx_qdnaobj@phenoData@data) <- c("name", "total.reads")
  # Calculate expected variance 
  expected.variance <- rep(NA_real_, 262)
  for (i in seq_len(262)) {
    expected.variance[i] <- (sum(QDNAseq:::binsToUse(wx_qdnaobj[,i])) / wx_qdnaobj[,i]$total.reads)
  }
  wx_qdnaobj@phenoData@data$expected.variance <- expected.variance
  metadata <-  matrix(NA_integer_, nrow= 3, ncol= 1,
                      dimnames=list(c('name', 'total.reads', 'expected.variance')))
  colnames(metadata) <- "labelDescription"
  wx_qdnaobj@phenoData@varMetadata <- as.data.frame(metadata)
  return(wx_qdnaobj)
}

