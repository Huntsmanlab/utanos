# This script is used to retrieve CN profiles from pipeline/computational objects/tables.
# And then convert these profiles into useable/helpful tables for wetlab researchers.

library(RMySQL)
library(GenomicRanges)

copy_number_segments_new <- function(copy_number) {
  
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

genHumanReadableCNprofile <- function(object, binsize) {
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
  save_dir <- '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/data/relativeCN_segs_and_cytoband_tables/'
  dir.create(file.path(save_dir, binsize))
  for (i in unique(collapsed_segs$sample)) {
    temp <- collapsed_segs[collapsed_segs$sample == i, ]
    file_name <- paste0(save_dir, '/', binsize, '/', i, '_', binsize, '_RsegsCytobandsTable.tsv')
    write.table(temp, file = file_name, sep = '\t', col.names = TRUE, row.names = FALSE)
  }
  dbDisconnect(connection)
  return(collapsed_segs)
}

# obj1 <- genHumanReadableCNprofile(cgh15kb, '15kb')
ccnvFILT_obj <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-12/Xchr_included/15kb_rCN_comCNVfilt.rds')
ccnvFILT_obj <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-12/Xchr_included/15kb_gl_rCN.rds')
obj1 <- genHumanReadableCNprofile(ccnvFILT_obj, '15kb')

# Scratch
test <- data.frame(chr = cgh15kb@featureData@data[["Chromosome"]], 
                   CC.VGH.0033a = cgh15kb@assayData[["segmented"]][,1], 
                   start = cgh15kb@featureData@data[["Start"]], 
                   end = cgh15kb@featureData@data[["End"]])



