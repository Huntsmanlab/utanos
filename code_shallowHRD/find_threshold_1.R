library(GenomicRanges)
source("helpers_threshold_functions.R")

#' TODO: WRITE DOCS
FindThreshold <- function(bam_ratios_frame, segments, num_simulations=100000, second_round) {
  segments = segments[which(segments$chr != 23),]
  #### Dropping the Feature column ####
  bam_ratios_frame = bam_ratios_frame[,-1]
  colnames(bam_ratios_frame) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  #### Creating GRanges objects ####
  granges_object_bam_ratios_frame = makeGRangesFromDataFrame(bam_ratios_frame, 
                                                             keep.extra.columns=TRUE,
                                                             ignore.strand=TRUE,
                                                             seqinfo=NULL,
                                                             seqnames.field="chr",
                                                             start.field="start",
                                                             end.field="end")
  for (i in 1:nrow(segments)) {
    gr = GRanges(seqnames=c(segments[i,1]),
                 ranges=IRanges(start=c(segments[i,3]), 
                                end=c(segments[i,4])),
                 strand=c("*"))
    subsetGRobject = subsetByOverlaps(granges_object_bam_ratios_frame, gr)
    segments[i,5] = median(subsetGRobject$ratio)
  }
  
  #### Getting segments with size > 3Mb ####
  print("Getting > 3mb data...")
  segments_3mb = segments[which(segments[,6] > 2999999),]
  print(dim(segments_3mb))
  N_3mb = dim(segments_3mb)[1]
  segments_3mb = data.matrix(segments_3mb)
  
  #### Getting the differences between all pairs of ratio_median's ####
  all_ratio_differences = c()
  for (i in 1:N_3mb) {
    v = i + 1
    for (j in v:N_3mb-1) {
      all_ratio_differences = c(all_ratio_differences, abs(segments_3mb[i,5] - segments_3mb[j, 5]))
    }
  }
  
  #### Simulations ####
  print("Entering simulations...")
  thr <- RunThresholdSimulations(num_simulations=num_simulations,
                                 all_ratio_differences=all_ratio_differences,
                                 second_round=second_round)
  #### final thr is min(max(0.025, thr), 0.045)
  if (thr > 0.45) {
    thr = 0.45
  }
  
  if (thr < 0.025) {
    thr = 0.025
  }
  
  thr
}
