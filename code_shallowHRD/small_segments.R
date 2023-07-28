library("GenomicRanges")
source("helpers_small_segments.R")

#' Returns the segments with 3Mb >= size >= 0.1Mb
#' 
#' @description
#' Filters the segments and returns only the small ones
#' 
#' @param segments A data frame containing the segment data. 
#' 
#' @export

GetSmallSegments <- function(segments) {
  result = segments[which(segments[,7] > 99999), ]
  result = result[which(result[,7] < 2999999), ]
  
  result
}

#' Deals with missing chromosome arm values in the given data frame of large segments
#' 
#' @description
#' TODO
#'
#' @param
#'
#' @export

LargeMissingChrArms <- function(large_segments, small_segments) {
  values = unique(small_segments[,3][!small_segments[,3] %in% large_segments[,3]])
  length_values = length(values)
  
  for (missing_chr_arm  in values) {
    i = 1
    N_large = dim(large_segments)[1]
    
    while (i < N_large+1) {
      if (missing_chr_arm > max(large_segments[,3])) {
        large_segments = rbind(large_segments,
                               small_segments[which(small_segments[,3] == missing_chr_arm),])
      } else if (missing_chr_arm < lage_segments[i,3]) {
        if (i == 1) {
          large_segments = rbind(small_segments[which(small_segments[,3] == missing_chr_arm),],
                                 large_segments[1:N_large,])
          i = N_large + 1
        } else {
          large_segments = rbind(large_segments[1:(i-1),],
                                 small_segments[which(small_segments[,3] == missing_chr_arm),],
                                 large_segments[i:N_large])
        }
        i = N_large + 1
      } else {
        i = i + 1
      }
    }
  }
  large_segments
}


#' Deals with missing chromosome arm values in the given data frame of small segments
#' 
#' @description
#' TODO
#'
#' @param
#'
#' @export

SmallMissingChrArms <- function(large_segments, small_segments) {
  values = unique(small_segments[,3][!small_segments[,3] %in% large_segments[,3]])
  length_values = length(values)
  
  for (missing_chr_arm in values) {
    small_segments <- small_segments[which(!(small_segments[,3] == missing_chr_arm)),]
  }
  
  small_segments
}

#' Adds small segments in `small_segments` that meet certain criteria (position wise w.r.t large segments)
#' into the `large_segments` data frame. Merges if necessary (i.e. if ratio differecne < THR).
#'
#' @description
#' TODO
#'
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param bam_ratios_frame A data frame: the cleaned-up version of the raw bam ratios file.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE.
#' 
#' @export

InsertSmallSegments <- function(large_segments, small_segments, bam_ratios_frame, threshold) {
  bam_ratios_frame = bam_ratios_frame[,-1]
  colnames(bam_ratios_frame) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  GRanges_object_ratio_file = makeGRangesFromDataFrame(bam_ratios_frame,
                                                       keep.extra.columns = TRUE,
                                                       ignore.strand = TRUE, 
                                                       seqinfo = NULL,
                                                       seqnames.field = "chr",
                                                       start.field = "start",
                                                       end.field = "end")
  
  output <- InitalizeSmallSegments(large_segments = large_segments,
                                   small_segments = small_segments,
                                   threshold = threshold,
                                   granges_obj = GRanges_object_ratio_file)
  large_segments <- output[1]
  small_segments <- output[2]
  
  large_segments <- FinalizeSmallSegments(large_segments = large_segments,
                                          small_segments = small_segments,
                                          threshold = threshold,
                                          granges_obj = GRanges_object_ratio_file)
  
  large_segments
}




