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

#' Returns a GRanges object used to obtain data from specified genomic regions.
#' 
#' @description
#' As above. Usually used to get data from a resulting merge, or from missing segments
#' in between segments we already have data for.
#'
#' @param bam_ratios_frame A data frame: the cleaned-up version of the raw bam ratios file.
#' 
#' @export

GetGRangesObject <- function(bam_ratios_frame) {
  bam_ratios_frame = bam_ratios_frame[,-1]
  colnames(bam_ratios_frame) <- c("chr", "start", "end", "ratio", "ratio_median")
  
  granges_obj = GenomicRanges::makeGRangesFromDataFrame(bam_ratios_frame, keep.extra.columns = TRUE, ignore.strand = TRUE, 
                                                        seqinfo = NULL, seqnames.field = "chr", start.field = "start",
                                                        end.field = "end")
  granges_obj
}
  

#' Deals with missing chromosome arm values in the given data frame of large segments.
#' 
#' @description
#' Finds the missing chromosome arms and reinserts them into `large_segments`. Specifically,
#' the arms that `large_segments` doesn't have, but `small_segments` does. For each of the missing
#' values, we iterate through the large segments, and re-insert the small segments wherever
#' appropriate. 
#' 
#' For example (first case in while loop), if the missing_chr_arm value is greater than the largest chr_arm
#' in the arm column (i.e., max(large_segments[,3])), then we should bind the small segments 
#' to the right of the large segments frame, since there are no more segments after
#' this largest chr_arm. Note that we only bind the small segments within this missing_chr_arm.
#' 
#' Otherwise, missing_chr_arm <= than the largest chr_arm. In this case we care about which large
#' segment we're looping over (the value of `i`). If i = 1, i.e. the first segment, then we bind
#' the small segments to the left of the large segments: this is trivial. Otherwise, we fit the
#' small segments in between `i-1` and `i`.
#'
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
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
      } else if (missing_chr_arm < large_segments[i,3]) {
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


#' Deals with missing chromosome arm values in the given data frame of small segments.
#' 
#' @description
#' Pretty much as in LargeMissingChrArms, except it removes the small segments
#' that are not in the missing_chr_arm. 
#'
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
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
#' into the `large_segments` data frame. Merges if necessary (i.e. if ratio_median difference < threshold).
#'
#' @description
#' Does a lot of things. In InitializeSmallSegments, we insert all the small segments that come before the first
#' large segment in `large_segments` (hence initializing). In FinalizeSmallSegments, we insert all the remaining small segments,
#' according to 6 different cases, which are listed in the FinalizeSmallSegments documentation. In the end, returns the updated
#' `large_segments` with small segments appropriately re-inserted.
#'
#' @param small_segments A data frame. Segments with 3Mb > size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= 3Mb.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE.
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' an overlap of segments and get the ratio_median of this overlap.
#' 
#' @export

InsertSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  output <- InitalizeSmallSegments(large_segments = large_segments,
                                   small_segments = small_segments,
                                   threshold = threshold,
                                   granges_obj = granges_obj)
  
  large_segments <- FinalizeSmallSegments(large_segments = as.data.frame(output[1]),
                                          small_segments = as.data.frame(output[2]),
                                          threshold = threshold,
                                          granges_obj = granges_obj)
  large_segments
}
