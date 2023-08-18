library("GenomicRanges")

#' Determines if absolute ratio_median difference is less than or equal to the given threshold.
#' 
#' @description 
#' Returns TRUE if large_segment and small_segment's ratio_median difference is <= threshold. FALSE otherwise.
#' 
#' @param large_segment An array: a single segment, doesn't necessarily have to be a 'large' segment. 
#' This name was chosen because most of the times we're comparing a large segment in `large_segments` with 
#' a small segment in `small_segments`
#' @param small_segment Same as above.
#' @param threshold A float: the threshold for ratio_median difference previously estimated.
#'
#' @export
RatioDiffBelowThreshold <- function(large_segment, small_segment, threshold) {
  output = (abs(small_segment[6] - large_segment[6])) <= threshold
  output
}


#' Helper for merging segments in the second case of FinalizeSmallSegments. Returns the updated `large_segments`
#' data frame with segments appropriately merged (or not... depending on parameters).
#' 
#' @description 
#' Given the `large_segments` data, a large_segment `j`, and a small_segment `small_segment`, we merge the large segment with
#' the small_segment according to various criteria. First, if their ratio_median difference is below the threshold, then
#' we typically merge large and small. If not, then small segment is separately added into the `large_segments` frame. Also, where
#' exactly in the `large_segments` frame we add this new data depends on `end_of_file`. For example, if it is TRUE, then there are no
#' other large segments after `j`, so we row-bind all of large_segments to the new data. If FALSE, then we 'squeeze' the new data between
#' j-1 and j+1. 
#' 
#' @param granges_obj A GRanges object: is used as reference to obtain true genomic data from specified regions.
#' @param small_segment An array: the given small segment that might be merged into the `large_segments` frame.
#' @param large_segments A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
#' small segments into it, smaller segments are progressively added.
#' @param j An integer/index: the position of the large segment we're comparing with small_segment, within `large_segments`.
#' @param below_thr A boolean: whether large_segment and small_segment ratio_median difference is below the threhsold.
#' @param end_of_file A boolean: whether there's any large segments after `j`.
#' 
#' @export
MergeSegmentsTwo <- function(granges_obj, small_segment, large_segments, j, below_thr, end_of_file) {
  N_large = dim(large_segments[1])
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)
  
  if (below_thr == TRUE) {
    gr = GRanges(seqnames = c(large_segment[2]),
                 ranges = IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(large_segment[1], 
                               large_segment[2], 
                               large_segment[3],
                               small_segment[4],
                               large_segment[5], 
                               median(subsetGRobject$ratio),
                               large_segment[5] - small_segment[4] + 1, 
                               large_segment[8]))
    } else {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(large_segment[1], 
                               large_segment[2], 
                               large_segment[3],
                               small_segment[4], 
                               large_segment[5], 
                               median(subsetGRobject$ratio),
                               large_segment[5] - small_segment[4] + 1, 
                               large_segment[8]),
                             large_segments[(j+1):N_large,])
    }
  } else {
    gr = GRanges(seqnames = c(small_segment[2]),
                 ranges = IRanges(start=c(small_segment[4]), end=c(large_segment[4])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               small_segment[4], 
                               large_segment[4], 
                               median(subsetGRobject$ratio),
                               large_segment[4] - small_segment[4] + 1, 
                               small_segment[8]),
                             large_segments[j:N_large,])
    } else {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               small_segment[4], 
                               large_segment[4], 
                               median(subsetGRobject$ratio),
                               large_segment[4] - small_segment[4] + 1, 
                               small_segment[8]))
    }
  }
  large_segments
}

#' Handles merging step for the fourth case of FinalizeSmallSegments.
#' 
#' @description
#' Pretty much the same as in MergeSegments
#'
#' @param granges_obj A GRanges object: is used as reference to obtain true genomic data from specified regions.
#' @param small_segment An array: the given small segment that might be merged into the `large_segments` frame.
#' @param large_segments A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
#' small segments into it, smaller segments are progressively added.
#' @param j An integer/index: the position of the large segment we're comparing with small_segment, within `large_segments`.
#' @param below_thr A boolean: whether large_segment and small_segment ratio_median difference is below the threhsold.
#' @param end_of_file A boolean: whether there's any large segments after `j`.
#' 
#' @export
MergeSegmentsFour <- function(granges_obj, small_segment, large_segments, j, below_thr, end_of_file) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)
  
  if (below_thr == TRUE) {
    gr = GRanges(seqnames = c(large_segment[2]),
                 ranges = IRanges(start=c(large_segment[4]), end=c(small_segment[5])),
                 strand = c("*"))
    granges_obj = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(large_segment[1], 
                               large_segment[2], 
                               large_segment[3],
                               large_segment[4],
                               small_segment[5], 
                               median(granges_obj$ratio),
                               small_segment[5] - large_segment[4] + 1, 
                               large_segment[8]))
    } else {
      large_segments = rbind(large_segments[1:(j-1),],
                             c(large_segment[1], 
                               large_segment[2], 
                               large_segment[3],
                               large_segment[4],
                               small_segment[5], 
                               median(granges_obj$ratio),
                               small_segment[5] - large_segment[4] + 1, 
                               large_segment[8]),
                             large_segments[(j+1):N_large,])
    }
  } else {
    gr = GRanges(seqnames = c(small_segment[2]),
                 ranges = IRanges(start=c(large_segment[5]), end=c(small_segment[5])),
                 strand = c("*"))
    granges_obj = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:j,],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               large_segment[5], 
                               small_segment[5], 
                               median(granges_obj$ratio),
                               small_segment[5] - large_segment[5] + 1, 
                               small_segment[8]))
    } else {
      large_segments = rbind(large_segments[1:j,],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               large_segment[5], 
                               small_segment[5], 
                               median(granges_obj$ratio),
                               small_segment[5] - large_segment[5] + 1, 
                               small_segment[8]),
                             large_segments[(j+1):N_large,])
    }
  }
  large_segments
}

MergeSegmentsFiveLargeSmall <- function(granges_obj, small_segment, large_segments, j) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)
  
  gr = GRanges(seqnames = c(large_segment[2]),
               ranges = IRanges(start=c(large_segment[4]), end=c(small_segment[5])),
               strand = c("*"))
  granges_obj = subsetByOverlaps(granges_obj, gr)
  
  large_segments = rbind(large_segments[1:(j-1),],
                         c(large_segment[1], 
                           large_segment[2], 
                           large_segment[3],
                           large_segment[4], 
                           small_segment[5], 
                           median(granges_obj$ratio),
                           small_segment[5] - large_segment[4] + 1, 
                           large_segment[8]),
                         large_segments[(j+1):N_large,])
  
  large_segments
}
MergeSegmentsFiveSmallNextLarge <- function(granges_obj, small_segment, large_segments, j, end_of_file) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)
  
  gr = GRanges(seqnames = c(large_segment[2]),
               ranges = IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
               strand = c("*"))
  granges_obj = subsetByOverlaps(granges_obj, gr)
  
  if (end_of_file == TRUE) {
    large_segments = rbind(large_segments[1:j-1,],
                           c(large_segment[1], 
                             large_segment[2], 
                             large_segment[3],
                             small_segment[4], 
                             large_segment[5], 
                             median(granges_obj$ratio),
                             large_segment[5] - small_segment[4] + 1, 
                             large_segment[8]))
  } else {
    large_segments = rbind(large_segments[1:j-1,],
                           c(large_segment[1], 
                             large_segment[2], 
                             large_segment[3],
                             small_segment[4], 
                             large_segment[5], 
                             median(granges_obj$ratio),
                             large_segment[5] - small_segment[4] + 1, 
                             large_segment[8]),
                           large_segments[(j+1):N_large,])
  }
  large_segments
}
