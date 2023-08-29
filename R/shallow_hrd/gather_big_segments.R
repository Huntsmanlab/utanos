#' Starting from the largest segments, merges segments if ratio_median 
#' between the largest segments and another segment is less than the threshold 
#'
#' @description 
#' 1. First, e find the largest segment (by size) in the copy_segments data frame and its corresponding index.
#' 2. Then, we find the segment that whose ratio is closest to this largest segment's ratio.
#' 3. 
#' 
#' @param segments Data frame with segment data. 
#' @param copy_segments A copy of segments.
#' @param threshold A float. The THR found previously.
#'
#' @export
GatherBigSegments <- function(segments, copy_segments, threshold) {
  while (dim(copy_segments)[1] > 1) {
    largest_segment = which.max(copy_seg)
  }
}