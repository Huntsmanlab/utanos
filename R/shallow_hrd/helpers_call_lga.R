library("GenomicRanges")
source("helpers_small_segments.R") # borrowing the helper that checks if segments ratio_median diff < THR

#' Re-assigns level column values for all segments.
#' 
#' @description 
#' Resets level to 0 for all segments, except for the first segment which is assigned 1.
#' Then, iterates through the rest:
#' For a segment 'i', if it is in the same chromosome arm as previous segment 'i-1' and
#' their ratio median difference is below the threshold, then segment 'i' gets 'i-1''s level.
#' Otherwise, it is 'i-1''s level + 1.
#' 
#' @param threshold A float: threshold for ratio_median difference previously estimated.
#' @param segments A data frame: segment data.
#' 
#' @export

GetSegmentID <- function(threshold, segments){
  # Level column
  segments[,8] <- 0
  segments[1,8] <- 1
  
  for (i in 2:dim(segments)[1]) {
    # Checking for two conditions:
    # 1) If segments are in the same chromosome arm
    # 2) If their ratio_median difference is below the threshold
    #    'large' and 'small' are irrelevant here, for the parameters of RatioDiffBelowThreshold,
    #     since we don't really care about segment size at this step.
    below_threshold = RatioDiffBelowThreshold(large_segment=segments[i,],
                                              small_segment=segments[i-1,],
                                              threshold=threshold)
    if (segments[i,3] == segments[i-1,3] & (below_threshold == TRUE)) {
      segments[i,8] <- segments[i-1,8]
    } else {
      segments[i,8] <- segments[i-1,8] + 1
    }
  }
  segments
}

#' (Kind of) merges segments that lie on the same level.
#' 
#' @description 
#' Iterates through segments in `segments`, and if two segments have the same level,
#' then they are essentially merged with help of `granges_obj`. Actual procedure described
#' in for-loop below.
#' Removes the segments that were not merged.
#' 
#' @param segments A data frame: segment data.
#' @param granges_obj A GRanges object.
#'
#' @export

ShrinkReprTMP <- function(segments, granges_obj) {
  for (i in 1:(dim(segments)[1]-1)) {
    # If segment (i) and next segment (i+1) have same level:
    # 1. Next segment's start position is set to segment's start position, essentially a merge.
    # 2. Then use GRanges to obtain genomic data from this merged segment, and sets 
    #    next segment's ratio_median to this genomic data's ratio_median.
    # 3. Then removes segment i.
    if (segments[i,8] == segments[i+1,8]) {
      segments[i+1,4] <- segments[i,4]
      
      gr = GRanges(seqnames=c(segments[i+1,2]),
                   ranges=IRanges(start=c(segments[i+1,4]), end=c(segments[i+1,5])),
                   strand=c("*"))
      subsetGRobject = subsetByOverlaps(granges_obj, gr)
      
      segments[i+1,6] <- median(subsetGRobject$ratio)
      segments[i,1] <- 0
    }
  }
  
  temp <- which(segments[,1] == 0)
  if (length(temp) != 0) {
    segments <- segments[-temp,]
  }
  
  segments
}

#' Determines if the given segments are LGAs of the given size
#' 
#' @descrption
#' Iterates through the segments in `segments` and firstly checks if the segment meets the definition
#' of an LGA: 
#' 1. Segment and Previous are in the same chromosome arm.
#' 2. Size of space between Segment and Previous is less than 3Mb.
#' Then, we check that:
#' 1. Size of Segment >= `size_lga`
#' 2. Size of Previous >= `size_lga`
#' Finally, we check that ratio_median difference between Segment and Previous > THR. If so, then mark segment as an LGA (via 1's in the output frame).
#'
#' @param threshold A float. Threshold in ratio_median difference previously estimated.
#' @param size_lga An integer: size of the LGA to look for.
#' @param segments A data frame: segment data.
#'
#' @export

DetermineNumberOfLGAs <- function(threshold, size_lga, segments) {
  segments[,1] <- seq(1,dim(segments)[1])
  segments_with_LGA <- segments[,6] * 0
  
  for (i in 2:dim(segments)[1]) {
    # Segment and Previous in same arm, size of space < 3Mb
    if (segments[i,3] == segments[i-1,3] & round((segments[i,4] - segments[i-1,5])/10^6,1) < 3) {
      # Segment & Previous have sizes >= size_lga
      if (round((segments[i,5] - segments[i,4]+1)/10^6,0) >= size_lga & round((segments[i-1,5] - segments[i-1,4]+1)/10^6,0) >= size_lga) {
        # If their median ratio difference > THR*coefficient (coefficient = 1 anyways so there's no param for it)
        below_threshold = RatioDiffBelowThreshold(large_segment=segments[i,],
                                                  small_segment=segments[i-1,],
                                                  threshold=threshold)
        if (!below_threshold) {
          segments_with_LGA[i] <- 1
        }
      }
    }
  }
  segments_with_LGA
}