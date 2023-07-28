library("GenomicRanges")

#' Initializes the process of re-inserting the small segments in between the large ones.
#'
#' @description
#' Essentially deals with the FIRST large segment. We add small segments from `small_segments`
#' segments to `large_segments` that:
#' 1. End before the first large segment and
#' 2. Overlap with the first large segment
#' We stop the process whenever we've changed chromosome arms 
#'
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE. Used to determine
#' whether we insert the small segment or not
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' an overlap of segments and get the ratio_median of this overlap 
#'
#' @export

InitalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_lage = dim(large_segments)[1]
  N_small = dim(small_segments)[1]
  
  halt = 0
  i = 1
  c = 1
  
  while (halt < 1) {
    # If small segment is in a greater chromosome arm than large segment, OR
    # If they are in the same chromosome arm but small segment start  >  large segment end
    # This means that we're no longer looking at the first large segments, so we proceed to the other large segments
    # with the other helper.
    
    # Update halt and stops the loop.
    if ((small_segments[i,3] > large_segments[c,3]) | (small_segments[i,3] == large_segments[c,3] && small_segments[i,4] > large_segments[c,5])) {
      halt = halt + 1
    }
    
    else {
      # If small segment ends before the start of the large segment
      if (small_segments[i,5] < large_segments[c,4]) {
        # If their ratio difference is > THR, then no merge. We add small segment before this large segment
        # in `large_segments` (i.e row bind before)
        if (abs(small_segments[i,6] - large_segments[c,6]) > threshold) {
          large_segments = rbind(c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                   small_segments[i,4], small_segments[i,5], small_segments[i,6],
                                   small_segments[i,5] - small_segments[i,4] + 1, small_segments[i,8]),
                                 large_segments[c:N_large,])
          i = i + 1 # move to next small segment
        } else {
          # Else: ratio difference <= threshold
          # Merged segment info:
          #   Start pos of small_segment, end pos of large_segment
          #   Ratio median calculated w help of GRanges objects
          #   Size is end pos of large segment - start pos of small segment
          
          # We add the data before c+1, not before c, meaning that this is essentially the new
          # segment c.
          gr = GRanges(seqnames = c(small_segments[i,2]),
                       ranges = IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                       strand = c("*"))
          subsetGRobject = subsetByOverlaps(granges_obj, gr)
          
          # Adding data into `large_segments`
          large_segments = rbind(c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                   small_segments[i,4], large_segments[c,5], median(subsetGRobject$ratio),
                                   large_segments[c,5] - small_segments[i,4] + 1, large_segments[c,8]),
                                 large_segments[(c+1):N_large,])
          i = i + 1
        }
      # If small segment start pos is before large segment start pos AND
      # If small segment end pos is after large segment start pos AND
      # If small small segment end pos is before large segment end pos
      # |----small----|
      #       |----large----|
      # *******
      } else if (small_segments[i,4] < large_segments[c,4] &&
                 small_segments[i,5] >= large_segments[c,4] &&
                 small_segments[i,5] <= large_segments[c,5]) {
        # If above, threshold, we add the segment in ***'s to `large_segments`
        if (abs(small_segments[i,6] - large_segments[c,6]) > threshold) {
          gr = GRanges(seqnames = c(small_segments[i,2]),
                       ranges = IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,4])),
                       strand = c("*"))
          subsetGRobject = subsetByOverlaps(granges_obj, gr)
          
          # Adding
          large_segments = rbind(c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                   small_segments[i,4], large_segments[c,4], median(subsetGRobject$ratio),
                                   large_segments[c,4] - small_segments[i,4] + 1, small_segments[c,8]),
                                 large_segments[c:N_large,])
          i = i + 1
        } else {
          # Else: below THR, we merge
          gr = GRanges(seqnames = c(large_segments[c,2]),
                       ranges = IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                       strand = c("*"))
          subsetGRobject = subsetByOverlaps(granges_obj, gr)
          
          # Adding
          large_segments = rbind(c(large_segments[i,1], large_segments[i,2], large_segments[i,3],
                                   small_segments[i,4], large_segments[c,5], median(subsetGRobject$ratio),
                                   large_segments[c,5] - small_segments[i,4] + 1, large_segments[c,8]),
                                 large_segments[(c+1):N_large,])
          i = i + 1
        }
      } else {
        # None of the 2 cases
        i = i +1
      }
    }
  }
  
  # Remove the rows from `small_segments` that we added into `large_segments`
  small_segments = small_segments[i:N_small,]
  
  output = c(large_segments, small_segments)
  output
}

#' Re-inserts small segments in between the large ones.
#'
#' @description
#' We add small segments from `small_segments` to `large_segments` in 6 cases, 
#' each depending on where the small segment fits in between the large ones. 
#' We merge segments if the ratio difference between the small and the large are <= THR. 
#' Otherwise, we keep small segment as its own segment in `large_segments`. 
#' 
#' Visually, these are the 6 cases:
#' 1. |--small-|    |-----large----|
#' 
#' 2. |----- previous large -----|    |----small----|
#'                                         |---------- large segment----------|
#'                                         
#' 3.   |------------- large ---------------|
#'          |----small----|
#'          
#' 4. |---------large----------|        |------------- next large--------|
#'               |--------small--------|
#'               
#' 5.   |-------large--------|                    |--------------next large-------------|
#'                            |----small----|   
#'                            
#' 6. chr N                                          |#| chr N+1
#' |-------- large --------|                         |#| |-------- next_large --------|
#'                          |------ small ------|    |#|        
#'                                                 
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE. Used to determine
#' whether we insert the small segment or not
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' an overlap of segments and get the ratio_median of this overlap 
#'
#' @export

FinalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_small = dim(small_segments)[1]
  i = 1
  
  while (i < N_small + 1) {
    c = 1
    N_large = dim(large_segments)[1]
    
    while  (c < N_large + 1) {
      # If in the same chromosome arm
      if (small_segments[i,3] == large_segments[c,3]) {
        below_threshold = RatioDiffBelowThreshold(large_segment = large_segments[c,],
                                                  small_segment = small_segments[i,],
                                                  threshold = threshold)
        #### First case ####
        if (CaseOne(large_segment = large_segments[c,], small_segment = small_segments[i,]) == TRUE) {
          # If below threshold: merge, and small+large become new segment c
          # So, new segment's data is:
          #   Start: small segment start pos
          #   End: large segment end pos
          #   Median: calculate w/ help of Granges
          #   Size: large segment end - small segment start 
          
          # Finally, we squeeze it in between segments :c-1 and c+1:
          if (below_threshold == TRUE) {
            gr = GRanges(seqnames = c(small_segments[i,2]),
                         ranges = IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                         strand = c("*"))
            subsetGRobject = subsetByOverlaps(granges_obj, gr)
            
            large_segments = rbind(large_segments[1:(c-1),],
                                   c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                     small_segments[i,4], large_segments[c,5], median(subsetGRobject$ratio), 
                                     large_segments[c,5] - small_segments[i,4] + 1, large_segments[c,8]),
                                   large_segments[(c+1):N_large,])
            i = i + 1
            c = N_large + 1
          } else {
            # Else: not below threshold: small_segment is its own separate thing
            # Add small_segment data, in between :c-1 and c:
            large_segments = rbind(large_segments[1:(c-1),],
                                   c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                     small_segments[i,4], small_segments[i,5], small_segments[i,6], 
                                     small_segments[i,5] - small_segments[i,4] + 1, small_segments[i,8]),
                                   large_segments[c:N_large,])
            i = i + 1
            c = N_large + 1
          }
          
        #### Second Case ####
        } else if (CaseTwo(prev_large_segment = large_segments[c-1,], small_segment = small_segments[i,], large_segment = large_segments[c,]) == TRUE) {
          # Ratio diff <= THR: merge
          if (below_threshold == TRUE) {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                 large_segments = large_segments, c = c, N_large = N_large,
                                                 below_thr = TRUE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                 large_segments = large_segments, c = c, N_large = N_large,
                                                 below_thr = TRUE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                 large_segments = large_segments, c = c, N_large = N_large,
                                                 below_thr = FALSE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                 large_segments = large_segments,c = c, N_large = N_large,
                                                 below_thr = FALSE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          }
        #### Third Case #### 
        } else if (CaseThree(large_segment = large_segments[c,], small_segment = small_segments[i,]) == TRUE) {
          i = i + 1
        #### Fourth Case ####
        } else if (CaseFour(large_segment = large_segments[c,], small_segment = small_segments[i,], next_large_segment = large_segments[c+1,]) == TRUE) {
          # Ratio diff <= THR: merge
          if (below_threshold == TRUE) {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                 large_segments = large_segments, c = c, N_large = N_large,
                                                 below_thr = TRUE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                  large_segments = large_segments, c = c, N_large = N_large,
                                                  below_thr = TRUE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segments[i,],
                                                  large_segments = large_segments,c = c, N_large = N_large,
                                                  below_thr = FALSE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj,
                                                  small_segment = small_segments[i,],
                                                  large_segments = large_segments,
                                                  c = c,
                                                  N_large = N_large,
                                                  below_thr = FALSE,
                                                  end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          }
        } else if (CaseFive(large_segment = large_segments[c,], small_segment = small_segments[i,], next_large_segment = large_segments[c+1,]) == TRUE) {
          small_nextlarge_below_threshold = RatioDiffBelowThreshold(large_segment = large_segments[c+1,],
                                                                    small_segment = small_segments[i,],
                                                                    threshold = threshold)
          # Ratio diff <= THR for large & small AND next_large & small:
          if (below_threshold == TRUE && small_nextlarge_below_threshold == TRUE) {
            if (abs(small_segments[i,6] - large_segments[c,6]) <= abs(small_segments[i,6] - large_segments[c+1,6])) {
              
            }
          }
        }
      }
    }
  }
}

RatioDiffBelowThreshold <- function(large_segment, small_segment, threshold) {
  output = (abs(small_segment[6] - large_segment[6])) <= threshold
  output
}

MergeSegmentsTwo <- function(granges_obj, small_segment, large_segments, c, N_large, below_thr, end_of_file) {
  if (below_thr == TRUE) {
    gr = GRanges(seqnames = c(large_segment[c,2]),
                 ranges = IRanges(start=c(small_segment[4]), end=c(large_segments[c,5])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(large_segments[c,1], 
                               large_segments[c,2], 
                               large_segments[c,3],
                               small_segment[4],
                               large_segments[c,5], 
                               median(subsetGRobject$ratio),
                               large_segments[c,5] - small_segment[4] + 1, 
                               large_segments[c,8]))
    } else {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(large_segments[c,1], 
                               large_segments[c,2], 
                               large_segments[c,3],
                               small_segment[4], 
                               large_segments[c,5], 
                               median(subsetGRobject$ratio),
                               large_segments[c,5] - small_segment[4] + 1, 
                               large_segments[c,8]),
                             large_segments[(c+1):N_large,])
    }
  } else {
    gr = GRanges(seqnames = c(small_segment[2]),
                 ranges = IRanges(start=c(small_segment[4]), end=c(large_segments[c,4])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               small_segment[4], 
                               large_segments[c,4], 
                               median(subsetGRobject$ratio),
                               large_segments[c,4] - small_segment[4] + 1, 
                               small_segment[8]),
                             large_segments[c:N_large,])
    } else {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(small_segment[1], 
                               small_segment[2], 
                               small_segment[3],
                               small_segment[4], 
                               large_segments[c,4], 
                               median(subsetGRobject$ratio),
                               large_segments[c,4] - small_segment[4] + 1, 
                               small_segment[8]))
    }
  }
  large_segments
}

MergeSegmentsFour <- function(granges_obj, small_segments,i, large_segments, c, N_large, below_thr, end_of_file) {
  if (below_thr == TRUE) {
    gr = GRanges(seqnames = c(large_segment[c,2]),
                 ranges = IRanges(start=c(large_segment[4]), end=c(small_segment[c,5])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(large_segments[c,1], 
                               large_segments[c,2], 
                               large_segments[c,3],
                               large_segments[c,4],
                               small_segments[i,5], 
                               median(subsetGRobject$ratio),
                               small_segments[i,5] - large_segments[c,4] + 1, 
                               large_segments[c,8]))
    } else {
      large_segments = rbind(large_segments[1:(c-1),],
                             c(large_segments[c,1], 
                               large_segments[c,2], 
                               large_segments[c,3],
                               large_segments[c,4],
                               small_segments[i,5], 
                               median(subsetGRobject$ratio),
                               small_segments[i,5] - large_segments[c,4] + 1, 
                               large_segments[c,8]),
                             large_segments[(c+1):N_large,])
    }
  } else {
    gr = GRanges(seqnames = c(small_segment[2]),
                 ranges = IRanges(start=c(large_segments[c,5]), end=c(small_segment[5])),
                 strand = c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    
    # Adding into `large_segments`. But this depends whether we're at the end of file or not.
    if (end_of_file == TRUE) {
      large_segments = rbind(large_segments[1:c,],
                             c(small_segments[i,1], 
                               small_segments[i,2], 
                               small_segments[i,3],
                               large_segments[c,5], 
                               small_segments[c, 5], 
                               median(subsetGRobject$ratio),
                               small_segments[c,5] - large_segments[c,5] + 1, 
                               small_segments[i,8]))
    } else {
      large_segments = rbind(large_segments[1:c,],
                             c(small_segments[i,1], 
                               small_segments[i,2], 
                               small_segments[i,3],
                               large_segments[c,5], 
                               small_segments[c, 5], 
                               median(subsetGRobject$ratio),
                               small_segments[c,5] - large_segments[c,5] + 1, 
                               small_segments[i,8]),
                             large_segments[(c+1):N_large,])
    }
  }
  large_segments
}
  

CaseOne <- function(large_segment, small_segment) {
  output = small_segment[5] < large_segment[4]
  output
}

CaseTwo <- function(prev_large_segment, small_segment, large_segment) {
  output =  (small_segment[4] < large_segment[4]) && 
    (small_segment[5] >= large_segment[4]) && 
    (small_segment[5] <= large_segment[5]) &&
    (small_segment[4] >= prev_large_segment[5])
  
  output
}

CaseThree <- function(large_segment, small_segment) {
  output = (small_segment[4] >= large_sement[4]) && (small_segment[5] <= large_segment[5])
  
  output
}

CaseFour <- function(large_segment, small_segment, next_large_segment) {
  output = (small_segment[4] >= large_segment[4]) &&
    (small_segment[4] <= large_segment[5]) &&
    (small_segment[5] > large_segment[5]) &&
    (small_segment[5] <= next_large_segment[5])
  
  output
}

CaseFive

CaseSix
