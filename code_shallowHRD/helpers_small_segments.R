library("GenomicRanges")

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

MergeSegmentsFour <- function(granges_obj, small_segments, i, large_segments, c, N_large, below_thr, end_of_file) {
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

MergeSegmentsFiveLargeSmall <- function(granges_obj, small_segment, large_segments, j, N_large) {
  large_segment = large_segments[j,]
  
  gr = GRanges(seqnames = c(large_segment[2]),
               ranges = IRanges(start=c(large_segment[4]), end=c(small_segment[5])),
               strand = c("*"))
  subsetGRobject = subsetByOverlaps(granges_obj, gr)
  
  large_segments = rbind(large_segments[1:(j-1),],
                         c(large_segment[1], large_segment[2], large_segment[3],
                           large_segment[4], small_segment[5], median(subsetGRobject$ratio),
                           small_segment[5] - large_segment[4] + 1, large_segment[8]),
                         large_segments[(j+1):N_large,])
  
  large_segments
}
MergeSegmentsFiveSmallNextLarge <- function(granges_obj, small_segment, large_segments, j, N_large, end_of_file) {
  large_segment = large_segments[j,]
  
  gr = GRanges(seqnames = c(large_segment[2]),
               ranges = IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
               strand = c("*"))
  subsetGRobject = subsetByOverlaps(granges_obj, gr)
  
  if (end_of_file == TRUE) {
    large_segments = rbind(large_segments[1:j-1,],
                           c(large_segment[1], large_segment[2], large_segment[3],
                             small_segment[4], large_segment[5], median(subsetGRobject$ratio),
                             large_segment[5] - small_segment[4] + 1, large_segment[8]))
  } else {
    large_segments = rbind(large_segments[1:j-1,],
                           c(large_segment[1], large_segment[2], large_segment[3],
                             small_segment[4], large_segment[5], median(subsetGRobject$ratio),
                             large_segment[5] - small_segment[4] + 1, large_segment[8]),
                           large_segments[(j+1):N_large,])
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

CaseFive <- function(large_segment, small_segment, next_large_segment) {
  output = (small_segment[4] >= large_segment[5]) &&
    (small_segment[5] <= next_large_segment[4]) &&
    (large_segment[3] == next_large_segment[3])
  
  output
}

CaseSix <- function(large_segment, small_segment, next_large_segment) {
  output = (small_segment[4] >= large_segment[5]) &&
    (large_segment[3] != next_large_segment[3]) |
    (is.null(next_large_segment[3]))
  
  output
}
