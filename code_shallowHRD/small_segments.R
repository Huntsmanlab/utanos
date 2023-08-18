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
  
  granges_obj = makeGRangesFromDataFrame(bam_ratios_frame, keep.extra.columns = TRUE, ignore.strand = TRUE, 
                                        seqinfo = NULL, seqnames.field = "chr", start.field = "start",
                                        end.field = "end")
  granges_obj
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


#' Deals with missing chromosome arm values in the given data frame of small segments.
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
#' into the `large_segments` data frame. Merges if necessary (i.e. if ratio difference < threshold).
#'
#' @description
#' Does a lot of things. In InitializeSmallSegments, we insert all the small segments that come before the first
#' large segment in `large_segments` (hence initializing). In FinalizeSmallSegments, we insert all the remaining small segments,
#' according to 6 different cases, which are listed in the FinalizeSmallSegments documentation. In the end, returns the updated
#' `large_segments` with small segments appropriately re-inserted.
#'
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= Mb.
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
  
  print("Exited InitializeSmallSegments")
  large_segments <- FinalizeSmallSegments(large_segments = as.data.frame(output[1]),
                                          small_segments = as.data.frame(output[2]),
                                          threshold = threshold,
                                          granges_obj = granges_obj)
  print("Exited FinalizeSmall")
  large_segments
}

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
#' whether we insert the small segment or not.
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' an overlap of segments and get the ratio_median of this overlap. 
#'
#' @export

InitalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_large = dim(large_segments)[1]
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
        # If their ratio difference is > threshold, then no merge. We add small segment before this large segment
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
      } else if (small_segments[i,4] < large_segments[c,4] && small_segments[i,5] >= large_segments[c,4] && small_segments[i,5] <= large_segments[c,5]) {
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
          # Else: below threshold, we merge
          gr = GRanges(seqnames = c(large_segments[c,2]),
                       ranges = IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                       strand = c("*"))
          subsetGRobject = subsetByOverlaps(granges_obj, gr)
          
          # Adding
          large_segments = rbind(c(large_segments[c,1], large_segments[c,2], large_segments[c,3],
                                   small_segments[i,4], large_segments[c,5], median(subsetGRobject$ratio),
                                   large_segments[c,5] - small_segments[i,4] + 1, large_segments[c,8]),
                                 large_segments[(c+1):N_large,])
          i = i + 1
        }
      } else {
        # None of the 2 cases
        i = i + 1
      }
    }
  }
  
  # Remove the rows from `small_segments` that we added into `large_segments`
  small_segments = small_segments[i:N_small,]
  
  output = list(large_segments, small_segments)
  output
}

#' Re-inserts small segments in between the large ones.
#'
#' @description
#' We add small segments from `small_segments` to `large_segments` in 6 cases, 
#' each depending on where the small segment fits in between the large ones. 
#' We merge segments if the ratio difference between the small and the large are <= threshold. 
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
#' @param large_segments A data frame. Segments with size >= 3Mb.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE. Used to determine
#' whether we insert the small segment or not.
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' an overlap of segments and get the ratio_median of this overlap. 
#'
#' @export

FinalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_small = dim(small_segments)[1]
  i = 1
  while (i < N_small + 1) {
    c = 1
    N_large = dim(large_segments)[1]
    
    while  (c < N_large + 1) {   
      print("checking...")
      print(c)
      print(N_large)
      # If in the same chromosome arm
      if (small_segments[i,3] == large_segments[c,3]) {
        
        # Getting the segments we'll work with
        large_segment = as.numeric(large_segments[c,])
        small_segment = as.numeric(small_segments[i,])
        
        # Ratio difference below threshold?
        below_threshold = RatioDiffBelowThreshold(large_segment = large_segment, small_segment = small_segment, threshold = threshold)
        
        #### First case ####
        if (small_segment[5] < large_segment[4]) {
          # If below threshold: merge, and small+large become new segment c
          # So, new segment's data is:
          #   Start: small segment start pos
          #   End: large segment end pos
          #   Median: calculate w/ help of Granges
          #   Size: large segment end - small segment start 
          
          # Finally, we squeeze it in between segments :c-1 and c+1:
          if (below_threshold == TRUE) {
            gr = GRanges(seqnames = c(small_segment[2]),
                         ranges = IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                         strand = c("*"))
            subsetGRobject = subsetByOverlaps(granges_obj, gr)
            
            large_segments = rbind(large_segments[1:(c-1),],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], large_segment[5], median(subsetGRobject$ratio), 
                                     large_segment[5] - small_segment[4] + 1, large_segment[8]),
                                   large_segments[(c+1):N_large,])
            i = i + 1
            c = N_large + 1
          } else {
            # Else: not below threshold: small_segment is its own separate thing
            # Add small_segment data, in between :c-1 and c:
            large_segments = rbind(large_segments[1:(c-1),],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6], 
                                     small_segment[5] - small_segment[4] + 1, small_segment[8]),
                                   large_segments[c:N_large,])
            i = i + 1
            c = N_large + 1
          }
          
          #### Second Case ####
        } else if (small_segment[4] < large_segment[4] && small_segment[5] >= large_segment[4] && small_segment[5] <= large_segment[5] && small_segment[4] >= large_segments[c-1,5]) {
          # Ratio diff <= THR: merge
          if (below_threshold == TRUE) {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, 
                                                 small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, 
                                                 small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
            # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, 
                                                 small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          }
          #### Third Case #### 
        } else if (small_segment[4] >= large_segment[4] && small_segment[5] <= large_segment[5]) {
          i = i + 1
          #### Fourth Case ####
        } else if (small_segment[4] >= large_segment[4] && small_segment[4] <= large_segment[5] && small_segment[5] > large_segment[5] && small_segment[5] <= large_segments[c+1,5]) {
          # Ratio diff <= THR: merge
          if (below_threshold == TRUE) {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
            # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
          }
          #### Case Five ####
        } else if (c == N_large && small_segment[4] >= large_segment[5]) {
          if (below_threshold == TRUE) {
            large_segments <- MergeSegmentsFiveLargeSmall(granges_obj = granges_obj, 
                                                          small_segment = small_segment,
                                                          large_segments = large_segments, 
                                                          j = c, 
                                                          N_large = N_large)
            i = i + 1
            c = N_large + 1
          } else {
            large_segments = rbind(large_segments[1:c,],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6],
                                     small_segment[7], small_segment[8]))
            i = i + 1
            c = N_large + 1
          }
        } else if (c+1 <= N_large && small_segment[4] >= large_segment[5] && small_segment[5] <= large_segments[c+1,4] && large_segment[3] == large_segments[c+1,3]) {
          small_nextlarge_below_threshold = RatioDiffBelowThreshold(large_segment = large_segments[c+1,],
                                                                    small_segment = small_segment,
                                                                    threshold = threshold)
          # Ratio diff <= THR for large & small AND next_large & small:
          if (below_threshold == TRUE && small_nextlarge_below_threshold == TRUE) {
            # Large segment & small segment ratio difference <= Next_large segment & small segment ratio difference
            # This essentially means that large segment is closer to small than next_large.
            
            # So, we merge large and small.
            if (abs(small_segment[6] - large_segment[6]) <= abs(small_segment[6] - large_segments[c+1,6])) {
              large_segments <- MergeSegmentsFiveLargeSmall(granges_obj = granges_obj, 
                                                            small_segment = small_segment,
                                                            large_segments = large_segments, 
                                                            j = c, 
                                                            N_large = N_large)
              i = i + 1
              c = N_large + 1
            } else {
              # Else: next_large and small are closer than large and small.
              # So, we merge small and next_large. Gotta check if c+2 even exists.
              # If c+2 is end of file:
              if (c + 2 > N_large) {
                large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj, 
                                                                  small_segment = small_segment,
                                                                  large_segments = large_segments, 
                                                                  j = c+1, 
                                                                  end_of_file = TRUE)
                i = i + 1
                c = N_large + 1
              } else {
                # Else: Not end of file  
                large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj, 
                                                                  small_segment = small_segment,
                                                                  large_segments = large_segments, 
                                                                  j = c+1, 
                                                                  end_of_file = FALSE)
                i = i + 1
                c = N_large + 1
              }
            }
            # Only large & small ratio diff <= THR
          } else if (below_threshold == TRUE) {
            large_segments <- MergeSegmentsFiveLargeSmall(granges_obj = granges_obj, 
                                                          small_segment = small_segment,
                                                          large_segments = large_segments, 
                                                          j = c)
            i = i + 1
            c = N_large + 1
            # Only small & next_large ratio diff <= THR
          } else if (small_nextlarge_below_threshold == TRUE) {
            # If end of file:
            if (c + 2 > N_large) {
              large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj, 
                                                                small_segment = small_segment,
                                                                large_segments = large_segments, 
                                                                j = c+1, 
                                                                end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            } else {
              # Else: not end of file
              # Checking if next small segment ends before next large segment (and same chr arm)
              # i.e. if there's a next_small segment in between small segment and next_large.
              if (small_segments[i+1,5] < large_segments[c+1,4] && small_segments[i+1,3] == large_segments[c+1,3]) {
                gr=GRanges(seqnames=c(small_segment[2]),
                           ranges=IRanges(start=c(small_segment[4]), end=c(small_segment[5])),
                           strand=c("*"))
                subsetGRobject = subsetByOverlaps(granges_obj, gr)
                
                large_segments=rbind(large_segments[1:c,], 
                                     c(small_segment[1], small_segment[2], small_segment[3], 
                                       small_segment[4], small_segment[5], median(subsetGRobject$ratio), 
                                       small_segment[5] - small_segment[4] + 1, small_segment[8]) , 
                                     large_segments[(c+1):N_large,])  
                i= i + 1
                c= N_large + 1
              } else {
                gr=GRanges(seqnames=c(small_segment[2]),
                           ranges=IRanges(start=c(small_segment[4]),end=c(large_segments[c+1,5])),
                           strand=c("*"))
                subsetGRobject = subsetByOverlaps(granges_obj, gr)
                
                large_segments=rbind(large_segments[1:c,], 
                                     c(large_segments[c+1,1], large_segments[c+1,2], large_segments[c+1,3], 
                                       small_segment[4], large_segments[c+1,5], median(subsetGRobject$ratio), 
                                       large_segments[c+1,5] - small_segment[4] + 1, large_segments[c+1,8]) , 
                                     large_segments[(c+2):N_large,])  
                i= i + 1
                c= N_large + 1
              }
            }
            # No segment ratio differences are <= THR: add small segment as its own thing
          } else {
            large_segments = rbind(large_segments[1:c,],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6],
                                     small_segment[7], small_segment[8]),
                                   large_segments[(c+1):N_large,])
            i = i + 1
            c = N_large + 1
          }
          #### Case Six ####
        } else if (small_segment[4] >= large_segment[5] && large_segment[3] != large_segments[c+1,3] | is.null(large_segments[c+1,3])) {
          # If large and small ratio diff <= THR
          if (below_threshold == TRUE) {
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, 
                                                  small_segment = small_segment,
                                                  large_segments = large_segments, 
                                                  j = c, 
                                                  below_thr = TRUE, 
                                                  end_of_file = TRUE)
              i = i + 1
              c = N_large + 1
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, 
                                                  small_segment = small_segment,
                                                  large_segments = large_segments, 
                                                  j = c, 
                                                  below_thr = TRUE, 
                                                  end_of_file = FALSE)
              i = i + 1
              c = N_large + 1
            }
            # Else: large and small ratio diff > THR
          } else {
            if (c + 1 > N_large) {
              large_segments = rbind(large_segments[1:c,],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], small_segment[5], small_segment[6],
                                       small_segment[7], small_segment[8]))
              i = i + 1
              c = N_large + 1
            } else {
              large_segments = rbind(large_segments[1:c,],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], small_segment[5], small_segment[6],
                                       small_segment[7], small_segment[8]),
                                     large_segments[(c+1):N_large,])
              i = i + 1
              c = N_large + 1
            }
          }
          # None of the 6 cases
        } else {
          c = c + 1
        }
        # Different chromosome arms
      } else {
        c = c + 1
      }
    }
  }
  rownames(large_segments) <- NULL
  large_segments
}

#' Adds 'leftsovers' to the segment frame after re-inserting small segments.
#'
#' @description
#' Iterates over the `segments` frame and if the space between adjacent segments is
#' greater than one, then we use GRanges to obtain genomic data from this missing region in the
#' data frame. We update the adjacent segment's start or end positions depending on 
#' which segment is closer to the missing segment.
#' 
#' Do this as long as adjacent segments are in the same chromosome arm. 
#'
#' @param segments A data frame. Segments with size >= 3Mb.
#' @param second_round A boolean. If true, some indices where we check in the frame
#' are different.
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' empty space between segments and get the ratio_median of this missing portion. 
#' 
#' @export
CorrectLeftovers <- function(segments, second_round, granges_obj) {
  c = 1 
  
  while (c < dim(segments)[1]) {
    # If in the same chromosome arm
    if (segments[c,3] == segments[c+1,3]) {
      # If space between them is > 1 (i.e. missing data)
      if (segments[c+1,4] - segments[c,5] > 1) {
        gr = GRanges(seqnames=c(segments[c,2]),
                     ranges=IRanges(start=c(segments[c,5]), end=c(segments[c+1,4])),
                     strand=c("*"))
        subsetGRobject = subsetByOverlaps(granges_obj, gr)
        
        # Figure out which adjacent segment's ratio_median is closest
        # to the subsetGRobject's ratio_medain
        segment_c_grobj = abs(segments[c,6] - median(subsetGRobject$ratio))
        segment_c1_grobj = abs(segments[c+1,6] - median(subsetGRobject$ratio))
        
        # If segment c is closer than segment c+1
        # Then segment c's end position is c+1's start position - 1
        # Otherwise:
        #     segment c+1's start position is c's end position + 1
        if (segment_c_grobj < segment_c1_grobj) {
          segments[c,5] = segments[c+1,4] - 1
        } else {
          segments[c+1,4] = segments[c,5] + 1
        }
      }
    }
    c = c + 1
  }
  
  rownames(segments) <- NULL
  colnames(segments) <= c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")
  
  # Updating sizes
  segments[,7] = segments[,5] - segments[,4] + 1
  segments
}




