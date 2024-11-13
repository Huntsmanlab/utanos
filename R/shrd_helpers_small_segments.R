#' Initializes the process of re-inserting the small segments in between the large ones.
#'
#' @description
#' Essentially deals with the FIRST large segment. We add small segments from `small_segments`
#' segments to `large_segments` that:
#' 1. End before the first large segment and
#' 2. Overlap with the first large segment
#' We stop the process whenever we've changed chromosome arms between the small segment and the large segment
#'
#' @param small_segments A data frame. Segments with 3Mb >= size >= 0.1Mb.
#' @param large_segments A data frame. Segments with size >= Mb.
#' @param threshold A float: the estimated threshold for ratio_median difference via KDE. Used to determine
#' whether we insert the small segment or not.
#' @param granges_obj A GRanges object: is used as genomic reference to check whenever we have
#' an overlap of segments and get the data (specifically the ratio_median) of this overlap.
#'

InitalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_large = dim(large_segments)[1]
  N_small = dim(small_segments)[1]

  # In some cases we have no small segments; therefore return right away
  if (N_small == 0) {
    return(list(large_segments, small_segments))
  }

  halt = 0 # whether we change chromosome arms, or small segment is after first large segment
  i = 1 # small segment
  c = 1 # large segment, always equal to 1, since we only care about the first one at this stage.

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
        # in `large_segments` (i.e row bind small segment to the left of the large_segments: before)
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
          gr = GenomicRanges::GRanges(seqnames = c(small_segments[i,2]),
                                      ranges = IRanges::IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                                      strand = c("*"))
          subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

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
          gr = GenomicRanges::GRanges(seqnames = c(small_segments[i,2]),
                                      ranges = IRanges::IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,4])),
                                      strand = c("*"))
          subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

          # Adding
          large_segments = rbind(c(small_segments[i,1], small_segments[i,2], small_segments[i,3],
                                   small_segments[i,4], large_segments[c,4], median(subsetGRobject$ratio),
                                   large_segments[c,4] - small_segments[i,4] + 1, small_segments[c,8]),
                                 large_segments[c:N_large,])
          i = i + 1
        } else {
          # Else: below threshold, we merge
          gr = GenomicRanges::GRanges(seqnames = c(large_segments[c,2]),
                                      ranges = IRanges::IRanges(start=c(small_segments[i,4]), end=c(large_segments[c,5])),
                                      strand = c("*"))
          subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

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
  # Ensure we keep the required column names after adding `small_segments` into `large_segments`; this prevents getting any empty data frame when indexing with column name
  large_segments = as.data.frame(large_segments)
  names(large_segments) = c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

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

FinalizeSmallSegments <- function(large_segments, small_segments, threshold, granges_obj) {
  N_small = dim(small_segments)[1]
  i = 1
  while (i < N_small + 1) {
    c = 1
    N_large = dim(large_segments)[1]
    while  (c < N_large + 1) {
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
            # End of file
            if (c + 1 > N_large) {
              gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                          ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                                          strand = c("*"))
              subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

              large_segments = rbind(large_segments[1:(c-1),],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], large_segment[5], median(subsetGRobject$ratio),
                                       large_segment[5] - small_segment[4] + 1, large_segment[8]))
              c = N_large + 1
              i = i + 1
            }
            else {
              gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                          ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                                          strand = c("*"))
              subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

              large_segments = rbind(large_segments[1:(c-1),],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], large_segment[5], median(subsetGRobject$ratio),
                                       large_segment[5] - small_segment[4] + 1, large_segment[8]),
                                     large_segments[(c+1):N_large,])
              c = N_large + 1
              i = i + 1
            }
          } else {
            # Else: not below threshold: small_segment is its own separate thing
            # Add small_segment data, in between :c-1 and c:
            large_segments = rbind(large_segments[1:(c-1),],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6],
                                     small_segment[5] - small_segment[4] + 1, small_segment[8]),
                                   large_segments[c:N_large,])
            c = N_large + 1
            i = i + 1
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
              c = N_large + 1
              i = i + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj,
                                                 small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
            }
            # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj,
                                                 small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = TRUE)
              c = N_large + 1
              i = i + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsTwo(granges_obj = granges_obj, small_segment = small_segment,
                                                 large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
            }
          }
          #### Third Case ####
        } else if (small_segment[4] >= large_segment[4] && small_segment[5] <= large_segment[5]) {
          i = i + 1
          c = N_large + 1
          #### Fourth Case ####
        } else if (small_segment[4] >= large_segment[4] && small_segment[4] <= large_segment[5] && small_segment[5] > large_segment[5] && small_segment[5] <= large_segments[c+1,5]) {
          # Ratio diff <= THR: merge
          if (below_threshold == TRUE) {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = TRUE)
              c = N_large + 1
              i = i + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = TRUE, end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
            }
            # Ratio diff > THR: no merge
          } else {
            # End of file
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = TRUE)
              c = N_large + 1
              i = i + 1
              # Not end of file
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj, small_segment = small_segment,
                                                  large_segments = large_segments, j = c, below_thr = FALSE, end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
            }
          }
          #### Case Five ####
        } else if (c == N_large && small_segment[4] >= large_segment[5]) {
          if (below_threshold == TRUE) {
            large_segments <- MergeSegmentsFiveLargeSmall(granges_obj = granges_obj,
                                                          small_segment = small_segment,
                                                          large_segments = large_segments,
                                                          j = c,
                                                          end_of_file = TRUE)
            c = N_large + 1
            i = i + 1
          } else {
            large_segments = rbind(large_segments[1:c,],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6],
                                     small_segment[7], small_segment[8]))
            c = N_large + 1
            i = i + 1
          }
        } else if (small_segment[4] >= large_segment[5] && small_segment[5] <= large_segments[c+1,4] && large_segment[3] == large_segments[c+1,3]) {
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
                                                            end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
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
                c = N_large + 1
                i = i + 1
              } else {
                # Else: Not end of file
                large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj,
                                                                  small_segment = small_segment,
                                                                  large_segments = large_segments,
                                                                  j = c+1,
                                                                  end_of_file = FALSE)
                c = N_large + 1
                i = i + 1
              }
            }
            # Only large & small ratio diff <= THR
          } else if (below_threshold == TRUE) {
            large_segments <- MergeSegmentsFiveLargeSmall(granges_obj = granges_obj,
                                                          small_segment = small_segment,
                                                          large_segments = large_segments,
                                                          j = c,
                                                          end_of_file = FALSE)
            c = N_large + 1
            i = i + 1
            # Only small & next_large ratio diff <= THR
          } else if (small_nextlarge_below_threshold == TRUE) {
            # If end of file:
            if (c + 2 > N_large) {
              large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj,
                                                                small_segment = small_segment,
                                                                large_segments = large_segments,
                                                                j = c+1,
                                                                end_of_file = TRUE)
              c = N_large + 1
              i = i + 1
            } else if (c + 2 <= N_large && i == N_small) {
              large_segments <- MergeSegmentsFiveSmallNextLarge(granges_obj = granges_obj,
                                                                small_segment = small_segment,
                                                                large_segments = large_segments,
                                                                j = c+1,
                                                                end_of_file = FALSE)
              c = N_large + 1
              i = i + 1
            } else {
              # Else: not end of file
              # Checking if next small segment ends before next large segment (and same chr arm)
              # i.e. if there's a next_small segment in between small segment and next_large.
              if (small_segments[i+1,5] < large_segments[c+1,4] && small_segments[i+1,3] == large_segments[c+1,3]) {
                gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                            ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(small_segment[5])),
                                            strand = c("*"))
                subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

                large_segments = rbind(large_segments[1:c,],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], small_segment[5], median(subsetGRobject$ratio),
                                       small_segment[5] - small_segment[4] + 1, small_segment[8]) ,
                                     large_segments[(c+1):N_large,])
                c= N_large + 1
                i = i + 1
              } else {
                gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                            ranges = IRanges::IRanges(start=c(small_segment[4]),end=c(large_segments[c+1,5])),
                                            strand = c("*"))
                subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

                large_segments = rbind(large_segments[1:c,],
                                     c(large_segments[c+1,1], large_segments[c+1,2], large_segments[c+1,3],
                                       small_segment[4], large_segments[c+1,5], median(subsetGRobject$ratio),
                                       large_segments[c+1,5] - small_segment[4] + 1, large_segments[c+1,8]) ,
                                     large_segments[(c+2):N_large,])
                c= N_large + 1
                i = i + 1
              }
            }
            # No segment ratio differences are <= THR: add small segment as its own thing
          } else {
            large_segments = rbind(large_segments[1:c,],
                                   c(small_segment[1], small_segment[2], small_segment[3],
                                     small_segment[4], small_segment[5], small_segment[6],
                                     small_segment[7], small_segment[8]),
                                   large_segments[(c+1):N_large,])
            c = N_large + 1
            i <- i + 1
          }
          #### Case Six ####
        } else if (small_segment[4] >= large_segment[5] && (large_segment[3] != large_segments[c+1,3] | is.null(large_segments[c+1,3]))) {
          # If large and small ratio diff <= THR
          if (below_threshold == TRUE) {
            if (c + 1 > N_large) {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj,
                                                  small_segment = small_segment,
                                                  large_segments = large_segments,
                                                  j = c,
                                                  below_thr = TRUE,
                                                  end_of_file = TRUE)
              c = N_large + 1
              i <- i + 1
            } else {
              large_segments <- MergeSegmentsFour(granges_obj = granges_obj,
                                                  small_segment = small_segment,
                                                  large_segments = large_segments,
                                                  j = c,
                                                  below_thr = TRUE,
                                                  end_of_file = FALSE)
              c = N_large + 1
              i <- i + 1
            }
            # Else: large and small ratio diff > THR
          } else {
            if (c + 1 > N_large) {
              large_segments = rbind(large_segments[1:c,],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], small_segment[5], small_segment[6],
                                       small_segment[7], small_segment[8]))
              c = N_large + 1
              i <- i + 1
            } else {
              large_segments = rbind(large_segments[1:c,],
                                     c(small_segment[1], small_segment[2], small_segment[3],
                                       small_segment[4], small_segment[5], small_segment[6],
                                       small_segment[7], small_segment[8]),
                                     large_segments[(c+1):N_large,])
              c = N_large + 1
              i <- i + 1
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
#' @param granges_obj A GRanges object: is used as reference to check whenever we have
#' empty space between segments and get the ratio_median of this missing portion.
#'

CorrectLeftovers <- function(segments, granges_obj) {
  c = 1

  while (c < dim(segments)[1]) {
    # If in the same chromosome arm
    if (segments[c,3] == segments[c+1,3]) {
      # If space between them is > 1 (i.e. missing data)
      if (segments[c+1,4] - segments[c,5] > 1) {
        gr = GenomicRanges::GRanges(seqnames = c(segments[c,2]),
                                    ranges = IRanges::IRanges(start=c(segments[c,5]), end=c(segments[c+1,4])),
                                    strand = c("*"))
        subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

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

MergeSegmentsTwo <- function(granges_obj, small_segment, large_segments, j, below_thr, end_of_file) {
  N_large = dim(large_segments[1])
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)

  if (below_thr == TRUE) {
    gr = GenomicRanges::GRanges(seqnames = c(large_segment[2]),
                                ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                                strand = c("*"))
    subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

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
    gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(large_segment[4])),
                                strand = c("*"))
    subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

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
                               small_segment[8]))
    } else {
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
    }
  }
  large_segments
}

#' Handles merging step for the fourth case of FinalizeSmallSegments.
#'
#' @description
#' Pretty much the same as in MergeSegmentsTwo. Both helpers difference in the genomic data we search for
#' through GRanges. as well as which segment (large or small) we add into the `large_segments` frame. For example,
#' start position of the large segment? Or the small segment?
#'
#' @param granges_obj A GRanges object: is used as reference to obtain true genomic data from specified regions.
#' @param small_segment An array: the given small segment that might be merged into the `large_segments` frame.
#' @param large_segments A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
#' small segments into it, smaller segments are progressively added.
#' @param j An integer/index: the position of the large segment we're comparing with small_segment, within `large_segments`.
#' @param below_thr A boolean: whether large_segment and small_segment ratio_median difference is below the threshold
#' @param end_of_file A boolean: whether there's any large segments after `j`.
#'

MergeSegmentsFour <- function(granges_obj, small_segment, large_segments, j, below_thr, end_of_file) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)

  if (below_thr == TRUE) {
    gr = GenomicRanges::GRanges(seqnames = c(large_segment[2]),
                                ranges = IRanges::IRanges(start=c(large_segment[4]), end=c(small_segment[5])),
                                strand = c("*"))
    granges_obj = IRanges::subsetByOverlaps(granges_obj, gr)

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
    gr = GenomicRanges::GRanges(seqnames = c(small_segment[2]),
                                ranges = IRanges::IRanges(start=c(large_segment[5]), end=c(small_segment[5])),
                                strand = c("*"))
    granges_obj = IRanges::subsetByOverlaps(granges_obj, gr)

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

#' Deals with merging large segment to small segment.
#'
#' @description
#' Procedure is overall similar to the previous cases. For case five, we have two options when merging:
#' - Merge large and small
#' - Merge small and next large
#' This helper deals with the first option.
#' Note that there is no 'below_thr' parameter for these two helpers, since we already know we're going to merge
#' when they're called in FinalizeSmallSegments.
#'
#' @param granges_obj A GRanges object: is used as reference to obtain true genomic data from specified regions.
#' @param small_segment An array: the given small segment that might be merged into the `large_segments` frame.
#' @param large_segments A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
#' small segments into it, smaller segments are progressively added.
#' @param j An integer/index: the position of the large segment we're comparing with small_segment, within `large_segments`.
#' @param end_of_file A boolean: whether there's any large segments after `j`.
#'

MergeSegmentsFiveLargeSmall <- function(granges_obj, small_segment, large_segments, j, end_of_file) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)

  gr = GenomicRanges::GRanges(seqnames = c(large_segment[2]),
                              ranges = IRanges::IRanges(start=c(large_segment[4]), end=c(small_segment[5])),
                              strand = c("*"))
  granges_obj = IRanges::subsetByOverlaps(granges_obj, gr)

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
  large_segments
}

#' Deals with merging small segment to next_large segment.
#'
#' @description
#' Deals with the second option.
#'
#' @param granges_obj A GRanges object: is used as reference to obtain true genomic data from specified regions.
#' @param small_segment An array: the given small segment that might be merged into the `large_segments` frame.
#' @param large_segments A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
#' small segments into it, smaller segments are progressively added.
#' @param j An integer/index: the position of the large segment we're comparing with small_segment, within `large_segments`.
#' @param end_of_file A boolean: whether there's any large segments after `j`.
#'

MergeSegmentsFiveSmallNextLarge <- function(granges_obj, small_segment, large_segments, j, end_of_file) {
  N_large = dim(large_segments)[1]
  large_segment = as.numeric(large_segments[j,])
  small_segment = as.numeric(small_segment)

  gr = GenomicRanges::GRanges(seqnames = c(large_segment[2]),
                              ranges = IRanges::IRanges(start=c(small_segment[4]), end=c(large_segment[5])),
                              strand = c("*"))
  granges_obj = IRanges::subsetByOverlaps(granges_obj, gr)

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
