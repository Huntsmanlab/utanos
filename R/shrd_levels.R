#' Calls some helpers to do something to gathered by ratio medians.... gotta check this
#' Also adds a Level and an Index column to the data frame
#'
#' @description
#' Adds couple of columns to the given segment dataframe + other processing things
#'
#' @param segments The segments gathered_by_ratio_medians
#' @param granges_obj A granges object containing the clean ratios file after CleanBamRatios()
#'

PrepForLevelsInitialization <- function(segments, granges_obj) {
  # SHRD recalculates the ratio_median, by finding the overlap with the original bin-wise ratios per segment and using their median, instead of the segmented value
  for (i in 1:nrow(segments)) {
    gr = GenomicRanges::GRanges(seqnames=c(segments[i,1]),
                                ranges=IRanges::IRanges(start=c(segments[i,3]),
                                                        end=c(segments[i,4])),
                                strand=c("*"))
    subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

    # Handles cases in which we have infinite negative values due to losses, i.e. ratios of 0.
    # We assume that the most extreme loss is 1 copy in 64 copies, as in a pure tumor bearing a single copy of a segment with an overall ploidy of 64n
    # So that log2(1/64) = -6
    new_ratio_medians <- as.vector(subsetGRobject$ratio)
    new_ratio_medians[which(new_ratio_medians == -Inf)] <- -6
    segments[i,5] = median(new_ratio_medians)
  }

  #### Adds index and level columns ####
  segments <- cbind(seq(1,dim(segments)[[1]]),
                    segments)
  colnames(segments)[1] <- c("index")

  segments <- cbind(segments,
                    rep(0, dim(segments)[1]))
  colnames(segments)[8] <- c("level")

  segments <- as.data.frame(segments)
  segments
}

#' Removes the short chromosome arms (21, 22, 24) & keeps segments with size >= 3Mb
#'
#' @description
#' Removes any segments that lie on a chromosome arm that:
#' - is greater than 23.6 (i.e chromosome arms 24 and beyond)
#' - is equal to 21 (i.e. first chromosome arm of chromosome 21)
#' - is equal to 22 (i.e. first chromosome arm of chromosome 22)
#' Then, updates sizes of the kept segments.
#'
#' @param segments The segments data frame, gathered by ratio medians.
#'

ExcludeShortArms <- function(segments) {
  #### Excludes the selected chromosome arms ####
  more_than_arm_24 <- which(segments[,3] > 23.6)
  if (length(more_than_arm_24) > 0) {
    segments <- segments[-more_than_arm_24,]
  }

  first_arm_chr21 <- which(segments[,3] == 21)
  if (length(first_arm_chr21) > 0) {
    segments <- segments[-first_arm_chr21,]
  }

  first_arm_chr22 <- which(segments[,3] == 22)
  if (length(first_arm_chr22) > 0) {
    segments <- segments[-first_arm_chr22,]
  }

  #### Updating sizes ####
  segments[,7] = segments[,5] - segments[,4] + 1
  segments
}

#' Returns segments with a size >= 3Mb.
#'
#' @description
#' Returns a data frame with segments whose size is larger than 3 mega bases.
#'
#' @param segments
#'

GetLargeSegments <- function(segments) {
  result = segments[which(segments[,7] > 2999999),]
  result
}


#' Gather big segments and assign levels iteratively
#'
#' @description
#' Iterates through the given segments_copy data frame. Does the following at each iteration:
#' 1. Finds the largest segment (by size). Then, finds all the segments 'close' to it: closeness
#'    here means that that ratio_median difference between them is less than threshold (`thr`).
#' 2. Once we identify these segments that are 'close' to each other (called closest_indices),
#'    we iterate through them, and set their level in `segments`.
#' 3. Finally, remove these closest_indices segments from `segments_copy`.
#'
#' @param segments  A data frame containing segment data.
#' @param thr A float: the threshold previously estimated via KDE. Used to determine whether
#' large segments are merged.
#'

AssignLevels <- function(segments, thr) {
  segments_copy <- segments

  #### Iterate through the rows in segments_copy, while there are still any segments left ####
  level = 1
  while (dim(segments_copy)[1] > 1) {
    # Starting from the largest segment, we find segments that are the closest to it.
    # We determine if a segment is closest to it if its ratio_median difference is less than the threshold.
    closest_indices = FindBigSegmentsIndices(segments_copy = segments_copy,
                                             thr = thr)
    # Now we assign `level` to all segments in `segments`.
    # We iterate through the segments in `closest_indices`, find the respective
    # row they belong to in `segments`, and set their level. Also remove them from segments_copy
    # so that we don't keep finding them at every iteration of the while loop.
    n = length(closest_indices)
    for (i in 1:n) {
      # Finding  segment i from `closest_indices` in `segments`
      segments[which(segments[,1] == closest_indices[i]), 8] = level

      # Getting the rest of the segments and keeping them in `segments_copy`.
      segments_copy <- segments_copy[which(!(segments_copy[,1] == closest_indices[i])), ]
    }
    level = level + 1
  }
  segments
}

#' Merges segments by level.
#'
#' @description
#' Iterates through given `segments` and merges by level or by chromosome arm. In other words,
#' we are keeping track of whenever we change chromosome arm or level.
#'
#' @param segments A data frame containing segment data. For the first pass, these segments
#' are those with size >= 3mb. Their last column should be their level, and should be non-zero
#' (i.e. already determined)
#' @param granges_obj A GRanges object: used to obtain genomic data from regions we're merging.
#'

GatherSegmentsByLevels <- function(segments, granges_obj) {
  #### Iterate through rows in segments ####
  N = dim(segments)[1]
  result = matrix(0, ncol=8, nrow=N)
  i = 1
  c = 1

  while (i < N+1) {
    # Last row: we keep the data from segment i as it is, since there is nothing
    # else to merge it to
    if (i == N) {
      result[c,] = c(segments[i,1], # index
                     segments[i,2], # chr_n
                     segments[i,3], # chr_arm
                     segments[i,4], # start
                     segments[i,5], # end
                     segments[i,6], # ratio_median
                     segments[i,5] - segments[i,4] + 1, # size
                     segments[i,8]) # level
      i = i + 1
      c = c + 1
    } else {
      # Else: not last row
      # If segments i and i+1 are in same chromosome arm
      if (segments[i,3] == segments[i+1,3]) {
        # If segments i and i+1 are in the same level, then we are going to iterate
        # through the segments starting at i + 1 (so the next one), until we no longer have same level or same arm.
        # As we iterate, we save the data from the segments we are merging.
        if (segments[i,8] == segments[i+1,8]) {
          n = 1
          vector_of_chr_numbers = c(segments[i,2])
          vector_of_seg_starts = c(segments[i,4])
          vector_of_seg_ends = c(segments[i,5])
          vector_strands = c("*")

          # While we are on the same level and same chr_arm, merge chr_n/start/end
          while (segments[i,8] == segments[i+n,8] && segments[i,3] == segments[i+n,3]) {
            vector_of_chr_numbers = c(vector_of_chr_numbers, segments[i+n,2])
            vector_of_seg_starts = c(vector_of_seg_starts, segments[i+n,4])
            vector_of_seg_ends = c(vector_of_seg_ends, segments[i+n,5])
            vector_strands = c(vector_strands, "*")

            n = n + 1

            if (i + n == N+1) {
              break
            }
          }

          # Merge the segments: GRanges helps us get ratio_median of the merge
          genomic_ranges = GenomicRanges::GRanges(seqnames = vector_of_chr_numbers,
                                                  ranges = IRanges::IRanges(start = vector_of_seg_starts, end = vector_of_seg_ends),
                                                  strand = vector_strands)
          subsetGRobject = IRanges::subsetByOverlaps(granges_obj, genomic_ranges)

          # Adding the merged segment into the result matrix
          result[c,] = c(segments[i,1],
                         segments[i,2],
                         segments[i,3],
                         segments[i,4],
                         segments[i+n-1,5], # i+n-1 is where we stopped being in same level/arm
                         median(subsetGRobject$ratio),
                         segments[i+n-1,5] - segments[i,4] + 1,
                         segments[i,8])
          i = i + n
          c = c + 1
        } else {
          # Else: not same level
          result[c,] = c(segments[i,1], # index
                         segments[i,2], # chr_n
                         segments[i,3], # chr_arm
                         segments[i,4], # start
                         segments[i,5], # end
                         segments[i,6], # ratio_median
                         segments[i,5] - segments[i,4] + 1, # size
                         segments[i,8]) # level
          i = i + 1
          c = c + 1
        }
      } else {
        # Else: not same chromosome arm
        result[c,] = c(segments[i,1], # index
                       segments[i,2], # chr_n
                       segments[i,3], # chr_arm
                       segments[i,4], # start
                       segments[i,5], # end
                       segments[i,6], # ratio_median
                       segments[i,5] - segments[i,4] + 1, # size
                       segments[i,8]) # level
        i = i + 1
        c = c + 1
      }
    }
  }

  #### Preparing result matrix for output ####
  result = subset(result, result[,1] != 0)
  rownames(result) <- NULL
  colnames(result) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")
  result
}


