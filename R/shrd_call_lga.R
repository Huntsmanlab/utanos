#' Preps `segments` for LGA call.
#'
#' @description
#' Does many things, details found within function body. First, assigns segments
#' that meet certain criteria an index of 0 (there's 4 cases + edge cases for assigning an
#' index of 0, inside the inner for-loop.) Then, removes the segments that were assigned this
#' index of 0. They're considered "small breaks".
#'
#' @param threshold A float. Threshold in ratio_median difference previously estimated.
#' @param segments A data frame: segment data that's (ideally) already been processed (merging segments,
#' small ones have been reinserted, and so on)
#' @param granges_obj A GRanges object.
#'
#' @export
BreakSmoothToLGA <- function(threshold, segments, granges_obj) {
  segments <- GetSegmentID(threshold=threshold,
                          segments=segments)
  segments <- ShrinkReprTMP(segments=segments,
                           granges_obj=granges_obj)

  small_breaks = 0
  pass = 1
  segments_lt3mb_and_pass_gt10000 = TRUE
  while (segments_lt3mb_and_pass_gt10000) {
    segments[,1] <- seq(1,dim(segments)[1])
    segments_lt3mb <- which(round((segments[,5] - segments[,4])/10^6,1) < 3)

    if (length(segments_lt3mb) > 0 & pass < 10000) {
      # Ordering segments (their indices, really) in segments_lt3mb by their size
      ordered_segments_lt3mb <- segments_lt3mb[order((segments[segments_lt3mb,5] - segments[segments_lt3mb,4])/10^6)]
      pass = pass + 1

      # Iterating through these small segments
      for (i in 1:length(ordered_segments_lt3mb)) {
        # If size of segment (i) is 1
        if (ordered_segments_lt3mb[i] == 1) {
          # If next segment (i+1) is in same chromosome as (i).
          if (segments[ordered_segments_lt3mb[i]+1, 3] == segments[ordered_segments_lt3mb[i], 3]) {
            # If segment (i+1)'s index is not 0
            if (segments[ordered_segments_lt3mb[i]+1, 1] != 0) {
              segments[ordered_segments_lt3mb[i], 1] <- 0
            }
          }
        # else: segment size of (i) is not 1
        } else {
          # If last segment
          if (ordered_segments_lt3mb[i] == dim(segments)[1]) {
            # If last segment and previous  segment in same chromosome arm AND
            # previous segment's index is not 0
            if (segments[ordered_segments_lt3mb[i]-1, 3] == segments[ordered_segments_lt3mb[i], 3] & segments[ordered_segments_lt3mb[i]-1, 1] != 0) {
                segments[ordered_segments_lt3mb[i], 1] <- 0
            }
          # Not last segment
          } else {
            # We have |-previous-| |-segment-| |-next-|
            previous_next_same_arm = segments[ordered_segments_lt3mb[i]-1, 3] == segments[ordered_segments_lt3mb[i]+1, 3]
            previous_segment_same_arm = segments[ordered_segments_lt3mb[i]-1, 3] == segments[ordered_segments_lt3mb[i], 3]
            segment_next_same_arm = segments[ordered_segments_lt3mb[i]+1, 3] == segments[ordered_segments_lt3mb[i], 3]

            # Case 1:
            # If Previous and Next in same chromosome arm AND
            # Previous index is not 0 AND
            # Next's index is not 0
            # Then set Segment's index to 0
            if (previous_next_same_arm &
                (segments[ordered_segments_lt3mb[i]-1, 1] != 0) &
                (segments[ordered_segments_lt3mb[i]+1, 1] != 0)) {
              segments[ordered_segments_lt3mb[i], 1] <- 0
            }

            # Case 2:
            # If Previous and Segment in same chromosome arm AND
            # Segment and Next are in different chromosome arms AND
            # Previous index is not 0
            # Then set Segment's index to 0
            if (previous_segment_same_arm &
                (segments[ordered_segments_lt3mb[i]+1, 3] != segments[ordered_segments_lt3mb[i], 3]) &
                (segments[ordered_segments_lt3mb[i]-1, 1] != 0)) {
              segments[ordered_segments_lt3mb[i], 1] <- 0
            }

            # Case 3:
            # If Segment and Next in same chromosome arm AND
            # Previous and Segment in different chromosome arms AND
            # Next index is not 0
            # Then set Segment's index to 0
            if (segment_next_same_arm &
                (segments[ordered_segments_lt3mb[i]-1, 3] != segments[ordered_segments_lt3mb[i], 3]) &
                (segments[ordered_segments_lt3mb[i]+1, 1] != 0)) {
              segments[ordered_segments_lt3mb[i], 1] <- 0
            }

            # Case 4:
            # If Segment and Next in different chromosome arms AND
            # Previous and Segment in different chromosome arms
            # Then set Segment's index to 0
            if (!segment_next_same_arm &
                (segments[ordered_segments_lt3mb[i]-1, 3] != segments[ordered_segments_lt3mb[i], 3])) {
              segments[ordered_segments_lt3mb[i], 1] <- 0
            }
          }
        }
      }

      zero_indexed <- which(segments[,1] == 0)
      if (length(zero_indexed) > 0) {
        segments <- segments[-zero_indexed,]
      }
      small_breaks = small_breaks + length(zero_indexed)

      for (i in 1:(dim(segments)[1]-1)) {
        # If space between Segment and Next is < 3mb
        if (round((segments[i+1,4] - segments[i,5])/10^6,1) < 3) {
          below_thr <- RatioDiffBelowThreshold(large_segment=segments[i,],
                                               small_segment=segments[i+1,],
                                               threshold=threshold)
          # If in same chromosome arm AND
          # ratio_median difference is below threshold
          if (segments[i,3] == segments[i+1,3] & below_thr) {
            segments[i+1,4] <- segments[i,4]

            gr = GenomicRanges::GRanges(seqnames = c(segments[i+1,2]),
                                        ranges = IRanges::IRanges(start=c(segments[i+1,4]), end=c(segments[i+1,5])),
                                        strand = c("*"))
            subsetGRobject = IRanges::subsetByOverlaps(granges_obj, gr)

            segments[i+1,6] <- median(subsetGRobject$ratio)
            segments[i,1] <- 0
          }
        }
      }

      zero_indexed <- which(segments[,1] == 0)
      if (length(zero_indexed) > 0) {
        segments <- segments[-zero_indexed,]
      }
    }
    else {
      segments_lt3mb_and_pass_gt10000 <- FALSE
    }
  }
  segments <- GetSegmentID(threshold=threshold,
                           segments=segments)
  segments <- ShrinkReprTMP(segments=segments,
                           granges_obj=granges_obj)
  segments[1,1] <- small_breaks
  segments
}

#' Final segmentation before LGA call
#'
#' @description
#' Merges the initial `bam_ratios_frame` by segments in the same chromosome,
#' and then assigns start/end pos to segments in `segments`, whenever
#' we change Chr_N in the latter.
#'
#' @param segments A data frame of segments with size > 3Mb. Ideally this
#' frame has already gone through merging, small segments reinsertion, etc.
#' @param bam_ratios_frame A data frame. For first pass: the cleaned-up version of the initial
#' ratio file tsv. For second pass: this clean-up segments but merged by ratio_median.
#' @param granges_obj A GRanges object: used to obtain genomic data from regions we're merging.
#' @param second_round A boolean: FALSE if first pass, TRUE if second pass.
#' @export
GetSegmentationBeforeLGACall <- function(segments, bam_ratios_frame, granges_obj, second_round) {
  if (second_round == FALSE) {
    # Setting up frame
    bam_ratios_frame = bam_ratios_frame[,-1]
    bam_ratios_frame = bam_ratios_frame[,-5]
    colnames(bam_ratios_frame) <- c("chr", "start", "end", "readcount")

    # Setting variables that are specific to the first pass
    attach(bam_ratios_frame)
    merged_bam = merge(aggregate(start ~ chr, bam_ratios_frame, min),
                       aggregate(end ~ chr, bam_ratios_frame, max))[seq(from=1, to=23, by=1),]
    detach(bam_ratios_frame)

    while_loop_index = 2 # chromosome number
    num_chr = 23
  } else {
    # Setting up frame
    colnames(bam_ratios_frame) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")

    # Setting variables that are specific to the second pass
    attach(bam_ratios_frame)
    merged_bam = merge(aggregate(start ~ chr_arm, bam_ratios_frame, min),
                       aggregate(end ~ chr_arm, bam_ratios_frame, max))[seq(from=1, to=41, by=1),]
    detach(bam_ratios_frame)

    while_loop_index = 3 # chromosome arm
    num_chr = 41
  }

  # First segment
  segments[1,4] = merged_bam[1,2]

  # We're going to iterate through the large segments.
  # If large_segment and next_large are in different Chr_N, then:
  #     1) large_segment's end pos is n's end pos
  #     2) next_large's start pos is n+1's start pos
  N_large = dim(segments)[1]
  c = 1
  n = 1
  while (c < N_large) {
    if (segments[c,while_loop_index] != segments[c+1,while_loop_index]) {
      segments[c,5] = merged_bam[n,3]
      n = n + 1
      segments[c+1,4] = merged_bam[n,2]
    }
    c = c + 1
  }

  segments[N_large,5] = merged_bam[num_chr,3]
  segments[,7] = segments[,5] - segments[,4] + 1

  # Again iterate and correct segment ratio medians with values from original cleaned data
  c = 1
  while(c < N_large + 1) {
    genomic_ranges <- GenomicRanges::GRanges(seqnames = c(segments[c, 2]),
                                             ranges = IRanges::IRanges(start = c(segments[c, 4]), end = c(segments[c, 5])),
                                             strand = c("*"))
    subsetGRobject = IRanges::subsetByOverlaps(granges_obj, genomic_ranges)
    segments[c, 6] = median(subsetGRobject$ratio)

    c = c + 1
  }


  segments
}

#' Determines number of large genomic alterations (LGAs) for sizes 3:11
#'
#' @description
#' Uses DetermineNumberOfLGAs to determine the number of LGAs for sizes 3 to 11 (Mb).
#'
#' @param threshold A float. Threshold in ratio_median difference previously estimated.
#' @param segments A data frame: segment data that's (ideally) already been processed (merging segments,
#' small ones have been reinserted, and so on)
#'
#' @export

CallLGA <- function(threshold, segments) {
  result <- as.data.frame(matrix(0, ncol=2, nrow=9))
  result[,1] <- c(3:11)
  colnames(result) <- c("Size_LGA", "Number_LGA")

  segments = segments[which(segments$chr != 23),]
  segments = as.matrix(segments)

  for (i in (3:11)) {
    segments_with_LGA <- DetermineNumberOfLGAs(threshold=threshold, size_lga=i, segments=segments)
    result[i-2,2] = sum(segments_with_LGA)
  }
  result
}

#' Returns the segments that are LGAs of the given size
#'
#' @description
#' With help of DetermineNumberOfLGAs, keeps segment data for LGA size of `size_lga`
#'
#' @param threshold A float. Threshold in ratio_median difference previously estimated.
#' @param size_lga An integer: size of the LGA to look for.
#' @param segments A data frame: segment data
#'
#' @export
GetLGAOfSize <- function(threshold, size_lga, segments) {
  segments_with_LGA <- DetermineNumberOfLGAs(threshold=threshold,
                                             size_lga=size_lga,
                                             segments=segments)
  result = cbind(segments, segments_with_LGA)
  colnames(result) <- c("index", "chr","chr_arm", "start", "end", "ratio_median", "size", "level", "WC")

  N = dim(result)[1]
  i = 1
  # Iterate through the segments, and remove the ones that aren't an LGA (their index in segments_with_LGA are not 1)
  # and whose next segment is also not an LGA.
  # Essentially keeping the segments that are LGAs (i.e their index in segments_with_LGA is 1)
  while (i + 1 < N) {
    N = dim(result)[1]
    if (result[i,9] == 0 & result[i+1,9] != 1) {
      result = result[-i,]
    } else {
      i = i + 1
    }
  }

  N_new = dim(result)[1]
  if (is.null(N_new) == FALSE) {
    # If last segment is not an LGA
    if (result[N_new,9] == 0) {
      # While the last segment is not an LGA, remove the last segment
      # Subtract 1 from N_new, and check the next last segment
      while (result[N_new,9] == 0) {
        result = result[-N_new,]
        N_new = N_new - 1

        # Only iterate if we have more than one non-LGA segment left
        if (N_new == 0) {
          break
        }
      }
    }
  }
  result
}

#' Returns the HRD status of a sample
#'
#' @description
#' Using the output of CallLGA, we determine the HRD status of a sample based on the number of LGAs >= 10 Mb in size.
#'
#' @param lga_calls The table of LGA counts by size, as returned by CallLGA.
#'
#' @returns The HRD status of the sample: TRUE if the number of LGAs >= 10 Mb exceeds 20; FALSE otherwise.
#'
GetHRDStatus <- function(lga_calls) {
  n_lga <- lga_calls[which(lga_calls$Size_LGA == 10), 2]

  if (n_lga >= 20) {
    hrd_status <- TRUE
  } else {
    hrd_status <- FALSE
  }

  return(hrd_status)
}
