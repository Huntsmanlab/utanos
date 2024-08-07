% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_helpers_small_segments.R
\name{MergeSegmentsFour}
\alias{MergeSegmentsFour}
\title{Handles merging step for the fourth case of FinalizeSmallSegments.}
\usage{
MergeSegmentsFour(
  granges_obj,
  small_segment,
  large_segments,
  j,
  below_thr,
  end_of_file
)
}
\arguments{
\item{granges_obj}{A GRanges object: is used as reference to obtain true genomic data from specified regions.}

\item{small_segment}{An array: the given small segment that might be merged into the \code{large_segments} frame.}

\item{large_segments}{A data frame: segment data. Initially contains only segments with size >= 3Mb, but as we merge
small segments into it, smaller segments are progressively added.}

\item{j}{An integer/index: the position of the large segment we're comparing with small_segment, within \code{large_segments}.}

\item{below_thr}{A boolean: whether large_segment and small_segment ratio_median difference is below the threshold}

\item{end_of_file}{A boolean: whether there's any large segments after \code{j}.}
}
\description{
Pretty much the same as in MergeSegmentsTwo. Both helpers difference in the genomic data we search for
through GRanges. as well as which segment (large or small) we add into the \code{large_segments} frame. For example,
start position of the large segment? Or the small segment?
}
