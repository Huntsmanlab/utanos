% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_helpers_small_segments.R
\name{MergeSegmentsTwo}
\alias{MergeSegmentsTwo}
\title{Helper for merging segments in the second case of FinalizeSmallSegments. Returns the updated \code{large_segments}
data frame with segments appropriately merged (or not... depending on parameters).}
\usage{
MergeSegmentsTwo(
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

\item{below_thr}{A boolean: whether large_segment and small_segment ratio_median difference is below the threhsold.}

\item{end_of_file}{A boolean: whether there's any large segments after \code{j}.}
}
\description{
Given the \code{large_segments} data, a large_segment \code{j}, and a small_segment \code{small_segment}, we merge the large segment with
the small_segment according to various criteria. First, if their ratio_median difference is below the threshold, then
we typically merge large and small. If not, then small segment is separately added into the \code{large_segments} frame. Also, where
exactly in the \code{large_segments} frame we add this new data depends on \code{end_of_file}. For example, if it is TRUE, then there are no
other large segments after \code{j}, so we row-bind all of large_segments to the new data. If FALSE, then we 'squeeze' the new data between
j-1 and j+1.
}
