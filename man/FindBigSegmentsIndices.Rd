% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_helpers_levels.R
\name{FindBigSegmentsIndices}
\alias{FindBigSegmentsIndices}
\title{Finds the indices of big segments that are close ratio_median-wise}
\usage{
FindBigSegmentsIndices(segments_copy, thr)
}
\arguments{
\item{segments_copy}{A data frame containing segment data. For the first pass, these segments
are those with size >= 3mb.}

\item{thr}{A float: the threshold previously estimated via KDE. Used to determine whether
big segments are merged.}
}
\description{
Starting from the largest segment, we find big segments if their ratio_median
difference is less than the the given threshold. We save their indices.
}
