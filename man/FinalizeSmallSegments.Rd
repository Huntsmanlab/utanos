% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_helpers_small_segments.R
\name{FinalizeSmallSegments}
\alias{FinalizeSmallSegments}
\title{Re-inserts small segments in between the large ones.}
\usage{
FinalizeSmallSegments(large_segments, small_segments, threshold, granges_obj)
}
\arguments{
\item{large_segments}{A data frame. Segments with size >= 3Mb.}

\item{small_segments}{A data frame. Segments with 3Mb >= size >= 0.1Mb.}

\item{threshold}{A float: the estimated threshold for ratio_median difference via KDE. Used to determine
whether we insert the small segment or not.}

\item{granges_obj}{A GRanges object: is used as reference to check whenever we have
an overlap of segments and get the ratio_median of this overlap.}
}
\description{
We add small segments from \code{small_segments} to \code{large_segments} in 6 cases,
each depending on where the small segment fits in between the large ones.
We merge segments if the ratio difference between the small and the large are <= threshold.
Otherwise, we keep small segment as its own segment in \code{large_segments}.

Visually, these are the 6 cases:
\enumerate{
\item |--small-|    |-----large----|
\item |----- previous large -----|    |----small----|
|---------- large segment----------|
\item |------------- large ---------------|
|----small----|
\item |---------large----------|        |------------- next large--------|
|--------small--------|
\item |-------large--------|                    |--------------next large-------------|
|----small----|
\item chr N                                          |#| chr N+1
|-------- large --------|                         |#| |-------- next_large --------|
|------ small ------|    |#|
}
}
