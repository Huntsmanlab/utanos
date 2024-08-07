% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_small_segments.R
\name{SmallMissingChrArms}
\alias{SmallMissingChrArms}
\title{Deals with missing chromosome arm values in the given data frame of small segments.}
\usage{
SmallMissingChrArms(large_segments, small_segments)
}
\arguments{
\item{large_segments}{A data frame. Segments with size >= Mb.}

\item{small_segments}{A data frame. Segments with 3Mb >= size >= 0.1Mb.}
}
\description{
Pretty much as in LargeMissingChrArms, except it removes the small segments
that are not in the missing_chr_arm.
}
