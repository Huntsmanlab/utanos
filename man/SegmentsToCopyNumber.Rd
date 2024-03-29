% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{SegmentsToCopyNumber}
\alias{SegmentsToCopyNumber}
\title{Transforms segment tables into per-bin copy-number tables.}
\usage{
SegmentsToCopyNumber(segs, bin_size, genome = "hg19", Xincluded = FALSE)
}
\arguments{
\item{segs}{A list of dataframes, where each dataframe is a dataframe of the segmented copy-numbers.}

\item{bin_size}{A natural number. The binsize used in the copy number object, ex. 30 -> 30kb, 100 -> 100kb.}

\item{genome}{A string. Refers to the reference genome, common reference genomes are: 'hg19', 'mm10', or 'hg38'.}

\item{Xincluded}{A boolean. Sex chromosomes as a Boolean (e.g. FALSE if not included).}
}
\value{
A dataframe of per-bin copy numbers (bounded for each sample).
}
\description{
SegmentsToCopyNumber() transforms segment tables into per-bin copy-number tables.
It is an expansion of the calls into per-bin style, where the bin size is user defined.
It is the inverse of the CopyNumberSegments function.
}
