% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrd_gather_segments_ratiomedian.R
\name{RemoveSpuriousRegions}
\alias{RemoveSpuriousRegions}
\title{Removes specific genomic regions (centromeres and telomeres) of all chromosomes}
\usage{
RemoveSpuriousRegions(bam_ratios_frame, include_chr_X)
}
\arguments{
\item{bam_ratios_frame}{A data frame: a cleaned up version of the original bam_ratios.txt file.}

\item{include_chr_X}{A boolean: whether we search for chromosome X or not.}
}
\description{
First, drops the 'Feature' and the 'Ratio' columns.
Spurious regions are unstable regions that are difficult to map, so we remove them. The \code{hg19_segments.R} file
contains positions of these regions. We find them in the \code{bam_ratios_frame} and, if found, they are removed.
}
