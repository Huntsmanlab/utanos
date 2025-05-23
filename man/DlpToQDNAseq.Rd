% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DlpToQDNAseq.R
\name{DlpToQDNAseq}
\alias{DlpToQDNAseq}
\title{Build a QDNAseq object from the DLP output files}
\usage{
DlpToQDNAseq(
  input_path,
  bin_size = 5e+05,
  genome = "hg19",
  Xincluded = FALSE,
  getfiltered = FALSE
)
}
\arguments{
\item{input_path}{(string) The path to the output files generated by DLP. A directory.}

\item{bin_size}{(integer) The bin size used during copy number calling in base pairs. (default: 500000)}

\item{Xincluded}{(logical) Whether or not the X chromosome was included in this analysis. (default: FALSE)}

\item{genome_used}{(string) The reference genome used for alignment. (default: hg19)}

\item{getfilter}{(logical) Grab the filtered read/segment files from the DLP output directories.}
}
\value{
A QDNAseq object.
}
\description{
In order to do 1-to-1 comparisons between QDNAseq and DLP outputs, it is useful to have outputs in the same format.
This function creates a QDNAseq object from the DLP output files and returns it.
}
