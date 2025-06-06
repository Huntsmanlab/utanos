% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainsignaturefunctions.R
\name{ExtractRelativeCopyNumberFeatures}
\alias{ExtractRelativeCopyNumberFeatures}
\title{Extract Relative Copy-Number Features}
\usage{
ExtractRelativeCopyNumberFeatures(
  CN_data,
  genome,
  cores = 1,
  log_features = FALSE,
  extra_features = FALSE
)
}
\arguments{
\item{CN_data}{List of datafames or S4 QDNAseq object. Segmented relative copy-number data for 1 or more samples.
If input is a list of dataframes, columns should be: \cr
\enumerate{
\item chromosome
\item start
\item end
\item segVal
}}

\item{genome}{Character string. The reference genome used for alignment. \cr
Options: 'hg19', 'hg38'}

\item{cores}{Integer. The number of cores to use for parallel processing.}

\item{log_features}{FALSE or char vector. If a vector of feature names is provided, take the log1p of these extracted CN-features.}

\item{extra_features}{Logical. If TRUE, extracts CN-feature data for two more features: nc50, and cdist.}
}
\value{
A list. Each list element contains feature data for a single feature.
}
\description{
Extract genome-wide copy-number features from either a list of dataframes or a QDNAseq S4 object for 1 or more samples.
This function is intended to be run on relative CN data.
}
\details{
This function is identical to the absolute calling equivalent other than for three features.
The \code{osCN}, \code{changepoint}, and \code{copynumber} features require slightly different modelling at the relative scale.

The extracted copy-number features are:
\enumerate{
\item Breakpoint count per 10MB - \code{bp10MB}
\item Copy-number value of each segment - \code{copynumber}
\item Copy-number difference between adjacent segments - \code{changepoint}
\item Breakpoint count per chromosome arm - \code{bpchrarm}
\item Lengths of oscillating CN segment chains - \code{osCN}
\item Size of copy-number segments in base-pairs - \code{segsize}
}

Extra features: \cr
7. Minimum number of chromosomes (a count) needed to account for 50\% of CN changes in a sample - \code{nc50} \cr
8. Distance in base pairs of each breakpoint to the centromere - \code{cdist}
}
