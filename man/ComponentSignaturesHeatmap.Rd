% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatureplotting.R
\name{ComponentSignaturesHeatmap}
\alias{ComponentSignaturesHeatmap}
\title{Create Heatmap of Component loadings per Signature}
\usage{
ComponentSignaturesHeatmap(
  signatures,
  transpose = FALSE,
  colour_scheme = "C",
  colour_dir = 1,
  cm_begin = 0
)
}
\arguments{
\item{signatures}{Dataframe. Expects a dataframe of signature exposures (rows) by samples (columns).}

\item{transpose}{(optional) Logical. If set to TRUE the function returns the order in which samples were plotted.}

\item{colour_scheme}{(optional) Character. Value passed to the \code{viridis::scale_fill_viridis} option parameter. ex. 'A' or 'B'}

\item{colour_dir}{(optional) Either 1 OR -1. Value passed to the \code{viridis::scale_fill_viridis} direction parameter.}

\item{cm_begin}{(optional) Float. Value within range \link{0,1}. Value passed to the \code{viridis::scale_fill_viridis} direction parameter.}
}
\value{
A ggplot2 object.
}
\description{
Converts Components-per-signature data to a heatmap, and returns the ggplot object.
Components and Signatures are sorted for display according to the order provided.
Heatmap is plotted using the viridis colour-schemes.
}
