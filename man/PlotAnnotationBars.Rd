% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{PlotAnnotationBars}
\alias{PlotAnnotationBars}
\title{Create annotation bars}
\usage{
PlotAnnotationBars(ann_df, cols, colors = LETTERS[1:length(cols)], order)
}
\arguments{
\item{ann_df}{A dataframe. Must contain at least \code{sample_id} and all columns in \code{cols} \cr}

\item{cols}{Character vector. Names of categorical columns to include as annotations, i.e. one column = one annotation bar \cr}

\item{colors}{Character vector. Letters A-H indicating the viridis color option to be used for each column. Must have same length as \code{cols} \cr
Example if \code{cols} lists two columns: \cr
colors <- c('A', 'B')}

\item{order}{vector containing all values of sample_id in desired order, controls order of tiles. \cr}
}
\value{
A list of ggplot2 objects
}
\description{
Make a list of coloured annotation bars for each sample provided and their corresponding categories
}
