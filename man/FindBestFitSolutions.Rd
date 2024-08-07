% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copynumberscaling.R
\name{FindBestFitSolutions}
\alias{FindBestFitSolutions}
\title{Find best fit solutions for fitting copy numbers through grid-based search}
\usage{
FindBestFitSolutions(
  relative_copy_numbers,
  weights = NULL,
  min_ploidy = 1.5,
  max_ploidy = 5.5,
  ploidy_step = 0.01,
  min_cellularity = 0.2,
  max_cellularity = 1,
  cellularity_step = 0.01,
  distance_function = c("MAD", "RMSD"),
  distance_filter_scale_factor = 1.25,
  max_proportion_zero = 0.05,
  min_proportion_close_to_whole_number = 0.5,
  max_distance_from_whole_number = 0.15,
  solution_proximity_threshold = 5,
  keep_all = FALSE
)
}
\arguments{
\item{relative_copy_numbers}{a numeric vector containing relative copy
numbers, i.e. ratios of copy numbers to the average copy number.}

\item{weights}{a numeric vector of weights to apply to each copy number value
(should be same length as relative_copy_numbers)}

\item{min_ploidy, max_ploidy}{the range of ploidies.}

\item{ploidy_step}{the stepwise increment of ploidies along the grid.}

\item{min_cellularity, max_cellularity}{the range of cellularities.}

\item{cellularity_step}{the stepwise increment of cellularities along the
grid.}

\item{distance_function}{the distance function to use, either "MAD" for the
mean absolute difference or "RMSD" for the root mean square difference, where
differences are between the fitted absolute copy number values and the
nearest whole number.}

\item{distance_filter_scale_factor}{the distance threshold above which
solutions will be discarded as a multiple of the solution with the smallest
distance.}

\item{max_proportion_zero}{the maximum proportion of fitted absolute copy
number values in the zero copy number state.}

\item{min_proportion_close_to_whole_number}{the minimum proportion of fitted
absolute copy number values sufficiently close to a whole number.}

\item{max_distance_from_whole_number}{the maximum distance from a whole
number that a fitted absolute copy number can be to be considered
sufficiently close.}

\item{solution_proximity_threshold}{how close two solutions can be before one
will be filtered; reduces the number of best fit solutions where there are
many minima in close proximity.}

\item{keep_all}{set to \code{TRUE} to return all solutions but with
additional \code{best_fit} column to indicate which are the local minima that
are acceptable solutions (may be useful to avoid computing the distance grid
twice)}
}
\value{
the distance in the fitted absolute copy numbers to whole numbers.
}
\description{
Find best fit solutions for a grid-based search through the given ranges of
ploidies and cellularities.
Copied from rascal.
}
\examples{
data(copy_number)
copy_number <- copy_number[copy_number$sample_id == "X17222", ]

segments <- CopyNumberSegments(copy_number)

solutions <- FindBestFitSolutions(
  segments$copy_number, segments$weight,
  min_ploidy = 1.5, max_ploidy = 5.5,
  distance_function = "MAD")
}
