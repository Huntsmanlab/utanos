#' Determines threshold for ratio_median difference via KDE.
#' 
#' @description 
#' The procedure is mostly similar for both rounds, except for a few changes. Cleaning/transformation leading
#' up to the first round is vastly different than second round. In second round we pretty much just use the data frame
#' we've worked up to Graph 5, whereas first round has to take different steps.
#' 
#' Estimating the threshold goes like:
#' 1. Prep the segment data, depending on `second_round`. First round works with only the large segments. Second round works
#'    with all (large + small) segments initially, then only with large in step 4.
#' 2. Get the differences between ratio_medians for all pairs of segments.
#' 3. Run the simulations: more description found in RunThresholdSimulations docs, but essentially samples a batch of the differences, and
#'    finds the local minima/maxima of this batch, via Kernel Density Estimation (KDE). Repeat `num_simulations` times, and in the end
#'    we get a list of a bunch of local minima. The threshold is the average of these values.
#' 4. If `second_round` == TRUE, then few extra steps (estimates threshold but only with large segments).
#' 5. For QC purposes, the final threshold is = min(max(0.025, thr), 0.045).
#' 
#' @param granges_obj A GRanges object: genomic data to obtain for reference.
#' @param segments A data frame: segment data. For the first round, these segments have been gathered by ratio_median.
#' For the second round, the segments have already been re-inserted with small ones, and LGAs have already been called.
#' @param num_simulations An integer: the number of simulations to run to estimate the critical points.
#' @param second_round A boolean: which route to take for the algorithm.

FindThreshold <- function(granges_obj, segments, num_simulations=100000, second_round) {
  #### Prepping data depending on first or second round ####
  if (second_round == FALSE) {
    segments_3mb = PrepFirstRound(granges_obj=granges_obj,
                                  segments=segments)
  } else {
    segments_3mb = PrepSecondRound(segments=segments) # Not actually segments with size > 3mb, it's all of them (large and small)
  }
  
  #### Getting the differences between all pairs of ratio_median's ####
  N_3mb = dim(segments_3mb)[1]
  all_ratio_differences = c()
  for (i in 1:N_3mb) {
    v = i + 1
    for (j in v:N_3mb-1) {
      all_ratio_differences = c(all_ratio_differences, abs(segments_3mb[i,5] - segments_3mb[j, 5]))
    }
  }
  
  #### Simulations ####
  thr <- RunThresholdSimulations(num_simulations=num_simulations,
                                 all_ratio_differences=all_ratio_differences,
                                 second_round=second_round)
  
  #### If second round, few extra steps: works with large segments size > 3Mb ####
  # Large here is also defined as size > (Q3 - Q1)/2
  if (second_round == TRUE) {
    segments_3mb = segments_3mb[which(segments_3mb[,6] > 2999999),]
    segments_3mb = segments_3mb[which(segments_3mb[,6] > ((quantile(segments_3mb[,6])[4] - quantile(segments_3mb[,6])[2])/2)),]
    segments_3mb = segments_3mb[which(segments_3mb$chr != 23),]
    
    N_3mb = dim(segments_3mb)[1]
    all_ratio_differences = c()
    for (i in 1:N_3mb) {
      v = i + 1
      for (j in v:N_3mb-1) {
        all_ratio_differences = c(all_ratio_differences, abs(segments_3mb[i,5] - segments_3mb[j, 5]))
      }
    }
    
    # First maxima
    maxima = GetCriticalPoints(all_ratio_differences, 'max')[1]    # grabbing the first maxima (y-pos: density)
    first_local_maxima = density(all_ratio_differences)$x[maxima] # getting x-pos: a ratio_median difference, of this maxima
    
    # Getting first minima: specifically the first minima with x-pos that is after the first_local_maxima
    minima = GetCriticalPoints(all_ratio_differences, 'min')
    first_local_minima = density(all_ratio_differences)$x[minima][which(density(all_ratio_differences)$x[minima] > first_local_maxima)][1]
    
    # Getting second maxima
    maxima = GetCriticalPoints(all_ratio_differences, 'max')[2]
    second_local_maxima = density(all_ratio_differences)$x[maxima]
    
    if (is.na(second_local_maxima) == FALSE) {
      if (abs((first_local_minima - first_local_maxima) - (second_local_maxima - first_local_minima)) < 0.10) {
        if (first_local_minima < thr) {
          thr = first_local_minima
        }
      }
    }
  }
  
  #### Final thr is min(max(0.025, thr), 0.045) ####
  if (thr > 0.45) {
    thr = 0.45
  }
  
  if (thr < 0.025) {
    thr = 0.025
  }
  
  thr
}
