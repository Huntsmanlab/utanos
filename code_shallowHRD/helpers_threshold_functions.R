#' Preps the segment data frame for use in the first round of threshold estimation.
#' 
#' @description
#' Only keeps segments whose size is >= 3Mb. First iterates through segments and assigns genomic
#' data from GRanges to this segment (specifically their ratio_median)
#' 
#' @param granges_obj A GRanges object: used to fine genomic data from the specified regions.
#' @param segmenta A data frame: segment data.
#' 
#' @export
PrepFirstRound <- function(granges_obj, segments) {
  segments = segments[which(segments$chr != 23),]
  for (i in 1:nrow(segments)) {
    gr = GRanges(seqnames=c(segments[i,1]),
                 ranges=IRanges(start=c(segments[i,3]), 
                                end=c(segments[i,4])),
                 strand=c("*"))
    subsetGRobject = subsetByOverlaps(granges_obj, gr)
    segments[i,5] = median(subsetGRobject$ratio)
  }
  
  #### Getting segments with size > 3Mb ####
  segments_3mb = segments[which(segments[,6] > 2999999),]
  segments_3mb = data.matrix(segments_3mb)
  
  segments_3mb
}

#' Preps the segment data drame for use in the second round of threshold estimation.
#'
#' @description 
#' Removes columns 1, 7, and 8 from `segments`, also readjusts segment sizes.
#'
#' @param segments A data frame: segment data
PrepSecondRound <- function(segments) {
  segments_3mb = segments[,-8]
  segments_3mb = segments[,-7]
  segments_3mb = segments[,-1]
  
  segments_3mb$size <- segments_3mb$end - segments_3mb$start + 1
  segments_3mb = segments_3mb[which(segments_3mb$chr != 23),]
  
  segments_3mb
}

#' Runs simulations to obtain a list of possible thresholds from the segment ratio differences
#'
#' @description 
#' At each simulation:
#' - Takes a random number of samples from all_ratio_differences (the batch)
#' - Passes this batch to GetCriticalPoints to estimate (via KDE) the batch's density minima/maxima
#' - Adds the minima as a possible threshold if it meets certain criteria. Most specifically, a local minima from a batch
#'   is considered a valid possible threshold if the difference between the distance between first maxima and the minima, and the
#'   second maxima and the minima, are 'close' to each other. If first round, this closeness is < 0.10 distance-wise. If second_round,
#'   this closeness is < 0.05 distance-wise. So, for the second round, we really care about the minima being very close to the two
#'   neighboring maxima. The point of these tolerances is that we do not get a 'low' and a 'flat' local minima, meaning that it is a ratio_median
#'   difference that's not very common. This would affect the rest of the algorithm as not many segments would meet this threshold.
#' Once simulations are done and we have a list of possible thresholds, we return the threshold (mean of list) 
#' 
#' @param num_simulations An integer: the number of simulations to run on all_ratio_differences. Set to 100,000
#' by default
#' @param all_ratio_differences An array with all the differences of ratio medians between all segment-pairs in the data set
#' @param second_round A boolean. False if first round of simulation, True if second round. This should be set to True once
#' the simulations have been ran before and we already have a possible THR that we want to re-calibrate with new differences
#' of ratio medians (i.e. with updated segments, such as post-merging). Set to False by default. 
#'
#' @export
#'
#'
RunThresholdSimulations <- function(num_simulations=100000, all_ratio_differences, second_round=FALSE) {
  list_of_possible_thresholds = c()
  tolerance = 0.10
  
  if (second_round == TRUE) {
    tolerance = 0.05
  }
  
  for (i in 1:num_simulations) {
    if (i %% 1000 == 0) {
      print(i)
    }
      
    size_of_batch = sample(2:length(all_ratio_differences), 1)
    sampled_batch = sample(all_ratio_differences, size_of_batch)
    
    # Getting first maxima
    maxima = GetCriticalPoints(sampled_batch, 'max')[1]    # grabbing the first maxima (y-pos: density)
    first_local_maxima = density(sampled_batch)$x[maxima] # getting x-pos: a ratio_median difference, of this maxima
    
    # Getting first minima: specifically the first minima with x-pos that is after the first_local_maxima
    minima = GetCriticalPoints(sampled_batch, 'min')
    first_local_minima = density(sampled_batch)$x[minima][which(density(sampled_batch)$x[minima] > first_local_maxima)][1]
    
    # Getting second maxima
    maxima = GetCriticalPoints(sampled_batch, 'max')[2]
    second_local_maxima = density(sampled_batch)$x[maxima]
    
    # We determine whether this first_local_minima is a possible Threshold if it is "close" to the two local_maxima.
    # Closeness is determined by tolerance: 0.010 if first round, 0.05 if second round
    if (is.na(second_local_maxima) == FALSE) {
      if (abs((first_local_minima - first_local_maxima) - (second_local_maxima - first_local_minima)) < tolerance) {
        list_of_possible_thresholds = c(list_of_possible_thresholds, first_local_minima)
      }
    }
  }
  
  # Checking if list is null
  if (is.null(list_of_possible_thresholds) == TRUE) {
    threshold <- NullListOfThresholds(all_ratio_differences)
  }
  
  # Getting first threshold:  mean of list_of_possible_thresholds
  if (is.na(mean(list_of_possible_thresholds)) == FALSE) {
    threshold = mean(list_of_possible_thresholds)
  }

  threshold
}

#' Finds critical points (minima or maxima) of the differences between
#' segments ratios.
#'
#' @description
#' GetCriticalPoints takes the estimated density values (via KDE) from a batch of samples
#' of segment ratio differences (S_i - S_j), a critical point (min/max), and returns a 
#' list of the densities where the critical point is found.
#' 
#' If crit_point == 'max', then densities > 0
#' If crit_point == 'min', then densities < 0
#' 
#' @param ratio_differences A Data Frame. The estimated density values for the given sampled batch
#' @param crit_point A string. The critical point to look for: either 'max' or 'min'.
#' 
#' @export
#'
GetCriticalPoints <- function(batch_ratio_differences, crit_point) {
  densities = density(batch_ratio_differences)$y
  
  if (crit_point == "min") {
    criticals <- diff(c(.Machine$integer.max, densities)) < 0L
  } else if (crit_point == "max") {
    criticals <- diff(c(-.Machine$integer.max, densities)) > 0L
  } else {
    stop("Crit_point must be either min or max.")
  }
  
  criticals <- cumsum(rle(criticals)$lengths)
  criticals <- criticals[seq.int(1L, length(criticals), 2L)]
  
  if (densities[[1]] == densities[[2]]) {
    criticals <- criticals[-1]
  }
  
  criticals
}

#' This function is ran in case that no possible thresholds were found. 
#'
#' @description 
#' NullListOfThresholds returns a threshold in the case where list_of_possible_threshold was null.
#' It finds the first local minima, but also finds the inflection points of the density function, and 
#' returns as threshold the min of these: min(local_minia, inflection_point)
#'
#' @param all_ratio_differenves An array: the ratio differences for all pairs of segments.
#' 
#' @export
NullListOfThresholds <- function(all_ratio_differences) {
  maxima = GetCriticalPoints(all_ratio_differences, 'max')[1]
  first_local_maxima = density(all_ratio_differences)$x[maxima]
  
  minima = GetCriticalPoints(all_ratio_differences, 'min')
  first_local_minima = density(all_ratio_differences)$x[minima][which(density(all_ratio_differences)$x[minima] > first_local_maxima)][1]
  
  h <- hns(all_ratio_differencesk, deriv.order=2)
  den <- kdde(all_ratio_differences, h=h, deriv.order=2)
  
  inflection_points <- c()
  for (i in 2:length(den$estimate)) {
    if (sign(den$estimate[i]) != sign(den$estimate[i-1])) {
      inflection_points <- c(inflection_points, i)
    }
  }
  
  xpos_all_inflection_points = den$x[inflection_points][which(den$x[inflection_points] > first_local_maxima)]
  threshold_changes_sign_second_derivative_2nd_value = sort(xpos_all_inflection_points)[1]
  
  threshold = min(first_local_minima, threshold_changes_sign_second_derivative_2nd_value)
  threshold
}
