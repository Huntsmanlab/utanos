#' Removes centromeres and telomeres from 22 chromosomes (23rd if specified)
#'
#' @description 
#' Removes each chromosome's centromere and telomeres. The positions of these centromeres/telomeres is specified in the
#' given arrays. So, the start positions of the 23 centromeres is in `centromere_starts`, the end positions of the 23 centromeres
#' is in `centromere_ends`, and so on. By 'telomere 2' we mean the second (bottom) telomere, while telomere 1 woud be the first (top) one.
#' Note that there are parameters for telomere 2 (`telomere_2_starts` and `telomere_2_ends`), but not for telomere 1. The reason for this
#' is that, in the algorithm, the start position and end position of telomere 1 is the same for ALl chromosomes (base 0 up to 10000).
#'
#' @param df A Data Frame with segment data. Requirements:
#' 1. Column 1 is chromosome number
#' 2. Column 2 is the start position of the segment
#' 3. column 3 is the end position of the segment. 
#' @param include_chr_X A boolean: True if you want its spurious regions removed. False otherwise. 
#' @param centromere_starts An array: contains the start positions of all 23 chromosome's centromeres IN ORDER from 1 to 23.
#' @param centromere_ends An array: contains the end positions of all 23 chromosome's centromeres IN ORDER from 1 to 23.
#' @param telomere_2_starts An array: contains the start positions of all 23 chromosome's second telomere IN ORDER from 1 to 23.
#' @param telomere_2_ends An array: contains the end positions of all 23 chromosome's second telomere IN ORDER from 1 to 23.
#'
#' @export
#' 
#' TODO: implement for hg38
#' TODO: can check that df's columns match what's needed (column 1= chr_n, etc.)
RemoveCentromereTelomeres <- function(df, include_chr_X, centromere_starts, centromere_ends, telomere_2_starts, telomere_2_ends) {
  if (include_chr_X == TRUE) {
    n_chr_to_check = 23 # so in for-loop  we look iterate 23 times: we check 23 chromosomes.
  } else {
    n_chr_to_check = 22 # otherwise skip it
  }
  for (i in 1:n_chr_to_check) {
    if (i == 17) {
      df <- df[which(!(df[,1] == i & df[,2] >= centromere_starts[i] & df[,3] <= centromere_ends[i])),]
    } else {
      df <- df[which(!(df[,1] == i & df[,2] >= 0 & df[,3] <= 10000)),] # telomere 1
      df <- df[which(!(df[,1] == i & df[,2] >= centromere_starts[i] & df[,3] <= centromere_ends[i])),] # centromere
      df <- df[which(!(df[,1] == i & df[,2] >= telomere_2_starts[i] & df[,3] <= telomere_2_ends[i])),] # telomere 2
    }
  }
  df
}

#' Adds a chromosome arm column to df at the respective positions for each chromosome.
#' 
#' @description
#' Determines the chromosome arm for each segment in each chromosome. This is determined by the positions given in `centromere_positions`.
#' Segments before centromere are in Chr_N (first arm), segments after centromere are in Chr_N.5 (second arm). '
#' 
#' The procedure goes like this. For each chromosome, the segments that are before its middle position have their chromosome arm set to 1.
#' The segments that are after the middle position have their chromosome arm set to 2. Then, the chromosome arm column is:
#'      chromosome_arm = chromosome_number + (chromosome_arm - 1)/2
#' So the segments whose arm was 1 have 0 added to their chromosome number, meaning that if a segment its chromosome arm is Chr_N.
#' But, the segments whose arm was 2 have 1 added to their chromosome number, meaning that its chromosome arm is Chr_N.5.
#' 
#' Edge case: segments that start before the middle pos but end after the middle pos. In that case, we split it into 2 and proceed like above.
#'
#' @param df The Data Frame to modify. Column numbers are same as in RemoveCentromereTelomeres. 
#' @param include_chr_X A boolean. Whether to set chromosome arms for segments in Chr X/23 or not.
#' @param centromere_positions An array: the middle positions for each chromosome IN ORDER from 1 to 23.
#'
#' @export
AddChromosomeArmHelper <- function(df, include_chr_X, centromere_positions) {
  if (include_chr_X == TRUE) {
    n_chr_to_check = 23 # so in for-loop next we look for the 23rd centromere
  } else {
    n_chr_to_check = 22 # otherwise skip it
  }
  
  for (i in 1:23) {
    pos = centromere_positions[i]
    
    # segments before centromere
    before_centromere <- df[,1] == i & df[,3] < pos & df[,4] < pos
    df[before_centromere, 2] <- 1
    
    # segments after centromere
    after_centromere <- df[,1] == i & df[,3] > pos & df[,4] > pos
    df[after_centromere, 2] <- 2
    
    # Dealing with segments that start before the centromere but end after the centromere
    if (i == 23) {
      new_start = pos
    } else {
      new_start = pos+1
    }
    
    before_then_after_centromere <- df[,1] == i & df[,3] < pos & df[,4] > pos
    df[before_then_after_centromere, 2] = 1
    
    end_pos = df[before_then_after_centromere, 4]
    df[before_then_after_centromere, 4] = pos
    
    if (length(which(before_then_after_centromere == TRUE)) == 0) {
      next
    }
    
    df = rbind(df[(1:which(before_then_after_centromere == TRUE)),],
                c(df[before_then_after_centromere,1], 2, new_start, end_pos, df[before_then_after_centromere,5]),
                df[-(1:which(before_then_after_centromere == TRUE)),])
    rownames(df) <- NULL
  }
  df[,2] <- df[,1]+(df[,2]-1)/2
  df
}

#' Merges segments that have the same ratio_median and are in the same chromosome arm.
#'
#' @description 
#' Takessegment data (df) and gathers/merges the segments according to their ratio_median
#' or chromosome arm. 
#' 
#' @param df The Data Frame with segment data: first column is chromosome number, second column is chromosome arm,
#' third column is start position, fourth is end position, and fifth is the ratio_median.
#'
#' @export
GatherByRatioMedian <- function(df) {
  N = dim(df)[1]
  gathered_by_ratio_median = matrix(0, ncol=6, nrow=N)
  
  i = 1
  c = 1
  
  while (i < N) {
    # if in the same chromosome arm
    if (df[i,2] == df[i+1,2]) {
      # if same ratio_median
      if (df[i,5] == df[i+1,5]) {
        j = 1
        # loop through all other segments until segment `i` and segment `j`
        # are no longer in the same chromosome arm and no longer have the same ratio_median
        
        # if `i+j` is the end of the file, then we merge data from `i` and `i+j`
        # if `i+j`` is not the end of the file, then we merge data from `i`` and `i+j-1`, because `j` 
        #                                     is the index at which we've already CHANGED chr_arm/ratio_median
        #                                     so `j-1` belongs to the same segment as `i`
        while (df[i,2] == df[i+j, 2] & df[i,5] == df[i+j,5] & i+j < N) {
          j = j + 1
        }
        
        # if we reached the end of the file
        if (i+j == N) {
          gathered_by_ratio_median[c,] = c(df[i,1], df[i,2], df[i,3], df[i+j, 4], df[i,5], df[i+j,4] - df[i,3] + 1)
          i = i + j # we move to the segment after i+j
          c = c + 1 # we move to the next row in gathered_by_ratio_median: next "merged" segment with unique ratio_median (or different chr_arm)
        } else {
          gathered_by_ratio_median[c,] = c(df[i,1], df[i,2], df[i,3], df[i+j-1, 4], df[i,5], df[i+j-1,4] - df[i,3] + 1)
          i = i + j
          c = c + 1
        }
      } else { # different ratio_median: save this segment `i` data, then move to next
        gathered_by_ratio_median[c,] = c(df[i,1], df[i,2], df[i,3], df[i,4], df[i,5], df[i,4] - df[i,3] + 1)
        i = i + 1
        c = c + 1
      }
    } else { # different chromosome arm: save this segment `i` data, then move to next
      gathered_by_ratio_median[c,] = c(df[i,1], df[i,2], df[i,3], df[i,4], df[i,5], df[i,4] - df[i,3] + 1)
      i = i + 1
      c = c + 1
    }
  }
  gathered_by_ratio_median
}
