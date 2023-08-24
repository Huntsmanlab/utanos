source("hg19_segments.R")
source("helpers_gather_segments_ratiomedian.R")

#' Performs some cleaning operations.
#' 
#' @description 
#' These pre-processing steps are:
#' 1. Replaces chromosome 'X' for '23'.
#' 2. Removes segments in chromosome 'Y'.
#' 3. Adds an extra column for the End Position of the segment.
#' 4. Removes segments whose ratio and ratio_median are -1.
#' 5. Log2 transforms the ratio and ratio_median columns.
#' 6. Removes segments with a ratio_median = -Inf.
#' 
#' @param raw_bam_ratios A data frame: the raw bam ratios file obtained from the bam_ratio.txt file. 
#' 
#' @export
CleanBamRatiosFrame <- function(raw_bam_ratios) {
  # Cleaning up raw_bam_ratios file
  raw_bam_ratios = raw_bam_ratios[,1:4]
  raw_bam_ratios[,1] <- gsub('X', '23', raw_bam_ratios[,1])
  raw_bam_ratios = raw_bam_ratios[which(!raw_bam_ratios[,1] == 'Y'),]
  raw_bam_ratios[,1] = as.numeric(as.character(raw_bam_ratios[,1]))
  raw_bam_ratios = raw_bam_ratios[order(raw_bam_ratios[,1]),]
  
  # Adding 'end' column to raw_bam_ratios
  size_window = raw_bam_ratios[2,2] - raw_bam_ratios[1,2]
  raw_bam_ratios = cbind(1:dim(raw_bam_ratios)[1], 
                     raw_bam_ratios[,1], # chr_N
                     raw_bam_ratios[,2], #start
                     raw_bam_ratios[,2] + size_window - 1, # end
                     raw_bam_ratios[,3],# ratio
                     raw_bam_ratios[,4]) # ratio_median
  raw_bam_ratios = as.data.frame(raw_bam_ratios)
  colnames(raw_bam_ratios) = c("feature", "chromosome", "start", "end", "ratio", "ratio_median")
  # More cleaning + log transform
  raw_bam_ratios = raw_bam_ratios[which(!raw_bam_ratios$ratio == -1),]
  raw_bam_ratios = raw_bam_ratios[which(!raw_bam_ratios$ratio_median == -1),]
  raw_bam_ratios[,5] = log2(raw_bam_ratios[,5])
  raw_bam_ratios[,6] = log2(raw_bam_ratios[,6])
  
  # Removing segments with a ratio_median = -Inf
  inf_ratio_median_indices = which(bam_ratios_frame$ratio_median == -Inf)
  if (length(inf_ratio_median_indices) >= 1) {
    bam_ratios_frame = bam_ratios_frame[-inf_ratio_median_indices,]
  }
  
  bam_ratios_frame <- data.frame(raw_bam_ratios)
  bam_ratios_frame 
}

#' Removes specific genomic regions (centromeres and telomeres) of all chromosomes
#' 
#' @description 
#' First, drops the 'Feature' and the 'Ratio' columns.
#' Spurious regions are unstable regions that are difficult to map, so we remove them. The `hg19_segments.R` file
#' contains positions of these regions. We find them in the `bam_ratios_frame` and, if found, they are removed.
#' 
#' @param bam_ratios_frame A data frame: a cleaned up version of the original bam_ratios.txt file.
#' @param include_chr_X A boolean: whether we search for chromosome X or not.
#' 
#' @export
RemoveSpuriousRegions <- function(bam_ratios_frame, include_chr_X) {
  # Dropping the 'feature' and the 'ratio' column.
  bam_ratios_frame <- bam_ratios_frame[,-1]
  bam_ratios_frame <- bam_ratios_frame[,-4]
  bam_ratios_frame <- RemoveCentromereTelomeres(df=bam_ratios_frame, 
                                                include_chr_X=include_chr_X,
                                                centromere_starts=hg19_centromere_starts,
                                                centromere_ends=hg19_centromere_ends,
                                                telomere_2_starts=hg19_telomere_2_starts,
                                                telomere_2_ends=hg19_telomere_2_ends)
  bam_ratios_frame
}

#' Adds a column for the chromosome arm in which each segment is in.
#'
#' @description 
#' More specifically, iterates through all the segments, and using genomic data for the
#' middle positions of all chromosomes, it determines whether the segment is to the left
#' or to the right of the middle position. For example, if a segment is in chromosome 5
#' and is to the left of the middle position, it is in chromosome arm 1. If it were to the right
#' instead, then its chromosome arm would be 1.5.
#' 
#' Note that most of this iterative process is done in the AddChromosomeArmHelper helper.
#' 
#' @param bam_ratios_frame A data frame: the cleaned up version of the bam_ratios.txt file. In the algorithm,
#' this segment data has already had its spurious regions removed. 
#' @param include_chr_X A boolean: True if we include chromosome X/23, False otherwise.
#' 
#' @export
AddChromosomeArm <- function(bam_ratios_frame, include_chr_X) {
  # Adding the new column and adding new names to the columns
  bam_ratios_frame$chr_arm <- rep(0, nrow(bam_ratios_frame))
  colnames(bam_ratios_frame) <- c("chr", "start", "end", "ratio_median", "chr_arm")
  bam_ratios_frame = bam_ratios_frame[c("chr", "chr_arm", "start", "end", "ratio_median")]
  
  # Setting them to numeric
  bam_ratios_frame[,3] = as.numeric(as.character(bam_ratios_frame[,3]))
  bam_ratios_frame[,4] = as.numeric(as.character(bam_ratios_frame[,4]))
  
  # Adding the chr_arm column
  bam_ratios_frame <- AddChromosomeArmHelper(df=bam_ratios_frame, 
                                             include_chr_X=include_chr_X,
                                             centromere_positions=hg19_middle_positions)
  bam_ratios_frame
}

#' Merges segments that are in the same chromosome arm AND have the same ratio_median.
#' 
#' @description
#' Most of the work is done by GatherSegmentsByRatioMedianHelper. Essentially, we iterate through the rows of 
#' `bam_ratios_frame`, and keep track of adjacent rows that have the same ratio_median and are in the same chromosome arm.
#' Whenever any of these 2 conditions is false, we merge the segments up until that point. More information is in the helper function.
#' 
#' Finally, adds a 'Size' column to `bam_ratios_frame`.
#' 
#' @param bam_ratios_frame A data frame: the cleaned up version of the bam_ratios.txt file. In the pipeline,
#' this segment data has already had its spurious regions removed, and a chromosome arm column has been added/determined.
#' 
#' @export
GatherSegmentsByRatioMedian <- function(bam_ratios_frame) {
  bam_ratios_frame = as.matrix(bam_ratios_frame)
  gathered_by_ratio_median = GatherSegmentsByRatioMedianHelper(df=bam_ratios_frame)
  
  gathered_by_ratio_median = subset(gathered_by_ratio_median, gathered_by_ratio_median[,1] != 0)
  rownames(gathered_by_ratio_median) <- NULL
  colnames(gathered_by_ratio_median) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")
  
  gathered_by_ratio_median <- as.data.frame(gathered_by_ratio_median)
  gathered_by_ratio_median
}