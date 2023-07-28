source("hg19_segments.R")
source("helpers_gather_segments_ratiomedian.R")

#' TODO: WRITE DOCUMENTATION
GetCleanBamRatiosFrame <- function(raw_bam_ratios) {
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
  bam_ratios_frame <- data.frame(raw_bam_ratios)
  bam_ratios_frame #bam_ratios_frame
}

#' TODO: WRITE DOCUMENTATION
GatherSegmentsByRatioMedian <- function(bam_ratios_frame, include_chr_X) {
  # Creating final data frames we'll actually work with
  copy_segment_ratios = bam_ratios_frame
  copy_segment_ratios = copy_segment_ratios[,-1]
  copy_segment_ratios = copy_segment_ratios[,-4]

  #### Removing spurious regions (cemtromeres and telomeres) ####
  copy_segment_ratios = RemoveCentromereTelomeres(df=copy_segment_ratios, 
                                                  include_chr_X=include_chr_X,
                                                  centromere_starts=hg19_centromere_starts,
                                                  centromere_ends=hg19_centromere_ends,
                                                  telomere_2_starts=hg19_telomere_2_starts,
                                                  telomere_2_ends=hg19_telomere_2_ends)
  
  print("Exited RemoveCentromereTelomeres")
  #### Adding Chr_Arm column ####
  # Adding the new column and adding new names to the columns
  copy_segment_ratios$chr_arm <- rep(0, nrow(copy_segment_ratios))
  colnames(copy_segment_ratios) <- c("chr", "start", "end", "ratio_median", "chr_arm")
  copy_segment_ratios = copy_segment_ratios[c("chr", "chr_arm", "start", "end", "ratio_median")]
  
  # Setting them to numeric
  copy_segment_ratios[,3] = as.numeric(as.character(copy_segment_ratios[,3]))
  copy_segment_ratios[,4] = as.numeric(as.character(copy_segment_ratios[,4]))
  
  # Adding the chr_arm column
  copy_segment_ratios <- AddChromosomeArmColumn(df=copy_segment_ratios, 
                                                include_chr_X=include_chr_X,
                                                centromere_positions=hg19_middle_positions)
  print("Exited AddChromosomeArmColumn")

  #### Removing segments with ratio_median=-Inf ####
  inf_ratio_median_indices = which(copy_segment_ratios$ratio_median == -Inf)
  if (length(inf_ratio_median_indices) >= 1) {
    copy_segment_ratios = copy_segment_ratios[-inf_ratio_median_indices,]
  }
  
  #### Gathering segments by ratio_median ####
  copy_segment_ratios = as.matrix(copy_segment_ratios)
  gathered_by_ratio_median = GatherByRatioMedian(df=copy_segment_ratios)
  print("Exited GatherByRatioMedian")
  gathered_by_ratio_median = subset(gathered_by_ratio_median, gathered_by_ratio_median[,1] != 0)
  rownames(gathered_by_ratio_median) <- NULL
  colnames(gathered_by_ratio_median) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")
  
  result = as.data.frame(gathered_by_ratio_median)
  result
}