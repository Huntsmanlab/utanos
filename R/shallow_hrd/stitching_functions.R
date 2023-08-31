source("gather_segments_ratiomedian.R")
source("find_threshold_1.R")
source("levels.R")
source("small_segments.R")
source("call_lga.R")
source("plot_segments.R")
source("hg19_segments.R")

#' Utanos version of the shallowHRD algorithm for detecting homologous recombination deficiency (HRD) on tumor samples from the number of 
#' large genomic alterations (LGAs). Original code found here (https://github.com/aeeckhou/shallowHRD), and original paper found here (https://academic.oup.com/bioinformatics/article/36/12/3888/5823300).
#' Essentially, the algorithm can be divided into 5 sequential steps. We perform these steps two times (I call them 'passes'), 
#' so the first pass is sort of a 'discovery' stage, with still good estimates, but the second pass improves on the results 
#' from the first pass. 
#'
#' The overall idea is not so complicated: we have some raw segment data and want to find the number of LGAs of size >= 10Mb. 
#' To do so, the 5 steps are (excluding the cleaning functions):
#' 1) Gather/merge segments by ratio_median: segments that have the same ratio_median might as well be considered as the same
#'    segment: imagine plotting the segments --segment 1-- --segment 2--, since they have the same ratio_median, then they're
#'    on the same 'level' in the y-axis, so they can also be: ---- a single segment----. We merge segments in the same ratio_median AND in the same
#'    chromosome arm. 
#'    
#'    This is what `GatherSegmentsByRatioMedian` does. The other functions in the section do other minor things: adding chromosome arm column,
#'    removing spurious/unstable regions from the chromosomes, and other cleaning functions. More description on their activites can be found in
#'    their respective documentation.
#'    
#'    NOTE: in `CleanBamRatiosFrame` you have the option to log-transform the ratio_median's. If some ratio_medians in your file are negative,
#'    then this will produce NANs, in which case there's no need to log-transform. 
#'    
#' 2) Finding the threshold: we use KDE to find the minima of the density of all the ratio_median differences (between all pairs of segments).
#'    The idea here is that by finding the local minima of these differences, then we are finding a 'tolerance' point that determines whether
#'    we merge two segments or not. For example, imagine a plot again, and we have two segments (different ratio_medians, so they're not on the same level
#'    in the y-axis) :
#'    
#'    ---- segment 1 ----   |                              <= threshold
#'                          |                           ------------------>          ---- segment 1 + segment 2 ----
#'                          |  ---- segment 2 ----
#'                          
#'    Then the threshold indicates whether their distance is meaningful enough such that we should consider them a single segment. So, if their ratio_median
#'    difference is <= threshold, then we consider them the same segment because we aren't completely sure that they are far apart enough. If their ratio_median
#'    difference is > threshold, then we don't merge and we treat each segment as separate segments. All of this is illustrated above. More details on the procedure
#'    of finding the threshold find in the `FindThreshold` + helpers documentation.
#' 3) Levels: while we are still not entirely sure the reasoning behind the section, we think that this section assigns segments a sort of 'confidence' in our merging. 
#'    The idea here is that the largest segment likely has the highest confidence in its ratio_median as it was easier to sequence. Then, we use it as our reference
#'    to determine how confident we are with the rest of the segments. We use the threshold to see how close the other segments are to this largest segment, and if they are
#'    close (i.e. ratio_median difference <= threshold), then we assign it the same level as the largest segment.
#'    `AssignLevels` takes care of this. Finally, `GatherSegmentsByLevels` then merges segments in the same chromosome arm and with the same level.
#' 4) Small segments. For the previous section, we were only working with 'large' segments, meaning their size is >= 3Mb. 


#### Importing data####
raw_ratios_file <- data.frame(read.table(file="./test_data/example_2_QDNAseq_final_chrX_hg19.bam_ratio.txt",header=TRUE))

#### Gathering segments by ratio_median and/or chromosome arm & removing centromeres, telomeres, etc.####
clean_ratios_file <- CleanBamRatiosFrame(raw_bam_ratios=raw_ratios_file,
                                         log_transform=FALSE) 
clean_ratios_file_copy <- clean_ratios_file # this is ratio_file_tsv
clean_ratios_file <- RemoveSpuriousRegions(bam_ratios_frame=clean_ratios_file, 
                                           include_chr_X=TRUE)
clean_ratios_file <- AddChromosomeArm(bam_ratios_frame=clean_ratios_file,
                                      include_chr_X=TRUE)
gathered_by_ratio_median <- GatherSegmentsByRatioMedian(bam_ratios_frame=clean_ratios_file) # this is ratio_median_gathered.txt

#### Finding the first threshold ####
granges_obj <- GetGRangesObject(bam_ratios_frame=clean_ratios_file_copy)
first_threshold <- FindThreshold(granges_obj=granges_obj, 
                                 segments=gathered_by_ratio_median, 
                                 second_round=FALSE, 
                                 num_simulations = 10000)

#### Levels ####
prepped_gathered_by_ratio_median <- PrepForLevelsInitialization(segments=gathered_by_ratio_median) # Graph 1 `tmp`
segments_wo_short_arms <- ExcludeShortArms(segments=prepped_gathered_by_ratio_median) # Graph 2
segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)
segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                               thr=first_threshold)
segments_gt3mb <- GatherSegmentsByLevels(segments=segments_gt3mb,
                                         granges_obj=granges_obj) # Graph 3

#### Re-inserting small segments into segments_gt3mb####'
segments_gt3mb <- LargeMissingChrArms(large_segments=segments_gt3mb,
                                      small_segments=segments_btw_0.1_3mb)
segments_btw_0.1_3mb <- SmallMissingChrArms(large_segments=segments_gt3mb,
                                            small_segments=segments_btw_0.1_3mb)
segments_gt3mb_small_reinserted <- InsertSmallSegments(large_segments=segments_gt3mb,
                                                       small_segments=segments_btw_0.1_3mb,
                                                       threshold=first_threshold,
                                                       granges_obj=granges_obj)
segments <- CorrectLeftovers(segments=segments_gt3mb_small_reinserted,
                             granges_obj=granges_obj) # Graph 4
#### Determining LGAs ####
segments <- BreakSmoothToLGA(threshold=first_threshold,
                             segments=segments,
                             granges_obj=granges_obj)
segments <- CorrectLeftovers(segments=segments,
                             granges_obj=granges_obj)
segments <- GetSegmentationBeforeLGACall(segments=segments,
                                         bam_ratios_frame=clean_ratios_file,
                                         second_round=FALSE) # Graph 5: final segmentation
initial_lga_res <- CallLGA(threshold=first_threshold,
                          segments=segments) # Could save this as a .txt
lga_segments <- GetLGAOfSize(threshold=first_threshold,
                             size_lga=10,
                             segments=segments) # Graph 6: to plot the LGAs 

#### v2 Finding the second threshold ####
second_threshold <- FindThreshold(granges_obj=granges_obj, 
                                  segments=segments, 
                                  second_round=TRUE, 
                                  num_simulations = 10000)

#### v2 Reading & Initialization ####
segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)

v2_segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                                  thr=second_threshold)
v2_segments_gt3mb <- GatherSegmentsByLevels(segments=v2_segments_gt3mb,
                                            granges_obj=granges_obj) # v2 Graph 3

#### v2 Re-inserting small segments ####
v2_segments_gt3mb <- LargeMissingChrArms(large_segments=v2_segments_gt3mb,
                                         small_segments=segments_btw_0.1_3mb)
v2_segments_btw_0.1_3mb <- SmallMissingChrArms(large_segments=v2_segments_gt3mb,
                                               small_segments=segments_btw_0.1_3mb)
v2_segments_gt3mb_small_reinserted <- InsertSmallSegments(large_segments=v2_segments_gt3mb,
                                                          small_segments=v2_segments_btw_0.1_3mb,
                                                          threshold=second_threshold,
                                                          granges_obj=granges_obj)
v2_segments <- CorrectLeftovers(segments=v2_segments_gt3mb_small_reinserted,
                                granges_obj=granges_obj) # v2 Graph 4
v2_graph_4_copy <- v2_segments

#### v2 Determining LGAs ####
v2_segments <- BreakSmoothToLGA(threshold=second_threshold,
                                segments=v2_segments,
                                granges_obj=granges_obj)
v2_segments <- LargeMissingChrArms(large_segments=v2_segments,
                                   small_segments=v2_graph_4_copy)
v2_segments <- CorrectLeftovers(segments=v2_segments,
                                granges_obj=granges_obj)
v2_segments <- GetSegmentationBeforeLGACall(segments=v2_segments,
                                            bam_ratios_frame=gathered_by_ratio_median,
                                            second_round=TRUE) # v2 Graph 5: final segmentation
final_lga_res <- CallLGA(threshold=second_threshold,
                         segments=v2_segments) # Could save this as a .txt
final_lga_segments <- GetLGAOfSize(threshold=second_threshold,
                                   size_lga=10,
                                   segments=v2_segments) # Graph 6: to plot the LGAs 
write.table(final_lga_res, "./2qdna_number_of_lgas.txt", sep="\t", row.names = FALSE)

#### Plotting ####
graph <- PlotSegments(gathered_by_ratio_median=prepped_gathered_by_ratio_median,
                      bam_ratios_frame=clean_ratios_file_copy,
                      segments=segments,
                      chr_mid_positions=hg19_middle_positions)
suppressWarnings(ggsave("./final_segmentation.jpeg",".jpeg", plot = graph, device = "jpeg", width = 23, height = 13))  
