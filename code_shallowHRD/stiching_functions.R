source("gather_segments_ratiomedian.R")
source("find_threshold_1.R")
source("reading_initialize.R")
source("small_segments.R")
source("call_lga.R")

#### Importing data####
raw_ratios_file <- data.frame(read.table(file="./test_data/example_1_controlfreec_final_hg19.bam_ratio.txt",header=TRUE))

#### Gathering segments by ratio_median and/or chromosome arm & removing centromeres, telomeres, etc.####
clean_ratios_file <- GetCleanBamRatiosFrame(raw_bam_ratios=raw_ratios_file) # this is ratio_file_tsv
gathered_by_ratio_median <- GatherSegmentsByRatioMedian(bam_ratios_frame=clean_ratios_file, 
                                                        include_chr_X=TRUE) # this is ratio_median_gathered.txt.

#### Finding the first threshold ####
granges_obj <- GetGRangesObject(bam_ratios_frame=clean_ratios_file)
first_threshold <- FindThreshold(granges_obj=granges_obj, 
                                 segments=gathered_by_ratio_median, 
                                 second_round=FALSE, 
                                 num_simulations = 10000)

#### Reading and Initialization ####
prepped_gathered_by_ratio_median <- InitializeReadingInitialization(segments=gathered_by_ratio_median) # Graph 1 `tmp`
segments_wo_short_arms <- ExcludeShortArms(segments=prepped_gathered_by_ratio_median) # Graph 2
segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)
segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                               segments_copy=segments_gt3mb,
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
                                         second_round=FALSE) # Graph 5
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
                                  segments_copy=segments_gt3mb,
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
write.table(v2_graph_4_copy, "./test_outputs/v2_small_segments/imported_v2_segments_small_inserted.txt", sep="\t", row.names = FALSE)

#### v2 Determining LGAs ####
v2_segments <- BreakSmoothToLGA(threshold=second_threshold,
                                segments=v2_segments,
                                granges_obj=granges_obj)
v2_segments <- LargeMissingChrArms(large_segments=v2_segments,
                                   small_segments=v2_graph_4_copy,
                                   second_round=TRUE)
v2_segments <- CorrectLeftovers(segments=v2_segments,
                                granges_obj=granges_obj)
v2_segments <- GetSegmentationBeforeLGACall(segments=v2_segments,
                                            bam_ratios_frame=gathered_by_ratio_median,
                                            second_round=TRUE) # v2 Graph 5
final_lga_res <- CallLGA(threshold=second_threshold,
                         segments=v2_segments) # Could save this as a .txt
write.table(final_lga_res, "./number_of_lgas.txt", sep="\t", row.names = FALSE)
final_lga_segments <- GetLGAOfSize(threshold=second_threshold,
                                   size_lga=10,
                                   segments=v2_segments) # Graph 6: to plot the LGAs 





