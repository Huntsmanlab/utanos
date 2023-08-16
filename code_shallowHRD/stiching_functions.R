source("gather_segments_ratiomedian.R")
source("find_threshold_1.R")
source("reading_initialize.R")
source("small_segments.R")
source("call_lga.R")

#### Importing data & cleaning it up ####
raw_ratios_file <- data.frame(read.table(file="./test_data/example_1_controlfreec_final_hg19.bam_ratio.txt",header=TRUE))
clean_ratios_file <- GetCleanBamRatiosFrame(raw_bam_ratios=raw_ratios_file) # this is ratio_file_tsv
 
#### Gathering segments by ratio_median and/or chromosome arm & removing centromeres, telomeres, etc.####
gathered_by_ratio_median <- GatherSegmentsByRatioMedian(bam_ratios_frame=clean_ratios_file, 
                                                        include_chr_X=TRUE) # this is ratio_median_gathered.txt.

#### Finding the first threshold ####
granges_obj <- GetGRangesObject(bam_ratios_frame=clean_ratios_file)
first_threshold <- FindThreshold(granges_obj=granges_obj, 
                                 segments=gathered_by_ratio_median, 
                                 second_round=FALSE, 
                                 num_simulations = 100000)

#### Reading and Initialization ####
prepped_gathered_by_ratio_median <- InitializeReadingInitialization(segments=gathered_by_ratio_median) # Graph 1 `tmp`
segments_wo_short_arms <- ExcludeShortArms(segments=prepped_gathered_by_ratio_median) # Graph 2

segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)

segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                               segments_copy=segments_gt3mb,
                               thr=first_threshold)
segments_gt3mb <- GatherSegmentsByLevels(segments=segments_gt3mb,
                                         bam_ratios_frame=clean_ratios_file) # Graph 3

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
                             threshold=first_threshold,
                             granges_obj=granges_obj) # Graph 4

#### Determining LGAs ####
segments <- BreakSmoothToLGA(threshold=first_threshold,
                             segments=segments,
                             granges_obj=granges_obj)
segments <- CorrectLeftovers(segments=segments,
                             threshold=first_threshold,
                             granges_obj=granges_obj)
segments <- GetSegmentationBeforeLGACall(segments=segments,
                                         bam_ratios_frame=clean_ratios_file) # Graph 5
intial_lga_res <- CallLGA(threshold=first_threshold,
                          segments=segments,
                          second_round=FALSE) # Could save this as a .txt
lga_segments <- GetLGAOfSize(threshold=first_threshold,
                             size_lga=10,
                             segments=segments) # Graph 6: to plot the LGAs 

#### v2 Finding the second threshold ####
second_threshold <- FindThreshold(granges_obj=granges_obj, 
                                  segments=segments, 
                                  second_round=TRUE, 
                                  num_simulations = 100000)





