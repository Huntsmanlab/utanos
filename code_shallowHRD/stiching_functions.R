source("gather_segments_ratiomedian.R")
source("find_threshold_1.R")
source("reading_initialize.R")
source("small_segments.R")

#### Importing data & cleaning it up####
raw_ratios_file <- data.frame(read.table(file="./test_data/example_1_controlfreec_final_hg19.bam_ratio.txt",header=TRUE))
clean_ratios_file <- GetCleanBamRatiosFrame(raw_bam_ratios=raw_ratios_file) # clean_ratios_file is ratio_file_tsv in shallowHRD code
 
#### Gathering segments by ratio_median and/or chromosome arm & removing centromeres, telomeres, etc.####
gathered_by_ratio_median <- GatherSegmentsByRatioMedian(bam_ratios_frame=clean_ratios_file, 
                                                        include_chr_X=TRUE) # this is ratio_median_gathered.txt.

#### Finding the first threshold ####
first_threshold <- FindThreshold(bam_ratios_frame=clean_ratios_file, 
                                 segments=gathered_by_ratio_median, 
                                 second_round=FALSE, 
                                 num_simulations = 100000)

#### Reading and Initialization ####
prepped_gathered_by_ratio_median <- InitializeReadingInitialization(segments=gathered_by_ratio_median) # Graph 1 `tmp`
segments_wo_short_arms <- ExcludeShortArms(segments=prepped_gathered_by_ratio_median) # Graph 2
segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                              segments_copy=segments_gt3mb,
                              thr=first_threshold)
segments_gt3mb <- GatherByLevels(segments=segments_gt3mb,
                                 bam_ratios_frame=clean_ratios_file) # Graph 3

#### Re-inserting small segments into segments_gt3mb####
segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)
segments_gt3mb <- LargeMissingChrArms(large_segments=segments_gt3mb,
                                      small_segments=segments_btw_0.1_3mb)
segments_btw_0.1_3mb <- SmallMissingChrArms(large_segments=segments_gt3mb,
                                            small_segments=segments_btw_0.1_3mb)
segments_gt3mb <- InsertSmallSegments(large_segments=segments_gt3mb,
                                      small_segments=segments_btw_0.1_3mb,
                                      bam_ratios_frame=clean_ratios_file,
                                      threshold=first_threshold)

