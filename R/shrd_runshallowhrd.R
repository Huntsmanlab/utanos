#' Run shallowHRD on a dataframe or file containing segmented relative copy-number data.
#'
#' @description
#' Utanos' version of the shallowHRD algorithm for detecting homologous recombination deficiency (HRD) on tumor samples from the number of
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
#'    removing spurious/unstable regions from the chromosomes, and other cleaning stuff. More description on their activities can be found in
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
#'    ---- segment 1 ----   |                             diff <= threshold
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
#' 4) Small segments. For the previous section, we were only working with 'large' segments, meaning their size is >= 3Mb. Now we have to figure out where the small
#'    segments fit in into our current segment data. `InsertSmallSegments` takes care of this, and it essentially iterates over the small segments and the large segments,
#'    and checks if we encounter any of the 6 cases it looks for. The cases depend on the positioning of the small segment with respect to the large segment, and they are all
#'    listed/visualized in the helper `FinalizeSmallSegments` documentation. Wherever appropriate, we merge small and large segments too (i.e. if their ratio_median diff <= threshold).
#'    By the end of the section, we have re-inserted small segments into our main segment data frame.
#' 5) LGAs: now that we have re-inserted the small segments, we can finally determine the number of LGAs. `CallLGA` does this, and it is basically iterating over the
#'    finalized segments (marked as Graph 5) and checking which of the segments meet the criteria for an LGA (these requirements are listed in the aforementioned function's documentation.)
#'    The function returns the number of LGAs for different LGA sizes, the paper mostly cares about LGAs of size 10Mb, but you'll get results for sizes 3 to 11.
#'    You can also call `GetLGAOfSize`, which returns the segments that were marked as an LGA for the given LGA size, in case you want to visualize them.
#'
#'    To determine HRD, check the number of LGAs for 10Mb, and if the number of LGAs >= 20, then HRD was detected.
#'
#' We repeat this procedure from the second step for the second pass. Note how the `FindThreshold` function doesn't receive the bam_ratio segments this time. Instead, we pass the
#' finalized segments found in the first pass (Graph 5) - this is what we meant by 'improving' our results, as we are going to try to find a better threshold for this second
#' round of function calling.
#'
#' To finalize, while some of these choices seem arbitrary (e.g. why 3Mb as 'large' segments, and not 5MB? Why is >= 20 LGAs equal to HRD?), perhaps the following paper by the same author can give
#' some intuition on their choices. Please refer to https://pubmed.ncbi.nlm.nih.gov/22933060/.

#' Calculates the number of large genomic alterations for various sizes.
#'
#' @export
RunShallowHRD <- function(raw_ratios_file, log_transform=TRUE, include_chr_X=FALSE, num_simulations=100000,
                          shrd_save_path=FALSE, sample=NULL, seed=1337, plot=FALSE) {
  #### Importing data####
  if(inherits(x = raw_ratios_file, what = "character")) {
    raw_ratios_file <- data.frame(read.table(file=raw_ratios_file, header=TRUE))
  } else if (inherits(x = raw_ratios_file, what = "data.frame")) {
    raw_ratios_file <- raw_ratios_file
  } else {
    stop("Supply either a dataframe or a path.")
  }

  #### Gathering segments by ratio_median and/or chromosome arm & removing centromeres, telomeres, etc ####
  clean_ratios_file <- CleanBamRatiosFrame(raw_bam_ratios=raw_ratios_file,
                                           log_transform=log_transform)
  clean_ratios_file_copy <- clean_ratios_file
  clean_ratios_file <- RemoveSpuriousRegions(bam_ratios_frame=clean_ratios_file,
                                             include_chr_X=include_chr_X)
  clean_ratios_file <- AddChromosomeArm(bam_ratios_frame=clean_ratios_file,
                                        include_chr_X=include_chr_X)
  gathered_by_ratio_median <- GatherSegmentsByRatioMedian(bam_ratios_frame=clean_ratios_file)

  #### Finding the first threshold ####
  granges_obj <- GetGRangesObject(bam_ratios_frame=clean_ratios_file_copy)

  first_threshold <- FindThreshold(granges_obj=granges_obj,
                                     segments=gathered_by_ratio_median,
                                     second_round=FALSE,
                                     num_simulations=num_simulations,
                                     seed=seed)

  #### Levels ####
  prepped_gathered_by_ratio_median <- PrepForLevelsInitialization(segments=gathered_by_ratio_median, granges_obj=granges_obj)
  segments_wo_short_arms <- ExcludeShortArms(segments=prepped_gathered_by_ratio_median)
  segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
  segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)
  segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                                 thr=first_threshold)
  segments_gt3mb <- GatherSegmentsByLevels(segments=segments_gt3mb,
                                           granges_obj=granges_obj)

  #### Re-inserting small segments into segments_gt3mb ####
  segments_gt3mb <- LargeMissingChrArms(large_segments=segments_gt3mb,
                                        small_segments=segments_btw_0.1_3mb)
  segments_btw_0.1_3mb <- SmallMissingChrArms(large_segments=segments_gt3mb,
                                              small_segments=segments_btw_0.1_3mb)
  segments_gt3mb_small_reinserted <- InsertSmallSegments(large_segments=segments_gt3mb,
                                                         small_segments=segments_btw_0.1_3mb,
                                                         threshold=first_threshold,
                                                         granges_obj=granges_obj)
  segments <- CorrectLeftovers(segments=segments_gt3mb_small_reinserted,
                               granges_obj=granges_obj)

  #### Determining LGAs ####
  segments <- BreakSmoothToLGA(threshold=first_threshold,
                               segments=segments,
                               granges_obj=granges_obj)
  segments <- CorrectLeftovers(segments=segments,
                               granges_obj=granges_obj)
  segments <- GetSegmentationBeforeLGACall(segments=segments,
                                           bam_ratios_frame=clean_ratios_file_copy,
                                           granges_obj=granges_obj,
                                           second_round=FALSE)
  initial_lga_res <- CallLGA(threshold=first_threshold,
                             segments=segments)
  lga_segments <- GetLGAOfSize(threshold=first_threshold,
                               size_lga=10,
                               segments=segments)

  #### v2 Finding the second threshold ####
    second_threshold <- FindThreshold(granges_obj=granges_obj,
                                      segments=segments,
                                      second_round=TRUE,
                                      num_simulations=num_simulations,
                                      seed=seed)

  #### v2 Reading & Initialization ####
  segments_gt3mb <- GetLargeSegments(segments=segments_wo_short_arms)
  segments_btw_0.1_3mb <- GetSmallSegments(segments=segments_wo_short_arms)
  v2_segments_gt3mb <- AssignLevels(segments=segments_gt3mb,
                                    thr=second_threshold)
  v2_segments_gt3mb <- GatherSegmentsByLevels(segments=v2_segments_gt3mb,
                                              granges_obj=granges_obj)

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
                                  granges_obj=granges_obj)
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
                                              granges_obj = granges_obj,
                                              second_round=TRUE)
  final_lga_res <- CallLGA(threshold=second_threshold,
                           segments=v2_segments)
  final_lga_segments <- GetLGAOfSize(threshold=second_threshold,
                                     size_lga=10,
                                     segments=v2_segments)

  if (shrd_save_path != FALSE) {
    write.table(x = final_lga_segments, file = paste0(shrd_save_path, sample, "_LGAs.tsv"), sep = "\t", row.names = FALSE)
    write.table(x = final_lga_res, file = paste0(shrd_save_path, sample, "_number_LGAs.tsv"), sep = "\t", row.names = FALSE)
  }

  #### Determine HRD Status ####
  hrd_status <- GetHRDStatus(lga_calls=final_lga_res)

  #### Plot the segments at the start, after the first pass, and at the end ####
  if(plot) {
    segments_plot <- PlotSegmentChanges(segments_initial = clean_ratios_file_copy,
                                        segments_pass_one = segments,
                                        segments_pass_two = v2_segments,
                                        sample = sample)

    return(list("hrd_status" = hrd_status, "n_lga" = final_lga_res, "lga" = final_lga_segments, "plot" = segments_plot))
  }

  return(list("hrd_status" = hrd_status, "n_lga" = final_lga_res, "lga" = final_lga_segments))
}

#' Run ShallowHRD on a QDNAseq object.
#'
#' @description
#' This function takes a QDNAseqCopyNumbers object, extracts the binned and segmented copy-number data, and runs ShallowHRD on each sample.
#'
#' @param qdna_obj An S4 object of type QDNAseqCopyNumbers. This object must contain a copynumber slot and a segmented slot.
#' @param include_chr_X Whether or not to include the X-chromosome when running ShallowHRD.
#' @param num_simulations The number of simulations to use during the FindThreshold() step.
#' @param shrd_save_path An optional path to save the ShallowHRD results to.
#' @param plot Whether or not to return a plot showing how segments are merged as the algorithm runs.
#' @param seed A seed to use for PRNG-dependent functions. Ensures reproduciblity between runs.
#' @param cores The number of cores to use for running samples in parallel. If set to 1, no parallel processing will be used.
#'
#' @returns A list of ShallowHRD results, with one entry per sample.
#'
#' @export
RunShallowHRDFromQDNA <- function(qdna_obj, include_chr_X=FALSE, num_simulations=100000, shrd_save_path=FALSE,
                                  plot=FALSE, seed = 1337, cores = 1) {
  # Note that ExportBinsQDNAObj does not log normalize
  bin_df <- ExportBinsQDNAObj(object = qdna_obj, type = "copynumber", filter = TRUE) %>%
    tidyr::pivot_longer(cols = !c("feature", "chromosome", "start", "end"), names_to = "sample", values_to = "ratio") %>%
    dplyr::select(!feature)

  seg_df <- ExportBinsQDNAObj(object = qdna_obj, type = "segments", filter = TRUE) %>%
    tidyr::pivot_longer(cols = !c("feature", "chromosome", "start", "end"), names_to = "sample", values_to = "ratio_median") %>%
    dplyr::select(!feature)

  df_group <- dplyr::inner_join(x = bin_df, y = seg_df, by = dplyr::join_by(sample, chromosome, start, end)) %>%
    dplyr::select(!end) %>%
    dplyr::group_by(sample)

  df_list <- dplyr::group_split(df_group, .keep = FALSE)
  names(df_list) <- dplyr::group_keys(df_group)$sample

  if (cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)

    shrd_result <- foreach(sample=names(df_list)) %dopar% {
      RunShallowHRD(raw_ratios_file = as.data.frame(df_list[[sample]]),
                    log_transform = TRUE,
                    include_chr_X = include_chr_X,
                    num_simulations = num_simulations,
                    shrd_save_path = shrd_save_path,
                    sample = sample,
                    plot = plot,
                    seed = seed)
    }

  }
  else {
    shrd_result <- lapply(X = names(df_list), FUN = function(sample) RunShallowHRD(raw_ratios_file = as.data.frame(df_list[[sample]]),
                                                                                   log_transform = TRUE,
                                                                                   include_chr_X = include_chr_X,
                                                                                   num_simulations = num_simulations,
                                                                                   shrd_save_path = shrd_save_path,
                                                                                   sample = sample,
                                                                                   plot = plot,
                                                                                   seed = seed))
  }

  names(shrd_result) <- names(df_list)

  return(shrd_result)
 }
