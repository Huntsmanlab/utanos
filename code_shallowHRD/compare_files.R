# clean_ratios_file
imported_clean_ratios_file = read.table(file="./test_outputs/gathered_by_ratio_median/clean_ratios_file.txt", sep="\t", header=TRUE)
shallow_clean_ratios_file = read.table(file="./test_outputs/gathered_by_ratio_median/example_1_controlfreec_final_hg19_ratio_file_tsv.txt", sep="\t", header=TRUE)

# gathered_by_ratio_median
imported_gathered_by_ratio_median = read.table(file="./test_outputs/gathered_by_ratio_median/gathered_by_ratio_median.txt", sep="\t", header=TRUE)
shallow_gathered_by_ratio_median = read.table(file="./test_outputs/gathered_by_ratio_median/example_1_controlfreec_final_hg19_ratio_median_gathered.txt", sep="\t", header=TRUE)

# reading/initialize
imported_wo_short_arms = read.table(file="./test_outputs/reading_init/segments_wo_short_arms.txt", sep="\t", header=TRUE)
shallow_wo_short_arms = read.table(file="./test_outputs/reading_init/gathered_by_ratio_medianshort_ams.txt", sep="\t", header=TRUE)

imported_graphIII = read.table(file="./test_outputs/reading_init/imported_graphIII.txt", sep="\t", header=TRUE)
shallow_graphIII = read.table(file="./test_outputs/reading_init/shallow_graphIII.txt", sep="\t", header=TRUE)

# small segments
imported_segments_gt3mb = read.table(file="./test_outputs/small_segments/imported_segments_gt3mb.txt", sep="\t", header=TRUE)
shallow_segments_gt3mb = read.table(file="./test_outputs/small_segments/shallow_segments_gt3mb.txt", sep="\t", header=TRUE)

imported_segments_btw_0.1_3mb = read.table(file="./test_outputs/small_segments/imported_segments_btw_0.1_3mb.txt", sep="\t", header=TRUE)
shallow_segments_btw_0.1_3mb = read.table(file="./test_outputs/small_segments/shallow_segments_btw_0.1_3mb.txt", sep="\t", header=TRUE)

imported_segments_reinserted = read.table(file="./test_outputs/small_segments/imported_segments_inserted.txt", sep="\t", header=TRUE)
shallow_segments_reinserted = read.table(file="./test_outputs/small_segments/shallow_segments_inserted.txt", sep="\t", header=TRUE)

imported_segments_leftovers = read.table(file="./test_outputs/small_segments/imported_segments_leftovers.txt", sep="\t", header=TRUE)
shallow_segments_leftovers = read.table(file="./test_outputs/small_segments/shallow_segments_leftovers.txt", sep="\t", header=TRUE)
