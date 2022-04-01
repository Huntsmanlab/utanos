# Read in Metadata
cellularity <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/cellularity_data.tsv', sep = '\t', header = TRUE)
samples <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/CC_p53abn_cases_sWGS_batch_456.tsv', sep = '\t', header = TRUE)
colnames(cellularity) <- c('study_id', 'cellularity', 'accession_number', 'block', 'type', 'date')
samples_new <- left_join(samples, cellularity, by = 'study_id') %>% select(sb_id, acc_num, study_id, cohort, hist, grade, cellularity)
write_tsv(samples_new, file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batches456.tsv', col_names = TRUE)

