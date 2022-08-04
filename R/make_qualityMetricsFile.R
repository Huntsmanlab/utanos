# make_qualityMetricsFile.R
suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

#####
#####
##### Read in shallowHRD output
s.hrd.files.path <- '~/Documents/projects/cn_sigs_swgs/shallowHRDoutput'
s.hrd.files <- list.files(path = s.hrd.files.path, full.names = TRUE)

s.hrd.out <- data.frame()
for (i in 1:length(s.hrd.files)) {
  df <- data.table::fread(file = s.hrd.files[i])
  out <- data.frame(cna_cutoff = df$V2[1], n_sims_cutoff = df$V2[2],
                    n_sims_succes = df$V2[3], MAX2 = df$V2[4],
                    cMAD = df$V2[5], n_lgas_per_10mb = df$V2[6],
                    hrd_status = df$V2[7])
  s.hrd.out <- rbind(s.hrd.out, out)
}

metadata <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/metadata/metadata.tsv')
sample_ids <- gsub(s.hrd.files, pattern = '/Users/mdouglas/Documents/projects/cn_sigs_swgs/shallowHRDoutput/', replacement = '')
s.hrd.out$sample_id <- gsub(sample_ids, pattern = '_summary_table.tsv', replacement = '')

# Begin quality file
metadata <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/metadata/metadata.tsv')
metadata <- metadata %>% dplyr::select(sample_id, tissue, status, cancer_type)
quality_file <- s.hrd.out %>% dplyr::select(sample_id, cMAD, MAX2) %>% 
                              dplyr::left_join(metadata, by = c('sample_id'))

#####
#####
##### Add coverage/depth metrics
# mean_prop_cov = on average, each base has been sequenced X number of times
coverage.files.path <- '~/Documents/projects/cn_sigs_swgs/data/coverage/se'
coverage.files <- list.files(path = coverage.files.path, full.names = TRUE)

coverage <- data.frame()
for (i in 1:length(coverage.files)) {
  df <- data.table::fread(file = coverage.files[i])
  df <- df[1:23,]
  out <- data.frame(percent_coverage = mean(df$coverage), 
                    mean_X_coverage = mean(df$meandepth),
                    nreads = sum(df$numreads))
  coverage <- rbind(coverage, out)
}
sample_ids <- gsub(coverage.files, pattern = '/Users/mdouglas/Documents/projects/cn_sigs_swgs/data/coverage/se/', replacement = '')
coverage$sample_id <- gsub(sample_ids, pattern = '.se.coverageTable.tsv', replacement = '')

quality_file <- quality_file %>% dplyr::left_join(coverage, by = c('sample_id'))
quality_file <- quality_file[,c(1,4,5,6,2,3,7,8,9)]


#####
#####
##### Add dispersion metrics
Xchr_30kb <- "~/Documents/projects/cn_sigs_swgs/copy_number_objects/Xchr_included/30kb_rCN_comCNVfilt.rds"
Xchr_30kb <- readRDS(file = Xchr_30kb)
dispersions <- data.frame(sample_id = sampleNames(Xchr_30kb), 
                       sd_dev = apply(Xchr_30kb@assayData[["copynumber"]], 
                                      MARGIN = 2, na.rm = TRUE, FUN = sd), 
                       mad = apply(Xchr_30kb@assayData[["copynumber"]], 
                                   MARGIN = 2, na.rm = TRUE, FUN = mad))
quality_file <- quality_file %>% dplyr::left_join(dispersions, by = c('sample_id'))


#####
#####
##### Add hand-annotated quality calls from Juliana and Maxwell
bio_met <- "~/Documents/projects/cn_sigs_swgs/metadata/swgs_sample_biologic_metadata.tsv"
bio_met <- data.table::fread(file = bio_met, sep = '\t')
bio_met <- bio_met %>% dplyr::select(other_aliases, quality_maxwell, quality_juliana)
quality_file <- quality_file %>% dplyr::left_join(bio_met, by = c("sample_id" = "other_aliases"))

quality_file <- quality_file[c(1,2,3,4,9,5,6,7,8,10,11,12,13)]
write.table(quality_file, file = "~/Documents/projects/cn_sigs_swgs/metadata/swgs_sample_quality_metrics.tsv", sep = '\t', row.names = FALSE)
