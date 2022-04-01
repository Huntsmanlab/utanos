library(rascal)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)


####################################
# Choose single solutions from many rascal batch runs across window sizes
####################################
metadata <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv", sep = '\t', header = TRUE)
metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")

rascal_batch_solutions_1 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_123/autosomes_only/sWGS_pilot_solutions_100kb_1.csv", sep = ',', header = TRUE)
rascal_batch_solutions_1$sample <- str_replace_all(rascal_batch_solutions_1$sample, "-", ".")
rascal_batch_solutions_2 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_123/autosomes_only/sWGS_pilot_solutions_500kb_1.csv", sep = ',', header = TRUE)
rascal_batch_solutions_2$sample <- str_replace_all(rascal_batch_solutions_2$sample, "-", ".")

rbs_xchr_15 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/15kb_rascal_solutions.csv", sep = ',', header = TRUE)
rbs_xchr_15$sample <- str_replace_all(rbs_xchr_15$sample, "-", ".")
rbs_xchr_30 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/30kb_rascal_solutions.csv", sep = ',', header = TRUE)
rbs_xchr_30$sample <- str_replace_all(rbs_xchr_30$sample, "-", ".")
rbs_xchr_50 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/50kb_rascal_solutions.csv", sep = ',', header = TRUE)
rbs_xchr_50$sample <- str_replace_all(rbs_xchr_50$sample, "-", ".")
rbs_xchr_100 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/100kb_rascal_solutions.csv", sep = ',', header = TRUE)
rbs_xchr_100$sample <- str_replace_all(rbs_xchr_100$sample, "-", ".")
rbs_xchr_500 <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/500kb_rascal_solutions.csv", sep = ',', header = TRUE)
rbs_xchr_500$sample <- str_replace_all(rbs_xchr_500$sample, "-", ".")

names(which(table(rbs_xchr_15$sample) == 1))

test1 <- rbs_xchr_15[!(rbs_xchr_15$sample %in% names(which(table(rbs_xchr_15$sample) == 1))),]
test2 <- rbs_xchr_30[!(rbs_xchr_30$sample %in% names(which(table(rbs_xchr_30$sample) == 1))),]
test3 <- rbs_xchr_50[!(rbs_xchr_50$sample %in% names(which(table(rbs_xchr_50$sample) == 1))),]
test4 <- rbs_xchr_100[!(rbs_xchr_100$sample %in% names(which(table(rbs_xchr_100$sample) == 1))),]
test5 <- rbs_xchr_500[!(rbs_xchr_500$sample %in% names(which(table(rbs_xchr_500$sample) == 1))),]
test1 <- as.data.frame(table(test1$sample))
test2 <- as.data.frame(table(test2$sample))
test3 <- as.data.frame(table(test3$sample))
test4 <- as.data.frame(table(test4$sample))
test5 <- as.data.frame(table(test5$sample))
test <- full_join(test1, test2, by="Var1")
test <- full_join(test, test3, by="Var1")
test <- full_join(test, test4, by="Var1")
test <- full_join(test, test5, by="Var1")
names(test) <- c('sample', '15kb', '30kb', '50kb', '100kb', '500kb')
test <- left_join(test, metadata, by=c('sample' = 'sample_id'))
test <- test[,c(1:7)]

# Filter out certain batches if wanted
test <- test %>% filter(batch == 3)
  
# Give the chart file a name.
png(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/NrascalSols_batch3_line_chart.png", width = 700, height = 700)
# Plot the bar chart.
plot(test$`15kb`,type = "o",col = "red", xlab = "Sample", ylab = "Solution Number", 
     main = "Number of rascal solutions for samples by bin size", ylim = c(0,18))
lines(test$`30kb`, type = "o", col = "orange")
lines(test$`50kb`, type = "o", col = "blue")
lines(test$`100kb`, type = "o", col = "green")
lines(test$`500kb`, type = "o", col = "black")
legend(15, 18, legend=c('15kb total: 49', '30kb total: 42', '50kb total: 40', '100kb total: 44', '500kb total: 42'),
       col=c("red", "orange", "blue", "green", "black"), lty=1, cex=1)
# Save the file.
dev.off()


#################################
### Weight rascal batch solutions by MAD Create Additional Parameter
#################################
rascal_batch_solutions <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/500kb_rascal_solutions.csv", sep = ',', header = TRUE)
rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")

test3 <- rascal_batch_solutions %>% 
  group_by(sample) %>% 
  summarise(avg = round(10*(((1-distance)-min(1-distance)+0.01)/sum((1-distance)-min(1-distance)+0.01))))
weighted_sols <- data.frame(sample = rep(rascal_batch_solutions$sample, test3$avg), copy_numbers = rep(rascal_batch_solutions$ploidy, test3$avg), cellularity = rep(rascal_batch_solutions$cellularity, test3$avg))
write_csv(weighted_sols, file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/500kb_rascal_multisolutions.csv', col_names = TRUE)

##############
## MAD optimal column addition and segs table creation
##############
solutions <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_copyNumbersSegmented_solutions.txt"
rcn_segs <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_copyNumbersSegmented.tsv"
acn_segs <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalMAD.rds"
getOptimalMADSolutions(solutions, acn_segs)

##############
## Add VAF column to rascal solutions and create corresponding segs table
##############
solutions <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_copyNumbersSegmented_solutions.txt"
rcn_segs <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_comCNVfilt.tsv"
acn_segs_save_path <- "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalVAF.rds"
calculate_vaf_acns(solutions, rcn_segs, acn_segs_save_path)

####################################
# Create CN-Collapsed Segment Tables
# Single chosen solutions
####################################



####################################
# Create Rascal CN-plots
####################################
# RDS read-in
qdnaseq_segs <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_copyNumbersSegmented.tsv'
qdnaseq_segs <- read.table(file = qdnaseq_segs, header = TRUE)
qdnaseq_cns <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_copyNumbersSmooth.tsv'
qdnaseq_cns <- read.table(file = qdnaseq_cns, header = TRUE)

rascal_batch_solutions <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_copyNumbersSegmented_solutions_mad_optimal.csv", sep = ',', header = TRUE)
rascal_batch_solutions$sample <- str_replace_all(rascal_batch_solutions$sample, "-", ".")

# Convert to long
qdnaseq_cns <- gather(qdnaseq_cns, sample, copy_number, `CC.CHM.1341`:`YW.EC052`, factor_key=TRUE)
qdnaseq_segments <- gather(qdnaseq_segs, sample, segmented, `CC.CHM.1341`:`YW.EC052`, factor_key=TRUE)

# Collapse to continuous segments
segments <- copy_number_segments(qdnaseq_segments)

# Have a combined dataframe
qdnaseq_cns <- qdnaseq_cns %>% mutate(segmented = qdnaseq_segments$segmented)

plot_abs_copy_number <- function(x, sample, sample_segments, sample_cns, output) {
  # browser()
  ploidy <- as.numeric(x["ploidy"])
  cellularity <- as.numeric(x["cellularity"])
  chr_order <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X')
  absolute_segments <- dplyr::mutate(sample_segments, copy_number = relative_to_absolute_copy_number(copy_number, ploidy, cellularity))
  absolute_copy_number <- dplyr::mutate(sample_cns, across(c(copy_number, segmented), relative_to_absolute_copy_number, ploidy, cellularity))
  chromosomes <- chromosome_offsets(absolute_copy_number)
  # chromosomes <- chromosomes %>% mutate(chromosome = factor(chromosome, levels = chr_order)) %>% arrange(chromosome)
  genomic_copy_number <- convert_to_genomic_coordinates(absolute_copy_number, "position", chromosomes)
  # genomic_copy_number <- genomic_copy_number %>% mutate(chromosome = factor(chromosome, levels = chr_order)) %>% arrange(chromosome)
  genomic_segments <- convert_to_genomic_coordinates(absolute_segments, c("start", "end"), chromosomes)
  # genomic_segments <- genomic_segments %>% mutate(chromosome = factor(chromosome, levels = chr_order)) %>% arrange(chromosome)
  p <- genome_copy_number_plot(genomic_copy_number, genomic_segments, chromosomes,
                          min_copy_number = 0, max_copy_number = 15,
                          copy_number_breaks = 0:15,
                          point_colour = "grey40",
                          ylabel = "absolute copy number") +
    ggplot2::ggtitle(paste0(sample, '  ploidy:', ploidy, '  cellularity:', cellularity)) +
    ggplot2::theme(plot.title = element_text(size = 20))
  ggsave(filename = paste0('~/Documents/projects/cn_signatures_shallowWGS/plotting/acn_rascal_plots/batch1-13_autosomes_30kb_rascal_plots/', sample, '.pl', ploidy, '.cel', cellularity, '.copynumberplot.png'),
         p, device = 'png', width = 16, height = 8)
}

for (i in unique(segments$sample)) {
  sample_segments <- dplyr::filter(segments, sample == i)
  sample_cns <- dplyr::filter(qdnaseq_cns, sample == i)
  sample_cns <- sample_cns %>% dplyr::mutate(position = round((start + end) / 2))
  print(i)
  rascal_sample_soln <- rascal_batch_solutions[rascal_batch_solutions$sample == i,]
  if (is.null(dim(rascal_sample_soln))) next
  if (dim(rascal_sample_soln)[1] == 0 ) next
  print(rascal_sample_soln)
  apply(rascal_sample_soln, 1, plot_abs_copy_number, sample = i, sample_segments, sample_cns, output = 'outputfile')
}


####################################
# Create CN-Collapsed Segment Tables
# Multi-solutions object creation
####################################

weighted_rascal_sols <- read.table(file = "~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/500kb_rascal_multisolutions.csv", sep = ',', header = TRUE)
weighted_rascal_sols$sample <- str_replace_all(weighted_rascal_sols$sample, "-", ".")

rows = function(tab) lapply(
  seq_len(nrow(tab)),
  function(i) unclass(tab[i,,drop=F])
)

chosenSegmentTablesList <- list()
j <- 1
for (i in unique(weighted_rascal_sols$sample)) {
  sample_segments <- filter(segments, sample == i)
  sols <- weighted_rascal_sols[weighted_rascal_sols$sample == i,]
  out<-c()
  for (k in rows(sols)) {
    absolute_segments <- mutate(sample_segments, copy_number = relative_to_absolute_copy_number(copy_number, k$copy_numbers, k$cellularity))
    absolute_segments <- absolute_segments %>% select(chromosome=chromosome, start=start, end=end, segVal=copy_number)
    out<-rbind(out,absolute_segments)
  }
  chosenSegmentTablesList[[j]] <- as.data.frame(out)
  j = j+1
}
names(chosenSegmentTablesList) <- unique(weighted_rascal_sols$sample)
saveRDS(chosenSegmentTablesList, file = "/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/500kb_rascal_multisolutions_CN_Collapsed_segments.rds")


####################################
# Examine samples individually
####################################

sample <- 'CC.VGH.1193.T'
ss_segs <- segments %>% filter(sample == sample)
ss_segs <- segments %>% filter(sample == 'CC.VGH.1193.N')
ss_cns <- qdnaseq_cns %>% filter(sample == sample)
absolute_copy_number_distance_heatmap(
  ss_segs$copy_number,
  ss_segs$weight,
  distance_function = "MAD"
)

test <- find_best_fit_solutions(
  ss_segs$copy_number, ss_segs$weight,
  min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
  min_cellularity = 0.1, max_cellularity = 1.0, cellularity_step = 0.01,
  distance_function = "MAD", keep_all = TRUE, max_distance_from_whole_number = 0.2
)

ploidy <- 4.19
cellularity <- 0.78
absolute_segments <- mutate(ss_segs, copy_number = relative_to_absolute_copy_number(copy_number, ploidy, cellularity))
ss_cns <- ss_cns %>% mutate(position = round((start + end) / 2))
absolute_copy_number <- mutate(ss_cns, across(c(copy_number, segmented), relative_to_absolute_copy_number, ploidy, cellularity))
chromosomes <- chromosome_offsets(absolute_copy_number)
genomic_copy_number <- convert_to_genomic_coordinates(absolute_copy_number, "position", chromosomes)
genomic_segments <- convert_to_genomic_coordinates(absolute_segments, c("start", "end"), chromosomes)
p <- genome_copy_number_plot(genomic_copy_number, genomic_segments, chromosomes,
                             min_copy_number = 0, max_copy_number = 10,
                             copy_number_breaks = 0:10,
                             point_colour = "grey40",
                             ylabel = "absolute copy number") +
  ggplot2::ggtitle(paste0(sample, '  ploidy:', ploidy, '  cellularity:', cellularity)) +
  ggplot2::theme(plot.title = element_text(size = 20))
p
