suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(data.table)
})

# Declare paths
input_path <- '~/Documents/projects/cn_signatures_shallowWGS/sequencing_metrics/coverage/se/'
metadata <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata_unix.tsv'

# Read-in data
metadata <- read.table(metadata, header = T, sep = '\t', stringsAsFactors = F)
combined <- data.frame()
allfiles <- list.files(input_path, full.names = TRUE)
for (i in allfiles) {
  cov_tab <- read.table(i, header = F, stringsAsFactors = F)
  cov_tab <- cov_tab %>% dplyr::slice(1:23)
  combined <- bind_rows(combined, cov_tab)
}

# Re-format dataframes
colnames(combined) <- c('chr', 'start', 'end', 'numreads', 'covbases', 'percent_coverage', 'mean_coverage', 'meanbaseq', 'meanmapq')
allfiles <- list.files(input_path)
allfiles <- word(allfiles, 1, sep = "\\.")
combined <- combined %>% mutate(sample = rep(allfiles, 1, each = 23)) %>% relocate(sample, .before = chr)

metadata <- metadata %>% arrange(sample_id)
combined <- left_join(combined, metadata, by = c("sample" = "sample_id"), keep = FALSE) %>%
            dplyr::select(sample, chr, start, end, numreads, batch, percent_coverage, mean_coverage, status, hist, cancer_type)

# Plot per sample
# mean_prop_cov = on average, each base has been sequenced X number of times
per_sample <- combined %>% 
              group_by(sample) %>% 
              summarise(mean_per_cov = mean(percent_coverage), mean_prop_cov = mean(mean_coverage), read_count = sum(numreads))
per_sample <- per_sample %>%
              mutate(batch = factor(metadata$batch, levels=unique(as.character(sort(metadata$batch, decreasing=FALSE)))),
                     status = metadata$status,
                     hist = metadata$hist,
                     cancer_type = metadata$cancer_type) %>%
              arrange(batch)

p <- ggplot(per_sample, aes(x=batch, y=mean_per_cov, color = batch)) +
  geom_point(size=4) +
  scale_colour_viridis_d() +
  # theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, .1))) +
  ggtitle("Average chromosomal percent coverage per sample")
p

p <- ggplot(per_sample, aes(x=batch, y=mean_prop_cov, color = batch)) +
  geom_point(size=4) +
  scale_colour_viridis_d() +
  # theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits = c(0, 2), expand = expansion(mult = c(0, .1))) +
  ggtitle("On average each base has been sequenced X number of times per sample")
p
# per_sample <- per_sample%>% arrange(status)
# p <- ggplot(per_sample, aes(x=sample, y=mean_per_cov, color = status)) +
#   geom_point(size=4) +
#   # theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, .1))) +
#   ggtitle("Average chromosomal percent coverage per sample")
# p
# 
# p <- ggplot(per_sample, aes(x=sample, y=mean_prop_cov, shape = as.factor(batch), color = as.factor(batch))) +
#   geom_point(size=4) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   scale_y_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, .1))) +
#   ggtitle("On average each base has been sequenced X number of times per sample")
# p
# 
# per_sample <- per_sample%>% arrange(status)
# p <- ggplot(per_sample, aes(x=sample, y=mean_prop_cov, color = status)) +
#   geom_point(size=4) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   scale_y_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, .1))) +
#   ggtitle("On average each base has been sequenced X number of times per sample")
# p
###########
# Read depth
# p <- ggplot(per_sample, aes(x=sample, y=read_count)) +
#   geom_col(aes(fill = status)) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   # scale_y_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, .1))) +
#   ggtitle("Read count per sample")
# p

p <- ggplot(per_sample, aes(x=batch, y=read_count, color=batch)) +
  geom_point(size=3) +
  scale_colour_viridis_d() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits = c(0, 4e7)) +
  ggtitle("Read count per sample")
p

######################
######################
# Targeted panel seq. depth - histogram version
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/per_sample/CC-LAV-0630/bedtools_coverage_3.txt'
input1 <- data.table::fread(file = input_path, sep = '\t', fill = TRUE)
hist_data <- input %>% dplyr::filter(V1 == 'all') %>% dplyr::select(!c(V6,V7,V8))
# Targeted panel seq. depth - mean depth version
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/per_sample/CC-LAV-0630/bedtools_coverage_2.txt'
input2 <- data.table::fread(file = input_path, sep = '\t', fill = TRUE)
hist_data <- input %>% dplyr::filter(V1 == 'all') %>% dplyr::select(!c(V6,V7,V8))
# Targeted panel seq. depth - multicov. depth version
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/bedtools_multicoverage.txt'
input3 <- data.table::fread(file = input_path, sep = '\t', fill = TRUE)
# Targeted panel seq. depth - multicov. depth version
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/bedtools_multicoverage_mapqFiltered.txt'
input4 <- data.table::fread(file = input_path, sep = '\t', fill = TRUE)

# Read-in data
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/coverage_files/'
combined <- data.frame()
table_list <- list()
allfiles <- list.files(input_path, full.names = TRUE)
sample_names <- list.files(input_path, full.names = FALSE)
sample_names <- str_replace_all(sample_names, '_bedtools_coverage_2.txt', '')
for (i in 1:length(allfiles)) {
  cov_tab <- read.table(allfiles[i], header = F, stringsAsFactors = F)
  table_list[[sample_names[i]]] <- cov_tab
  png(filename = paste0("~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/coverage_plots/", sample_names[i], "_depth_histogram.png"), 
      width=2000, height=1800, pointsize = 30, units = 'px')
  hist(cov_tab$V5, breaks = 2000, main = paste0("Average count depth per region for " , 
                                                sample_names[i], " \n Axis cut-off at 3000 \n", 
                                                "Proportion of regions > 250: ", as.character(signif(sum(cov_tab$V5 > 250)/dim(cov_tab)[1], digits = 3)), " ; ", 
                                                " > 500: ", as.character(signif(sum(cov_tab$V5 > 500)/dim(cov_tab)[1], digits = 3)), " ; ",
                                                " > 1000: ", as.character(signif(sum(cov_tab$V5 > 1000)/dim(cov_tab)[1], digits = 3))), 
       xlim = c(0,3000), xlab = c('Depth'))
  dev.off()
}

# Filtered coverage histogram plots
input_path <- '~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/bedtools_multicoverage_mapqFiltered.txt'
cov_tab <- read.table(input_path, header = F, stringsAsFactors = F)
# sample_names <- c('CC-RJH0119', 'CC-RJH0132', 'CC-RJH0175', 'CC-RJH0177', 'CC-HAM-0374', 'CC-HAM-0454')
sample_names <- c('CC-HAM-0369', 'CC-HAM-0374', 'CC-HAM-0379', 'CC-HAM-0383', 'CC-HAM-0385', 
                  'CC-HAM-0387', 'CC-HAM-0392', 'CC-HAM-0394', 'CC-HAM-0396', 'CC-HAM-0402', 
                  'CC-HAM-0405', 'CC-HAM-0406', 'CC-HAM-0407', 'CC-HAM-0421', 'CC-HAM-0422', 
                  'CC-HAM-0423', 'CC-HAM-0424', 'CC-HAM-0427', 'CC-HAM-0429', 'CC-HAM-0432', 
                  'CC-HAM-0435', 'CC-HAM-0438', 'CC-HAM-0442', 'CC-HAM-0446', 'CC-HAM-0448', 
                  'CC-HAM-0449', 'CC-HAM-0450', 'CC-HAM-0453', 'CC-HAM-0454', 'CC-HAM-0458', 
                  'CC-HAM-0459', 'CC-HAM-0461', 'CC-JGH-0469', 'CC-JGH-0487', 'CC-JGH-0492', 
                  'CC-JGH-0499', 'CC-JGH-0500', 'CC-JGH-0506', 'CC-JGH-0508', 'CC-JGH-0510', 
                  'CC-JGH-0511', 'CC-JGH-0515', 'CC-JGH-0523', 'CC-JGH-0525', 'CC-LAV-0570', 
                  'CC-LAV-0577', 'CC-LAV-0578', 'CC-LAV-0586', 'CC-LAV-0587', 'CC-LAV-0590', 
                  'CC-LAV-0598', 'CC-LAV-0600', 'CC-LAV-0615', 'CC-LAV-0616', 'CC-LAV-0630', 
                  'CC-LAV-0646', 'CC-LAV-0649', 'CC-LAV-0660', 'CC-LAV-0664', 'CC-LAV-0671', 
                  'CC-LAV-0681', 'CC-LAV-0687', 'CC-LAV-0698', 'CC-LAV-0699', 'CC-LAV-0719', 
                  'CC-LAV-0725', 'CC-LAV-0726', 'CC-LGH-0979', 'CC-NSH-0243', 'CC-NSH-0253', 
                  'CC-NSH-0254', 'CC-NSH-0269', 'CC-NSH-0279', 'CC-NSH-0285', 'CC-NSH-0309', 
                  'CC-NSH-0337', 'CC-NSH-0338', 'CC-NSH-0351', 'CC-NSH-0356', 'CC-RJH0119', 
                  'CC-RJH0132', 'CC-RJH0175', 'CC-RJH0177', 'CC-RJH0178', 'CC-RJH0182', 
                  'CC-RJH0201', 'CC-SMH-0751')
for (i in 5:dim(cov_tab)[2]) {
  png(filename = paste0("~/Documents/projects/cn_sigs_swgs/targeted_panel_seq/coverage_plots/", sample_names[i-4], "_depth_histogram_filtered.png"), 
      width=2000, height=1800, pointsize = 30, units = 'px')
  hist(cov_tab[,i], breaks = 2000, main = paste0("Average count depth per region for " , 
                                                sample_names[i-4], " \n Axis cut-off at 3000 \n", 
                                                "Proportion of regions > 250: ", as.character(signif(sum(cov_tab[,i] > 250)/dim(cov_tab)[1], digits = 3)), " ; ", 
                                                " > 500: ", as.character(signif(sum(cov_tab[,i] > 500)/dim(cov_tab)[1], digits = 3)), " ; ",
                                                " > 1000: ", as.character(signif(sum(cov_tab[,i] > 1000)/dim(cov_tab)[1], digits = 3))), 
       xlim = c(0,3000), xlab = c('Depth'))
  dev.off()
}
