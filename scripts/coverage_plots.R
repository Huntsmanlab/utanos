suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
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
