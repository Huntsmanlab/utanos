# This is a script containing plotting functions for shallow WGS analysis
# January 18th, 2022

# Packages
suppressPackageStartupMessages({ 
  library(dplyr)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggalt)
  library(viridis)
  library(caret)
  library(ggpubr)
})

# Create Heatmap of Signature Exposures
#####
# Converts signature-per-sample data to a heatmap, saves it as a png, and returns the ggplot. 
# Samples are sorted for display based on copy number similarities in the state column.
# Supports addition of alternative state number for comparison purposes.
# Parameters:
#   (DF)      signatures:           csv file where columns are samples and rows are normalized signature exposures
#																		Must have at least 2 samples and 2 signatures
#   (char)    save_file:            Path where to save the plot. Otherwise 'FALSE' to not save.
#   (Boolean) order:                
#   (Boolean) transpose:            
#
# Return: ggplot
###
plot_signature_exposures <- function (signatures, save_file = FALSE, 
                                      order = FALSE, transpose = FALSE) {
  long_data <- gather(signatures)
  long_data$max_sig <- rep(apply(output1, 2, function(x) which.max(x)), times = 1, each = nsigs)
  long_data$sigs <- rep(1:nsigs,dim(output1)[2])
  colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
  long_data <- long_data %>% arrange(max_sig)
  long_data$X <- factor(long_data$X, levels = unique(long_data$X))
  
  if (order != FALSE) {
    slice$sample_id <- factor(slice$sample_id, level = order)
  }
  
  g <- ggplot(long_data, aes(X, Y, fill=Z)) +  
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    theme(plot.title = element_text(size = 20),
          axis.ticks.x = element_blank(),
          # axis.text.x = element_text(size = 15, angle = 75, vjust = 0.5, hjust=0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18), 
          legend.title = element_text(size = 18)) + 
    labs(fill = 'Signature \nExposure', x = "Samples", y = " ") + 
    scale_y_discrete(limits = c('S1', 'S2', 'S3', 'S4', 'S5', 'S6')) +
    ggtitle("Signature exposures called per sample")
  g
  ggsave(paste0(save_file,".png"), plot = g, width = 15, height = 10)
  
  rownames(output1) <- c(1:7)
  write.csv(output1, file = "~/repositories/cnsignatures/data/pancan_UCEC_output/signature_exposures_absolute_autosomesOnly_pancanUCEC_vanSignatures.csv")
  
  return(g)
}


# Create Heatmaps of CNV calls
#####
# Converts provided sample data to a Copy Number heatmap and saves it as a png. Samples are sorted for display based on copy number similarities in the state column.
# Supports addition of alternative state number for comparison purposes.
# Parameters:
#   (DF)      reads:                Reads data frame with an additional column for experimental condition
#																		Required fields: (chr, start, cell_id, state)
#   (String)  save_path:            Path where the PNG file will be saved (Default: cnv_heatmap.png)
#   (Boolean) low_map_mask:         Specify whether or not to filter regions with low mappability from the heatmap
#		(String)	view:									Specify alternative views for data (Default: "")
#																			altState: Show copy numbers specified by the alternative state provided
#																			compare:	Show copy number differences between the state and altenative state (state - altState)
# Return: Void
###
swgs_cnv_heatmaps <- function(reads = data.frame(), save_path = "cnv_heatmap", 
                              obj_name = "gyne_cancer_obj", 
                              low_map_mask = FALSE,
                              order = FALSE) {
  
  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19
  
  # Determine if experimental condition is included in the data and build corresponding dataframe
  exp = FALSE
  if("condition" %in% colnames(reads)) {
    exp_cond = distinct(reads[, c("cell_id", "condition")]) 
    exp = TRUE
  }
  
  # Overwrite copy number state if region has low mappability
  if (low_map_mask) { reads$state[reads$is_low_mappability == TRUE] <- 'low mappability'}
  
  reads <- reads[, c("chromosome", "start", "sample_id", "state")]
  
  # Generate standard CNV heatmap
  reads$state[reads$state < 0] <- 0
  slice <- sort_heatmap(reads)
  
  slice$state <- round(slice$state)
  slice$state[slice$state >= 15] <- 15
  slice$state <- as.character(slice$state)
  slice$state[slice$state == '15'] <- '15+'                                     # Assign all copy-numbers greater than 15 to the same value
  slice$state <- as.factor(slice$state)
  
  # Set colours for the heatmap
  cols <- c("#3182BD", "#9ECAE1", "#CCCCCC", "#FDD49E", "#FE9F4F", "#EA8251",  
            "#E1571A", "#B33015", "#972913", "#621408", "#430227", "#730343", 
            "#A61568", "#CE4191", "#CE8FB4", "#F5B1D9")                            # #FFDBF0
  names(cols) <- c(0:14, "15+")
  
  if (order != FALSE) {
    slice$sample_id <- factor(slice$sample_id, level = order)
  }
  
  g <- ggplot(slice, aes(start, sample_id, fill = state)) + geom_tile() + 
    scale_fill_manual(na.value = "white", values = cols, name = 'Copy-Number \nColour Code') + 
    facet_grid(~chromosome, scales = "free", space = "free", switch = "x") + 
    scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
    theme(panel.spacing = unit(0.1, "lines"), 
          plot.title = element_text(size = 28, hjust = 0.5), 
          axis.text.y=element_blank(), 
          axis.title = element_text(size = 22), 
          legend.title = element_text(size = 22), 
          legend.text = element_text(size = 18)) +
    labs(x = "Chromosomes", y = "Samples", fill = "Copy Number") + 
    ggtitle("Absolute copy number calls across the genome")
  
  ggsave(paste0(save_path, "cnv_heatmap_", obj_name,".png"), plot = g, width = 24, height = 24)
  return(g)
}

###
# Custom palette scale colour function to go with ggplot2 
# scm = function(palette=cols) {
#   scale_color_manual(values=palette, na.value="#000000")
# }
# 
# ggplot(mydata, aes(x, y)) + geom_line() + scm()




### Add viral presence pseudo-chromosome to CNV matrix
# Adds an extra chromosome that represents the binary presence of a viral insert
# Parameters:
#   data: Dataframe containing observations describing DLP data   
#   viral_presence: 
# Return: data parameter modified to include a viral presence indicator
add_viral_presence <- function(data, viral_presence) {
  viral_range <- seq(1,50000001, 500000)
  viral_presence <- viral_presence[viral_presence$cell_id %in% unique(data$cell_id), ]
  l <- length(viral_range)
  w <- dim(viral_presence)[1]
  viral_state <- data.frame(chr = rep('V', l*w), 
                            start = rep(viral_range, w), 
                            cell_id = rep(viral_presence$cell_id, 1, each = l), 
                            state = rep('TP53', l*w))
  data <- rbind(viral_state, data)
  return(data)
}


### Sort the CN Heatmap using simple hierarchical clustering
# DESCRIPTION
# Parameters:
#   
sort_heatmap <- function(slice) {
  slice$chromosome <- factor(slice$chromosome, levels = c("V", 1:22, "X", "Y"))
  slice <- dplyr::mutate(slice, pos = paste0(chromosome, ":", start))
  wide <- spread(slice[, c("sample_id", "pos", "state")], pos, state)
  rownames(wide) <- wide$sample_id
  wide$sample_id <- NULL
  cluster <- hclust(dist(wide), method = "ward.D")
  slice$sample_id <- factor(slice$sample_id, level = rownames(wide)[cluster$ord])
  return(slice)
}

### Plot cellularity and VAF values
# DESCRIPTION
# Parameters: 
# slice: As input we want the order of the samples in the heatmap
# 
plot_cellularity_and_maxvaf <- function(reads = data.frame(), save_path, obj_name) {
  
  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19
  reads <- reads[, c("chromosome", "start", "sample_id", "state")]
  reads$state[reads$state < 0] <- 0
  slice <- sort_heatmap(reads)
  
  vafscels <- data.table::fread(file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/vafs_and_cellularities.tsv', sep = '\t', header = TRUE)
  vafscels$sample_id <- str_replace_all(vafscels$sample_id, "-", ".")
  vafscels <- vafscels %>% dplyr::filter(sample_id %in% unique(slice$sample_id))
  # vafscels <- vafscels[3:171,]
  
  # Cellularity ETL section
  celstab <- data.frame(sample_id = levels(slice$sample_id))
  celstab <- celstab %>% left_join(vafscels, by = c('sample_id'))
  celstab$sample_id <- factor(levels(slice$sample_id), level = levels(slice$sample_id))
  preproc2 <- preProcess(as.data.frame(celstab$cellularity), method=c("range"))
  minmaxed <- predict(preproc2, as.data.frame(celstab$cellularity))
  celstab$minmaxed_cel <- minmaxed$`celstab$cellularity`
  
  # Plot cellularity colour bar
  p1 <- ggplot(celstab, aes(sample_id, sample_type)) + 
    geom_tile(aes(fill = minmaxed_cel)) + 
    scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1) ) + 
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
  
  # VAFs ETL section
  # To simplify matters lets sum all unique vafs for a given gene/sample and treat them as a single vaf
  vafscels[,8:13] <- as.data.frame(apply(vafscels[,8:13], 
                                         MARGIN = c(1,2), 
                                         function(x) sum_delimited_elements(x, ';')))
  celstab$max_vaf <- NA
  for (i in 1:dim(vafscels)[1]) {
    vafs <- vafscels[i,]
    if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN
    vafs <- vafs %>% dplyr::select(-c(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity))
    vaf_idx <- which.max(as.double(vafs[1,]))                                   # Choose max vaf
    if (length(vaf_idx) == 0) next                                              # If there isn't a vaf for the sample skip finding an ACN
    celstab[celstab$sample_id == vafscels$sample_id[i],]$max_vaf <- as.double(vafs[1, ..vaf_idx])
    vaf_name <- names(vafs[1, ..vaf_idx])                                         
    vaf_name <- toupper(str_split(vaf_name, pattern = '\\.')[[1]][1])           # vaf name
  }
  celstab$max_vaf[celstab$max_vaf > 100] <- 100
  preproc2 <- preProcess(as.data.frame(celstab$max_vaf), method=c("range"))
  minmaxed <- predict(preproc2, as.data.frame(celstab$max_vaf))
  celstab$minmaxed_vaf <- minmaxed$`celstab$max_vaf`
  
  # Plot vaf colour bar
  p2 <- ggplot(celstab, aes(sample_id, sample_type)) + 
    geom_tile(aes(fill = minmaxed_vaf)) + 
    scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1) ) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none")
  
  # Plot both cellularity and max. vaf together
  p <- ggpubr::ggarrange(p1 + coord_flip(), p2 + coord_flip(), ncol = 2)
  ggsave(paste0(save_path, "colourbars_vafscellularities_clustered_", obj_name, ".png"), plot = p, width = 20, height = 24)
  
  # Plot both cellularity and max. vaf together but as a line plot
  lineplotdata <- data.frame(sample_id = as.character(rep(celstab$sample_id, 2)), 
                             groups = c(rep("cellularities", dim(celstab)[1]), rep("vafs", dim(celstab)[1])), 
                             vals = c(celstab$minmaxed_cel, celstab$minmaxed_vaf))
  lineplotdata <- lineplotdata %>% dplyr::arrange(desc(vals))
  lineplotdata$sample_id <- factor(lineplotdata$sample_id, levels = (lineplotdata %>% dplyr::filter(groups == 'cellularities'))$sample_id)
  p3 <- ggplot(data=lineplotdata, aes(x=sample_id, y=vals, group=groups)) +
    geom_line(aes(color=groups))+
    geom_point(aes(color=groups)) +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),)
  
  ggsave(paste0(save_path, "lineplot_vafscellularities_", obj_name, ".png"), plot = p3, width = 16, height = 4)
  
  message("The Pearson correlation of the cellularity and max. vaf for this subset is:")
  cor(celstab$minmaxed_cel, celstab$minmaxed_vaf, method = 'p')
}

plot_BABAM <- function(reads = data.frame(), save_path, obj_name) {
  
  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19
  reads <- reads[, c("chromosome", "start", "sample_id", "state")]
  reads$state[reads$state < 0] <- 0
  slice <- sort_heatmap(reads)
  
  qual <- data.table::fread(file = '~/Downloads/swgs_sample_biologic_metadata.csv', sep = ',', header = TRUE)
  qual <- qual %>% dplyr::filter(status == 'p53_abn')
  qual$other_aliases <- str_replace_all(qual$other_aliases, "-", ".")
  qual <- qual %>% dplyr::filter(!is.na(qual$`BABAM IHC H Score`))
  qual <- qual %>% dplyr::filter(other_aliases %in% unique(slice$sample_id))
  
  # Cellularity ETL section
  celstab <- data.frame(sample_id = levels(slice$sample_id))
  celstab <- celstab %>% left_join(qual, by = c('sample_id' = 'other_aliases'))
  celstab$sample_id <- factor(levels(slice$sample_id), level = levels(slice$sample_id))
  
  # Plot BABAM status
  preproc2 <- preProcess(as.data.frame(celstab$`BABAM IHC H Score`), method=c("range"))
  minmaxed <- predict(preproc2, as.data.frame(celstab$`BABAM IHC H Score`))
  celstab$minmaxed_babam <- minmaxed[,1]
  p1 <- ggplot(celstab, aes(sample_id, sample_type)) + 
    geom_tile(aes(fill = minmaxed_babam)) + 
    scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1) ) + 
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="none") + coord_flip()
  
  ggsave(paste0(save_path, "colourbars_bioIndicator_clustered_", obj_name, ".png"), plot = p1, width = 20, height = 24)
}

