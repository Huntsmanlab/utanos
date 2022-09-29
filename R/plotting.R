# This is a script containing plotting functions for shallow WGS analysis
# January 18th, 2022

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
#' @export
plot_signature_exposures <- function (signatures, save_path = FALSE,
                                      obj_name = 'sig_exposures_obj',
                                      order = FALSE, transpose = FALSE) {
  long_data <- gather(signatures)
  long_data$max_sig <- rep(apply(output1, 2, function(x) which.max(x)),
                           times = 1,
                           each = nsigs)
  long_data$sigs <- rep(1:nsigs,dim(output1)[2])
  colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
  long_data <- long_data %>% arrange(max_sig)
  long_data$X <- factor(long_data$X, levels = unique(long_data$X))

  if (order != FALSE) {
    long_data$X <- factor(long_data$X, levels = order)
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
    scale_y_discrete(limits = paste0('S', 1:dim(signatures)[1])) +
    ggtitle("Signature exposures called per sample")
  g

  if (transpose != FALSE) {
    g <- g + theme(plot.title = element_text(size = 20),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_text(size = 18),
                   axis.title.x = element_text(size = 18),
                   axis.text.y = element_blank(),
                   legend.title = element_text(size = 18)) + coord_flip()
  }
  if (save_path != FALSE) {
    if (transpose != FALSE) {
      ggsave(paste0(save_path, "signatures_heatmap_", obj_name,".png"), plot = g, width = 10, height = 15)
    } else {
      ggsave(paste0(save_path, "signatures_heatmap_", obj_name,".png"), plot = g, width = 15, height = 10)
    }
  }

  output <- list(plot = g, ordering = levels(long_data$X))

  return(output)
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
#		(FALSE or char vector)	order:	If false used the sorted ordering of the heatmap, if true use new ordering.
#
# Return: list of ggplot and the ordering of the heatmap in sample names.
###
#' @export
swgs_cnv_heatmaps <- function(reads = data.frame(), save_path = FALSE,
                              obj_name = "gyne_cancer_obj",
                              order = FALSE, ret_order = FALSE) {

  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19

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

  if (save_path != FALSE) {
    ggsave(paste0(save_path, "cnv_heatmap_", obj_name,".png"), plot = g, width = 24, height = 24)
  }
  if (ret_order != FALSE) {
    output <- list(plot = g, ordering = levels(slice$sample_id))
  }

  return(output)
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
#' @export
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
#' @export
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

#' @export
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


### Plot metadata values in heatmap form (essentially make annotation bars)
# DESCRIPTION
# Parameters:
#   (DF)      met_basic: Basic metadata dataframe
#   (DF)      met_bio: Biologic metadata
#   (DF)      met_quality: Quality metrics metadata dataframe
#   (DF)      met_vafs: VAFs and Cellularities per sample
#
# Return:
#   (ggplot)  output: A cowplot combination of several heatmap bars
###
#' @export
plot_metadata_heatmaps <- function(met_basic, met_bio = FALSE, met_quality = FALSE, met_vafs = FALSE) {


  stopifnot("The first argument must be a dataframe and should contain
            basic metadata pertaining to your samples." = is.data.frame(met_basic),
            "The basic metadata DF must contain a column for the names
            of the samples titled 'sample_id'." = 'sample_id' %in% colnames(test))


  # Massage data into a single dataframe

  metadata <- data.frame(sample_id = met_basic$sample_id)

  if (is.data.frame(met_bio)) {
    met_bio$ihc_data_present <- is.na(met_bio$`CCNE1 H score`) & is.na(met_bio$`GRB7 IHC H score`)
    met_bio <- met_bio %>% dplyr::rename(sample_id = other_aliases,
                                         CCNE1_H = `CCNE1 H score`,
                                         GRB7_H = `GRB7 IHC H score`,
                                         BABAM_H = `BABAM IHC H Score`) %>%
                           dplyr::mutate(HER2_H = suppressWarnings(as.numeric(`HER2 IHC int`) * as.numeric(`HER2 % pos`)))
    met_bio$normality_classification <- 'abnormal'
    met_bio$normality_classification[met_bio$notes == 'looks pretty normal'] <- 'looks normal'
    met_bio$normality_classification[met_bio$notes == 'normal looking but with odd half-ploidy changes'] <- 'near normal'
    met_bio$CCNE1_H <- minmax_column(met_bio$CCNE1_H)
    met_bio$GRB7_H <- minmax_column(met_bio$GRB7_H)
    met_bio$BABAM_H <- minmax_column(met_bio$BABAM_H)
    met_bio$HER2_H <- minmax_column(met_bio$HER2_H)

    met_bio <- met_bio %>% dplyr::select(sample_id, ihc_data_present, BABAM_H,
                                         normality_classification,
                                         CCNE1_H, GRB7_H, HER2_H)
    metadata <- metadata %>% dplyr::left_join(met_bio, by = c('sample_id'))
  }

  if (is.data.frame(met_quality)) {
    met_quality <- met_quality %>% dplyr::select(sample_id, nreads, percent_coverage,
                                         mean_X_coverage, std_dev, mad, contains('quality_'))
    metadata <- metadata %>% dplyr::left_join(met_quality, by = c('sample_id'))
  }

  if (is.data.frame(met_vafs)) {
    celstab <- met_vafs

    # To simplify matters lets sum all unique vafs for a given gene/sample and treat them as a single vaf
    met_vafs[,8:13] <- as.data.frame(apply(met_vafs[,8:13],
                                           MARGIN = c(1,2),
                                           function(x) sum_delimited_elements(x, ';')))
    # Select maximum vaf
    celstab$max_vaf <- as.numeric(NA, rep(dim(celstab)[1]))
    for (i in 1:dim(met_vafs)[1]) {
      vafs <- met_vafs[i,]
      if (is.na(vafs$sample_id[1])) next                                          # If there isn't a vaf for the sample skip finding an ACN
      vafs <- vafs %>% dplyr::select(-c(sample_id, batch, tissue, sample_type, cancer_type, status, cellularity))
      vaf_idx <- which.max(as.double(vafs[1,]))                                   # Choose max vaf
      if (length(vaf_idx) == 0) next                                              # If there isn't a vaf for the sample skip finding an ACN
      celstab[celstab$sample_id == met_vafs$sample_id[i],]$max_vaf <- as.double(vafs[1, ..vaf_idx])
    }
    celstab$max_vaf[celstab$max_vaf > 100] <- 100
    celstab$minmaxed_cel <- minmax_column(celstab$cellularity)
    celstab$minmaxed_vaf <- minmax_column(celstab$max_vaf)

    browser()
    # test
    xx <- met_vafs %>% mutate(maxvaf = apply(X = met_vafs, MARGIN = 1, function(x) select_max_element(x)))

    celstab <- celstab %>% dplyr::select(sample_id, minmaxed_cel, minmaxed_vaf, max_vaf)
    metadata <- metadata %>% dplyr::left_join(celstab, by = c('sample_id'))
  }

  browser()
  metadata$sample_id <- factor(metadata$sample_id, levels = metadata$sample_id)
  metadata$swgs <- c('swgs_sample')
  # # ETL
  # p1 <- ggplot(metadata, aes(sample_id, swgs)) +
  #   geom_tile(aes(fill = sd_dev)) +
  #   scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1) ) +
  #   theme(axis.text.x = element_blank(),
  #         # axis.text.y = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         legend.position="none") +
  #   labs(title = 'std.dev') + coord_flip()
  # p2 <- ggplot(metadata, aes(sample_id, swgs)) +
  #   geom_tile(aes(fill = mad)) +
  #   scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1) ) +
  #   theme(axis.text.x = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         legend.position="none") +
  #   labs(title = 'med.abs.dev') + coord_flip()
  # p3 <- ggplot(metadata, aes(sample_id, swgs)) +
  #   geom_tile(aes(fill = max_vaf)) +
  #   scale_fill_gradientn(colors = c('blue', 'red'), values = c(0, 0.6, 1), na.value = "white") +
  #   theme(axis.text.x = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         legend.position="none") +
  #   labs(title = 'max. VAF') + coord_flip()
  #
  # p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 1, rel_widths = c(1.5,1,1,1,1,1,1,1))
  # ggsave2(plot = p, filename = '~/Downloads/quality_heatmaps.pdf', width = 10, height = 25)
}

select_max_element <- function (element_vector) {
  element_vector <- element_vector[grepl('vaf', names(element_vector), fixed = TRUE)]
  max_value <- max(as.numeric(element_vector), na.rm = TRUE)
  return(max_value)
}

# Adapted from CGHcall or QDNAseq package
# What role does this function have to play?
#
#' @export
testSummaryPlot <- function (x, main='Summary Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  nclass <-3
  if (!is.null(probamp(x))) nclass <- nclass+1
  if (!is.null(probdloss(x))) nclass <- nclass+1

  chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)

  if(nclass==3) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))}
  if(nclass==4) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}
  if(nclass==5) {loss.freq <- rowMeans(probloss(x))+rowMeans(probdloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}

  # remove probabilities of bins that fall below 0.2
  loss.freq[loss.freq < 0.2] <- 0.001
  gain.freq[gain.freq < 0.2] <- 0.001
  browser()
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='mean probability', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
}

#'
#' @description Function to plot the copy number profile from relative copy numbers 
#' @details Vosualizing relative copy number calls from WiseCondorX and QDNASeq as heatmaps 
#' @param x *dataframe* containing the copy number segments in long format: columns should be start, sample id (named sample) and segmented
#' @param order_by optional *character vector* containing the order the samples in the heatmap should be in 
#' @param order_by optional *character vector* to plot a subset of samples
#' @param breaks *vector* containing the limits of your colour gradient. Passed to scale_fill_gradientn
#' @param limits *vector* containing the limits of your colour gradient. Passed to scale_fill_gradientn
#' @return Heatmap of relative copy number calls
#' 
#' 
rCNplotProfile <- function(x, order_by = NULL, cluster = TRUE, subset = NULL, limits = c(-1, 1), breaks = c(-2.5, -1, -0.75, 0, 0.75, 1)) {
  if(isTRUE(cluster)) {
    sort_heatmap(x)
  }
  if(sum(is.na(x) > 0)) {
    x <- na.omit(x)
    warning("Removed NA values")
  }
  if(!is.empty(order_by)) {
    x$sample<- factor(x$sample, level = order)
  }
  if(!is.empty(subset)) {
    x <- x %>%
      filter(sample_id %in% subset)
  }
  
  p <- ggplot(x, aes(start, sample, 
                     fill = segmented)) + geom_tile() + 
    facet_grid(~ chromosome, scales = "free", space = "free", switch = "x") +
    scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
    theme(panel.spacing = unit(0.1, "lines")) + scale_fill_gradientn(colours = c( "blue", "grey", "red"), 
                                                                     breaks = breaks, limits = limits + theme_bw() + 
                                                                       guides(fill = guide_legend(override.aes = list(size = 5))) + labs(fill = "Relative Copy Number")
                                                                     return(p)
                                                                     
}

#' @description Updated QDNASeq unction to plot the segments for relative copy number calls from WiseCondorX
setMethod("plot", signature(x="QDNAseqSignals", y="missing"),
          function (x, y, main=NULL, includeReadCounts=TRUE,
                    logTransform=TRUE, scale=TRUE, sdFUN="sdDiffTrim",
                    delcol=getOption("QDNAseq::delcol", "darkred"),
                    losscol=getOption("QDNAseq::losscol", "red"),
                    gaincol=getOption("QDNAseq::gaincol", "blue"),
                    ampcol=getOption("QDNAseq::ampcol", "darkblue"),
                    pointcol=getOption("QDNAseq::pointcol", "black"),
                    segcol=getOption("QDNAseq::segcol", "chocolate"),
                    misscol=getOption("QDNAseq::misscol", NA),
                    pointpch=getOption("QDNAseq::pointpch", 1L),
                    pointcex=getOption("QDNAseq::pointcex", 0.1),
                    xlab=NULL, ylab=NULL, ylim=NULL, xaxt="s", yaxp=NULL,
                    showDataPoints=TRUE, showSD=TRUE, doSegments=TRUE, doCalls=TRUE, ...,
                    verbose=getOption("QDNAseq::verbose", TRUE)) {
            
            oopts <- options("QDNAseq::verbose"=verbose)
            on.exit(options(oopts))
            
            ## Import private functions
            ns <- asNamespace("CGHbase")
            .getChromosomeLengths <- get(".getChromosomeLengths", envir=ns, mode="function")
            .makeSegments <- get(".makeSegments", envir=ns, mode="function")
            
            if (inherits(x, c("QDNAseqCopyNumbers", "QDNAseqReadCounts"))) {
              condition <- QDNAseq:::binsToUse(x)
            } else {
              condition <- rep(TRUE, times=nrow(x))
            }
            baseLine <- NA_real_
            doCalls <- "calls" %in% assayDataElementNames(x) & doCalls
            doSegments <- "segmented" %in% assayDataElementNames(x) & doSegments
            if (doCalls) {
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-5, 5)
                } else {
                  ylim <- c(-2, 4)
                }
            }
            if ("copynumber" %in% assayDataElementNames(x)) {
              copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~ratio), "ratio")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(-3, 5)
                } else {
                  ylim <- c(0, 4)
                }
              if (is.null(yaxp))
                yaxp <- c(ylim[1], ylim[2], ylim[2]-ylim[1])
              baseLine <- ifelse(logTransform, 0, 1)
            } else {
              copynumber <- assayDataElement(x, "counts")[condition, , drop=FALSE]
              if (is.null(ylab))
                ylab <- ifelse(logTransform, expression(log[2]~read~count),
                               "read count")
              if (is.null(ylim))
                if (logTransform) {
                  ylim <- c(0, max(QDNAseq:::log2adhoc(copynumber)))
                } else {
                  ylim <- range(copynumber)
                }
            }
            if (is.null(main))
              main <- sampleNames(x)
            if (includeReadCounts && "total.reads" %in% names(pData(x)))
              main <- paste(main, " (",
                            format(x$total.reads, trim=TRUE, big.mark=","), " reads)", sep="")
            if (length(ylab) == 1)
              ylab <- rep(ylab, times=ncol(x))
            all.chrom <- QDNAseq:::chromosomes(x)
            if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
              all.chrom <- as.character(all.chrom)
            chrom <- all.chrom[condition]
            uni.chrom <- unique(chrom)
            uni.chrom <- as.character(str_sort(uni.chrom, numeric = TRUE))
            chrom.num <- as.integer(factor(chrom, levels= uni.chrom, ordered=TRUE))
            uni.chrom.num <- unique(chrom.num)
            uni.chrom.num <- sort( uni.chrom.num )
            if (!scale) {
              pos <- pos2 <- 1:sum(condition)
              chrom.ends <- aggregate(pos,
                                      by=list(chromosome=chrom), FUN=max)$x
            } else {
              if (inherits(x, c("cghRaw", "cghSeg", "cghCall"))) {
                chrom.lengths <- .getChromosomeLengths("GRCh37")
              } else {
                all.chrom.lengths <- aggregate(bpend(x),
                                               by=list(chromosome=all.chrom), FUN=max)
                chrom.lengths <- all.chrom.lengths$x
                names(chrom.lengths) <- all.chrom.lengths$chromosome
              }
              pos <- as.numeric(bpstart(x)[condition])
              pos2 <- as.numeric(bpend(x)[condition])
              chrom.lengths <- chrom.lengths[uni.chrom]
              chrom.ends <- integer()
              cumul <- 0
              for (i in seq_along(uni.chrom)) {
                pos[chrom.num > uni.chrom.num[i]] <-
                  pos[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                pos2[chrom.num > uni.chrom.num[i]] <-
                  pos2[chrom.num > uni.chrom.num[i]] +
                  chrom.lengths[uni.chrom[i]]
                cumul <- cumul + chrom.lengths[uni.chrom[i]]
                chrom.ends <- c(chrom.ends, cumul)
              }
              names(chrom.ends) <- names(chrom.lengths)
            }
            if (length(uni.chrom) == 1) {
              xax <- pretty(c(0, chrom.lengths[uni.chrom]))
              xaxlab <- xax / 1e6L
              if (is.null(xlab))
                xlab <- paste0("chromosome ", uni.chrom, ", Mbp")
            } else {
              xax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
              xaxlab <- uni.chrom
              if (is.null(xlab))
                xlab <- "chromosome"
            }
            if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
              copynumber <- log2adhoc(copynumber, inv=TRUE)
            if (is.character(sdFUN) && sdFUN == "sdDiffTrim") {
              symbol <- quote(hat(sigma)[Delta^"*"])
            } else if (is.character(sdFUN) && length(grep("Diff", sdFUN)) == 1) {
              symbol <- quote(hat(sigma)[Delta])
            } else {
              symbol <- quote(hat(sigma))
            }
            sdFUN <- QDNAseq:::sdDiffTrim
            noise <- apply(copynumber, MARGIN=2L, FUN=sdFUN, na.rm=TRUE)
            if (logTransform)
              copynumber <- QDNAseq:::log2adhoc(copynumber)
            for (i in seq_len(ncol(x))) {
              QDNAseq:::vmsg("Plotting sample ", main[i], " (", i, " of ", ncol(x), ") ...",
                             appendLF=FALSE)
              cn <- copynumber[, i]
              if (doSegments) {
                segmented <- assayDataElement(x, "segmented")[condition, i]
                if (inherits(x, c("cghRaw", "cghSeg", "cghCall")))
                  segmented <- QDNAseq:::log2adhoc(segmented, inv=TRUE)
                if (logTransform)
                  segmented <- QDNAseq:::log2adhoc(segmented)
                segment <- .makeSegments(segmented, chrom)
              }
              if (doCalls) {
                losses <- probloss(x)[condition, i]
                gains <- probgain(x)[condition, i]
                if (!is.null(probdloss(x)))
                  losses <- losses + probdloss(x)[condition, i]
                if (!is.null(probamp(x)))
                  gains <- gains + probamp(x)[condition, i]
                par(mar=c(5, 4, 4, 4) + 0.2)
                plot(NA, main=main[i], xlab=NA, ylab=NA, las=1,
                     xlim=c(0, max(chrom.ends)), ylim=ylim, xaxs="i", xaxt="n",
                     yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs="i")
                lim <- par("usr")
                lim[3:4] <- c(0, 1)
                par(usr=lim)
                dticks <- seq(0, 1, by=0.2)
                axis(4, at=dticks, labels=NA, srt=270, las=1, cex.axis=1,
                     cex.lab=1, tck=-0.015)
                axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1,
                     cex.lab=1, line=-0.4, lwd=0)
                mtext("probability", side=4, line=2, cex=par("cex"))
                if (!is.na(misscol)) {
                  rect(0, -1, max(pos2), 1, col=misscol, border=NA)
                  rect(pos, -1, pos2, 1, col="white", border=NA)
                }
                rect(pos[segment[,2]], 0, pos2[segment[,3]], losses[segment[,2]],
                     col=losscol, border=losscol)
                if (!is.null(probdloss(x)))
                  rect(pos[segment[,2]], 0, pos2[segment[,3]],
                       probdloss(x)[condition, i][segment[,2]],
                       col=delcol, border=delcol)
                rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
                     col=gaincol, border=gaincol)
                if (!is.null(probamp(x)))
                  rect(pos[segment[,2]], 1, pos2[segment[,3]],
                       1-probamp(x)[condition, i][segment[,2]],
                       col=ampcol, border=ampcol)
                axis(3, at=pos[which(probamp(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                axis(1, at=pos[which(probdloss(x)[condition,i] >= 0.5)],
                     labels=FALSE, col=delcol, col.axis="black", srt=270, las=1,
                     cex.axis=1, cex.lab=1)
                box()
                lim[3:4] <- ylim
                par(usr=lim)
                points(pos, cn, cex=pointcex, col=pointcol, pch=pointpch)
              } else {
                plot(pos, cn, cex=pointcex, col=pointcol, main=main[i],
                     xlab=NA, ylab=NA, ylim=ylim, xaxt="n", xaxs="i", yaxs="i",
                     yaxp=yaxp, tck=-0.015, las=1, pch=pointpch)
              }
              mtext(text=xlab, side=1, line=2, cex=par("cex"))
              mtext(text=ylab[i], side=2, line=2, cex=par("cex"))
              abline(h=baseLine)
              abline(v=chrom.ends[-length(chrom.ends)], lty="dashed")
              if (!is.na(xaxt) && xaxt != "n") {
                axis(side=1, at=xax, labels=NA, cex=.2, lwd=.5, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015)
                axis(side=1, at=xax, labels=xaxlab, cex=.2, lwd=0, las=1,
                     cex.axis=1, cex.lab=1, tck=-0.015, line=-0.4)
              }
              if (doSegments) {
                for (jjj in seq_len(nrow(segment))) {
                  segments(pos[segment[jjj,2]], segment[jjj,1],
                           pos[segment[jjj,3]], segment[jjj,1], col=segcol, lwd=3)
                }
              }
              par(xpd=TRUE)
              amps <- cn
              amps[amps <= ylim[2]] <- NA_real_
              amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
              dels <- cn
              dels[dels >= ylim[1]] <- NA_real_
              dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
              points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
              points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
              if (doSegments) {
                amps <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  amps <- QDNAseq:::log2adhoc(amps)
                amps[amps <= ylim[2]] <- NA_real_
                amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
                dels <- assayDataElement(x, "segmented")[condition, i]
                if (logTransform)
                  dels <- QDNAseq:::log2adhoc(dels)
                dels[dels >= ylim[1]] <- NA_real_
                dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
                points(pos, amps, pch=24, col=segcol, bg=segcol, cex=0.5)
                points(pos, dels, pch=25, col=segcol, bg=segcol, cex=0.5)
              }
              par(xpd=FALSE)
              ### estimate for standard deviation
              if (showSD) {
                if (!is.na(x$expected.variance[i])) {
                  sdexp <- substitute(paste(E~sigma==e, ", ", symbol==sd),
                                      list(e=sprintf("%.3g", sqrt(x$expected.variance[i])),
                                           symbol=symbol, sd=sprintf("%.3g", noise[i])))
                } else {
                  sdexp <- substitute(symbol==sd,
                                      list(symbol=symbol, sd=sprintf("%.3g", noise[i])))
                }
                mtext(sdexp, side=3, line=0, adj=1, cex=par("cex"))
              }
              ### number of data points
              if (showDataPoints) {
                str <- paste(round(sum(condition) / 1000), "k x ", sep="")
                probe <- median(bpend(x)-bpstart(x)+1)
                if (probe < 1000) {
                  str <- paste(str, probe, " bp", sep="")
                } else {
                  str <- paste(str, round(probe / 1000), " kbp", sep="")
                }
                if (doSegments)
                  str <- paste(str, ", ", nrow(segment), " segments", sep="")
                mtext(str, side=3, line=0, adj=0, cex=par("cex"))
              }
              QDNAseq:::vmsg()
            }
            options("QDNAseq::plotLogTransform"=logTransform)
            options("QDNAseq::plotScale"=scale)
          })


