# This is a script containing plotting functions (and helpers) for shallow WGS analysis

# List of functions:
# PlotSignatureExposures
# SwgsCnvHeatmaps
# AddViralPresence <- does this still serve a purpose? Is it needed?
# SortHeatmap
# PlotCellularityAndMaxVaf <- too dataset specific, don't see a clear generalization, should likely be removed...
# PlotBABAM <- too dataset specific, don't see a clear generalization, should likely be removed...
# PlotMetadataHeatmaps <- too dataset specific, don't see a clear generalization, should likely be removed...
# SelectMaxElement
# SummaryCNPlot <- Should be converted to a ggplot whenever convenient
# RCNDiversityPlot
# RelativeCNSegmentsPlot <- made redundant? Can be removed soon, just test CNSegmentsPlot on WisecondorX output to make sure
# CNSegmentsPlot
# AddGenesToPlot
# QualityPlot

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
PlotSignatureExposures <- function (signatures, save_path = FALSE,
                                      obj_name = 'sig_exposures_obj',
                                      order = FALSE, transpose = FALSE) {
  nsigs <- nrow(signatures)
  long_data <- tidyr::gather(signatures)
  long_data$max_sig <- rep(apply(signatures, 2, function(x) which.max(x)),
                           times = 1,
                           each = nsigs)
  long_data$sigs <- rep(1:nsigs,dim(signatures)[2])
  colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
  long_data <- long_data %>% dplyr::arrange(max_sig)
  long_data$X <- factor(long_data$X, levels = unique(long_data$X))

  if (!isFALSE(order)) {
    long_data$X <- factor(long_data$X, levels = order)
  }

  g <- ggplot2::ggplot(long_data, ggplot2::aes(X, Y, fill=Z)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(discrete=FALSE) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
          axis.ticks.x = ggplot2::element_blank(),
          # axis.text.x = element_text(size = 15, angle = 75, vjust = 0.5, hjust=0.5),
          axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 18),
          axis.text.y = ggplot2::element_text(size = 18),
          legend.title = ggplot2::element_text(size = 18)) +
    ggplot2::labs(fill = 'Signature \nExposure', x = "Samples", y = " ") +
    ggplot2::scale_y_discrete(limits = paste0('S', 1:dim(signatures)[1])) +
    ggplot2::ggtitle("Signature exposures called per sample")
  g

  if (transpose != FALSE) {
    g <- g + ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 18),
                   axis.title.x = ggplot2::element_text(size = 18),
                   axis.text.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = 18)) +
      ggplot2::coord_flip()
  }
  if (save_path != FALSE) {
    if (transpose != FALSE) {
      ggplot2::ggsave(paste0(save_path, "/signatures_heatmap_", obj_name,".png"), plot = g, width = 10, height = 15)
    } else {
      ggplot2::ggsave(paste0(save_path, "/signatures_heatmap_", obj_name,".png"), plot = g, width = 15, height = 10)
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
SwgsCnvHeatmaps <- function(reads = data.frame(), save_path = FALSE,
                              obj_name = "gyne_cancer_obj",
                              order = FALSE, ret_order = FALSE) {

  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = RemoveBlacklist(reads)                                                # remove blacklist regions from hg19

  reads <- reads[, c("chromosome", "start", "sample_id", "state")]

  # Generate standard CNV heatmap data
  reads$state[reads$state < 0] <- 0
  slice <- SortHeatmap(reads)
  slice$state <- round(slice$state)
  slice$state[slice$state >= 15] <- 15
  slice$state <- as.character(slice$state)
  slice$state[slice$state == '15'] <- '15+'                                     # Assign all copy-numbers greater than 15 to the same value
  slice <- slice %>% dplyr::mutate(state = factor(state,
                                           levels = c(0:14, "15+"))) %>%
                     dplyr::arrange(chromosome)

  # Set colours for the heatmap
  cols <- c("#3182BD", "#9ECAE1", "#CCCCCC", "#FDD49E", "#FE9F4F", "#EA8251",
            "#E1571A", "#B33015", "#972913", "#621408", "#430227", "#730343",
            "#A61568", "#CE4191", "#CE8FB4", "#F5B1D9")                            # #FFDBF0
  names(cols) <- c(0:14, "15+")

  if (!isFALSE(order)) {
    slice$sample_id <- factor(slice$sample_id, level = order)
  }

  g <- ggplot2::ggplot(slice, ggplot2::aes(start, sample_id, fill = state)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(na.value = "white", values = cols, name = 'Copy-Number \nColour Code') +
    ggplot2::facet_grid(~chromosome, scales = "free", space = "free", switch = "x") +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.1, "lines"),
          plot.title = ggplot2::element_text(size = 28, hjust = 0.5),
          axis.text.y = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(size = 22),
          legend.title = ggplot2::element_text(size = 22),
          legend.text = ggplot2::element_text(size = 18)) +
    ggplot2::labs(x = "Chromosomes", y = "Samples", fill = "Copy Number") +
    ggplot2::ggtitle("Absolute copy number calls across the genome")

  output <- g

  if (save_path != FALSE) {
    # ggplot2::ggsave(paste0(save_path, "/cnv_diversity_heatmap_", obj_name,".png"), plot = g, width = 24, height = 24)
    png(paste0(save_path, "/cnv_diversity_heatmap_", obj_name,".png"), width=20, height=20, units = 'in', res = 400, type = 'cairo-png')
    print(g)
    dev.off()
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
AddViralPresence <- function(data, viral_presence) {
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
SortHeatmap <- function(slice) {
  slice$chromosome <- factor(slice$chromosome, levels = c("V", 1:22, "X", "Y"))
  slice <- dplyr::mutate(slice, pos = paste0(chromosome, ":", start))
  wide <- tidyr::spread(slice[, c("sample_id", "pos", "state")], pos, state)
  wide_ <- wide[,-which(colnames(wide) == "sample_id")]
  rownames(wide_) <- wide$sample_id
  cluster <- hclust(dist(wide_), method = "ward.D")
  slice$sample_id <- factor(slice$sample_id, level = rownames(wide_)[cluster$ord])
  return(slice)
}


### Plot cellularity and VAF values
# DESCRIPTION
# Parameters:
# slice: As input we want the order of the samples in the heatmap
#
#' @export
PlotCellularityAndMaxVaf <- function(reads = data.frame(), save_path, obj_name) {

  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19
  reads <- reads[, c("chromosome", "start", "sample_id", "state")]
  reads$state[reads$state < 0] <- 0
  slice <- SortHeatmap(reads)

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
PlotBABAM <- function(reads = data.frame(), save_path, obj_name) {

  if(nrow(reads) == 0) { stop("No reads data provided.") }                      # If reads DF is empty, return
  reads = removeBlacklist(reads)                                                # remove blacklist regions from hg19
  reads <- reads[, c("chromosome", "start", "sample_id", "state")]
  reads$state[reads$state < 0] <- 0
  slice <- SortHeatmap(reads)

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
PlotMetadataHeatmaps <- function(met_basic, met_bio = FALSE, met_quality = FALSE, met_vafs = FALSE) {


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
    xx <- met_vafs %>% mutate(maxvaf = apply(X = met_vafs, MARGIN = 1, function(x) SelectMaxElement(x)))

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

SelectMaxElement <- function (element_vector) {
  element_vector <- element_vector[grepl('vaf', names(element_vector), fixed = TRUE)]
  max_value <- max(as.numeric(element_vector), na.rm = TRUE)
  return(max_value)
}

#' Generate Summary Copy-Number Aberrations Plot
#'
#' @description
#' This function was copied and modified from another package (CGHbase).
#' The original code can be found here on github: https://github.com/tgac-vumc/CGHbase
#' The bioconductor page: https://bioconductor.org/packages/release/bioc/html/CGHbase.html
#'
#' This function creates a summary plot of the copy-number changes of the samples in the provided QDNAseq object.
#'
#' @param x An S4 object of type QDNAseqCopyNumbers containing multiple samples.
#' @param main (optional) String. Plot Title.
#' @param summarytype (optional) String. One of either 'probability' or 'frequency'. \cr
#' 'probability' - The vertical bars represent the average probability that the positions along the chromosome they cover are gained (blue bars) or lost (red bars) across all samples. \cr
#' 'frequency' - The vertical bars represent the frequency that the positions along the chromosome they cover are gained (blue bars) or lost (red bars) across all samples. \cr
#' @param maskprob (optional) Numeric. Gain or loss probabilities that fall below this number are masked. Used to exclude noise.
#' @param maskaberr (optional) Numeric. Copy-Number gains or losses that fall below this number are masked. Used to exclude noise.
#'
#' @export
SummaryCNPlot <- function (x, main='Relative Copy-Number Summary Plot',
                           summarytype = 'probability',
                           maskprob = 0.2, maskaberr = 0.1,
                           gaincol='blue', losscol='red', misscol=NA,
                           build='GRCh37',... ) {
  chrom <- QDNAseq::chromosomes(x)
  pos <- QDNAseq::bpstart(x)
  pos2 <- QDNAseq::bpend(x)
  uni.chrom <- unique(chrom)
  cns <- log2(CGHbase::segmented(x))
  com_cns <- as.data.frame(cns) %>% dplyr::mutate(mean = rowMeans(dplyr::across(where(is.numeric))))
  yaxis_label <- 'mean probability'

  nclass <-3
  if (!is.null(CGHbase::probamp(x))) nclass <- nclass+1
  if (!is.null(CGHbase::probdloss(x))) nclass <- nclass+1

  chrom.lengths <- GetChromosomeLengths(build)[as.character(uni.chrom)]

  # Convert the segment positions from a relative, per-chromosome basis to absolute positions in the genome
  abs_seg_pos <- RelToAbsSegPos(chromosomes = chrom, rel_start_pos = pos, rel_end_pos = pos2, build = build)
  pos <- abs_seg_pos$abs_start_pos
  pos2 <- abs_seg_pos$abs_end_pos
  chrom.ends <- abs_seg_pos$chrom_ends
  names(chrom.ends) <- names(chrom.lengths)

  if (nclass==3) {
    loss.freq <- rowMeans(CGHbase::probloss(x))
    gain.freq <- rowMeans(CGHbase::probgain(x))
  }
  if (nclass==4) {
    loss.freq <- rowMeans(CGHbase::probloss(x))
    gain.freq <- rowMeans(CGHbase::probgain(x)) + rowMeans(CGHbase::probamp(x))
  }
  if (nclass==5) {
    loss.freq <- rowMeans(CGHbase::probloss(x)) + rowMeans(CGHbase::probdloss(x))
    gain.freq <- rowMeans(CGHbase::probgain(x)) + rowMeans(CGHbase::probamp(x))
  }

  # Make Frequency plot option instead (rather than probability)
  if (summarytype == 'frequency') {
    calls <- CGHbase::calls(x)
    loss.freq <- rowMeans(calls < 0)
    gain.freq <- rowMeans(calls > 0)
    yaxis_label <- summarytype
  }

  # Mask gain/loss probability bins corresponding to CN aberrations that fall between maskaberr and zero
  loss.freq[abs(com_cns$mean) < maskaberr] <- 0.001
  gain.freq[abs(com_cns$mean) < maskaberr] <- 0.001

  # remove probabilities of bins that fall below maskprob
  loss.freq[loss.freq < maskprob] <- 0.001
  gain.freq[gain.freq < maskprob] <- 0.001

  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes',
       ylab=yaxis_label, xaxs='i', xaxt='n', yaxs='i', yaxt='n',
       main=main,...)
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
  uni.chrom.labels <-uni.chrom
  uni.chrom.labels[seq(2, length(uni.chrom.labels), 2)] <- "" # Make every second x-axis label empty, in order to fit labels neatly on the plot
  axis(side=1,at=ax,labels=uni.chrom.labels,
      lwd=.5,las=1,cex.axis=0.7,cex.lab=1)
  axis(side=2, at=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
       labels=c('100 %', '75 %', '50 %', '25 %', '0 %', '25 %', '50 %', '75 %', '100 %'),
       las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(QDNAseq::bpend(x)-QDNAseq::bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
  ### number of samples
  nsamps <- paste0('n=', as.character(dim(x)[2]))
  masks <- paste0('masks: prob=', as.character(maskprob), ', aberr=',
                  as.character(maskaberr))
  mtext(paste0(masks, '  |  ', nsamps), side=3, line=0, adj=1)
}


#' Plot relative CN diversity heatmap
#'
#' @description Function to plot the copy number profile from relative copy numbers
#' @details Visualizing relative copy number calls from WiseCondorX and QDNASeq as heatmaps
#' @param x *dataframe* containing the copy number segments in long format: columns should be start, sample id (named sample) and segmented
#' @param order_by optional *character vector* containing the order the samples in the heatmap should be in
#' @param cluster optional *logical* Do you want to do a row-wise sort of the  heatmap (sort the samples)
#' @param subset optional *character vector* A subset of samples that you'd like to plot.
#' @param breaks *vector* containing the limits of your colour gradient. Passed to scale_fill_gradientn
#' @param limits *vector* containing the limits of your colour gradient. Passed to scale_fill_gradientn
#' @return Heatmap of relative copy number calls
#'
#' @export
RCNDiversityPlot <- function(qdnaseq_obj, order_by = NULL, cluster = TRUE,
                             subset = NULL,
                             Xchr = FALSE,
                             limits = c(-1.5, 1.5),
                             breaks = c( 2.5, 1, 0.5, 0, -0.5, -1, -2.5)) {

  # Extract and re-shape data
  wide_binwise_segments <- ExportBinsQDNAObj(qdnaseq_obj, type = 'segments',
                                             filter = FALSE)
  long_segs <- wide_binwise_segments %>%
                        tidyr::gather(sample, segmented,
                                      5:dim(.)[2], factor_key=TRUE)
  segments <- CopyNumberSegments(long_segs) %>%
                    dplyr::select(sample, chromosome, start, end, copy_number) %>%
                    dplyr::rename(segVal = copy_number)
  long_segs <- SegmentsToCopyNumber(segments, 1000000,
                                    genome = 'hg19',
                                    Xincluded = Xchr)

  if (isTRUE(cluster)) {
    colnames(long_segs) <- c("chromosome", "start", "end", "state", "sample_id")
    long_segs <- SortHeatmap(long_segs)
    colnames(long_segs) <- c("chromosome", "start", "end", "segmented", "sample_id", "pos")
  }
  if (sum(is.na(long_segs) > 0)) {
    long_segs <- na.omit(long_segs)
    warning("Bins with NA values removed. May or may not be an issue. Up to the user.")
  }
  if (!is.null(order_by)) {
    long_segs$sample_id <- factor(long_segs$sample_id,
                                           levels = order_by)
  }
  if (!is.null(subset)) {
    long_segs <- long_segs %>% filter(sample_id %in% subset)
  }

  long_segs$segmented <- log2(long_segs$segmented)
  long_segs <- long_segs[!(long_segs$chromosome == 'Y'), ]
  long_segs$chromosome <- factor(long_segs$chromosome,
                                 level = c('1', '2', '3', '4', '5', '6',
                                           '7', '8', '9', '10', '11', '12',
                                           '13', '14', '15', '16', '17',
                                           '18', '19', '20', '21', '22',
                                           'X'))

  p <- ggplot2::ggplot(long_segs, ggplot2::aes(start, sample_id, fill = segmented)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(~ chromosome, scales = "free",
                     space = "free", switch = "x") +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::theme(panel.spacing = grid::unit(0.05, "lines")) +
    ggplot2::scale_fill_gradientn(colours = c( "blue", "grey", "red"),
                               breaks = breaks,
                               limits = limits) +
    ggplot2::theme_bw() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::labs(fill = "Relative Copy Number") +
    ggplot2::xlab('chromosomes') +
    ggplot2::ylab('samples')

  return(p)
}

#' Plot Relative Copy Numbers
#'
#' @description
#' Just an updated version of the \pkg{QDNAseq} package's plot() method.
#' This version just allows for relative copy number calls from WiseCondorX too.
#' All else remains the same.
#' See the <[`QDNAseq plotting`][QDNAseq::plot]> method for details.
#' @details
#' Supply a QDNAseq object.
#'
#' @export
RelativeCNSegmentsPlot <- function (x, main=NULL, includeReadCounts=TRUE,
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
                            showDataPoints=TRUE, showSD=TRUE, doSegments=TRUE,
                            doCalls=TRUE, ...,
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
  uni.chrom <- as.character(stringr::str_sort(uni.chrom, numeric = TRUE))
  chrom.num <- as.integer(factor(chrom, levels = uni.chrom, ordered=TRUE))
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
      all.chrom.lengths <- aggregate(QDNAseq::bpend(x),
                                     by=list(chromosome=all.chrom), FUN=max)
      chrom.lengths <- all.chrom.lengths$x
      names(chrom.lengths) <- all.chrom.lengths$chromosome
    }
    pos <- as.numeric(QDNAseq::bpstart(x)[condition])
    pos2 <- as.numeric(QDNAseq::bpend(x)[condition])
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
      losses <- CGHbase::probloss(x)[condition, i]
      gains <- CGHbase::probgain(x)[condition, i]
      if (!is.null(CGHbase::probdloss(x)))
        losses <- losses + CGHbase::probdloss(x)[condition, i]
      if (!is.null(CGHbase::probamp(x)))
        gains <- gains + CGHbase::probamp(x)[condition, i]
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
             CGHbase::probdloss(x)[condition, i][segment[,2]],
             col=delcol, border=delcol)
      rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
           col=gaincol, border=gaincol)
      if (!is.null(probamp(x)))
        rect(pos[segment[,2]], 1, pos2[segment[,3]],
             1-CGHbase::probamp(x)[condition, i][segment[,2]],
             col=ampcol, border=ampcol)
      axis(3, at=pos[which(CGHbase::probamp(x)[condition,i] >= 0.5)],
           labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1,
           cex.axis=1, cex.lab=1)
      axis(1, at=pos[which(CGHbase::probdloss(x)[condition,i] >= 0.5)],
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
      probe <- median(QDNAseq::bpend(x)-QDNAseq::bpstart(x)+1)
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
}


#' Segment Copy Number Plot across the Genome
#'
#' Whole genome copy number plot for one sample with partitioning lines between chromosomes.
#'
#' The \code{copynumber} values passed into this function can be on the relative
#' or absolute scales, or be log2 ratios. However, the same scale should be used
#' consistently in each slot of \code{cnobj} (and the \code{copy_number_steps} dataframe as the case may be).
#'
#' @param cnobj A *QDNAseq* object contains \code{copynumber}, and \code{segmented} assay slots. \cr
#' This function expects copy-number data that has been segmented. May contain multiple samples.
#' @param sample A *character string* indicating the sample. \cr
#' Required if the object passed to \code{cnobj} contains data for multiple samples. \cr
#' Ex. \code{c('sample1')} \cr
#' Must match one of the sample names exactly.
#' @param copy_number_steps a data frame containing \code{absolute_copy_number}
#' and \code{copy_number} columns.
#' @param max_points_to_display maximum number of copy number points to display
#' (downsampling may be carried out if there are more copy number values than
#' this number).
#' @param highlight_masks A *character string* indicating any feature data slots
#' containing logical masks over the bins. Each slot will be given its own colour. \cr
#' Keep in mind that NA values will be filtered out and should masks have overlap
#' the last one will take precedence. \cr
#' Ex. \code{c('comCNV.mask', 'centro.telo.mask')}
#' @param min_copy_number The minimum \code{copy_number} to display.
#' @param max_copy_number The maximum \code{copy_number} to display.
#' @param copy_number_breaks Breaks at which grid lines will be displayed.
#' @param def_point_colour A *character string* indicating the colour for the
#' default copy number points.
#' @param point_alpha The transparency of the copy number points.
#' @param point_size The size of the copy number points.
#' @param point_shape The point shape. See the ggplot2::geom_point() docs and
#' the \code{shape} parameter for details.
#' @param segment_colour The colour of the copy number segments.
#' @param segment_alpha The transparency of the copy number segments.
#' @param segment_line_size The size of the lines for copy number segments.
#' @param copy_number_step_colour The colour of the copy number step lines.
#' @param copy_number_step_alpha The transparency of the copy number step lines.
#' @param copy_number_step_line_size The size of the lines for the copy number.
#' steps.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return A \code{ggplot} object.
#' @examples
#'
#' @export
CNSegmentsPlot <- function(cnobj,
                           sample = NULL,
                           copy_number_steps = NULL,
                           max_points_to_display = Inf,
                           highlight_masks = NULL,
                           min_copy_number = NULL, max_copy_number = NULL,
                           ylims = NULL,
                           copy_number_breaks = NULL,
                           dolog2 = NULL,
                           def_point_colour = "black", point_alpha = 0.6,
                           point_size = 1, point_shape = 16,
                           segment_colour = "red", segment_alpha = 1,
                           segment_line_size = 1.5,
                           copy_number_step_colour = "blue", copy_number_step_alpha = 0.35,
                           copy_number_step_line_size = 0.75,
                           xlabel = "chromosome", ylabel = "copy number") {

  stopifnot(is(cnobj, "QDNAseqCopyNumbers"))
  stopifnot("chromosome" %in% names(cnobj@featureData@data))
  stopifnot("start" %in% names(cnobj@featureData@data), is.numeric(cnobj@featureData@data$start))
  stopifnot("end" %in% names(cnobj@featureData@data), is.numeric(cnobj@featureData@data$end))
  stopifnot("copynumber" %in% names(cnobj@assayData), is.numeric(cnobj@assayData$copynumber))
  stopifnot("segmented" %in% names(cnobj@assayData), is.numeric(cnobj@assayData$segmented))

  # compute offsets and genome coordinates for each chromosome
  chr_ord <- c(as.character(1:22), 'X', 'Y')
  chromosome_lengths <- cnobj@featureData@data %>%
    dplyr::group_by(chromosome) %>%
    dplyr::summarise(length = max(end)) %>%
    dplyr::mutate(chromosome = factor(chromosome,
                                      levels = chr_ord)) %>%
    dplyr::arrange(chromosome)
  chromosomes <- chromosome_lengths %>%
    dplyr::mutate(offset = dplyr::lag(cumsum(as.numeric(length)), default = 0)) %>%
    dplyr::mutate(start = offset + 1, end = offset + length) %>%
    dplyr::mutate(mid = offset + round(length / 2))
  offsets <- dplyr::select(chromosomes, chromosome, offset)

  # pull data out into data frames and filter
  wide_cns <- ExportBinsQDNAObj(cnobj, type = "copynumber", filter = FALSE)
  wide_cns$colouring <- 'normal.bins'
  wide_cns <- wide_cns %>% dplyr::relocate(colouring, .after = feature)
  if (!is.null(highlight_masks)) {
    for (i in highlight_masks) {
      wide_cns$colouring[!cnobj@featureData@data[[i]]] <- i
    }
  }
  long_cns <- wide_cns %>% tidyr::gather(sample, copy_number,
                                         6:dim(.)[2], factor_key=TRUE)
  wide_segs <- ExportBinsQDNAObj(cnobj, type = "segments", filter = FALSE)
  long_segs <- wide_segs %>% tidyr::gather(sample, segmented,
                                           5:dim(.)[2], factor_key=TRUE)

  if(!is.null(dolog2)) {
    long_segs$segmented[long_segs$segmented < 0] <- 0
    long_segs$segmented <- log2(long_segs$segmented)
    long_cns$copy_number[long_cns$copy_number < 0] <- 0
    long_cns$copy_number <- log2(long_cns$copy_number)
  }

  if (!is.null(sample)) {
    stopifnot("sample" %in% names(long_cns))
    stopifnot("sample" %in% names(long_segs))
    copy_number <- dplyr::filter(long_cns,
                                 sample == !!sample, is.finite(copy_number))
    segments <- dplyr::filter(long_segs,
                              sample == !!sample, is.finite(segmented))
  }

  # collapse bin-wise segment values down to proper segments
  segments <- CopyNumberSegments(segments) %>%
    dplyr::select(sample, chromosome, start, end, copy_number)

  # compute mid-point position for the copy number bins
  copy_number <- dplyr::mutate(copy_number, position = (start + end) / 2)

  # convert to genome coordinates
  copy_number <- copy_number %>%
    dplyr::left_join(offsets, by = "chromosome") %>%
    dplyr::mutate(position = position + as.numeric(offset)) %>%
    dplyr::select(-offset)

  segments <- segments %>%
    dplyr::left_join(offsets, by = "chromosome") %>%
    dplyr::mutate(across(c(start, end), ~ . + as.numeric(offset))) %>%
    dplyr::select(-offset)

  # filter bins outside intended plotting range
  if (is.null(min_copy_number)) {
    min_copy_number <- min(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(min_copy_number),
              length(min_copy_number) == 1,
              !is.na(min_copy_number))
    copy_number <- dplyr::filter(copy_number, copy_number >= min_copy_number)
    segments <- dplyr::filter(segments, copy_number >= min_copy_number)
  }

  if (is.null(max_copy_number)) {
    max_copy_number <- max(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(max_copy_number),
              length(max_copy_number) == 1,
              !is.na(max_copy_number))
    copy_number <- dplyr::filter(copy_number, copy_number <= max_copy_number)
    segments <- dplyr::filter(segments, copy_number <= max_copy_number)
  }

  # optional down-sampling of bins to improve speed
  stopifnot(is.numeric(max_points_to_display),
            length(max_points_to_display) == 1,
            !is.na(max_points_to_display))
  if (max_points_to_display < nrow(copy_number))
    copy_number <- dplyr::sample_n(copy_number, max_points_to_display)

  segment_lines <- segments %>%
    dplyr::mutate(segment_number = dplyr::row_number()) %>%
    dplyr::select(segment_number, start, end, copy_number) %>%
    tidyr::pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
    dplyr::arrange(segment_number)

  # Only plot the chromosomes needed (in case any have been filtered out/removed/excluded)
  chromosomes <- chromosomes[chromosomes$chromosome %in%
                               unique(copy_number$chromosome), ]
  xmin <- min(chromosomes$start)
  xmax <- max(chromosomes$end)

  # sort out point colouring if mask highlighting is requested
  if (!is.null(highlight_masks)) {
    viridis_colors <- viridis::viridis(length(highlight_masks), option = "turbo")
    custom_colors <- c("normal.bins" = def_point_colour, viridis_colors)
    names(custom_colors) <- c("normal.bins", highlight_masks)
  } else {
    custom_colors = c(def_point_colour)
  }

  # Create the plots
  plot <- ggplot2::ggplot(data = copy_number,
                          mapping = ggplot2::aes(x = position, y = copy_number)) +
    ggplot2::geom_vline(xintercept = chromosomes$end, colour = "grey90")

  if (!is.null(copy_number_steps)) {
    xmax <- xmin + (xmax - xmin) * 1.04

    copy_number_steps <- copy_number_steps %>%
      dplyr::filter(copy_number >= min_copy_number,
                    copy_number <= max_copy_number) %>%
      dplyr::arrange(dplyr::desc(absolute_copy_number))

    if (nrow(copy_number_steps) > 0) {
      plot <- plot +
        ggplot2::geom_hline(yintercept = copy_number_steps$copy_number,
                            colour = copy_number_step_colour,
                            alpha = copy_number_step_alpha,
                            size = copy_number_step_line_size) +
        ggplot2::geom_label(data = copy_number_steps,
                            mapping = ggplot2::aes(x = xmin + 0.98 * (xmax - xmin),
                                                   y = copy_number,
                                                   label = absolute_copy_number)) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
    }
  }

  if (is.null(copy_number_breaks)) { copy_number_breaks = ggplot2::waiver() }

  plot <- plot +
    ggplot2::geom_point(ggplot2::aes(colour = colouring),
                        alpha = point_alpha,
                        size = point_size,
                        pch = point_shape) +
    ggplot2::geom_line(data = segment_lines,
                       mapping = ggplot2::aes(x = position,
                                              y = copy_number,
                                              group = segment_number),
                       colour = segment_colour,
                       alpha = segment_alpha,
                       size = segment_line_size) +
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_x_continuous(limits = c(xmin, xmax),
                                expand = ggplot2::expansion(mult = 0),
                                breaks = chromosomes$mid,
                                labels = chromosomes$chromosome) +
    ggplot2::scale_y_continuous(limits = c(min_copy_number,
                                           max_copy_number),
                                breaks = copy_number_breaks,
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::labs(x = xlabel, y = ylabel) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(size = 0.2),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )

  if (!is.null(ylims)) {
    plot <- plot + ggplot2::ylim(ylims)
  }

  if (is.null(highlight_masks)) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(copy_number_steps) && nrow(copy_number_steps) > 0) {
    plot <- plot +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  return(plot)
}


### Given chromosomal segment positions, convert the relative (to a given chromosome) start and end positions of each segment to absolute positions in the genome.
# DESCRIPTION
# Parameters:
#   (character vector)  chromosomes: A vector of length(segments) indicating the chromosome each segment is on, as returned by QDNAseq::chromosomes
#   (numeric vector)    rel_start_pos: A vector of length(segments) indicating the relative start position of each segment, as returned by QDNAseq::bpstart
#   (numeric vector)    rel_end_pos: A vector of length(segments) indicating the relative end position of each segment, as returned by QDNAseq::bpend
#   (string)            build: The reference genome build to use when retrieving chromosome lengths.
#
# Returns:
#   (list)              abs: A list containing the absolute start/end positions and the total chromosome length.
#   (numeric vector)    abs_start_pos: A vector of length(segments) indicating the absolute start position of each segment
#   (numeric vector)    abs_end_pos: A vector of length(segments) indicating the absolute end position of each segment
#   (integer)           tot_length: The total length of all the chromosomes
#   (numeric vector)    chrom_ends: The absolute end position of each chromosome
###
#' @export
RelToAbsSegPos <- function(chromosomes, rel_start_pos, rel_end_pos, build = "GRCh37") {
  chrom_lengths <- GetChromosomeLengths(build)[as.character(unique(chromosomes))]

  chromosomes <- replace(chromosomes, chromosomes == "X", "23")
  chromosomes <- as.integer(replace(chromosomes, chromosomes == "Y", "24")) # Map all chromosomes to integers so we can iterate in numeric order
  uni_chrom <- unique(chromosomes)

  abs_start_pos <- rel_start_pos
  abs_end_pos <- rel_end_pos
  tot_length <- 0
  chrom_ends <- integer()

  for (chrom in uni_chrom) {
    # Since a one-to-one correspondence exists between the vector indicating the chromosome each segment is on, and the vectors of segment start/end positions, we can index the
    # start/end positions based on the chromosomes which are numbered higher than the current chromosome. We then iteratively add the length of lower numbered chromosomes
    # to the relative start/end positions, obtaining absolute positions.
    abs_start_pos[chromosomes > chrom] <- abs_start_pos[chromosomes > chrom] + chrom_lengths[chrom]
    abs_end_pos[chromosomes > chrom] <- abs_end_pos[chromosomes > chrom] + chrom_lengths[chrom]

    tot_length <- tot_length + chrom_lengths[chrom]
    chrom_ends <- c(chrom_ends, tot_length)
  }

  names(chrom_ends) <- as.character(chrom_lengths)
  abs <- list(abs_start_pos = abs_start_pos, abs_end_pos = abs_end_pos, tot_length = tot_length, chrom_ends = chrom_ends)
  return(abs)
}

#' Label genes of interest on an absolute copy-number plot.
#'
#' @description
#' This function takes an existing absolute copy-number plot, a vector of gene names, and optionally an Ensembl DB to use for looking up gene annotations.
#' Ensure that you use an Ensembl DB which annotates the same genome version used for your existing plot.
#' A plot with the genes of interest labeled is returned.
#'
#' @param plot A ggplot object, corresponding to the existing absolute copy-number plot.
#' @param genes A vector of gene symbols you'd like to label on the plot -- use upper case.
#' @param edb The name of the local Ensembl DB package to use for looking up gene information. Defaults to [EnsDb.Hsapiens.v75], which annotates hg19/GRCh37.
#' @param ... Additional options controlling the visual display of the labels; passed to ggrepel::geom_label_repel(), and supports any of its arguments.
#'
#' @returns A ggplot object with the gene labels added.
#'
#' @seealso [ggrepel::geom_label_repel()] for supported arguments that alter the visual display of the labels.
#'
#' @export
AddGenesToPlot <- function(plot, genes, edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, ...) {
  # Get the locations of the genes we are interested in.
  genes_ens <- as.data.frame(ensembldb::genes(edb, filter = ~ gene_name %in% genes, return.type = "DataFrame"))

  # We will bin the genes by their midpoint.
  genes_ens <- genes_ens %>%
    dplyr::mutate(midpoint = ((gene_seq_start + gene_seq_end) / 2))

  # Put each gene in its corresponding bin. This lets us borrow the existing absolute bin position for graphing.
  plot$data <- plot$data %>%
    dplyr::left_join(y = genes_ens, by = dplyr::join_by(chromosome == seq_name, start < midpoint, end > midpoint))

  plot <- plot +
    ggrepel::geom_label_repel(aes(label = symbol), ...)

  return(plot)
}


#' Plot sample quality according to a specified metric.
#'
#' @param sample_quality_df The sample quality data-frame returned by GetSampleQualityDecision().
#' @param metric The quality metric passed to GetSampleQualityDecision(). One of quanitle, or DecisionTree.
#'
#' @return A ggplot object of the quality plot.
#' @export
QualityPlot <- function(sample_quality_df, metric = "quantile") {

  stopifnot(is.character(metric))

  if (metric == "quantile") {
    quality_plot <- ggplot2::ggplot(data = sample_quality_df, mapping = ggplot2::aes(x = seg_counts, y = median_sd, colour = decision)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Segment count", y = "Median segment SD", title = "Sample quality via quantile") +
      ggplot2::scale_colour_discrete(name = "Quality decision") +
      ggplot2::geom_vline(xintercept = sample_quality_df$seg_cut[1], linetype = "dashed", colour = "red", size = 0.5) +
      ggplot2::geom_hline(yintercept = sample_quality_df$med_cut[1], linetype = "dashed", colour = "red", size = 0.5) +
      ggplot2::theme_bw()
  }

  else if (metric == "DecisionTree") {
    quality_plot <- ggplot2::ggplot(data = sample_quality_df, mapping = ggplot2::aes(x = seg_counts, y = median_sd, colour = decision)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Segment count", y = "Median segment SD", title = "Sample quality via decision tree") +
      ggplot2::scale_colour_discrete(name = "Quality decision") +
      ggplot2::theme_bw()
  }

  return(quality_plot)
}

