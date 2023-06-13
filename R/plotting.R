# This is a script containing plotting functions for shallow WGS analysis
# January 18th, 2022


PlotAbsCopyNumber <- function(x, sample, sample_segments, sample_cns, output) {
  ploidy <- as.numeric(x["ploidy"])
  cellularity <- as.numeric(x["cellularity"])
  chr_order <- c(as.character(1:22), 'X')
  absolute_segments <- dplyr::mutate(sample_segments,
                                     copy_number = relative_to_absolute_copy_number(copy_number, ploidy, cellularity))
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

# Adapted from CGHcall or QDNAseq package
# What role does this function have to play?
#
#' @export
TestSummaryPlot <- function (x, main='Summary Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
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

