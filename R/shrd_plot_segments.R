#suppressMessages(library("DescTools"))
#suppressMessages(library("ggplot2"))
#suppressMessages(library("gridExtra"))
#suppressMessages(library("GenomicRanges"))
#suppressMessages(library("ks"))
#suppressMessages(library("ggrepel"))

#' Plots segments in the given data frame.
#'
#' @description
#' Plots segments in their respective chromosomes with the genome's middle positions
#' from each chromosome, obtained in `hg19_segments.R`
#'
#' @param gathered_by_ratio_median A data frame: the segments merged by the same
#' ratio_median and/or chromosome arm. Most importantly, the frame must have already been
#' prepped for initialization of Levels in `PrepForInitialization()`. This is also
#' marked as 'Graph 1' in `stitching_functions.R`. The frame is used to find the max and min
#' values of the ratio_median, which corresponds to the y-axis of the graph.
#' @param bam_ratios_frame A data frame: the cleaned-up version of the bam_ratios.txt file.
#' @param segments A data frame: segment data. Usually is Graph 5 in the stitching script.
#' @param chr_mid_positions An array: the middle positions for each chromosomes from 1 to 23, in that order.
#'
#' @export

PlotSegments <- function(gathered_by_ratio_median, bam_ratios_frame, segments, chr_mid_positions) {
  print(dim(gathered_by_ratio_median))
  print(dim(bam_ratios_frame))
  print(dim(segments))
  gathered_by_ratio_median <- gathered_by_ratio_median[which(gathered_by_ratio_median$chr != 23),]

  #### Setting the y-axis limits ####
  graph_lower_limit = -2
  graph_higher_limit = -2

  if (max(abs(gathered_by_ratio_median$ratio_median)) > 2) {
    graph_lower_limit = -max(abs(gathered_by_ratio_median$ratio_median))
    graph_higher_limit = max(abs(gathered_by_ratio_median$ratio_median))
  }

  #### Setting up merged_bam frame ####
  copy_bam_ratios_frame = bam_ratios_frame
  copy_bam_ratios_frame = copy_bam_ratios_frame[,-1]
  copy_bam_ratios_frame = copy_bam_ratios_frame[,-5]
  colnames(copy_bam_ratios_frame) <- c("chr", "start", "end", "readcount")
  attach(copy_bam_ratios_frame)
  merged_bam = merge(aggregate(start ~ chr, copy_bam_ratios_frame, min),
                     aggregate(end ~ chr, copy_bam_ratios_frame, max))[seq(from=1,to=22,by=1),]
  merged_bam <- as.data.frame(merged_bam) #df
  detach(copy_bam_ratios_frame)

  #### Creating centromere/mid-chromosome data frame ####
  chr_mid_positions <- chr_mid_positions[-c(13,14,15,21,22,23)]
  chr = c(1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,20)
  mid_pos_frame = data.frame(chr_mid_positions, chr)
  colnames(mid_pos_frame) <- c("start_centromere", "chr") #table_adding_centromere

  #### Setting up copy_bam_ratios_frame and copy_segments ####
  segments$size = segments$end - segments$start + 1

  copy_bam_ratios_frame = bam_ratios_frame # B
  copy_bam_ratios_frame <- copy_bam_ratios_frame[,-1]
  colnames(copy_bam_ratios_frame) <- c("chr", "start", "end", "ratio", "ratio_median")
  copy_bam_ratios_frame <- copy_bam_ratios_frame[which(copy_bam_ratios_frame$chr != 23),]

  copy_segments = segments
  copy_segments <- copy_segments[which(copy_segments$chr != 23),] # C

  closest_higlight = Closest(copy_bam_ratios_frame$start, 30302805)[1]

  higlight_CCNE1 = copy_bam_ratios_frame[copy_bam_ratios_frame$chr == 19 & copy_bam_ratios_frame$start == closest_higlight,]

  data.segm = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.4, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
  data.text = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.5, chr = 19, label = "CCNE1")

  if (higlight_CCNE1[1,5] >= 0.5){
    data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.85, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)
    data.text = data.frame(x=higlight_CCNE1[1,2], y=1.95, chr = 19, label = "CCNE1")
  }  else if(higlight_CCNE1[1,5] >= 1.7) {
    data.segm = data.frame(x=higlight_CCNE1[1,2], y=1.1, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] - 0.05, chr = 19)
    data.text = data.frame(x=higlight_CCNE1[1,2], y=1, chr = 19, label = "CCNE1")
  }

  graph <- ggplot() +
    geom_rect(data=merged_bam,
              aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill= chr %% 2 == 0)) +
    geom_point(data=copy_bam_ratios_frame,
               aes(x=start, y=ratio),
               size=0.1,
               color="grey60",
               na.rm=TRUE) +
    geom_segment(data=copy_segments,
                 aes(x=start, xend=end, y=ratio_median, yend=ratio_median),
                 size=3,
                 color="#990000",
                 na.rm=TRUE) +
    geom_segment(data=data.segm,
                 mapping=aes(x=x, y=y, xend=xend, yend=yend),
                 arrow=arrow(unit(0.30,"cm"), angle = 20),
                 size=0.6,
                 color="black",
                 inherit.aes=FALSE,
                 na.rm=TRUE) +
    ylim(graph_lower_limit, graph_higher_limit) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle("Segmentation") +
    geom_vline(data=mid_pos_frame,
               mapping=aes(xintercept=start_centromere),
               color="black",
               linetype="dotted") +
    geom_text(data=data.text,
              mapping=aes(x=x, y=y, label= label),
              size=3.3,
              inherit.aes=TRUE,
              na.rm=TRUE) +
    geom_point(data=higlight_CCNE1,
               aes(x=start, y=ratio_median),
               color="orange",
               size=2,
               na.rm=TRUE) +
    geom_text(aes(x=start,y=graph_higher_limit, hjust = "left", vjust = "top", label= chr),
              data=merged_bam,
              fontface="bold",
              size=4.5) +
    scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "white")) +
    facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")


  graph <- graph + theme(plot.title = element_text(hjust = 0.5),
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y = element_text(size=15),
                         panel.spacing = unit(0, "lines"),
                         strip.text.x = element_blank(),
                         line = element_blank(),
                         legend.position = "none",
                         panel.background = element_blank())

  graph = ggplotGrob(x = graph)
  graph$layout$clip = "off"
  suppressWarnings(ggsave("./final_segmentation.jpeg", plot = graph, device = "jpeg", width = 23, height = 13))
}
