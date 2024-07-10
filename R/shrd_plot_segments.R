#' A plot illustrating the merging of segments by the SHRD algorithm.
#'
#' @param segments_initial A data frame of the initial gathered segments.
#' @param segments_pass_one A data frame of the gathered segments, after processing and merging in the first pass.
#' @param segments_pass_two A data frame of the gathered segments, after processing and merging in the second pass.
#' @param sample The sample name.
#' @param segment_line_size The line width for the initial segments. Defaults to 4; first-pass & second-pass segments will be 1/2 and 1/4 this width, respectively.
#'
#' @return A plot illustrating the initial, first-pass, and second-pass (final) segments.
#' @export
#'
PlotSegmentChanges <- function(segments_initial, segments_pass_one, segments_pass_two, sample, segment_line_size = 4) {

  segments_initial <- segments_initial %>%
    dplyr::select(chromosome, start, end, ratio, ratio_median) %>%
    dplyr::rename("copynumber" = "ratio", "segmented" = "ratio_median") %>%
    dplyr::rename_with(.fn = ~ paste0(sample, "_", .), .cols = c("copynumber", "segmented")) %>%
    dplyr::mutate(chromosome = dplyr::case_when(chromosome == "23" ~ "X",
                                                chromosome == "24" ~ "Y",
                                                .default = as.character(chromosome)))

  chromosome_lengths <- segments_initial %>%
    dplyr::group_by(chromosome) %>%
    dplyr::summarise(length = max(end)) %>%
    dplyr::mutate(chromosome = factor(chromosome,
                                      levels = unique(segments_initial$chromosome))) %>%
    dplyr::arrange(chromosome)

  chromosomes <- chromosome_lengths %>%
    dplyr::mutate(offset = dplyr::lag(cumsum(as.numeric(length)), default = 0)) %>%
    dplyr::mutate(start = offset + 1, end = offset + length) %>%
    dplyr::mutate(mid = offset + round(length / 2))

  offsets <- dplyr::select(chromosomes, chromosome, offset)

  segments_pass_one <- segments_pass_one %>%
    dplyr::select(chr, start, end, ratio_median) %>%
    dplyr::rename("chromosome" = "chr", "segmented" = "ratio_median") %>%
    dplyr::mutate(chromosome = dplyr::case_when(chromosome == "23" ~ "X",
                                                chromosome == "24" ~ "Y",
                                                .default = as.character(chromosome))) %>%
    dplyr::left_join(offsets, by = "chromosome") %>%
    dplyr::mutate(across(c(start, end), ~ . + as.numeric(offset))) %>%
    dplyr::select(-offset) %>%
    dplyr::mutate(segment_number = dplyr::row_number()) %>%
    dplyr::select(segment_number, start, end, segmented) %>%
    tidyr::pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
    dplyr::arrange(segment_number)

  segments_pass_two <- segments_pass_two %>%
    dplyr::select(chr, start, end, ratio_median) %>%
    dplyr::rename("chromosome" = "chr", "segmented" = "ratio_median") %>%
    dplyr::mutate(chromosome = dplyr::case_when(chromosome == "23" ~ "X",
                                                chromosome == "24" ~ "Y",
                                                .default = as.character(chromosome))) %>%
    dplyr::left_join(offsets, by = "chromosome") %>%
    dplyr::mutate(across(c(start, end), ~ . + as.numeric(offset))) %>%
    dplyr::select(-offset) %>%
    dplyr::mutate(segment_number = dplyr::row_number()) %>%
    dplyr::select(segment_number, start, end, segmented) %>%
    tidyr::pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
    dplyr::arrange(segment_number)

  segments_initial_qdna <- DfToQDNAseq(df = segments_initial)
  initial_plot <- CNSegmentsPlot(cnobj = segments_initial_qdna, sample = sample,
                                 segment_colour = "darkgreen",
                                 segment_alpha = 1,
                                 segment_line_size = segment_line_size)

  final_plot <- initial_plot +
    ggplot2::geom_line(data = segments_pass_one,
                       mapping = ggplot2::aes(x = position, y = segmented, group = segment_number),
                       colour = "turquoise", alpha = 1, size = segment_line_size / 2) +
    ggplot2::geom_line(data = segments_pass_two,
                       mapping = ggplot2::aes(x = position, y = segmented, group = segment_number),
                       colour = "green", alpha = 1, size = segment_line_size / 4) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Segment"))

  return(final_plot)
}
