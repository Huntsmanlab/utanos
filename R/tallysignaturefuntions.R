tally_components_helper <- function(min_value, max_value, feature_label, feature_data) {
  if (feature_label == "point") {
    feature_data %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value == min_value, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else if (feature_label == "range") {
    feature_data %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(as.numeric(.data$value) > as.numeric(min_value) & as.numeric(.data$value) <= as.numeric(max_value), na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else {
    stop()
  }
}

tally_components <- function(df, feature, cn_feature_setting) {
  unique_samples <- unique(df$ID)
  df <- df %>% dplyr::as_tibble()
  df$ID <- factor(df$ID, levels = unique_samples)
  df$value <- as.numeric(df$value)

  if (feature == "segsize") {
    df$value <- log10(df$value)
  }

  features_subset <- cn_feature_setting[cn_feature_setting$feature == feature,]

  features_subset %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$component) %>%
    dplyr::summarize(
      count = paste(tally_components_helper(.data$min, .data$max, .data$label, df), collapse = ","),
    ) %>%
    tidyr::separate(.data$count, unique_samples, sep = ",", convert = TRUE) %>%
    data.table::as.data.table()
}

#' @export
GetCNTally <- function(CN_data, genome = 'hg19', relative=FALSE, extra_features= FALSE, subset_features = NULL){
  if(relative){
    if(extra_features){
      feature_setting <- utanos::feature_relative_extra
    }
    feature_setting <- utanos::feature_relative_default
  }
  else{
    if(extra_features){
      feature_setting <- utanos::feature_absolute_extra
    }
    feature_setting <- utanos::feature_absolute_default

  }

  cn_components <- purrr::map2(CN_data, names(CN_data),
                               tally_components,
                               cn_feature_setting = feature_setting
  )

  cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
    dplyr::as_tibble() %>%
    tibble::column_to_rownames(var = "component") %>%
    as.matrix()

  # Order the matrix as feature_setting
  cn_matrix <- cn_matrix[feature_setting$component, ] %>%
    t()

  if (!is.null(subset_features)) {
    pattern <- paste(subset_features, collapse = "|")

    cn_matrix <- cn_matrix[, grepl(pattern, colnames(cn_matrix))]
  }

  return(cn_matrix)
}

#' @export
PlotCNTallySig <- function(signatures, extra_features = FALSE){
  sig_matrix <- NMF::basis(signatures)
  sig_matrix <- as.data.frame(sig_matrix)
  sig_matrix <- sig_matrix[grep("\\[|\\]", rownames(sig_matrix)),]

  colnames(sig_matrix) <- paste("Signature", 1:ncol(sig_matrix))

  prefixes <- gsub("\\[.*\\]", "", rownames(sig_matrix))
  suffixes <- gsub(".*\\[(.*)\\]", "\\1", rownames(sig_matrix))

  # Add a new column for row name prefixes
  sig_matrix$prefix <- prefixes
  sig_matrix$suffix <- suffixes

  sig_long <- tidyr::gather(sig_matrix, key = "signature", value = "value", -prefix, -suffix)

  suffix_order <- c(
    "0", "1", "<=2", "2", "3", "4",
    ">4 & <=8", "5", ">5",
    "6", "7", "8", "9", "10",
    ">10 & <=20", ">20 & <=30", ">30",
    ">4 & <=10", ">10",
    ">2 & <=3", ">3 & <=4",
    ">4 & <=5", ">5 & <=6",
    ">6 & <=7", ">7 & <=8",
    ">7", ">8",
    "<= -1",
    ">-1 & <=-0.9",
    ">-0.9 & <=-0.8", ">-0.8 & <=-0.7",
    ">-0.7 & <=-0.6", ">-0.6 & <=-0.5",
    ">-0.5 & <=-0.4", ">-0.4 & <=-0.3",
    ">-0.3 & <=-0.2", ">-0.2 & <=-0.1",
    ">-0.1 & <=0", "<=0",
    ">0 & <=0.1", ">0.1 & <=0.2",
    ">0.2 & <=0.3", ">0.3 & <=0.4",
    ">0.4 & <=0.5", ">0.5 & <=0.6",
    ">0.6 & <=0.7", ">0.7 & <=0.8",
    ">0.8 & <=0.9", ">0.9 & <=1",
    ">1", "> 1"
  )

  sig_long$suffix <- factor(sig_long$suffix, levels = suffix_order)

  sig_long <- sig_long %>%
    dplyr::group_by(signature, prefix) %>%
    dplyr::mutate(value = value / sum(value)) %>%
    dplyr::ungroup()

  ggplot2::ggplot(sig_long, ggplot2::aes(x = suffix, y = value, fill = prefix)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(x = "Components", y = "Contributions", fill = "Signature") +
    ggplot2::facet_grid(signature ~ prefix, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(breaks = unique(sig_long$suffix)) +
    ggplot2::theme(
      legend.position = "none",
      axis.line = ggplot2::element_line(size = 0.3, colour = "black"),
      panel.spacing.x = ggplot2::unit(0.1, "line"),
      strip.background.y = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(
        color = "black",
        face = "bold",
        size = 12
      ),
      strip.text.y = ggplot2::element_text(
        size = 12,
        vjust = 1,
        color = "black",
        face = "bold",
        angle = 0
      ),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10)  # Use prefix color for facet labels
    )
}
