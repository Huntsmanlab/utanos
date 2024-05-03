# File contains plotting functions related to CN-Signature visualization

###########################
### Functions
###########################
# TwoFeatureScatterPlot
# MixtureModelPlots
# GaussiansMixturePlot
# PoissonsMixturePlot
# WassDistancePlot



#' Make a two-way scatterplot of two CN-features
#'
#' @description
#'
#' This function makes a two-way scatterplot of two CN-features.
#' It then draws contour lines for density.
#'
#'
#'
#' @export
TwoFeatureScatterPlot <- function(featA, featB) {

  pp <- 'test'

  return(pp)
}


#' Plot Mixture Model Components by Signature
#'
#' @description
#'
#' This function (MixtureModelPlots), for an indicated CN-Signature (S), returns a list of six ggplots.
#' Each plot visualizes the mixture models used for each feature.
#' All mixture components are plotted, but those elevated for S are shaded.
#'
#' @param signatures A dataframe with components along the y-axis and signatures along the x.
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param sig_of_interest (optional) A single integer corresponding the signature for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @returns A list of ggplots.
#'
#'
#' @export
MixtureModelPlots <- function(signatures, components, sig_of_interest = 1) {

  # Validate input
  stopifnot("Components argument must be a list of flexmix objects." =
              typeof(components) == 'list')
  stopifnot("Components argument must be a list of flexmix objects." =
              typeof(components[[1]]) == 'S4')
  stopifnot("The indicated signature (sig_of_interest parameter) must be an
            integer less than or equal to the # of columns in 'signatures'." =
              dim(signatures) >= sig_of_interest)

  all_plots <- list()
  all_plots[["segmentsize"]] <- GaussiansMixturePlot(signatures, components,
                                                     sig_of_interest,
                                                     component = 'segsize')
  all_plots[["breakpoint10MB"]] <- PoissonsMixturePlot(signatures, components,
                                                       sig_of_interest,
                                                       component = 'bp10MB')
  all_plots[["oscillating"]] <- PoissonsMixturePlot(signatures, components,
                                                    sig_of_interest,
                                                    component = 'osCN')
  all_plots[["changepoint"]] <- GaussiansMixturePlot(signatures, components,
                                                     sig_of_interest,
                                                     component = 'changepoint')
  all_plots[["copynumber"]] <- GaussiansMixturePlot(signatures, components,
                                                    sig_of_interest,
                                                    component = 'copynumber')
  all_plots[["breakpointsarm"]] <- PoissonsMixturePlot(signatures, components,
                                                       sig_of_interest,
                                                       component = 'bpchrarm')

  return(all_plots)
}


#' Plot a Mixture of Gaussians
#'
#' @description
#'
#' This function plots the mixture models for a component that is composed of gaussians.
#' All mixture components are plotted, but those elevated for the indicated signature are shaded.
#' \cr
#' Additionally, an inlay plot is created out of the 'important' components if... \cr
#' 1. Some components have been 'squashed' down \cr
#' 2. There are any. - i.e. at least 1 weight > 0.05 \cr
#'
#' @param signatures A dataframe with components along the y-axis and signatures along the x.
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param sig_of_interest (optional) A single integer or a vector of integers corresponding the signatures for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @param component (optional) Which component to draw the gaussian curves for. (options: segsize, changepoint, copynumber)
#' @param inlay_flag (optional) Boolean flag to turn on/off the inlay plot.
#' @returns A ggplot.
#'
#'
#' @export
GaussiansMixturePlot <- function(signatures, components,
                                 sig_of_interest = 1, component = 'segsize',
                                 inlay_flag = TRUE) {

  # Palette for plotting
  cbPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                 "#A6761D", "#666666", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "black")

  norm_const <- apply(signatures, 1, sum)
  sig_mat_norm <- data.frame(apply(signatures,
                                   2,
                                   function(x){x/norm_const}))
  weights <- as.numeric(sig_mat_norm[,sig_of_interest])

  # plotting prep work
  mask <- grepl(component, rownames(sig_mat_norm), fixed=TRUE)
  plotparam <- flexmix::parameters(components[[component]])
  # CalculateSumOfPosteriors reorders its matrix on output...
  # So the signature components are correspondingly in ascending order of their means...
  # So this reordering is then necessary to correspond to those weights
  plotparam <- plotparam[,order(plotparam[1,])]
  segpalette <- cbPalette[1:sum(mask)]
  shading <- weights[mask]
  max_comp <- plotparam[,which(plotparam[1,] == max(plotparam[1,]))]
  xmax <- max_comp[1] + (max_comp[2]*1.75)
  digits <- nchar(as.character(round(xmax/2)))
  plotbreaks <- c(10^(digits-1), (10^digits)/2, 10^digits)
  plotbreaks <- plotbreaks[xmax > plotbreaks]

  # Make the main plot
  main_plot <- ggplot2::ggplot(data = data.frame(x = c(1,xmax)), aes(x)) +
    ylab("") +
    xlab(paste0(component, " mixture components")) +
    ggplot2::theme_bw() +
    theme(axis.text = ggplot2::element_text(size = 10),
          axis.title = ggplot2::element_text(size = 12),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank()) +
    scale_x_continuous(breaks = plotbreaks)

  for (i in 1:ncol(plotparam)) {
    main_plot <- main_plot +
                    geom_area(stat = "function",
                              fun = dnorm,
                              args = list(mean = plotparam[1,i],
                                          sd = plotparam[2,i]),
                              color = "black",
                              size = 0.05,
                              fill = segpalette[i],
                              alpha = shading[i])
  }
  final_plot <- main_plot

  # Make an inlay plot of the 'important' components if...
  # 1. Some components have been 'squashed' down
  # 2. There are any. - i.e. at least 1 weight > 0.05
  if ((plotparam[1,dim(plotparam)[2]]/plotparam[1,1] > 20) &&
      (sum(shading > 0.05) >= 1)) {
    mask <- shading > 0.05
    plotparam2 <- plotparam[,mask, drop = FALSE]
    shading2 <- shading[mask]
    segpalette <- segpalette[mask]
    min_comp <- plotparam2[,which(plotparam2[1,] == min(plotparam2[1,]))]
    max_comp <- plotparam2[,which(plotparam2[1,] == max(plotparam2[1,]))]
    xmin <- min_comp[1] - (min_comp[2]*1.75)
    xmax <- max_comp[1] + (max_comp[2]*1.75)
    digits <- nchar(as.character(round(xmax/2)))
    plotbreaks <- c(10^(digits-1), (10^digits)/2, 10^digits)
    plotbreaks <- plotbreaks[xmax > plotbreaks]
    inlay_plot <- ggplot2::ggplot(data = data.frame(x = c(1, xmax)), aes(x)) +
                      ylab("") + xlab("") +
                      ggplot2::theme_bw() +
                      theme(axis.text = ggplot2::element_text(size = 8),
                            axis.title = ggplot2::element_text(size = 8),
                            panel.grid.minor = ggplot2::element_blank(),
                            panel.grid.major = ggplot2::element_blank()) +
                      scale_x_continuous(breaks = plotbreaks)

    for (i in 1:ncol(plotparam2)) {
      inlay_plot <- inlay_plot + geom_area(stat = "function",
                                           fun = dnorm,
                                           args = list(mean = plotparam2[1,i],
                                                       sd = plotparam2[2,i]),
                                           color = "black",
                                           size = 0.05,
                                           fill = segpalette[i],
                                           alpha = shading2[i])
    }
    # Overwrite earlier final plot with inset plot version
    final_plot <-
      cowplot::ggdraw() +
      cowplot::draw_plot(main_plot) +
      cowplot::draw_plot(inlay_plot, x = 0.5, y = .5, width = .45, height = .45)
  }

  return(final_plot)
}


#' Plot a Mixture of Poissons
#'
#' @description
#'
#' This function plots the mixture models for a component that is composed of poisson distributions.
#' All mixture components are plotted, but those elevated for the specified signature are shaded.
#'
#' @param signatures A dataframe with components along the y-axis and signatures along the x.
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param sig_of_interest (optional) A single integer or a vector of integers corresponding the signatures for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @param component (optional) Which component to draw the poisson 'curves' for. (options: bp10MB, osCN, bpchrarm)
#' @returns A ggplot.
#'
#'
#' @export
PoissonsMixturePlot <- function(signatures, components,
                                sig_of_interest = 1, component = 'bp10MB') {

  # Palette for plotting
  cbPalette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                 "#A6761D", "#666666", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "black")

  norm_const <- apply(signatures, 1, sum)
  sig_mat_norm <- data.frame(apply(signatures,
                                   2,
                                   function(x){x/norm_const}))
  weights <- as.numeric(sig_mat_norm[,sig_of_interest])
  mask <- grepl(component, rownames(sig_mat_norm), fixed=TRUE)
  plotparam <- flexmix::parameters(components[[component]])

  # CalculateSumOfPosteriors reorders its matrix on output...
  # So the signature components are correspondingly in ascending order of their means...
  # So this reordering is then necessary to correspond to those weights
  plotparam <- plotparam[order(plotparam)]
  shading <- weights[mask]
  max_comp <- plotparam[which(plotparam == max(plotparam))]
  xmax <- max_comp[1] + (sqrt(max_comp[1])*3)

  if (component == 'bp10MB') {
    xlabel <- 'breakpoints per 10MB'
  } else if (component == 'osCN') {
    xlabel <- 'oscillating CN regions'
  } else if (component == 'bpchrarm') {
    xlabel <- 'breakpoints per chr arm'
  }

  main_plot <- ggplot2::ggplot(data = data.frame(x = c(0,round(xmax))), aes(x = x)) +
              labs(y = "", x = paste0("Number of ", xlabel)) +
              ggplot2::theme_bw() +
              theme(axis.text = ggplot2::element_text(size = 10),
                    axis.title = ggplot2::element_text(size = 12),
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.grid.major = ggplot2::element_blank())
  for (i in 1:length(plotparam)) {
    main_plot <- main_plot +
                    stat_function(geom = "area", n = round(xmax) + 1,
                                  fun = dpois, args = list(lambda = plotparam[i]),
                                  color = "black", size = 0.05,
                                  fill = cbPalette[i], alpha = shading[i])
  }
  return(main_plot)
}


#' Plot a Heatmap of Component-wise Wasserstein Distances
#'
#' @description
#'
#' This function calculates and builds a heatmap of component-wise wasserstein distances between two mixture models. \cr
#' Also known as calculating the 'Earth Mover's Distance'.
#' A CN-feature is indicated by the user using the "component" parameter for this function.
#' The mixture models for this feature in each set of component models are then pulled out and compared.
#' In short, for each component in each mixture model indicated, a distance will be calculated.
#'
#' @param cm_a A list of S4 objects each of class 'flexmix'. Likely the output from the FitMixtureModels function.
#' @param cm_b A list of S4 objects each of class 'flexmix'. Likely the output from the FitMixtureModels function.
#' @param component Which CN-feature models to compare to one another. (options: segsize, copynumber, bp10MB, osCN, bpchrarm)
#' @returns A *list*. One ggplot and its corresponding matrix.
#'
#'
#' @export
WassDistancePlot <- function(cm_a, cm_b, component) {

  # Pull component parameters out of the mixture model
  x <- flexmix::parameters(cm_a[[component]])
  y <- flexmix::parameters(cm_b[[component]])

  # Calculate Wasserstein distances
  if (cm_a[[component]]@call$model == "flexmix::FLXMCnorm1()" &&
      cm_b[[component]]@call$model == "flexmix::FLXMCnorm1()") {

    x <- as.data.frame(x)
    x <- x[, order(apply(x, 2, function(a) a[1]))]
    y <- as.data.frame(y)
    y <- y[, order(apply(y, 2, function(a) a[1]))]
    dist_matrix <- matrix(0, nrow = ncol(x), ncol = ncol(y))
    for (i in 1:ncol(x)) {
      for (j in 1:ncol(y)) {
        x_n <- rnorm(10000,x[1,i], x[2,i])
        y_n <- rnorm(10000,y[1,j], y[2,j])
        dist_matrix[i, j] <- waddR::wasserstein_metric(x_n, y_n)
      }
    }
    dist_df <- expand.grid(x = 1:dim(x)[2], y = 1:dim(y)[2])

  } else if (cm_a[[component]]@call$model == "flexmix::FLXMCmvpois()" &&
             cm_b[[component]]@call$model == "flexmix::FLXMCmvpois()") {

    x <- as.data.frame(t(x))
    x <- x[, order(apply(x, 2, function(a) a[1]))]
    y <- as.data.frame(t(y))
    y <- y[, order(apply(y, 2, function(a) a[1]))]
    dist_matrix <- matrix(0, nrow = length(x), ncol = length(y))
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        x_n <- rpois(10000,x[1,i])
        y_n <- rpois(10000,y[1,j])
        dist_matrix[i, j] <- waddR::wasserstein_metric(x_n, y_n)
      }
    }
    dist_df <- expand.grid(x = 1:length(x), y = 1:length(y))

  } else {
    stop("Unsuported mixture model(s) used.
       Please use either normals OR poissons ('flexmix::FLXMCnorm1()' or 'flexmix::FLXMCmvpois()')")
  }

  # Build Heatmap
  dist_df$value <- as.vector(dist_matrix)
  dist_df$value <- (dist_df$value - min(dist_df$value)) / (max(dist_df$value) - min(dist_df$value))
  dist_df$value <- -(dist_df$value - 0.5) + 0.5
  dist_df$x <- paste0("c.", as.character(dist_df$x))
  dist_df$y <- paste0("c.", as.character(dist_df$y))
  dist_df$x <- factor(dist_df$x, levels = unique(dist_df$x))
  dist_df$y <- factor(dist_df$y, levels = unique(dist_df$y))

  p1 <- ggplot2::ggplot(dist_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(option = "B",
                                name = "Similarity",
                                breaks = seq(0, 1, by = 0.2),
                                labels = c("", "Low", "", "", "High", "")) +
    ggplot2::labs(title = paste0("Heatmap of ",
                                 component,
                                 " component-wise wasserstein distance"),
                  x = "A - components", y = "B - components")

  return(list(plot = p1,
              matrix = dist_matrix))
}
