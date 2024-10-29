# File contains plotting functions related to CN-Signature visualization

###########################
### Functions
###########################
# StackedExposuresPlot
# TwoFeatureScatterPlot
# MixtureModelPlots
# GaussiansMixturePlot
# PoissonsMixturePlot
# WassDistancePlot


#' Create a Stacked Bar-plot of Signature Exposures
#'
#' @description
#' Converts signature-per-sample data into a ggplot2 stacked bar plot and returns that object.
#' Samples are sorted by signature exposure (default is the first row).
#' This can be changed to a different row, or a char vector or sample names can be supplied for a custom ordering.
#' The turbo viridis colour-scheme is used.
#'
#' @param input_matrix A dataframe/matrix/datatable with signatures in the rows and samples in the columns..
#' @param user_colours (optional) A vector of colours for the signatures, one per signature.
#' @param orderbyrow (optional) An integer. Single integer corresponding to the row by which to sort.
#' @param orderbyname (optional) NULL or a character vector of sample names.
#' @param do_transpose (optional) Logical. If TRUE then flip bars horizontal.
#' @returns A ggplot2 object.
#'
#' @export
StackedExposuresPlot <- function (input_matrix,
                                  user_colours = NULL,
                                  orderbyrow = 1,
                                  orderbyname = NULL,
                                  do_transpose = FALSE) {

  # Convert input matrix to data.table
  dt <- as.data.table(input_matrix, keep.rownames = TRUE)
  signames <- factor(dt$rn, levels = dt$rn)
  dt[, rn := NULL]

  if (orderbyrow > nrow(input_matrix) || orderbyrow < 1) {
    stop("The order_row must be between 1 and the number of rows in the matrix.")
  }

  # Order samples
  if (is.null(orderbyname)) {
    column_order <- order(input_matrix[orderbyrow, ], decreasing = TRUE)
  } else {
    if (!all(orderbyname %in% colnames(input_matrix))) {
      stop("All elements of orderbyname must match the column names of the input matrix.")
    }
    column_order <- match(orderbyname, colnames(input_matrix))
  }
  dt <- dt[, ..column_order]
  colnames(dt) <- colnames(input_matrix)[column_order]

  # Convert to long
  dt[, row_id := signames]
  dt_melted <- melt(dt, variable.name = "Category", value.name = "Value", id.vars = "row_id")
  dt_melted[, normalized_value := Value / sum(Value), by = Category]

  # Assign colours - use user-supplied colours if provided
  if (is.null(user_colours)) {
    user_colours <- viridis::turbo(nrow(input_matrix))
  } else if (length(user_colours) != nrow(input_matrix)) {
    stop("The length of user_colors must match the number of rows in the matrix")
  }

  # Create the plot
  p <- ggplot2::ggplot(dt_melted, ggplot2::aes(x = Category,
                                               y = normalized_value,
                                               fill = row_id)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 1) +
    ggplot2::scale_fill_manual(values = user_colours) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Samples", y = "Exposures", fill = "Signatures")
  p
  if (do_transpose) {
    p <- p + ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_line())
  }

  return(p)
}



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
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param signatures (optional) A dataframe with components along the x-axis and signatures along the y.
#' @param sig_of_interest (optional) A single integer corresponding the signature for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @returns A list of ggplots.
#'
#'
#' @export
MixtureModelPlots <- function(components,
                              signatures,
                              sig_of_interest = 1,
                              threshold = 0.3) {

  # Validate input
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              inherits(components, 'list'))
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              (inherits(components[[1]], 'flexmix') |
                 inherits(components[[1]], 'data.frame')))
  stopifnot("Signatures argument must be tabular in the form of a matrix/dataframe/datatable." =
              (inherits(signatures, 'matrix') |
                 inherits(signatures, 'data.frame') |
                 inherits(signatures, 'data.table')))
  stopifnot("The indicated signature (sig_of_interest parameter) must be an
            integer less than or equal to the # of columns in 'signatures'." =
              dim(signatures) >= sig_of_interest)

  all_plots <- list()
  all_plots[["segmentsize"]] <- GaussiansMixturePlot(components, signatures,
                                                     sig_of_interest,
                                                     component = 'segsize',
                                                     threshold = threshold)
  all_plots[["breakpoint10MB"]] <- PoissonsMixturePlot(components, signatures,
                                                       sig_of_interest,
                                                       component = 'bp10MB',
                                                       threshold = threshold)
  all_plots[["oscillating"]] <- PoissonsMixturePlot(components, signatures,
                                                    sig_of_interest,
                                                    component = 'osCN',
                                                    threshold = threshold)
  all_plots[["changepoint"]] <- GaussiansMixturePlot(components, signatures,
                                                     sig_of_interest,
                                                     component = 'changepoint',
                                                     threshold = threshold)
  all_plots[["copynumber"]] <- GaussiansMixturePlot(components, signatures,
                                                    sig_of_interest,
                                                    component = 'copynumber',
                                                    threshold = threshold)
  all_plots[["breakpointsarm"]] <- PoissonsMixturePlot(components, signatures,
                                                       sig_of_interest,
                                                       component = 'bpchrarm',
                                                       threshold = threshold)
  if("nc50" %in% names(components)){
    all_plots[["nc50"]] <- PoissonsMixturePlot(components, signatures,
                                               sig_of_interest,
                                               component = 'nc50',
                                               threshold = threshold)
  }
  if("cdist" %in% names(components)){
    all_plots[["cdist"]] <- GaussiansMixturePlot(components, signatures,
                                                 sig_of_interest,
                                                 component = 'cdist',
                                                 threshold = threshold)
  }


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
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param signatures (optional) A dataframe with components in the rows and signatures in the columns.
#' @param sig_of_interest (optional) A single integer or a vector of integers corresponding the signatures for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @param component (optional) Which component to draw the gaussian curves for. (options: segsize, changepoint, copynumber)
#' @param inlay_flag (optional) Boolean flag to turn on/off the inlay plot.
#' @returns A ggplot.
#'
#'
#' @export
GaussiansMixturePlot <- function(components,
                                 signatures = NULL,
                                 sig_of_interest = 1,
                                 log_flag = FALSE,
                                 component = 'segsize',
                                 inlay_flag = TRUE,
                                 threshold = 0.3) {

  # Validate input
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              inherits(components, 'list'))
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              (inherits(components[[1]], 'flexmix') |
                 inherits(components[[1]], 'data.frame') |
                 inherits(components[[1]], 'matrix')))

  # Extract components for a cn-feature and sort in ascending order
  if (inherits(components[[1]], "data.frame") | inherits(components[[1]], "matrix")) {
    plotparam <- components[[component]]
  } else if (inherits(components[[1]], "flexmix")) {
    plotparam <- flexmix::parameters(components[[component]])
  } else {
    stop("Input format not recognized for component models, please provide a dataframes/matrix/flexmix object. \n
         See function docs.")
  }
  plotparam <- as.data.frame(plotparam[,order(plotparam[1,])])

  # Log-transform if needed
  if (log_flag) {
    plotparam[2,] <- plotparam[2,]/plotparam[1,]
    plotparam[1,] <- log(plotparam[1,])
  }

  # If signatures provided, then shade by weight, otherwise set to uniform opacity
  if (!is.null(signatures)) {
    stopifnot("Signatures argument must be tabular in the form of a matrix/dataframe/datatable." =
                (inherits(signatures, 'matrix') |
                   inherits(signatures, 'data.frame') |
                   inherits(signatures, 'data.table')))
    norm_const <- apply(signatures, 2, sum)
    sig_mat_norm <- sweep(signatures, 2, norm_const, "/")
    weights <- as.numeric(sig_mat_norm[,sig_of_interest])
    mask <- grepl(component, rownames(sig_mat_norm), fixed=TRUE)
    shading <- weights[mask]
    segpalette <- viridis::turbo(ncol(plotparam))

  } else {
    shading <- rep(1,ncol(plotparam))
    segpalette <- viridis::turbo(ncol(plotparam))
  }

  # Define plot layout
  max_comp <- plotparam[,which(plotparam[1,] == max(plotparam[1,]))]
  min_comp <- plotparam[,which(plotparam[1,] == min(plotparam[1,]))]
  xmax <- max_comp[1] + (max_comp[2]*3)
  xmin <- pmax(min_comp[1] - (min_comp[2]*3), 0)
  digits <- nchar(as.character(round(xmax/2)))
  plotbreaks <- c(10^(digits-1), (10^digits)/2, 10^digits)
  plotbreaks <- plotbreaks[xmax > plotbreaks]
  plotparam_df <- as.data.frame(t(plotparam))
  colnames(plotparam_df) <- c("mean", "sd")
  plotparam_df$component <- factor(1:nrow(plotparam_df))
  x_vals <- seq(xmin, xmax, length.out = 1000)  # Increase the number of points for a smoother line

  # Make the main plot
  main_plot <- ggplot2::ggplot(data = data.frame(x = c(1,xmax)),
                               ggplot2::aes(x)) +
    ggplot2::ylab("") +
    ggplot2::xlab(paste0(component, " mixture components")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10),
                   axis.title = ggplot2::element_text(size = 12),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks = plotbreaks)

  for (i in 1:ncol(plotparam)) {
    if (shading[i] < threshold) {
      linecolour <- "grey"
      linealpha <- 0.4
    } else {
      linecolour <- segpalette[i]
      linealpha <- ifelse(shading[i] < 0.5, 0.5, shading[i])
    }
    y_vals <- dnorm(x_vals, mean = plotparam[1, i], sd = plotparam[2, i])

    main_plot <- main_plot +
      ggplot2::geom_line(data = data.frame(x = x_vals, y = y_vals),
                         ggplot2::aes(x = x, y = y),
                         size = 1,
                         color = linecolour,
                         alpha = linealpha) +
      ggplot2::geom_vline(xintercept = plotparam[1, i],
                          linetype = "dashed",
                          color = linecolour,
                          alpha = linealpha) +
      ggplot2::geom_text(data = plotparam_df[i, , drop = FALSE],
                         ggplot2::aes(x = mean, y = -1,
                                      label = if (mean > 10) {
                                        round(mean)
                                        } else { round(mean, 3)} ),
                         angle = 90,
                         vjust = if (i == 1) {
                           -0.5
                         } else if (round(plotparam_df[i, ][[1]],3) -
                                    round(plotparam_df[i-1, ][[1]], 3) < 0.01) {
                           1.5
                         } else {
                           -0.5
                         },
                         color = linecolour,
                         alpha = linealpha)
  }
  final_plot <- main_plot

  # Make an inlay plot of the 'important' components if..
  # 1. Some components have been 'squashed' down
  # 2. There are any. - i.e. at least 1 weight > 0.05
  if ((plotparam[1, dim(plotparam)[2]]/plotparam[1,1] > 20) &&
      (sum(shading > 0.05) >= 1) && inlay_flag) {
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
    inlay_plot <- ggplot2::ggplot(data = data.frame(x = c(0, xmax)),
                                  ggplot2::aes(x)) +
      ggplot2::ylab("") + ggplot2::xlab("") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 8),
                     axis.title = ggplot2::element_text(size = 8),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(breaks = plotbreaks)

    plotparam_df <- as.data.frame(t(plotparam2))
    colnames(plotparam_df) <- c("mean", "sd")
    plotparam_df$component <- factor(1:nrow(plotparam_df))

    x_vals <- seq(0, xmax, length.out = 1000)  # Increase the number of points for a smoother line
    for (i in 1:ncol(plotparam2)) {
      y_vals <- dnorm(x_vals, mean = plotparam[1, i], sd = plotparam[2, i])
      inlay_plot <- inlay_plot +
        ggplot2::geom_line(data = data.frame(x = x_vals, y = y_vals),
                           ggplot2::aes(x = x, y = y),
                           size = 1,
                           color = ifelse(shading2[i] < threshold, "grey", segpalette[i]),
                           alpha = ifelse(shading2[i] < threshold, 0.5, shading2[i]))+
        ggplot2::geom_vline(xintercept = plotparam2[1, i],
                            linetype = "dashed",
                            color = ifelse(shading2[i] < threshold, "grey", segpalette[i]),
                            alpha = ifelse(shading2[i] < threshold, 0.5, shading2[i]))+
        ggplot2::geom_text(data = plotparam_df[i, , drop = FALSE],
                           ggplot2::aes(x = mean, y = -2,
                                        label = if (mean > 10) {
                                          round(mean)
                                        } else { round(mean, 3)} ),
                           angle = 90,
                           vjust = if (i == 1) {
                             -0.5
                           } else if (round(plotparam_df[i, ][[1]],3) - round(plotparam_df[i-1, ][[1]],3)<0.05) {
                             1.5
                           } else {
                             -0.5
                           },
                           color = ifelse(shading2[i] < threshold, "grey", segpalette[i]),
                           alpha = ifelse(shading2[i] < threshold, 0.5, shading2[i]))
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
#' @param components A list of S4 objects belonging to the class 'flexmix'.
#' @param signatures (optional) A dataframe with components in the rows and signatures in the columns.
#' @param sig_of_interest (optional) A single integer or a vector of integers corresponding the signatures for which to make plots. \cr
#' In the absence of any value passed in for this parameter, plots for just the first signature will be returned.
#' @param component (optional) Which component to draw the poisson 'curves' for. (options: bp10MB, osCN, bpchrarm)
#' @returns A ggplot.
#'
#'
#' @export
PoissonsMixturePlot <- function(components,
                                signatures = NULL,
                                sig_of_interest = 1,
                                log_flag = FALSE,
                                component = 'bp10MB',
                                threshold = 0.3) {

  # Validate input
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              inherits(components, 'list'))
  stopifnot("Components argument must be a list of flexmix objects or tables." =
              (inherits(components[[1]], 'flexmix') |
                 inherits(components[[1]], 'data.frame') |
                 inherits(components[[1]], 'matrix')))

  # Extract components for a cn-feature and sort in ascending order
  if (inherits(components[[1]], "data.frame") |
      inherits(components[[1]], "matrix")) {
    plotparam <- components[[component]]
  } else if (inherits(components[[1]], "flexmix")) {
    plotparam <- flexmix::parameters(t(components[[component]]))
  } else {
    stop("Input format not recognized for component models, please provide a dataframes/matrix/flexmix object. \n
         See function docs.")
  }
  plotparam <- plotparam[order(plotparam)]

  # Log-transform if needed
  if (log_flag) {
    plotparam <- log(plotparam)
  }

  # If signatures provided, then shade by weight, otherwise set to uniform opacity
  if (!is.null(signatures)) {
    stopifnot("Signatures argument must be tabular in the form of a matrix/dataframe/datatable." =
                (inherits(signatures, 'matrix') |
                   inherits(signatures, 'data.frame') |
                   inherits(signatures, 'data.table')))
    norm_const <- apply(signatures, 2, sum)
    sig_mat_norm <- sweep(signatures, 2, norm_const, "/")
    weights <- as.numeric(sig_mat_norm[,sig_of_interest])
    mask <- grepl(component, rownames(sig_mat_norm), fixed=TRUE)
    shading <- weights[mask]
    segpalette <- viridis::turbo(length(plotparam))

  } else {
    shading <- rep(1, length(plotparam))
    segpalette <- viridis::turbo(length(plotparam))
  }

  # Define plot layout
  max_comp <- max(plotparam)
  min_comp <- min(plotparam)
  xmax <- max_comp + 3*sqrt(max_comp)
  xmin <- pmax(min_comp - 3*sqrt(min_comp), 0)
  plotparam_df <- as.data.frame(plotparam)
  colnames(plotparam_df) <- c("lambda")
  x_vals <- round(seq(0, xmax, length.out = 1000))  # Increase the number of points for a smoother line

  # Make the main plot
  main_plot <- ggplot2::ggplot(data = data.frame(x = c(0,round(xmax))),
                               ggplot2::aes(x = x)) +
    ggplot2::labs(y = "", x = paste0(component, " mixture components")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10),
                   axis.title = ggplot2::element_text(size = 12),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())

  for (i in 1:length(plotparam)) {
    if (shading[i] < threshold) {
      linecolour <- "grey"
      linealpha <- 0.4
    } else {
      linecolour <- segpalette[i]
      linealpha <- ifelse(shading[i] < 0.5, 0.5, shading[i])
    }
    y_vals <- dpois(x = x_vals, lambda = plotparam[i])
    main_plot <- main_plot +
      ggplot2::geom_line(data = data.frame(x = x_vals, y = y_vals),
                         ggplot2::aes(x = x, y = y),
                         size = 1,
                         color = linecolour,
                         alpha = linealpha) +
      ggplot2::geom_vline(xintercept = plotparam[i],
                          linetype = "dashed",
                          color = linecolour,
                          alpha = linealpha) +
      ggplot2::geom_text(data = plotparam_df[i, , drop = FALSE],  # Ensure it's treated as a data frame
                         ggplot2::aes(x = lambda, y = 0,
                                      label = if (lambda > 10) {
                                        round(lambda)
                                      } else { round(lambda, 3)} ),
                         angle = 90,
                         vjust = if (i == 1) {
                           -0.5
                         } else if (round(plotparam_df$lambda[i], 3) -
                                    round(plotparam_df$lambda[i - 1], 3) < 0.01) {
                           1.5
                         } else {
                           -0.5
                         },
                         hjust = 1,
                         color = linecolour,
                         alpha = linealpha)
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

  # Get model params
  params <- list()
  mmodels <- list(cm_a = cm_a, cm_b = cm_b)
  for (i in c(1:length(mmodels))) {
    if (inherits(mmodels[[i]][[1]], "data.frame") | inherits(mmodels[[i]][[1]], "matrix")) {
      params[[names(mmodels[i])]] <- as.data.frame(mmodels[[i]][[component]])
    } else if (inherits(mmodels[[i]][[1]], "flexmix")) {
      if (grepl("FLXMCnorm1",
                as.character(mmodels[[i]][[component]]@call$model),
                fixed = TRUE)) {
        params[[names(mmodels[i])]] <- as.data.frame(flexmix::parameters(mmodels[[i]][[component]]))
      } else if (grepl("FLXMCmvpois",
                  as.character(mmodels[[i]][[component]]@call$model),
                  fixed = TRUE)) {
        params[[names(mmodels[i])]] <- as.data.frame(t(flexmix::parameters(mmodels[[i]][[component]])))
      } else {
        stop("Unsuported flexmixmixture model(s) used.
             Please use either normals OR poissons ('flexmix::FLXMCnorm1()' or 'flexmix::FLXMCmvpois()')")
      }
    } else {
      stop("Input format not recognized, please provide dataframes/matrices/flexmix objects. \n
           See function docs.")
    }
  }

  # Calculate Wasserstein distances
  if (dim(params$cm_a)[1] == 2) {
    x <- params$cm_a[, order(apply(params$cm_a, 2, function(a) a[1]))]
    y <- params$cm_b[, order(apply(params$cm_b, 2, function(a) a[1]))]
    dist_matrix <- matrix(0, nrow = ncol(x), ncol = ncol(y))
    for (i in 1:ncol(x)) {
      for (j in 1:ncol(y)) {
        x_n <- rnorm(10000,x[1,i], x[2,i])
        y_n <- rnorm(10000,y[1,j], y[2,j])
        dist_matrix[i, j] <- waddR::wasserstein_metric(x_n, y_n)
      }
    }
    dist_df <- expand.grid(x = 1:dim(x)[2], y = 1:dim(y)[2])

  } else if (dim(params$cm_a)[1] == 1) {
    x <- params$cm_a[, order(apply(params$cm_a, 2, function(a) a[1]))]
    y <- params$cm_b[, order(apply(params$cm_b, 2, function(a) a[1]))]
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
    stop("Unsuported mixture model parameter format provided.
         Double check input dataframes.")
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
    ggplot2::labs(title = paste0(component,
                                 " component-wise wasserstein distance"),
                  x = "A - components", y = "B - components")

  return(list(plot = p1,
              matrix = dist_matrix))
}


#' Make an Alluvial Plot based on maximum signature exposure
#'
#' @description
#'
#' SEAlluvialPlot builds a ggplot using the easyalluvial R package.
#' For two sets of exposures on the same set of samples it visualizes how the maximum exposure for one set of signatures is related to the other.
#'
#' @param exposA A matrix of signature exposures. Signatures in rows, samples in columns. Likely the output from `CallSignatureExposures()`.
#' @param exposB A matrix of signature exposures. Signatures in rows, samples in columns.
#' @param orderA (Optional) A character vector of signature names that will enforce an order. \cr
#' Signature set A is assigned names S.A1 -> S.A#, and B S.B1 -> S.B# \cr
#' ex. c('S.B2', 'S.B5', 'S.B3', 'S.B4', 'S.B1', 'S.B6')
#' @param orderB (Optional) A character vector of signature names that will enforce an order. \cr
#' Signature set A is assigned names S.A1 -> S.A#, and B S.B1 -> S.B# \cr
#' ex. c('S.B2', 'S.B5', 'S.B3', 'S.B4', 'S.B1', 'S.B6')
#' Must be for the same set of same set of samples as provided to `exposA`. Likely the output from `CallSignatureExposures()`.
#' @returns A *list*. A ggplot and its corresponding data matrix.
#'
#'
#' @export
SEAlluvialPlot <- function (exposA, exposB,
                            orderA = NULL,
                            orderB = NULL) {

  exposA <- as.data.frame(t(exposA))
  colnames(exposA) <- paste0('S.A', as.character(1:dim(exposA)[2]))
  exposA$max_sigA <- apply(exposA, 1, function(x) which.max(x))
  exposA$max_sigA <- paste0('S.A', exposA$max_sigA)
  if (!is.null(orderA)) { factor(exposA$max_sigA, levels = orderA) }
  exposA$sample_ids <- rownames(exposA)

  exposB <- as.data.frame(t(exposB))
  colnames(exposB) <- paste0('S.B', as.character(1:dim(exposB)[2]))
  exposB$max_sigB <- apply(exposB, 1, function(x) which.max(x))
  exposB$max_sigB <- paste0('S.B', exposB$max_sigB)
  if (!is.null(orderB)) { factor(exposB$max_sigB, levels = orderB) }
  exposB$sample_ids <- rownames(exposB)

  plotting_data <- exposA %>% dplyr::left_join(exposB, by = c('sample_ids'))
  rownames(plotting_data) <- plotting_data$sample_ids
  plotting_data$sample_ids <- NULL

  ggp <- easyalluvial::alluvial_wide(dplyr::select(plotting_data,
                                                    max_sigA, max_sigB),
                                     fill_by = 'first_variable',
                                     stratum_label_size = 3.5) +
    ggplot2::labs(y = "Samples", caption = '') +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 7),
                   axis.title.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_blank(),
                   title = ggplot2::element_blank())

  return(list(plot = ggp, exposures_matrix = plotting_data))
}
