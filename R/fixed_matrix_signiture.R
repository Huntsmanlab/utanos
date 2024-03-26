#' Estimate Signature Number
#'
#' Use **NMF** package to evaluate the optimal number of signatures.
#' This is used along with [sig_extract].
#' Users should `library(NMF)` firstly. If NMF objects are returned,
#' the result can be further visualized by NMF plot methods like
#' `NMF::consensusmap()` and `NMF::basismap()`.
#'
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' More custom features please directly use [NMF::nmfEstimateRank].
#'
#' @name sig_estimate
#' @param nmf_matrix a `matrix` used for NMF decomposition with rows indicate samples and columns indicate components.
#' @param range a `numeric` vector containing the ranks of factorization to try. Note that duplicates are removed
#' and values are sorted in increasing order. The results are notably returned in this order.
#' @param keep_nmfObj default is `FALSE`, if `TRUE`, keep NMF objects from runs, and the result may be huge.
#' @param nrun a `numeric` giving the number of run to perform for each value in `range`, `nrun` set to 30~50 is
#' enough to achieve robust result.
#' @param what a character vector whose elements partially match one of the following item, which correspond to
#' the measures computed by summary on each multi-run NMF result: ‘all’, ‘cophenetic’, ‘rss’, ‘residuals’,
#'  ‘dispersion’, ‘evar’, ‘silhouette’ (and more specific .coef, .basis, .consensus), ‘sparseness’
#'  (and more specific .coef, .basis). It specifies which measure must be plotted
#'  (what='all' plots all the measures).
#'
#' @param cores number of cpu cores to run NMF.
#' @param seed specification of the starting point or seeding method, which will compute a starting point,
#'  usually using data from the target matrix in order to provide a good guess.
#' @param use_random Should generate random data from input to test measurements. Default is `TRUE`.
#' @param save_plots if `TRUE`, save signature number survey plot to local machine.
#' @param plot_basename when save plots, set custom basename for file path.
#' @param method specification of the NMF algorithm. Use 'brunet' as default.
#' Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
#' @param verbose if `TRUE`, print extra message.
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @return - sig_estimate: a `list` contains information of NMF run and rank survey.
#' @export
#' @seealso [sig_extract] for extracting signatures using **NMF** package, [sig_auto_extract] for
#' extracting signatures using automatic relevance determination technique.
sig_estimate <-
  function(nmf_matrix,
           range = 2:5,
           nrun = 10,
           use_random = FALSE,
           method = "brunet",
           seed = 123456,
           cores = 1,
           keep_nmfObj = FALSE,
           save_plots = FALSE,
           plot_basename = file.path(tempdir(), "nmf"),
           what = "all",
           verbose = FALSE) {
    if (nrow(nmf_matrix) < max(range)) {
      stop("The 'range' should not greater than ", nrow(nmf_matrix), " in your case.")
    }

    eval(parse(text = "suppressMessages(library('NMF'))"))
    if (cores > 1) cores <- min(cores, future::availableCores())
    mat <- t(nmf_matrix)

    ii <- colSums(mat) < 0.01
    if (any(ii)) {
      message(
        "The follow samples dropped due to null catalogue:\n\t",
        paste0(colnames(mat)[ii], collapse = ", ")
      )
      mat <- mat[, !ii, drop = FALSE]
    }

    # To avoid error due to NMF
    mat <- mat + 1e-12

    if (cores > 1) {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed,
          .opt = paste0("p", cores)
        )
    } else {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed
        )
    }
    print("hello")

    nmf.sum <- NMF::summary(estim.r) # Get summary of estimates
    if (verbose) {
      message("Estimation of rank based on observed data.")
      print(nmf.sum)
    }

    if (use_random) {
      if (verbose) message("Generating random matrix and run NMF...")
      V.random <- NMF::randomize(mat)

      if (cores > 1) {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed,
            .opt = paste0("p", cores)
          )
      } else {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed
          )
      }

      nmf.sum.random <- NMF::summary(estim.r.random) # Get summary of estimates
      if (verbose) {
        message("Estimation of rank based on random data.")
        print(nmf.sum.random)
      }
    } else {
      estim.r.random <- NULL
      nmf.sum.random <- NULL
    }

    if (save_plots) {
      if (use_random) {
        p <- NMF::plot(
          estim.r,
          estim.r.random,
          what = what,
          xname = "Observed",
          yname = "Randomised",
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      } else {
        p <- NMF::plot(
          estim.r,
          what = what,
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      }

      destdir <- dirname(plot_basename)
      if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
      pdf(
        paste0(plot_basename, "_survey.pdf"),
        bg = "white",
        pointsize = 9,
        width = 6 + abs((nrow(nmf.sum) - 6) / 3),
        height = 6,
        paper = "special"
      )
      print(p)
      dev.off()
      if (verbose) message("Created ", paste0(plot_basename, "_survey.pdf"))
    }

    if (keep_nmfObj) {
      res <- list(
        nmfEstimate = estim.r,
        nmfEstimate.random = estim.r.random,
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    } else {
      res <- list(
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    }

    class(res) <- "Survey"
    res
  }

#' Show Comprehensive Signature Number Survey
#'
#' `show_sig_number_survey2()` is modified from **NMF** package to
#' better help users to explore survey of signature number.
#'
#' @rdname sig_estimate
#' @param x a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @param y for random simulation,
#' a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @param what a character vector whose elements partially match one of the following item,
#' which correspond to the measures computed by `summary()` on each – multi-run – NMF result:
#' 'all', 'cophenetic', 'rss', 'residuals', 'dispersion', 'evar', 'silhouette'
#' (and more specific `*.coef`, `*.basis`, `*.consensus`), 'sparseness'
#' (and more specific `*.coef`, `*.basis`).
#' It specifies which measure must be plotted (what='all' plots all the measures).
#' @inheritParams NMF::nmfEstimateRank
#'
#' @return - show_sig_number_survey2: a `ggplot` object
#' @export
show_sig_number_survey <- function(x, y = NULL, what = c(
  "all", "cophenetic", "rss", "residuals",
  "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef",
  "silhouette", "silhouette.coef", "silhouette.basis", "silhouette.consensus"
),
na.rm = FALSE, xlab = "Total signatures",
ylab = "", main = "Signature number survey using NMF package") {

  # Useless, just store it in case I need
  # to modify in the future
  xname <- "x"
  yname <- "y"

  if (is.character(y) && missing(what)) {
    what <- y
    y <- NULL
  }
  what <- match.arg(what, several.ok = TRUE)
  if ("all" %in% what) {
    what <- c(
      "cophenetic", "rss", "residuals", "dispersion",
      "evar", "sparseness", "silhouette"
    )
  }
  .getvals <- function(x, xname) {
    measures <- x
    iwhat <- unlist(lapply(
      paste("^", what, sep = ""), grep,
      colnames(measures)
    ))
    if (na.rm) {
      measures <- measures[apply(measures, 1, function(row) !any(is.na(row[iwhat]))), ]
    }
    vals <- measures[, iwhat, drop = FALSE]
    x <- as.numeric(measures$rank)
    xlim <- range(x)
    measure.type <- setNames(
      rep("Best fit", ncol(measures)),
      colnames(measures)
    )
    cons.measures <- c(
      "silhouette.consensus", "cophenetic",
      "cpu.all"
    )
    measure.type[match(cons.measures, names(measure.type))] <- "Consensus"
    measure.type[grep("\\.coef$", names(measure.type))] <- "Coefficients"
    measure.type[grep("\\.basis$", names(measure.type))] <- "Basis"
    measure.type <- factor(measure.type)
    pdata <- tidyr::pivot_longer(cbind(rank = x, vals),
                                 cols = colnames(vals),
                                 names_to = "variable"
    )
    pdata$Type <- measure.type[as.character(pdata$variable)]
    pdata$Measure <- gsub("^([^.]+).*", "\\1", pdata$variable)
    pdata$Data <- xname
    pdata
  }

  if (inherits(x, "NMF.rank")) {
    x <- x$measure
    pdata <- .getvals(x, xname)
  } else {
    pdata <- .getvals(x, xname)
  }

  if (!is.null(y)) {
    if (inherits(y, "NMF.rank")) {
      y <- y$measure
      pdata.y <- .getvals(y, yname)
    } else {
      pdata.y <- .getvals(y, yname)
    }

    pdata <- rbind(pdata, pdata.y)
  }

  p <- ggplot2::ggplot(pdata, ggplot2::aes_string(x = "rank", y = "value")) +
    ggplot2::geom_line(ggplot2::aes_string(linetype = "Data", colour = "Type")) +
    ggplot2::geom_point(size = 2,ggplot2:: aes_string(shape = "Data", colour = "Type")) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(xlab, breaks = unique(pdata$rank)) +
    ggplot2::scale_y_continuous(ylab) +
    ggplot2::ggtitle(main)
  if (!is(y, "NMF.rank")) {
    p <- p + ggplot2::scale_shape(guide = "none") + ggplot2::scale_linetype(guide = "none")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    send_stop("Please install 'RColorBrewer' package firstly.")
  }
  myColors <- RColorBrewer::brewer.pal(5, "Set1")
  names(myColors) <- levels(pdata$Type)
  p <- p + ggplot2::scale_colour_manual(name = "Measure type", values = myColors)
  p <- p + ggplot2::facet_wrap(~Measure, scales = "free")
  p
}

#' Extract Signatures through NMF
#'
#' Do NMF de-composition and then extract signatures.
#'
#' @param n_sig number of signature. Please run [sig_estimate] to select a suitable value.
#' @param optimize if `TRUE`, then refit the denovo signatures with QP method, see [sig_fit].
#' @param pynmf if `TRUE`, use Python NMF driver [Nimfa](http://nimfa.biolab.si/index.html).
#' The seed currently is not used by this implementation.
#' @param ... other arguments passed to [NMF::nmf()].
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @return a `list` with `Signature` class.
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' res <- sig_extract(cn_tally_W$nmf_matrix, 2, nrun = 1)
#' }
#' @seealso [sig_tally] for getting variation matrix,
#' [sig_estimate] for estimating signature number for [sig_extract]
sig_extract <- function(nmf_matrix,
                        n_sig,
                        nrun = 10,
                        cores = 1,
                        method = "brunet",
                        seed = 123456, ...) {
  eval(parse(text = "suppressMessages(library('NMF'))"))
  if (cores > 1) cores <- min(cores, future::availableCores())
  # transpose matrix
  mat <- t(nmf_matrix)

  ii <- colSums(mat) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(colnames(mat)[ii], collapse = ", ")
    )
    mat <- mat[, !ii, drop = FALSE]
  }


  # To avoid error due to NMF
  mat <- mat + 1e-12

  nmf.res <- NMF::nmf(
    mat,
    n_sig,
    seed = seed,
    nrun = nrun,
    method = method,
    .opt = paste0("vp", cores),
    ...
  )

  # Signature loading
  W <- NMF::basis(nmf.res)
  # Exposure loading
  H <- NMF::coef(nmf.res)
  # Signature number
  K <- ncol(W)

  Signature <- W
  Exposure <- H

  for (j in seq_len(K)) {
    Signature[, j] <- Signature[, j] * rowSums(Exposure)[j]
    Exposure[j, ] <- Exposure[j, ] * colSums(Signature)[j]
  }

  scal_res <- list(Signature = Signature, Exposure = Exposure)
  Signature <- scal_res$Signature
  Exposure <- scal_res$Exposure

  Signature.norm <- apply(Signature, 2, function(x) x / sum(x, na.rm = TRUE))
  Exposure.norm <- apply(Exposure, 2, function(x) x / sum(x, na.rm = TRUE))


  # When only one signature
  if (!is.matrix(Exposure.norm)) {
    Exposure.norm <- matrix(Exposure.norm, nrow = 1, dimnames = list(NULL, names(Exposure.norm)))
  }

  if (ncol(Signature) > 1) {
    # Get orders
    sig_orders <- Exposure %>%
      rowSums() %>%
      order(decreasing = TRUE)

    Signature <- Signature[, sig_orders]
    Signature.norm <- Signature.norm[, sig_orders]
    Exposure <- Exposure[sig_orders, ]
    Exposure.norm <- Exposure.norm[sig_orders, ]

    W <- W[, sig_orders]
    H <- H[sig_orders, ]
  }

  sig_names <- paste0("Sig", seq_len(K))
  colnames(W) <- colnames(Signature) <- colnames(Signature.norm) <- sig_names
  rownames(H) <- rownames(Exposure) <- rownames(Exposure.norm) <- sig_names

  res <- list(
    Signature = Signature,
    Signature.norm = Signature.norm,
    Exposure = Exposure,
    Exposure.norm = Exposure.norm,
    K = K,
    Raw = list(
      nmf_obj = if (exists("nmf.res")) nmf.res else NULL,
      W = W,
      H = H
    )
  )
  class(res) <- "Signature"
  attr(res, "nrun") <- nrun
  attr(res, "method") <- method
  attr(res, "seed") <- seed
  attr(res, "call_method") <- "NMF"

  invisible(res)
}

count_components <- function(min, max, label, feature) {
  # Count all samples
  if (label == "point") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value == min, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else if (label == "range") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(as.numeric(.data$value) > as.numeric(min) & as.numeric(.data$value) <= as.numeric(max), na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else {
    stop("Bad labels for feature setting, can only be 'point' and 'range'!")
  }
}


count_components_wrapper <- function(feature_df, f_name, cn_feature_setting) {
  samps <- unique(feature_df$ID)
  feature_df <- feature_df %>% dplyr::as_tibble()
  # Make sure sample order is consistent
  feature_df$ID <- factor(feature_df$ID, levels = samps)
  feature_df$value <- as.numeric(feature_df$value)
  if (f_name == "segsize"){
    feature_df$value <- log10(feature_df$value)
  }

  specific_f <- cn_feature_setting[cn_feature_setting$feature == f_name,]

  specific_f %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$component) %>%
    dplyr::summarize(
      count = paste(count_components(.data$min, .data$max, .data$label, feature_df), collapse = ","),
    ) %>%
    tidyr::separate(.data$count, samps, sep = ",", convert = TRUE) %>%
    data.table::as.data.table()
}

#' @export
GetCNTally <- function(CN_data, genome = 'hg19', relative=FALSE, extra_features = FALSE){
  if(relative){
    CN_features <- ExtractRelativeCopyNumberFeatures(CN_data, genome = genome, extra_features = extra_features)
    feature_setting <- utanos::feature_setting_default
    if(extra_features){
      feature_setting <- utanos::feature_setting_extra
    }
  }
  else{
    CN_features <- ExtractCopyNumberFeatures(CN_data, genome = genome, extra_features = extra_features)
    feature_setting <- utanos::feature_setting_default_absoulute
    if(extra_features){
      feature_setting <- utanos::feature_setting_extra_absoulute
    }
  }

  cn_components <- purrr::map2(CN_features, names(CN_features),
                               count_components_wrapper,
                               cn_feature_setting = feature_setting
  )

  cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
    dplyr::as_tibble() %>%
    tibble::column_to_rownames(var = "component") %>%
    as.matrix()

  # Order the matrix as feature_setting
  cn_matrix <- cn_matrix[feature_setting$component, ] %>%
    t()
}

#' @export
PlotCNTallySig <- function(sig_matrix, extra_features = FALSE){

  sig_matrix <- as.data.frame(sig_matrix)

  prefixes <- gsub("\\[.*\\]", "", rownames(sig_matrix))
  suffixes <- gsub(".*\\[(.*)\\]", "\\1", rownames(sig_matrix))

  # Add a new column for row name prefixes
  sig_matrix$prefix <- prefixes
  sig_matrix$suffix <- suffixes

  sig_long <-  tidyr::gather(sig_matrix, key = "signature", value = "value", -prefix, -suffix)

  suffix_order <- c("0", "1", "<=2", "2", "3", "4", "5", ">5", "6", "7", "8", "9", "10",
                    ">10 & <=20", ">20 & <=30", ">30", ">4 & <=10", ">10",
                    ">2 & <=3", ">3 & <=4", ">4 & <=5", ">5 & <=6", ">6 & <=7",
                    ">7 & <=8", ">8", ">7", "<= -1", ">-1 & <=-0.9", ">-0.9 & <=-0.8", ">-0.8 & <=-0.7",
                    ">-0.7 & <=-0.6", ">-0.6 & <=-0.5", ">-0.5 & <=-0.4", ">-0.4 & <=-0.3",
                    ">-0.3 & <=-0.2", ">-0.2 & <=-0.1", ">-0.1 & <=0", "<=0", ">0 & <=0.1", ">0.1 & <=0.2",
                    ">0.2 & <=0.3", ">0.3 & <=0.4", ">0.4 & <=0.5", ">0.5 & <=0.6",
                    ">0.6 & <=0.7", ">0.7 & <=0.8", ">0.8 & <=0.9", ">0.9 & <=1",
                    ">1", "> 1")

  sig_long$suffix <- factor(sig_long$suffix, levels = suffix_order)

  sig_long <- sig_long %>%
    dplyr::group_by(signature, prefix) %>%
    dplyr::mutate(value = value / sum(value)) %>%
    dplyr::ungroup()

  ggplot2::ggplot(sig_long, ggplot2::aes(x = suffix, y = value, fill = prefix)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(x = "Components", y = "Contributions", fill = "Signature") +
    ggplot2::facet_grid(signature~prefix,scales = "free_x", space = "free_x")+
    ggplot2::scale_x_discrete(breaks = unique(sig_long$suffix)) +
    ggplot2::theme(
      legend.position = "none",
      axis.line = ggplot2::element_line(size = 0.3, colour = "black"),
      panel.spacing.x =  ggplot2::unit(0.1, "line"),
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
        angle = 0),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10))  # Use prefix color for facet labels)
}
