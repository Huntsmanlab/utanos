# mainsignaturefunctions.R
# The functions in this file were copied and modified from elsewhere.
# The original code can be found here: https://bitbucket.org/britroc/cnsignatures/src/master/
# The original code was published with this paper: https://doi.org/10.1038%2Fs41588-018-0179-8

#' Quantify Signature Exposures
#'
#' @description
#'
#' Given a set of \cr
#' 1. signatures and \cr
#' 2. posteriors for each extracted CN-feature across samples to a set of component models \cr
#' calculates signature exposures.  \cr
#'
#' Said more concretely, given a sample-by-component matrix and a component-by-signature matrix, finds the linear combination decomposition (LCD).
#' See the YAPSA::LCD function for details.
#'
#' @param sample_by_component Matrix. Provide a matrix of posteriors for each sample to each component in the model used.
#' @param component_by_signature String or object of type `NMFfitX1`. If a string, the function expects a path to CN-Signatures (an `NMFfitX1` object) stored as an rds.
#' @returns A matrix of the exposures for each sample to each signature.
#'
#' @export
QuantifySignatures <- function(sample_by_component,
                               component_by_signature = NULL) {

  # Validate input
  stopifnot("No signatures provided. Please provide a matrix, NMFfitX1 object, or filepath of signatures." =
              !is.null(component_by_signature))

  if ('character' %in% class(component_by_signature)) {
    if (file.exists(component_by_signature)) {
      component_by_signature <- NMF::basis(readRDS(file = component_by_signature))
    } else {
      stop("Unable to resolve path in component_by_signature variable to a file.")
    }
  } else if (grepl("NMFfitX1", class(component_by_signature), fixed = TRUE)) {
    component_by_signature <- NMF::basis(component_by_signature)
  }

  signature_by_sample <- YAPSA::LCD(t(sample_by_component),
                                    YAPSA:::normalize_df_per_dim(component_by_signature, 2))
  signature_by_sample <- NormaliseMatrix(signature_by_sample)
  signature_by_sample
}


#' Create CN-Signatures
#'
#' @description
#' Create copy-number signatures by decomposing a sample-by-component matrix.
#'
#' @details
#' Given a sample-by-component sum-of-posteriors matrix and a matrix rank, decompose into two distinct matrices.
#' One matrix will be a sample-by-signature matrix and the other a signature-by-component matrix.
#' In this package, we refer to the signature-by-component matrix as a set of 'copy-number signatures'.
#'
#' @param sample_by_component Matrix. Expects a sum-of-posteriors probability matrix.
#' @param nsig Integer. The number of signatures. Likely decided by examining the plot from `ChooseNumberSignatures()`.
#' @param seed Integer. (flexmix param) The random seed to use while modelling.
#' @param nmfalg Character string. (NMF param) The specification of the NMF algorithm.
#' @param cores Integer. (NMF param) The number of cores to use when running NMF.
#'
#' @returns An NMFfit object. Contains the two matrices post decomposition.
#'
#'
#' @export
GenerateSignatures <- function(sample_by_component, nsig,
                               seed = 77777, nmfalg="brunet",
                               nruns = 1000, cores = 4) {
    NMF::nmf(t(sample_by_component), nsig, seed = seed, nrun = nruns,
             method = nmfalg,.opt = paste0("p", cores) )
}


#' Choose an Optimal Number of Signatures
#'
#' @description
#' Uses the NMF package to estimate the rank of the matrix provided.
#' In the context of signature finding, this is equivalent to determining the optimal number of signatures.
#'
#' @details
#' The NMF Package (Brunet algorithm specification) is used to deconvolute the sample-by-component sum-of-posteriors matrix into a sample-by-signature matrix and a signature-by-component matrix.
#' By default a search interval of 3-12 is used for the rank.
#' To compare/for context the input matrix is also shuffled and the same procedure performed. The results are co-plotted to give a baseline.
#'
#' @param sample_by_component Matrix. Expects a sum-of-posteriors probability matrix.
#' Likely the output from `GenerateSampleByComponentMatrix()`.
#' @param min_sig Integer. Minimum number of signatures to consider.
#' @param max_sig Integer. Maximum number of signatures to consider.
#' @param iter Integer. The number of runs of NMF to perform at each rank with a random seed.
#' @param cores Integer. (NMF param) The number of cores to use when running NMF.
#' @returns A figure. A figure plotting the cophenetic distance, the dispersion, sparseness and silhouette score.
#'
#' @export
ChooseNumberSignatures <- function(sample_by_component,
                                   min_sig = 3, max_sig = 12,
                                   iter = 100, cores = 4,
                                   nmfalg = "brunet",
                                   seed = 77777) {

  estim.r <- NMF::nmfEstimateRank(t(sample_by_component), min_sig:max_sig,
                                  seed = seed, nrun = iter, verbose = FALSE,
                                  method = nmfalg,
                                  .opt = list(shared.memory = FALSE,
                                              paste0("p", cores)))
  V.random <- NMF::randomize(t(sample_by_component))
  estim.r.random <- NMF::nmfEstimateRank(V.random, min_sig:max_sig,
                                         seed = seed, nrun = iter,
                                         verbose = FALSE, method = nmfalg,
                                         .opt = list(shared.memory = FALSE,
                                                     paste0("p", cores)))

  output <- list(estrank = estim.r, randomized_estrank = estim.r.random)
  return(output)
}


#' Extract Absolute Copy-Number Features
#'
#' @description
#' Extract genome-wide copy-number features data from either a list of dataframes or a QDNAseq S4 object for 1 or more samples.
#' This function is intended to be run on absolute CN data.
#'
#' @details
#' The extracted copy-number features are:
#' 1. Breakpoint count per 10MB - `bp10MB`
#' 2. Copy-number value of each segment - `copynumber`
#' 3. Copy-number difference between adjacent segments - `changepoint`
#' 4. Breakpoint count per chromosome arm - `bpchrarm`
#' 5. Lengths of oscillating CN segment chains - `osCN`
#' 6. Size of copy-number segments in base-pairs - `segsize`
#'
#' Extra features: \cr
#' 7. Minimum number of chromosomes (a count) needed to account for 50% of CN changes in a sample - `nc50` \cr
#' 8. Distance in base pairs of each breakpoint to the centromere - `cdist`
#'
#' @param CN_data List of dataframes or S4 QDNAseq object. Segmented copy-number data for for 1 or more samples.
#' If input is a list of dataframes, columns should be: \cr
#' 1. chromosome
#' 2. start
#' 3. end
#' 4. segVal
#'
#' @param genome Character string. The reference genome used for alignment. \cr
#' Options: 'hg19', 'hg38'
#' @param cores Integer. The number of cores to use for parallel processing.
#' @param log_gaussians Logical. If TRUE, take the log1p of the extracted CN-features for those being modeled by gaussians. (segsize, changepoint, copynumber)
#' @param extra_features Logical. If TRUE, extracts CN-feature data for two more features: nc50, and cdist.
#' @returns A list. Each list element contains feature data for a single feature.
#'
#' @export
ExtractCopyNumberFeatures <- function(CN_data,
                                      genome,
                                      cores = 1,
                                      log_gaussians = FALSE,
                                      extra_features = FALSE) {

  # Get chromosome and centromere locations
  if (genome == 'hg19') {
    chrlen <- as.data.frame(hg19.chrom.sizes[1:24,])
    centromeres <- gaps.hg19[gaps.hg19$type == "centromere",]
  } else if (genome == 'hg38') {
    chrlen <- as.data.frame(hg38.chrom.sizes[1:24,])
    centromeres <- gaps.hg38[gaps.hg38$type == "centromere",]
  }

  # Extract CN-Features
  if (cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)

    temp_list = foreach::foreach(i=1:6) %dopar% {
        if (i == 1) {
          list(segsize = GetSegSize(CN_data) )
        } else if (i == 2) {
          list(bp10MB = GetBPNum(CN_data,chrlen) )
        } else if (i == 3) {
          list(osCN = GetOscilation(CN_data,chrlen) )
        } else if (i == 4) {
          list(bpchrarm = GetBPChromArmCounts(CN_data,centromeres,chrlen) )
        } else if (i == 5) {
          list(changepoint = GetChangePointCN(CN_data) )
        } else {
          list(copynumber = GetCN(CN_data))
        }
    }
    unlist( temp_list, recursive = FALSE )
  } else {
    segsize <- GetSegSize(CN_data)
    bp10MB <- GetBPNum(CN_data, chrlen)
    osCN <- GetOscilation(CN_data, chrlen)
    bpchrarm <- GetBPChromArmCounts(CN_data, centromeres, chrlen)
    changepoint <- GetChangePointCN(CN_data)
    copynumber <- GetCN(CN_data)

    if (isTRUE(log_gaussians)) {
      segsize$value = log1p(as.numeric(segsize$value))
      changepoint$value = log1p(as.numeric(changepoint$value))
      copynumber$value = log1p(as.numeric(copynumber$value))
    }
    features <- list(segsize = segsize,
                     bp10MB = bp10MB,
                     osCN = osCN,
                     bpchrarm = bpchrarm,
                     changepoint = changepoint,
                     copynumber = copynumber)

    if (extra_features) {
      nc50 <- GetNC50(CN_data)
      cdist <- GetDistsFromCentromere(CN_data, centromeres, chrlen)
      features[["nc50"]] <- nc50
      features[["cdist"]] <- cdist
    }
  }
  return(features)
}


#' Extract Relative Copy-Number Features
#'
#' @description
#' Extract genome-wide copy-number features from either a list of dataframes or a QDNAseq S4 object for 1 or more samples.
#' This function is intended to be run on relative CN data.
#'
#' @details
#' This function is identical to the absolute calling equivalent other than for three features.
#' The `osCN`, `changepoint`, and `copynumber` features require slightly different modelling at the relative scale.
#'
#' The extracted copy-number features are:
#' 1. Breakpoint count per 10MB - `bp10MB`
#' 2. Copy-number value of each segment - `copynumber`
#' 3. Copy-number difference between adjacent segments - `changepoint`
#' 4. Breakpoint count per chromosome arm - `bpchrarm`
#' 5. Lengths of oscillating CN segment chains - `osCN`
#' 6. Size of copy-number segments in base-pairs - `segsize`
#'
#' Extra features: \cr
#' 7. Minimum number of chromosomes (a count) needed to account for 50% of CN changes in a sample - `nc50` \cr
#' 8. Distance in base pairs of each breakpoint to the centromere - `cdist`
#'
#' @param CN_data List of datafames or S4 QDNAseq object. Segmented relative copy-number data for 1 or more samples.
#' If input is a list of dataframes, columns should be: \cr
#' 1. chromosome
#' 2. start
#' 3. end
#' 4. segVal
#'
#' @param genome Character string. The reference genome used for alignment. \cr
#' Options: 'hg19', 'hg38'
#' @param cores Integer. The number of cores to use for parallel processing.
#' @param log_gaussians Logical. If TRUE, take the log1p of the extracted CN-features for those being modeled by gaussians. (segsize, changepoint, copynumber)
#' @param extra_features Logical. If TRUE, extracts CN-feature data for two more features: nc50, and cdist.
#' @returns A list. Each list element contains feature data for a single feature.
#'
#' @export
ExtractRelativeCopyNumberFeatures <- function(CN_data,
                                              genome,
                                              cores = 1,
                                              log_gaussians = FALSE,
                                              extra_features = FALSE) {

  # Get chromosome and centromere locations
  if (genome == 'hg19') {
    chrlen <- as.data.frame(hg19.chrom.sizes[1:24,])
    centromeres <- gaps.hg19[gaps.hg19$type == "centromere",]
  } else if (genome == 'hg38') {
    chrlen <- as.data.frame(hg38.chrom.sizes[1:24,])
    centromeres <- gaps.hg38[gaps.hg38$type == "centromere",]
  }

  # Extract CN-Features
  if (cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)

    temp_list = foreach::foreach(i=1:6) %dopar% {
      if (i == 1) {
        list(segsize = GetSegSize(CN_data) )
      } else if (i == 2) {
        list(bp10MB = GetBPNum(CN_data,chrlen) )
      } else if (i == 3) {
        list(osCN = GetRelativeOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        list(bpchrarm = GetBPChromArmCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        list(changepoint = GetRelativeChangePointCN(CN_data) )
      } else {
        list(copynumber = GetRelativeCN(CN_data))
      }
    }
    unlist( temp_list, recursive = FALSE )
  } else {
    segsize <- GetSegSize(CN_data)
    bp10MB <- GetBPNum(CN_data, chrlen)
    osCN <- GetRelativeOscilation(CN_data, chrlen)
    bpchrarm <- GetBPChromArmCounts(CN_data, centromeres, chrlen)
    changepoint <- GetRelativeChangePointCN(CN_data)
    copynumber <- GetRelativeCN(CN_data)

    if (isTRUE(log_gaussians)) {
      segsize$value = log1p(as.numeric(segsize$value))
      changepoint$value = log1p(as.numeric(changepoint$value))
      copynumber$value = log1p(as.numeric(copynumber$value))
    }
    features <- list(segsize = segsize,
                     bp10MB = bp10MB,
                     osCN = osCN,
                     bpchrarm = bpchrarm,
                     changepoint = changepoint,
                     copynumber = copynumber)

    if (extra_features) {
      nc50 <- GetNC50(CN_data)
      cdist <- GetDistsFromCentromere(CN_data, centromeres, chrlen)
      features[["nc50"]] <- nc50
      features[["cdist"]] <- cdist
    }
  }
  return(features)
}


#' Fit Mixture Models for each CN-Feature
#'
#' @description
#' Perform mixture modelling on CN-features using either a mixture of gaussians or poissons.
#'
#' @details
#' The segment size, changepoint copy number, and segment copy-number value CN-features are modelled with a mixture of Gaussians.
#' For the breakpoint count per 10MB, length of segments with oscillating copy-number, and breakpoint count per chromosome a mixture of Poissons is used instead.
#' Mixture modelling is done using the FlexMix package.
#'
#' @param CN_features A list. The output from either the ExtractRelativeCopyNumberFeatures or ExtractCopyNumberFeatures functions.
#' @param seed Integer. (flexmix param) The random seed to use while modelling.
#' @param min_comp Integer. (flexmix param) The minimum number of components for each CN-feature to consider.
#' @param max_comp Integer. (flexmix param) The maximum number of components for each CN-feature to consider.
#' @param min_prior Numeric. (flexmix param) Minimum prior probability of clusters, components falling below this threshold are removed during the iteration.
#' @param model_selection Integer or character. (flexmix param) Which model to get. Choose by number or name of the information criterion.
#' @param nrep Integer. (flexmix param) The number of times flexmix is run for each k (number of components).
#' @param niter Integer. (flexmix param) The maximum number of iterations for the EM-algorithm.
#' @param cores Integer. The number of cores to use for parallel processing.
#' @param featsToFit Integer vector. The CN-features to fit.
#'
#' @returns A list of flexmix objects. One for each CN-feature.
#'
#' @export
FitMixtureModels <- function(CN_features, seed = 77777, min_comp = 2,
                             max_comp = 10, min_prior = 0.001, model_selection = "BIC",
                             nrep = 1, niter = 1000, cores = 1, featsToFit = seq(1, 6)) {

  if (cores > 1) {

    require(foreach)
    doMC::registerDoMC(cores)
    temp_list = foreach(i=1:6) %dopar% {

          if (i == 1 & i %in% featsToFit ) {

              dat<-as.numeric(CN_features[["segsize"]][,2])
              list( segsize = FitComponent(dat,seed=seed,model_selection=model_selection,
                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

          } else if (i == 2 & i %in% featsToFit ) {

              dat<-as.numeric(CN_features[["bp10MB"]][,2])
              list( bp10MB = FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

          } else if (i == 3 & i %in% featsToFit ) {

              dat<-as.numeric(CN_features[["osCN"]][,2])
              list( osCN = FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

          } else if (i == 4 & i %in% featsToFit ) {

              dat<-as.numeric(CN_features[["bpchrarm"]][,2])
              list( bpchrarm = FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

          } else if (i == 5 & i %in% featsToFit ) {

              dat<-as.numeric(CN_features[["changepoint"]][,2])
              list( changepoint = FitComponent(dat,seed=seed,model_selection=model_selection,
                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp) )

          } else if (i == 6 & i %in% featsToFit) {

              dat<-as.numeric(CN_features[["copynumber"]][,2])
              list( copynumber = FitComponent(dat,seed=seed,model_selection=model_selection,
                  nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000) )
          }
      }
      unlist( temp_list, recursive = FALSE )
  } else {
      dat<-as.numeric(CN_features[["segsize"]][,2])
      segsize_mm<-FitComponent(dat,seed=seed,model_selection=model_selection,
                       min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

      dat<-as.numeric(CN_features[["bp10MB"]][,2])
      bp10MB_mm<-FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                              min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

      dat<-as.numeric(CN_features[["osCN"]][,2])
      osCN_mm<-FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                            min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

      dat<-as.numeric(CN_features[["bpchrarm"]][,2])
      bpchrarm_mm<-FitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

      dat<-as.numeric(CN_features[["changepoint"]][,2])
      changepoint_mm<-FitComponent(dat,seed=seed,model_selection=model_selection,
                                   min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp)

      dat<-as.numeric(CN_features[["copynumber"]][,2])
      copynumber_mm<-FitComponent(dat,seed=seed,model_selection=model_selection,
                              nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.005,niter=2000)

      output <- list(segsize = segsize_mm, bp10MB = bp10MB_mm, osCN = osCN_mm,
                     bpchrarm = bpchrarm_mm, changepoint = changepoint_mm,
                     copynumber = copynumber_mm)

      if (exists("nc50", CN_features)) {
        dat <- as.numeric(CN_features[["nc50"]][,2])
        nc50_mm <- FitComponent(dat, seed = seed,
                                model_selection = model_selection,
                                min_prior = min_prior, niter = niter,
                                nrep = nrep, min_comp = min_comp,
                                max_comp = max_comp)
        output[["nc50"]] <- nc50_mm
      }

      if (exists("cdist", CN_features)) {
        dat <- as.numeric(CN_features[["cdist"]][,2])
        cdist_mm <- FitComponent(dat, dist="pois", seed = seed,
                                model_selection = model_selection,
                                min_prior = min_prior, niter = niter,
                                nrep = nrep, min_comp = min_comp,
                                max_comp = max_comp)
        output[["cdist"]] <- cdist_mm
      }

  }
  return(output)
}


#' Generate sum-of-posteriors probability matrix
#'
#' @description
#' Given a set of extracted copy-number features and mixture models for each feature, generate a sum-of-posteriors matrix.
#'
#' @details
#' For each copy-number event for each sample the posterior probability of belonging to a component is computed.
#' These posterior event vectors are then summed resulting in a sum-of-posterior probabilities vector.
#' All sum-of-posterior vectors are combined into a single patient-by-component sum-of-posterior probabilities matrix.
#'
#' @param CN_features A list. The output from either the `ExtractRelativeCopyNumberFeatures()` or `ExtractCopyNumberFeatures()` functions.
#' @param all_components A list of flexmix objects. One for each CN-feature.
#' Likely the output from `FitMixtureModels()`.
#' @param cores Integer. The number of cores to use for parallel processing.
#' @param rowIter Integer. Number of rows to ingest per iteration if using multiple cores. Otherwise ignored.
#' @param subcores Integer. In parallel mode each CN-feature will be run separately, this is the number for cores to allocate to each of these sub-jobs.
#' @returns A sum-of-posteriors probability matrix.
#'
#' @export
GenerateSampleByComponentMatrix <- function (CN_features, all_components = NULL,
                                             cores = 1, rowIter = 1000, subcores = 2) {

  if ((class(all_components) == 'character') && (file.exists(all_components))) {
    all_components<-readRDS(file = all_components)
  } else if (class(all_components) == 'character') {
    stop(paste0('Component models path not valid. Please fix this path: ', all_components))
  } else if ((class(all_components) == 'list') && (class(all_components[[1]]) != 'flexmix')) {
    stop('Component models object not valid. Expecting list of flexmix S4 objects.')
  }
  if (cores > 1) {
    require(foreach)

    feats = c( "segsize", "bp10MB", "osCN",
               "changepoint", "copynumber", "bpchrarm", "nc50" )
    doMC::registerDoMC(cores)
    full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
      CalculateSumOfPosteriors(CN_features[[feat]],
                               all_components[[feat]],
                               feat, rowIter = rowIter, cores = subcores)
    }
  } else {
    full_mat <- cbind(
      CalculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
      CalculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
      CalculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
      CalculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
      CalculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
      CalculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm")
    )
    if (exists("nc50", CN_features)) {
      sumofpost <- CalculateSumOfPosteriors(CN_features[["nc50"]],all_components[["nc50"]],"nc50")
      full_mat <- cbind(full_mat, sumofpost)
    }
    if (exists("dist", CN_features)) {
      sumofpost <- CalculateSumOfPosteriors(CN_features[["dist"]],all_components[["dist"]],"dist")
      full_mat <- cbind(full_mat, sumofpost)
    }
  }
  rownames(full_mat) <- unique(CN_features[["segsize"]][,1])
  full_mat[is.na(full_mat)] <- 0

  return(full_mat)
}
