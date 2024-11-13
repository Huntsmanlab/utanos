# helpersignaturefunctions.R
# The functions in this file were copied from another code repository and in several cases modified.
# The original code can be found here: https://bitbucket.org/britroc/cnsignatures/src/master/
# The original code was a part of a publication in Nat. Gen. (August 2018)

#' @export
FitComponent <- function(dat, dist = "norm", seed = 77777,
                         model_selection = "BIC",
                         min_prior = 0.001,
                         niter = 1000,
                         nrep = 1,
                         min_comp = 2, max_comp = 10) {

  control <- new("FLXcontrol")
  control@minprior <- min_prior
  control@iter.max <- niter
  set.seed(seed)
  if (dist == "norm") {
    if (min_comp == max_comp) {
      fit <- flexmix::flexmix(dat ~ 1,
                              model = flexmix::FLXMCnorm1(),
                              k = min_comp, control = control)
    } else {
      fit <- flexmix::stepFlexmix(dat ~ 1,
                                  model = flexmix::FLXMCnorm1(),
                                  k = min_comp:max_comp,
                                  nrep = nrep, control = control)
      fit <- flexmix::getModel(fit, which = model_selection)
    }
  } else if (dist == "pois") {
    if (min_comp == max_comp) {
      fit <- flexmix::flexmix(dat ~ 1,
                              model = flexmix::FLXMCmvpois(),
                              k = min_comp, control = control)
    } else {
      fit <- flexmix::stepFlexmix(dat ~ 1,
                                  model = flexmix::FLXMCmvpois(),
                                  k = min_comp:max_comp,
                                  nrep = nrep, control = control)
      fit <- flexmix::getModel(fit, which = model_selection)
    }
  }
  fit
}

MultiSeedFitComponent<-function(dat,dist="norm",num_seed=100,model_selection="BIC",
                                min_prior=0.001,niter=1000,nrep=1,min_comp=2,max_comp=10, cores = 1) {
  control<-new("FLXcontrol")
  control@minprior<-min_prior
  control@iter.max<-niter
  if (cores > 1) {
    require(foreach)
    require(doMC)

    registerDoMC(cores)

    all_bic = foreach(i = 1:num_seed) %dopar% {
      set.seed(i)

      # Use tryCatch to handle errors
      fit <- tryCatch({
        if (min_comp == max_comp) {
          if (dist == "pois") {
            flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCmvpois(), k = min_comp, control = control)
          } else {
            flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCnorm1(), k = min_comp, control = control)
          }
        } else {
          if (dist == "pois") {
            flexmix::stepFlexmix(dat ~ 1, model = flexmix::FLXMCmvpois(), k = min_comp:max_comp, nrep = nrep, control = control)
          } else {
            flexmix::stepFlexmix(dat ~ 1, model = flexmix::FLXMCnorm1(), k = min_comp:max_comp, nrep = nrep, control = control)
          }
        }
      }, error = function(e) {
        message(paste("Error occurred during fitting for seed", i, ":", e$message))
        return(NULL)  # Return NULL if an error occurs
      })

      # Check if fitting was successful
      if (!is.null(fit)) {
        bic <- t(as.data.frame(BIC(fit)))
        rownames(bic) <- as.character(i)
        bic
      } else {
        bic <- NA
        rownames(bic) <- as.character(i)
        bic
      }
    }
  } else { #only one core
    all_bic <- list()
    for (i in 1:num_seed) {
      set.seed(i)
      fit <- tryCatch({
        if (min_comp == max_comp) {
          if (dist == "pois") {
            flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCmvpois(), k = min_comp, control = control)
          } else {
            flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCnorm1(), k = min_comp, control = control)
          }
        } else {
          if (dist == "pois") {
            flexmix::stepFlexmix(dat ~ 1, model = flexmix::FLXMCmvpois(), k = min_comp:max_comp, nrep = nrep, control = control)
          } else {
            flexmix::stepFlexmix(dat ~ 1, model = flexmix::FLXMCnorm1(), k = min_comp:max_comp, nrep = nrep, control = control)
          }
        }
      }, error = function(e) {
        message(paste("Error occurred during fitting for seed", i, ":", e$message))
        return(NULL)  # Return NULL or an alternative value if an error occurs
      })

      # Check if fitting was successful
      if (!is.null(fit)) {
        bic <- t(as.data.frame(BIC(fit)))
        rownames(bic) <- as.character(i)
        all_bic[[i]] <- bic
      }
      else {
        all_bic[[i]] <- NA  # Or another placeholder value if fitting failed
      }
    }
  }
  if (min_comp!=max_comp) {
    all_bic <- all_bic[sapply(all_bic, length) == max_comp-min_comp+1]
  }
  all_bic <- do.call(rbind, all_bic)
  #choose best seed and number of comp
  if (min_comp == max_comp) {
    all_bic_df <- as.data.frame(all_bic)
    min_bic_index <- which.min(all_bic_df[, 1])
    best_seed <- rownames(all_bic_df)[min_bic_index]
    best_comp <- as.integer(best_comp)
    best_seed <- as.integer(best_seed)
    set.seed(best_seed)
    if (dist == "pois") {
      fit <-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCmvpois(),k=min_comp,control=control)
    } else {
      fit<-flexmix::flexmix(dat ~ 1,model=flexmix::FLXMCnorm1(),k=min_comp,control=control)
    }

  } else {
    all_bic_df <- as.data.frame(all_bic)
    min_bic_per_seed <- apply(all_bic_df, 1, min)
    min_bic_index <- which.min(min_bic_per_seed)
    best_seed <- rownames(all_bic_df)[min_bic_index]
    best_comp <- names(all_bic_df)[which.min(all_bic_df[min_bic_index, ])]
    best_comp <- as.integer(best_comp)
    best_seed <- as.integer(best_seed)
    set.seed(best_seed)
    if (dist == "pois") {
      fit <- flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCmvpois(), k = best_comp, control = control)
    } else {
      fit <- flexmix::flexmix(dat ~ 1, model = flexmix::FLXMCnorm1(), k = best_comp, control = control)
    }
  }
  colnames(all_bic_df) <- paste0("component_", colnames(all_bic_df))
  rownames(all_bic_df) <- paste0("seed_", rownames(all_bic_df))

  return(list(fit, all_bic_df))
}

CalculateSumOfPosteriors <- function(CN_feature, components,name,
                                     rowIter = 1000, cores = 1) {
  if (cores > 1) {
    require(foreach)
    require(doMC)

    len = dim(CN_feature)[1]
    iters = floor( len / rowIter ) - 1
    lastiter = iters + 1

    registerDoMC(cores)
    curr_posterior = foreach(i = 0:iters,
                             .combine = rbind) %dopar% {
      start = i * rowIter + 1
      if (i != lastiter) { end = (i + 1) * rowIter } else { end = len }
      flexmix::posterior(components,
                         data.frame(dat = as.numeric(CN_feature[start:end,2])))
    }

    # Handle the last chunk if it exists and len is not perfectly divisible by rowiter
    if (len %% rowIter != 0) {
      start <- lastiter * rowIter + 1
      end <- len

      last_chunk <- flexmix::posterior(components, data.frame(dat = as.numeric(CN_feature[start:end, 2])))
      curr_posterior <- rbind(curr_posterior, last_chunk)
    }
  } else {
    curr_posterior <- flexmix::posterior(components,
                                         data.frame(dat = as.numeric(CN_feature[,2])))
  }
  mat <- cbind(CN_feature, curr_posterior)
  posterior_sum <- c()

  # note - foreach and parallelization doesn't make the following code faster.
  for (i in unique(mat$ID)) {
      posterior_sum <- rbind(posterior_sum,
                             colSums(mat[mat$ID == i, c(-1, -2)]))
  }
  params <- flexmix::parameters(components)

  if (!is.null(nrow(params))) {
    posterior_sum <- posterior_sum[, order(params[1,])]
  } else {
    posterior_sum <- posterior_sum[, order(params)]
  }

  colnames(posterior_sum) <- paste0(name, 1:ncol(posterior_sum))
  rownames(posterior_sum) <- rownames(unique(mat$ID))
  posterior_sum
}

GetSegSize <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
    seglen <- (as.numeric(segTab$end) - as.numeric(segTab$start))
    seglen <- seglen[seglen > 0]
    out <- rbind(out, cbind(ID = rep(i, length(seglen)),
                            value = seglen))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

GetBPNum <- function(cn_profiles, chrlen) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    chrs <- unique(segTab$chromosome)
    allBPnum <- c()
    for (c in chrs) {
      currseg <- segTab[segTab$chromosome == c,]
      intervals <- seq(1,
                       chrlen[chrlen[,1] == paste0("chr", c), 2] + 10000000,
                       10000000)
      res <- hist(as.numeric(currseg$end[-nrow(currseg)]),
                  breaks = intervals, plot = FALSE)$counts
      allBPnum <- c(allBPnum, res)
    }
    out <- rbind(out, cbind(ID = rep(i, length(allBPnum)),
                              value = allBPnum))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

GetOscilation <- function(cn_profiles, chrlen) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    chrs <- unique(segTab$chromosome)
    segTab <- segTab[!is.na(segTab$segVal),]
    oscCounts <- c()
    for (c in chrs) {
      currseg <- segTab[segTab$chromosome == c,"segVal"]
      currseg <- round(as.numeric(currseg))                          # Change to currseg$segVal from currseg in order to work with datatables (but that usually screws up more stuff)
      if (length(currseg) > 3) {
        prevval <- currseg[1]
        count = 0
        for (j in 3:length(currseg)) {
          if ((currseg[j] == prevval) & (currseg[j] != currseg[j-1])) {
            count <- count + 1
          } else {
            oscCounts <- c(oscCounts, count)
            count <- 0
          }
          prevval <- currseg[j-1]
        }
      }
    }
    out <- rbind(out, cbind(ID = rep(i, length(oscCounts)),
                            value = oscCounts))
    if (length(oscCounts) == 0) {
      out <- rbind(out, cbind(ID = i, value = 0))
    }
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}


GetRelativeOscilation <- function(cn_profiles, chrlen) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    chrs <- unique(segTab$chromosome)
    segTab <- segTab[!is.na(segTab$segVal),]
    # Rather than being rounded to integers...
    # In the relative case we round to the nearest 0.1 to determine an oscillation
    # Change to currseg$segVal from currseg in order to work with datatables (but that usually screws up more stuff)
    segTab$segVal <- RoundToNearest(log(as.numeric(segTab$segVal)), 0.1)
    oscCounts <- c()
    for (c in chrs) {
      currseg <- segTab[segTab$chromosome == c, "segVal"]
      if (length(currseg) > 3) {
        prevval <- currseg[1]
        count <- 0
        for (j in 3:length(currseg)) {
          # Adapt the if statement to handle the decimals -> fudge factor of 0.1 on oscillation detection
          if (dplyr::between(currseg[j], prevval - 0.1, prevval + 0.1) &
              (currseg[j] != currseg[j - 1]) &
              (prevval != currseg[j - 1])) {
            count <- count + 1
          } else {
            oscCounts <- c(oscCounts, count)
            count <- 0
          }
          prevval<-currseg[j-1]
        }
      }
    }
    out <- rbind(out, cbind(ID = rep(i, length(oscCounts)), value = oscCounts))
    if (length(oscCounts) == 0) {
      out <- rbind(out, cbind(ID = i, value = 0))
    }
  }
  rownames(out) <- NULL
  data.frame(out,stringsAsFactors = F)
}


GetBPChromArmCounts <- function (cn_profiles, centromeres, chrlen) {
  out <- c()
  samps <- GetSampNames(cn_profiles)

  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[,which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    chrs <- unique(segTab$chromosome)
    all_dists <- c()

    for (c in chrs) {
      if (nrow(segTab) > 1) {
        starts <- as.numeric(segTab[segTab$chromosome == c, 2])[-1]
        segstart <- as.numeric(segTab[segTab$chromosome == c, 2])[1]
        ends <- as.numeric(segTab[segTab$chromosome == c, 3])
        segend <- ends[length(ends)]
        ends <- ends[-length(ends)]
        centstart <- as.numeric(centromeres$start[substr(centromeres$chromosome,4,5)==c])
        centend <- as.numeric(centromeres$end[substr(centromeres$chromosome,4,5)==c])
        chrend <- chrlen[substr(chrlen[,1],4,5)==c,2]
        ndist <- cbind(rep(NA,length(starts)),rep(NA,length(starts)))
        ndist[starts<=centstart,1] <- (centstart-starts[starts<=centstart])/(centstart-segstart)*-1
        ndist[starts>=centend,1] <- (starts[starts>=centend]-centend)/(segend-centend)
        ndist[ends<=centstart,2] <- (centstart-ends[ends<=centstart])/(centstart-segstart)*-1
        ndist[ends>=centend,2] <- (ends[ends>=centend]-centend)/(segend-centend)
        ndist <- apply(ndist,1,min)

        # If any segments spill over into a centromere an NA is thrown.
        # We don't want to count these seg. ends as breakpoints nor...
        # have later functions error out due to the presence of NA values.
        # So we toss the NA values.
        ndist <- ndist[!is.na(ndist)]

        all_dists <- rbind(all_dists, sum(ndist > 0))
        all_dists <- rbind(all_dists, sum(ndist <= 0))
      }
    }
    if (nrow(all_dists) > 0) {
      out <- rbind(out, cbind(ID = i, value = all_dists[,1]))
    }
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}


GetChangePointCN <- function (cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)

  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segTab <- segTab[!is.na(segTab$segVal),]
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
    chrs <- unique(segTab$chromosome)
    allcp <- c()

    for (c in chrs) {
      currseg <- as.numeric(segTab[segTab$chromosome == c, "segVal"])
      allcp <- c(allcp, abs(currseg[-1] - currseg[-length(currseg)]))
    }
    if (length(allcp) == 0) {
      allcp <- 0 # if there are no changepoints
    }
    out <- rbind(out, cbind(ID = rep(i, length(allcp)), value = allcp))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}


# The minimal number of chromosome with 50% CNV
GetNC50 <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    # Count the breaks per chrom, arrange, then identify the min # needed for 50%
    breaks_per_chrom <- as.data.frame(table(segTab$chromosome)) %>%
                            dplyr::arrange(desc(Freq))
    breaks_per_chrom$cumsum <- cumsum(breaks_per_chrom$Freq)
    nc50_count <- data.frame(
      ID = i,
      value = sum(!breaks_per_chrom$cumsum > max(breaks_per_chrom$cumsum)/2) + 1
      )
    out <- rbind(out, nc50_count)
  }
  return(out)
}


GetRelativeChangePointCN <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[,which(colnames(cn_profiles) == i)])
    }
    else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segTab <- segTab[!is.na(segTab$segVal),]
    segTab$segVal[as.numeric(segTab$segVal) <= 0] <- 0.00001
    segTab$segVal <- log(as.numeric((segTab$segVal))) #log transform segVal for relative CN
    chrs <- unique(segTab$chromosome)
    allcp <- c()
    for (c in chrs) {
      currseg <- as.numeric(segTab[segTab$chromosome == c, "segVal"])
      allcp <- c(allcp, abs(currseg[-1] - currseg[-length(currseg)]))
    }
    if (length(allcp) == 0) {
      allcp <- 0 #if there are no changepoints
    }
    out <- rbind(out, cbind(ID = rep(i, length(allcp)), value = allcp))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}


GetCN <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[,which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segTab <- segTab[!is.na(segTab$segVal),]
    segTab$segVal[as.numeric(segTab$segVal) < 0] <- 0
    cn <- as.numeric(segTab$segVal)
    out <- rbind(out, cbind(ID = rep(i, length(cn)), value = cn))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

# Return relative copy-numbers
# Note: In the relative case we return log transformed segVal
GetRelativeCN <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[,which(colnames(cn_profiles) == i)])
    }
    else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segTab <- segTab[!is.na(segTab$segVal),]
    segTab$segVal[as.numeric(segTab$segVal) <= 0] <- 0.00001
    segTab$segVal <- log(as.numeric((segTab$segVal))) # log transform segVal for relative CN
    relative_cn <- as.numeric(segTab$segVal)
    out <- rbind(out,cbind(ID = rep(i, length(relative_cn)), value = relative_cn))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

GetSampNames <- function(cn_profiles) {
  if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
    samps <- colnames(cn_profiles)
  } else {
    samps <- names(cn_profiles)
  }
  samps
}

GetSegTable <- function(x) {
  dat <- x
  sn <- Biobase::assayDataElement(dat, "segmented")
  fd <- Biobase::fData(dat)
  fd$use -> use
  fdfiltfull <- fd[use,]
  sn <- sn[use,]
  segTable <- c()

  for (c in unique(fdfiltfull$chromosome)) {
    snfilt <- sn[fdfiltfull$chromosome == c]
    fdfilt <- fdfiltfull[fdfiltfull$chromosome == c,]
    sn.rle <- rle(snfilt)
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- fdfilt$start[starts[s]]
      to <- fdfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp),
                                     ncol = 4,
                                     byrow = T),
                              stringsAsFactors = F)
    segTable <- rbind(segTable, segTableRaw)
  }

  colnames(segTable) <- c("chromosome", "start", "end", "segVal")
  segTable
}


GetPloidy <- function(cn_profiles) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[, which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    segLen <- (as.numeric(segTab$end) - as.numeric(segTab$start))
    ploidy <- sum((segLen/sum(segLen)) * as.numeric(segTab$segVal))
    out <- c(out, ploidy)
  }
  data.frame(out, stringsAsFactors = F)
}


NormaliseMatrix <- function(signature_by_sample,
                            sig_thresh = 0.01) {
  norm_const <- colSums(signature_by_sample)
  sample_by_signature <- apply(signature_by_sample, 1, function(x){x/norm_const})
  sample_by_signature <- apply(sample_by_signature, 1, LowerNorm, sig_thresh)
  signature_by_sample <- t(sample_by_signature)
  norm_const <- apply(signature_by_sample, 1, sum)
  sample_by_signature <- apply(signature_by_sample, 2, function(x){x/norm_const})
  signature_by_sample <- t(sample_by_signature)
  signature_by_sample
}

LowerNorm <- function(x, sig_thresh = 0.01) {
  new_x <- x
  for (i in 1:length(x)) {
    if (x[i] < sig_thresh) {
      new_x[i] <- 0
    }
  }
  new_x
}

GetDistsFromCentromere <- function(cn_profiles, centromeres, chrlen) {
  out <- c()
  samps <- GetSampNames(cn_profiles)
  for (i in samps) {
    if (inherits(cn_profiles, "QDNAseqCopyNumbers")) {
      segTab <- GetSegTable(cn_profiles[,which(colnames(cn_profiles) == i)])
    } else {
      segTab <- cn_profiles[[i]]
      colnames(segTab)[4] <- "segVal"
    }
    chrs <- unique(segTab$chromosome)
    all_dists <- c()
    for (c in chrs) {
      if (nrow(segTab) > 1) {
        starts <- as.numeric(segTab[segTab$chromosome == c, 2])[-1]
        segstart <- as.numeric(segTab[segTab$chromosome == c, 2])[1]
        ends <- as.numeric(segTab[segTab$chromosome == c, 3])
        segend <- ends[length(ends)]
        ends <- ends[-length(ends)]
        centstart <- as.numeric(centromeres$start[substr(centromeres$chromosome, 4, 5) == c])
        centend <- as.numeric(centromeres$end[substr(centromeres$chromosome, 4, 5) == c])
        chrend <- chrlen[substr(chrlen[,1], 4, 5) == c, 2]
        ndist <- cbind(rep(NA, length(starts)), rep(NA, length(starts)))
        ndist[starts <= centstart,1] <- (centstart - starts[starts <= centstart])/(centstart - segstart) * -1
        ndist[starts >= centend,1] <- (starts[starts >= centend] - centend)/(segend - centend)
        ndist[ends <= centstart,2] <- (centstart - ends[ends <= centstart])/(centstart - segstart) * -1
        ndist[ends >= centend,2] <- (ends[ends >= centend] - centend)/(segend - centend)
        ndist <- apply(ndist, 1, min)

        # If any segments spill over into a centromere an NA is thrown.
        # We don't want to count these seg. ends as breakpoints nor...
        # have later functions error out due to the presence of NA values.
        # So we toss the NA values.
        ndist <- ndist[!is.na(ndist)]
        ndist <- abs(ndist)
        ndist <- data.frame(ndist)

        all_dists <- rbind(all_dists, ndist)
      }
    }
    if (nrow(all_dists) > 0) {
      out <- rbind(out, cbind(ID = i, value = all_dists$ndist))
    }
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}
