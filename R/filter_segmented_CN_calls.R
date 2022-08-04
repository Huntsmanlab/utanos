# Split genome into commonly CN-variable regions in endometrial carcinoma and others.
# Filter out commonly CN-variable regions across most human populations
# Input: Segmented QDNAseq Object of samples.
# Output: Two distinct Segmented QDNAseq objects of samples with fewer bins.
#             1) Genome regions that are commonly CN variable in EC
#             2) The rest

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(QDNAseq)
  library(CGHcall)
  library(EnsDb.Hsapiens.v75)
  library(readr)
})


###########################
### Functions
###########################

gene_cr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    if (sum(toupper(queryset$Gene) %in% toupper(line)) > 0) {
      queryset_matches <- c(queryset_matches, i)
    }
  }
  return(targetset[unique(queryset_matches),])
}

annotation_cr <- function(queryset, targetset) {
  queryset_matches <- c()
  for (i in 1:dim(targetset)[1]) {
    line <- str_split(targetset$all[i], pattern = ',')[[1]]
    line <- unique(line[line != ""])
    queryset_matches <- c(queryset_matches, which(toupper(queryset$gene_name) %in% toupper(line)))
  }
  return(queryset[queryset_matches,])
}

segments_to_copy_number <- function(segs, bin_size, genome = 'autosomes') {
  
  # Stop execution if we don't have the required input
  stopifnot(is.list(segs))
  stopifnot(length(segs) != 0)
  stopifnot("chromosome" %in% names(segs[[1]]))
  stopifnot("start" %in% names(segs[[1]]), is.numeric(segs[[1]]$start))
  stopifnot("end" %in% names(segs[[1]]), is.numeric(segs[[1]]$end))
  stopifnot("segVal" %in% names(segs[[1]]), is.numeric(segs[[1]]$segVal))
  # Create template binned genome 
  genome_chrs <- 22
  if (genome == 'Xchr_included') { genome_chrs <- 23}
  chroms <- getChromInfoFromUCSC("hg19") %>% 
                head(genome_chrs) %>%
                dplyr::mutate(nbins = ceiling(size/bin_size))
  chroms$chrom <- sub('chr', '', chroms$chrom)
  genome_template <- data.frame(chromosome = rep(chroms$chrom, chroms$nbins), 
                                start = rep(rep(1,dim(chroms)[1]), chroms$nbins),
                                end = rep(chroms$size, chroms$nbins)) %>% 
                dplyr::group_by(chromosome,start,end) %>% 
                dplyr::mutate(id = row_number()) %>% 
                dplyr::mutate(start = start + ((id-1) * bin_size), 
                              end = start + bin_size - 1) %>% 
                dplyr::select(-id) %>% 
                dplyr::ungroup()
  # Build dataframe of per-bin copy numbers
  out <- c()
  for (name in names(segs)) {
    sample <- segs[[name]] %>% dplyr::mutate(size = end - start + 1) %>% 
              dplyr::mutate(nbins = size/bin_size)
    sample <- as.data.frame(lapply(sample, rep, sample$nbins)) %>% 
              group_by(start,end,segVal) %>% 
              dplyr::mutate(id = row_number()) %>% 
              dplyr::mutate(chromosome = as.character(chromosome),
                            start = start + ((id-1) * bin_size), 
                            end = start + bin_size - 1) %>% 
              dplyr::select(-c(size,nbins,id)) %>% 
              dplyr::ungroup()
    per_bin_copy_numbers <- left_join(genome_template, sample, by = c('chromosome', 'start', 'end'))
    per_bin_copy_numbers$sample_id <- name
    out <- rbind(out, per_bin_copy_numbers)
  }
  names(out) <- c("chromosome", "start", "end", "state", "sample_id")
  return(out)
}

StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

GetGRangesFromEnsDb <- function(
  ensdb,
  standard.chromosomes = TRUE,
  biotypes = c(listGenebiotypes(ensdb)),
  verbose = TRUE
) {
  # Code taken from Signac
  # https://github.com/timoast/signac/blob/master/R/utilities.R
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }
  
  # extract genes from each chromosome
  if (verbose) {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
    })
  }
  
  # combine
  tx <- do.call(what = c, args = tx)
  tx <- tx[tx$gene_biotype %in% biotypes]
  return(tx)
}


###########################
### Function calls
###########################

#####
##### Get input CN calls
#####
input_obj <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_bins_rCN.rds'
input_obj <- readRDS(input_obj)
# Get commonly CN variable regions in hg19
com_CNVs <- '~/Documents/projects/cn_signatures_shallowWGS/data/external_datasets/GRCh37_hg19_common_variants_20200225.tsv'
com_CNVs <- read.table(file = com_CNVs, header = TRUE, sep = "\t")
# Get commonly CN variable regions in EC
EC_com_CNVs <- '~/Documents/projects/cn_signatures_shallowWGS/data/external_datasets/common_hg19_EC_CNA_genes_TCGA.txt'
EC_com_CNVs <- data.table::fread(file = EC_com_CNVs, header = TRUE, sep = "\t", fill = TRUE)
# Cross-referencing dataset
cr_ref <- '~/Documents/projects/cn_signatures_shallowWGS/data/external_datasets/gene_names_symbols_crossreference.txt'
cr_ref <- read_tsv(file = cr_ref, col_names = TRUE)

#####
##### Input transformations
#####
EC_com_CNVs$Freq <- as.numeric(str_sub(EC_com_CNVs$Freq, 1, nchar(EC_com_CNVs$Freq)-1))
EC_com_CNVs <- EC_com_CNVs %>% dplyr::filter(Freq > 2) %>% dplyr::select(Gene, CNA, Freq)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
# restrict cross-referencing dataset to applicable columns
colnames(cr_ref) <- c('hgnc_id', 'approved_symbol', 'previous_symbols', 'alias_symbols', 'Chromosome', 'LNCipedia')
cr_ref <- as.data.frame(cr_ref) %>% 
            dplyr::select(approved_symbol, previous_symbols, alias_symbols, LNCipedia) %>% 
            mutate(all = paste(approved_symbol,
                               replace_na(previous_symbols, ''),
                               replace_na(alias_symbols, ''),
                               replace_na(LNCipedia, ''),
                               sep = ','))
# Filter com_CNVs to useful cols
com_CNVs <- com_CNVs %>% 
            dplyr::select(chr, start, end, samplesize, observedlosses, observedgains, genes) %>% 
            dplyr::mutate(size = end-start) %>% 
            dplyr::filter(((observedlosses + observedgains)/samplesize > 0.01) & (observedgains+observedlosses > 9) & (size > 10000))

#####
##### Generate new filtered Segmented CN-call objects
#####
found_genes <- gene_cr(EC_com_CNVs, cr_ref)
gene_annos <- annotation_cr(annotations, found_genes)
sample_ranges <- StringToGRanges(featureNames(input_obj), sep = c(':', '-'))
seqlevelsStyle(sample_ranges) <- 'UCSC'
genome(sample_ranges) <- "hg19"
hits <- findOverlaps(sample_ranges, gene_annos)
cEC_CNV_obj <- input_obj[unique(hits@from),]
new_obj <- input_obj[setdiff(1:length(sample_ranges), unique(hits@from)),]

gcom_CNVs <- makeGRangesFromDataFrame(com_CNVs)
seqlevelsStyle(gcom_CNVs) <- 'UCSC'
genome(gcom_CNVs) <- "hg19"
hits <- findOverlaps(sample_ranges, gcom_CNVs, minoverlap = 5000)
cCNV_obj <- input_obj[unique(hits@from),]
ccnvFILT_obj <- input_obj[setdiff(1:length(sample_ranges), unique(hits@from)),]

#####
##### Save output
#####
saveRDS(ccnvFILT_obj, file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_comCNVfilt.rds')
saveRDS(cEC_CNV_obj, file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_eccnvs.rds')
saveRDS(new_obj, file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_NOeccnvs.rds')
