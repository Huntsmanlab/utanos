suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# test if there is one argument: if not, return an error
if (length(args)==0) {
  stop("An input path, output path, and bin size must be provided.", call.=FALSE)
} else if (length(args)==3) {
  listofbams <- args[1]
  outpath <- args[2]
  binsize <- args[3]
} else {
  stop("An input path, output path, and bin size must be provided. Nothing else", call.=FALSE)
}

# Declare paths
bin_annos <- paste0('/projects/molonc/huntsman_lab/madouglas/PE150_bin_annotation/qdnaseq_docs/1000genomes_sWGS_bins/hg19_bins_150bp_', binsize, '_SE.rds')

# Read-in bams list
bams <- read.table(listofbams, header = F, stringsAsFactors = F)
bams <- c(bams$V1)

# Load Bin Annotations
bins <- readRDS(file = bin_annos)

# Read in bams
sample_names <- word(bams,9,sep = "\\/")
readCounts <- binReadCounts(bins, bamfiles=bams, bamnames = word(sample_names,1,sep = "\\."))

# QDNAseq processing and CN calling
readCountsFiltered <- applyFilters(readCounts, chromosomes="Y", residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# QDNAseq Segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

# Fix sample names
copyNumbersSegmented@phenoData@data[["name"]] <- word(copyNumbersSegmented@phenoData@data[["name"]],1,sep = "\\.")

# Save results
saveRDS(copyNumbersSegmented, file = file.path(outpath,paste0(binsize, "_copyNumbersSegmented.rds")))
exportBins(copyNumbersSmooth, file = file.path(outpath,paste0(binsize, "_copyNumbersSmooth.tsv")), format="tsv", logTransform=FALSE)
exportBins(copyNumbersSegmented, file = file.path(outpath,paste0(binsize, "_copyNumbersSegmented.tsv")), type = 'segments', format="tsv", logTransform=FALSE)
