suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(pheatmap)
  library(viridis)
  library(QDNAseq)
  library(CGHcall)
  library(cowplot)
})


# Inputs
segs15kb <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/15kb_copyNumbersSegmented.rds'
segs30kb <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rCN_comCNVfilt.rds'
segs30kb <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_gl_rCN.rds'
segs50kb <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-11/Xchr_included/50kb_gl_rCN.rds'
segs100kb <- '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-11/Xchr_included/100kb_gl_rCN.rds'
segs15kb <- readRDS(file = segs15kb)
segs30kb <- readRDS(file = segs30kb)
segs50kb <- readRDS(file = segs50kb)
segs100kb <- readRDS(file = segs100kb)

metadata <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv', sep = '\t', header = TRUE)
# metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
outcomes <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/nsmp_outcomes.tsv', sep = '\t', header = TRUE)
# outcomes$sample_id <- str_replace_all(outcomes$sample_id, "-", ".")
metadata <- metadata[,c(1,2,4:8)] %>% left_join(outcomes, by = c('sample_id' = 'patient_id'))
metadata <- metadata %>% dplyr::filter(dawn_tfri == 1)
# Computing
# glBins15kb <- callBins(segs15kb, ncpus=4)
# glBins30kb <- callBins(segs30kb, ncpus=4)
# glBins50kb <- callBins(segs50kb, ncpus=4)
# glBins100kb <- callBins(segs100kb, ncpus=4)
cgh15kb <- makeCgh(glBins15kb)
cgh30kb <- makeCgh(segs30kb)
cgh50kb <- makeCgh(glBins50kb)
cgh100kb <- makeCgh(glBins100kb)

bad_nsmp <- metadata %>% dplyr::filter(outcome == 'bad')
good_nsmp <- metadata %>% dplyr::filter(outcome == 'good')
endo_samples <- metadata[metadata$tissue == 'endometrium']
endo_samples <- endo_samples[endo_samples$sample_type == 'primary_tumour']
nsmp_samples <- endo_samples %>% dplyr::filter(status %in% c('NSMP'))
p53abn_samples <- endo_samples %>% filter(status %in% c('p53_abn'))
pathwt <- nsmp_samples %>% dplyr::filter(path_AI_pred == 'p53wt')
pathabn <- nsmp_samples %>% dplyr::filter(path_AI_pred == 'p53_abn')

plot(cgh100kb[,NSMP_samples$sample_id])
summaryPlot(cgh30kb[,good_nsmp$sample_id])
summaryPlot(cgh30kb[,bad_nsmp$sample_id])
frequencyPlotCalls(cgh100kb[,NSMP_samples$sample_id])

p <- plot_grid(plot(cgh100kb[,NSMP_samples$sample_id]), nrow = 12, ncol = 2)

library(gridGraphics)
library(grid)
plot(cgh100kb[,"CC-VGH-1239"])
grid.echo()
a <- grid.grab()
plot(cgh100kb[,"CC-VGH-0033a"])
grid.echo()
b <- grid.grab()
figs <- grid.arrange(a, b, nrow = 2, ncol = 1)

par(mfrow = c(2, 1)) 
summaryPlot(segs30kb[,pathabn$sample_id])
summaryPlot(segs30kb[,pathwt$sample_id])
par(mfrow = c(2, 1)) 
summaryPlot(segs15kb[,good_nsmp$sample_id])
summaryPlot(segs15kb[,bad_nsmp$sample_id])
par(mfrow = c(2, 1))
summaryPlot(cgh50kb[,NSMP_samples$sample_id])
summaryPlot(cgh100kb[,NSMP_samples$sample_id])
# frequencyPlotCalls(cgh15kb[,NSMP_samples$sample_id])

par(mfrow = c(4, 1))
plot(cgh100kb[,NSMP_samples$sample_id[5:8]])

par(mfrow = c(12, 2))
plot(cgh100kb[,NSMP_samples$sample_id])

dir_path <- '~/Documents/projects/cn_signatures_shallowWGS/plotting/30kb_unfiltered_rel_segs_plots/'
for (i in 1:dim(segs_obj)[2]) {
  saveFile <- paste0(dir_path, segs_obj@phenoData@data[["name"]][i], '_segs_plot.png')
  png(file = saveFile, width = 2000, height = 1200, pointsize = 22, units = 'px')
  par(mfrow = c(1, 1))
  plot(segs_obj[,i])
  dev.off()
}



testSummaryPlot <- function (x, main='Summary Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh37',... )
{
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  nclass <-3
  if (!is.null(probamp(x))) nclass <- nclass+1 
  if (!is.null(probdloss(x))) nclass <- nclass+1 
  
  chrom.lengths <- .getChromosomeLengths(build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  
  if(nclass==3) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))}
  if(nclass==4) {loss.freq <- rowMeans(probloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}
  if(nclass==5) {loss.freq <- rowMeans(probloss(x))+rowMeans(probdloss(x)); gain.freq <- rowMeans(probgain(x))+rowMeans(probamp(x))}
  
  # remove probabilities of bins that fall below 0.2
  loss.freq[loss.freq < 0.2] <- 0.001
  gain.freq[gain.freq < 0.2] <- 0.001
  browser()
  plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='mean probability', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  if (length(chrom.ends) > 1)
    for (j in names(chrom.ends)[-length(chrom.ends)])
      abline(v=chrom.ends[j], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
}

.getChromosomeLengths <- function(build) {
  build <- as.integer(gsub('[^0-9]', '', build))
  if (build == 34 || build == 16) {
    chromosome.lengths <- c(246127941, 243615958, 199344050, 191731959, 181034922, 170914576, 158545518, 146308819, 136372045, 135037215, 134482954, 132078379, 113042980, 105311216, 100256656, 90041932, 81860266, 76115139, 63811651, 63741868, 46976097, 49396972, 153692391, 50286555)
  } else if (build == 35 || build == 17) {
    chromosome.lengths <- c(245522847, 243018229, 199505740, 191411218, 180857866, 170975699, 158628139, 146274826, 138429268, 135413628, 134452384, 132449811, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49554710, 154824264, 57701691)
  } else if (build == 36 || build == 18) {
    chromosome.lengths <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else {
    chromosome.lengths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  }
  names(chromosome.lengths) <- 1:24
  chromosome.lengths
}
