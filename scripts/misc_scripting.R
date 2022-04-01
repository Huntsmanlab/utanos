library(tidyr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)

# Combine survival data
survival <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/survival_data.tsv'
survival <- read.table(file = survival, sep = '\t', header = TRUE)
all_survival <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/survival_data.tsv'
all_survival <- read_tsv(file = all_survival)
metD <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv'
metD <- read.table(file = metD, sep = '\t', header = TRUE)
exp_anno <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/annotation_pool1-12.tsv'
exp_anno <- data.table::fread(file = exp_anno, sep = '\t', header = TRUE, fill = TRUE)
tp53_cel <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/CellularityET_CrossCan/TP53abn-Table 1.tsv'
tp53_cel <- read.table(file = tp53_cel, sep = '\t', header = TRUE)
nsmp_cel <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/CellularityET_CrossCan/NSMP-Table 1.tsv'
nsmp_cel <- read.table(file = nsmp_cel, sep = '\t', header = TRUE)
# met_5to12 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/samples_batches5to12.tsv'
# met_5to12 <- read_tsv(file = met_5to12)

all_survival <- all_survival %>% dplyr::select(study_id, eclass, eclass2)
colnames(all_survival) <- c('sample_id', 'status', 'eclass2')

exp_anno <- exp_anno %>% dplyr::mutate(sample_id = metD$sample_id) %>% 
                         dplyr::select(sample_id, cancer_type, status, batch)
# exp_anno <- exp_anno %>% left_join(all_survival, by=c('sample_id'))

# met_5to12 <- met_5to12 %>% select(study_id, hist, `GSC submission`)
# colnames(met_5to12) <- c('sample_id', 'cancer_type', 'batch')
# met_5to12 <- met_5to12 %>% left_join(all_survival, by=c('sample_id'))
# met_5to12$tissue <- c('endometrium')
# met_5to12$sample_type <- c('primary_tumour')
# met_5to12$batch[met_5to12$batch == 'HRD-pool4'] <- 5
# met_5to12$batch[met_5to12$batch == 'HRD-pool5'] <- 6
# met_5to12$batch[met_5to12$batch == 'HRD-pool6'] <- 7
# met_5to12$batch[met_5to12$batch == 'HRD-pool7'] <- 8
# met_5to12$batch[met_5to12$batch == 'HRD-pool8'] <- 9
# met_5to12$batch[met_5to12$batch == 'HRD-pool9'] <- 10
# met_5to12$batch[met_5to12$batch == 'HRD-pool10'] <- 11
# met_5to12$batch[met_5to12$batch == 'HRD-pool11'] <- 12
# met_5to12 <- met_5to12 %>% transform(batch = as.integer(batch))

tp53_cel <- tp53_cel %>% dplyr::select(Study.ID, Cellularity.ET...., VAF, VAF.1, VAF.2, VAF.3, VAF.4, VAF.5) %>% 
                          transform(VAF.4 = as.character(VAF.4), VAF.5 = as.character(VAF.5))
nsmp_cel <- nsmp_cel %>% dplyr::select(Study.ID, Cellularity.ET...., VAF, VAF.1, VAF.2, VAF.3, VAF.4, VAF.5) %>% 
                          transform(VAF = as.character(VAF), VAF.3 = as.character(VAF.3), VAF.5 = as.character(VAF.5))
cellularity <- bind_rows(tp53_cel, nsmp_cel)
colnames(cellularity) <- c('sample_id', 'cellularity', 'pole.vaf', 'tp53.vaf', 'pik3ca.vaf', 'ctnnb1.vaf', 'pten.vaf', 'kras.vaf')
# cellularity <- cellularity %>% left_join(all_survival, by=c('sample_id'))
  
# Cleaning entries
exp_anno <- exp_anno %>% dplyr::filter(str_sub(sample_id,-1,-1) != 'N')
exp_anno$sample_id[str_sub(exp_anno$sample_id,-2,-1) == '-T'] <- str_sub(exp_anno$sample_id,1,11)[str_sub(exp_anno$sample_id,-2,-1) == '-T']
exp_anno$sample_id[str_sub(exp_anno$sample_id,-2,-1) == '-T'] <- str_sub(exp_anno$sample_id,1,8)[str_sub(exp_anno$sample_id,-2,-1) == '-T']
exp_anno$sample_id[str_sub(exp_anno$sample_id,-1,-1) == 'T'] <- str_sub(exp_anno$sample_id,1,11)[str_sub(exp_anno$sample_id,-1,-1) == 'T']

# metD$sample_id[str_sub(metD$sample_id,-2,-1) == '-T'] <- str_sub(metD$sample_id,1,11)[str_sub(metD$sample_id,-2,-1) == '-T']
# metD$sample_id[str_sub(metD$sample_id,-2,-1) == '-T'] <- str_sub(metD$sample_id,1,8)[str_sub(metD$sample_id,-2,-1) == '-T']
# metD$sample_id[str_sub(metD$sample_id,-1,-1) == 'T'] <- str_sub(metD$sample_id,1,11)[str_sub(metD$sample_id,-1,-1) == 'T']
# metadata_all <- bind_rows(metD, met_5to12)

# Save the highlighting missing cellularities tsv file
exp_anno <- exp_anno %>% left_join(cellularity, by = c('sample_id'))
# metadata_all$cancer_type[metadata_all$cancer_type == 'endometrioid'] <- 'endometrioid endometrial carcinoma'
# metadata_all$cancer_type[metadata_all$cancer_type == 'serous'] <- 'serous endometrial carcinoma'
# metadata_all$cancer_type[metadata_all$cancer_type == 'clear cell'] <- 'clear cell endometrial carcinoma'
# metadata_all$cancer_type[metadata_all$cancer_type == 'SEIC'] <- 'serous endometrial intraepithelial carcinoma'
# metadata_all$cancer_type[metadata_all$cancer_type == 'undifferentiated'] <- 'undifferentiated endometrial carcinoma'
write.table(exp_anno, file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/status_eclass_fix.tsv', sep = '\t', col.names = TRUE, row.names = FALSE)





# match sequencing indices from gsc file to janine's records
batch1_2 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch1&2_annotation.tsv'
batch3 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch3_annotation.tsv'
batch4 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch4_annotation.tsv'
batch5 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch5_annotation.tsv'
batch6 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch6_annotation.tsv'
batch7 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch7_annotation.tsv'
batch8 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch8_annotation.tsv'
batch9 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch9_annotation.tsv'
batch10 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch10_annotation.tsv'
batch11 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch11_annotation.tsv'
batch12 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch12_annotation.tsv'
batch13 <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/batch13_annotation.tsv'

annos <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/annotation_pool1-12.tsv'
annos <- read_tsv(file = annos)

fastqs <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/fastqs_batches_1to13.tsv'
fastqs <- read.table(file = fastqs, sep = '\t', header = TRUE)
fastqs <- fastqs %>% 
            mutate(variable = rep(c("forward", "backward"), nrow(fastqs) / 2), 
              key = rep(1:(nrow(fastqs) / 2), each = 2)) %>%
            pivot_wider(id_cols = key, names_from = variable, values_from = paths) %>% 
            select(-key)
fastqs <- cbind(as.data.frame(str_split_fixed(fastqs$forward, "/", 9)), fastqs$forward, fastqs$backward)
fastqs[,c(1:6, 8:9)] <- NULL
colnames(fastqs) <- c('library', 'fastq_path_1', 'fastq_path_2')

strReverse <- function(x) { 
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") 
}
add_rev_col <- function(batch_anno) {
  batch_anno <- read.table(file = batch_anno, sep = '\t', header = TRUE)
  batch_anno$FLIPPED_INDEX <- paste0(str_split_fixed(batch_anno$INDEX, "-", 2)[,1], '-', chartr("ATGC","TACG",strReverse(str_split_fixed(batch_anno$INDEX, "-", 2)[,2])))
  batch_anno <- batch_anno %>% select(SOW, LIBRARY, INDEX, FLIPPED_INDEX, FLOWCELL, EXTERNAL_IDENTIFIER)
  colnames(batch_anno) <- c('sow', 'library', 'index', 'flipped_index', 'flowcell', 'external_identifier')
  return(batch_anno)
}
batch1_2 <- add_rev_col(batch1_2)
batch3 <- add_rev_col(batch3)
batch4 <- add_rev_col(batch4)
batch5 <- add_rev_col(batch5)
batch6 <- add_rev_col(batch6)
batch7 <- add_rev_col(batch7)
batch8 <- add_rev_col(batch8)
batch9 <- add_rev_col(batch9)
batch10 <- add_rev_col(batch10)
batch11 <- add_rev_col(batch11)
batch12 <- add_rev_col(batch12)
batch13 <- add_rev_col(batch13)
batches <- rbind(batch1_2, batch3, batch4, batch5, batch6, batch7, batch8, batch9, batch10, batch11, batch12, batch13)

colnames(annos) <- c('sample_id', 'sb_id', 'acc_num', 'study_id', 'cohort', 'hist', 'cancer_type', 'status', 'grade', 'sequence', 'gsc_submission', 'batch')
annos <- left_join(annos, batches, by = c("gsc_submission" = "external_identifier", "sequence" = "flipped_index"))
annos <- left_join(annos, fastqs, by = c("library"))
write.table(annos, file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata_new_2.tsv', sep = '\t', col.names = TRUE, row.names = FALSE)


# Are the cellularities agreeing with the VAFs?
multiple_vafs_to_long <- function(df) {
  multi_cells <- c()
  for (col in c(8:12)) {
    temp <- which(do.call(rbind, lapply(str_split(metD[,col], '/|;'), function(x) length(x))) > 1, arr.ind = TRUE)
    temp[,2] <- col
    multi_cells <- rbind(multi_cells, temp)
  }
  for (i in c(1:dim(multi_cells)[1])) {
    entry <- str_split(metD[multi_cells[1,][[1]], multi_cells[1,][[2]]], '/|;')[[1]]
    line <- metD[multi_cells[i,][[1]],]
    metD[multi_cells[i,][[1]], multi_cells[i,][[2]]] <- str_split(metD[multi_cells[1,][[1]], multi_cells[1,][[2]]], '/|;')[[1]][1]
    metD <- rbind(metD, )
  }
  which(df == "horse", arr.ind = TRUE)
}

p <- ggplot(data=metD, aes(x=dose, y=len, group=supp)) +
  geom_line()+
  geom_point()


######################
### plotting
######################
# Using aheatmap built into NMF package
pdf(file = paste0('~/Desktop/colv_consensus_2tracks.pdf'), 7, 7, onefile = FALSE)
coefmap(mat2, Colv="consensus", tracks=c("basis:", "consensus:"), main="Patient x Signature matrix")
dev.off()

pdf(file = paste0('~/Desktop/colv_consensus_basistrack.pdf'), 7, 7, onefile = FALSE)
coefmap(mat2, Colv="consensus", tracks=c("basis:"), main="Patient x Signature matrix")
dev.off()

pdf(file = paste0('~/Desktop/colv_basis_basistrack_abs.pdf'), 7, 7, onefile = FALSE)
testmat1@fit@W <- mat1@fit@W[,c(2,6,5,4,7,3,1)]
testmat1@fit@H <- mat1@fit@H[c(2,6,5,4,7,3,1),]
coefmap(testmat1, Colv="basis", tracks=c("basis:"), main="Patient x Signature matrix")
dev.off()

pdf(file = paste0('~/Desktop/basismatrix_abs.pdf'), 7, 7, onefile = FALSE)
basismap(testmat1, Rowv = NA, main = "Signature x Component matrix")
dev.off()

pdf(file = paste0('~/Desktop/colv_basis_basistrack.pdf'), 7, 7, onefile = FALSE)
testmat2@fit@W <- mat2@fit@W[,c(2,3,4,5,6,7,1)]
testmat2@fit@H <- mat2@fit@H[c(2,3,4,5,6,7,1),]
coefmap(testmat2, Colv="basis", tracks=c("basis:"), main="Patient x Signature matrix")
dev.off()

pdf(file = paste0('~/Desktop/basismatrix.pdf'), 7, 7, onefile = FALSE)
basismap(testmat2, Rowv = NA, main = "Signature x Component matrix")
dev.off()

pdf(file = paste0('plots/', datasrc, "_plots/", signatures_type, "/CS_heatmap_", cn_caller, "_", window_size, "_", datasrc, ".pdf"), 7, 7, onefile = FALSE)
basismap(component_by_signature, Rowv = NA, main = "Signature x Component matrix")
dev.off()

pdf(file = paste0('~/Desktop/consensusmatrix_britroc_abs.pdf'), 7, 7, onefile = FALSE)
consensusmap(testmat1)
dev.off()

pdf(file = paste0('~/Desktop/consensusmatrix_rbritroc6.pdf'), 7, 7, onefile = FALSE)
consensusmap(component_by_signature_6)
dev.off()

#########
#########
#########
# Plot multisample Heatmaps of CNs across the genome while colour-coding variations.
source("~/Documents/projects/cn_sigs_swgs/scripts/utils.R")

segs <- readRDS("~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalMAD.rds")
segs <- readRDS("~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalVAF.rds")
save_path <- '~/Documents/projects/cn_sigs_swgs/plotting/cnv_heatmaps/'
obj_name <- 'pancan_UCEC_endometrioid_100kb_absolute'
  
# Just p53abn
qual <- data.table::fread(file = '~/Downloads/swgs_sample_biologic_metadata.csv', sep = ',', header = TRUE)
vafcel <- data.table::fread(file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/vafs_and_cellularities.tsv', sep = '\t', header = TRUE)
metadata <- data.table::fread(file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv', sep = '\t', header = TRUE)
metadata <- metadata %>% dplyr::filter(status == 'p53_abn')
metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
p53_segs <- segs[names(segs) %in% metadata$sample_id]

segs_1 <- absolute_called_segments_tcga_pancan_374[names(absolute_called_segments_tcga_pancan_374) %in% serous_samples$bcr_id]

# convert segs to long format table (per-bin CNs)
copy_numbers <- segments_to_copy_number(segs_1, 100000)

# Run 
setDT(copy_numbers)
swgs_cnv_heatmaps(copy_numbers, save_path, obj_name)
plot_cellularity_and_maxvaf(copy_numbers, save_path, obj_name)
plot_BABAM(copy_numbers, save_path, obj_name)

# Retrieve comparable bin CNs & annotations 
obj1 <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_copyNumbersSegmented.rds')
obj2 <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_gl_rCN.rds')
obj3 <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/Xchr_included/30kb_rCN_comCNVfilt.rds')
objs <- list(obj1,obj2,obj3)
bin_area <- 2010000
sample <- 'CC-CHM-1341'

outlist <- compare_bin_CNs(objs, sample, bin_area)


########
compare_mat <- data.frame(snames = colnames(testmat1@fit@H),
                          max_basis_1 = apply(testmat1@fit@H, 2, function(x) which.max(x)),
                          max_basis_2 = apply(testmat2@fit@H, 2, function(x) which.max(x)),
                          correlation = mapply(cor, as.data.frame(testmat1@fit@H), as.data.frame(testmat2@fit@H)))
compare_mat <- compare_mat %>% arrange(max_basis_1)
compare_mat$snames <- factor(compare_mat$snames, levels = compare_mat$snames)
g <- ggplot(compare_mat, aes(snames, c('1'), fill=correlation)) +  
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme(plot.title = element_text(size=20), 
        axis.text.x = element_text(size = 10, angle = 75, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 10)) + 
  labs(x = NULL, y = NULL) +
  ggtitle("Per sample correlations", )
g


############################
############################
### NSMP Metadata Assemblage
############################
############################
# Path_AI metadata from Ali Bashashati
path_ai <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/pathAI_classes_survival.tsv'
path_ai <- read_tsv(file = path_ai)
met <- '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv'
met <- read.table(file = met, sep = '\t', header = TRUE)

todel <- seq(1, nrow(path_ai), 2)
path_ai <- path_ai[-todel,]
met <- met %>% dplyr::filter(status == 'NSMP') %>% dplyr::select(!c(quality,sequence:fastq_path_2))

met <- left_join(met, path_ai, by = c('sample_id' = 'study_id'))
