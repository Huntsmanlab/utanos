suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(purrr)
  library(NMF)
  library(flexmix)
  library(YAPSA)
  library(pheatmap)
  library(viridis)
  library(hrbrthemes)
})
source(file = 'scripts/main_functions.R')
source(file = 'scripts/helper_functions.R')


# ichorCNA
# sample_list <- list()
# j <- 1
# for (i in metadata$sample_id) {
#   print(i)
#   sample <- data.table::fread(file = paste0('/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/se_ichorCNA_outs_lowNormFrac_noSubCloneFracMax_agctNormPanel/', i, '.seg'), sep = '\t', header = TRUE)
#   sample <- as.data.frame(sample)
#   sample <- sample[ , c('chr', 'start', 'end', 'copy.number')]
#   colnames(sample) <- c('chromosome', 'start', 'end', 'segVal')
#   sample <- rbind(c(1,1,1000000,2), sample)
#   sample[is.na(sample)] <- 2
#   sample_list[[i]] <- sample
#   j = j+1
# }
# sample_list <- sample_list[7:22]
# CN_features <- extractCopynumberFeatures(sample_list)


# Function that creates plots
make_diagnostics_cnsignatures_plots <- function(signatures, metadata, component_by_signature, cn_caller, window_size, signatures_type, datasrc) {
  
  # Order samples by type
  signatures <- as.data.frame(signatures)
  snames <- str_split(colnames(signatures), '.se')
  snames <- map(snames, as.data.table)
  snames <- as.data.frame(t(as.data.frame(snames)))
  colnames(signatures) <- snames$V1
  rownames(signatures) = paste0("signature ", rownames(signatures))
  metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
  plotting_metadata <- metadata %>% select(sample_id, batch, tissue, status, cancer_type)
  annotation <- left_join(snames, plotting_metadata, by = c("V1" = "sample_id"))
  annotation.col <- data.frame(Status = paste0(annotation$status,'_',str_replace_all(annotation$cancer_type, " ", ".")))
  rownames(annotation.col) <- annotation$V1
  
  # Using aheatmap built into NMF package
  pdf(file = paste0('plots/', datasrc, "_plots/", signatures_type, "/PS_heatmap_", cn_caller, "_", window_size, "_", datasrc, ".pdf"), 7, 7, onefile = FALSE)
  coefmap(component_by_signature, Colv="consensus", annCol = annotation.col, tracks=c("basis:"), main="Patient x Signature matrix")
  dev.off()
  
  pdf(file = paste0('plots/', datasrc, "_plots/", signatures_type, "/CS_heatmap_", cn_caller, "_", window_size, "_", datasrc, ".pdf"), 7, 7, onefile = FALSE)
  basismap(component_by_signature, Rowv = NA, main = "Signature x Component matrix")
  dev.off()
  
  return(signatures)
}

# Implement function that does the signature calling for both custom signatures and using the britroc-based ones
# Gen plots automatically too
gen_signatures_for_custom_and_builtin <- function (copy_numbers_input, metadata, cn_caller, window_size, signatures_type, nsigs, 
                                                   datasrc, all_components, component_by_signature = NULL, multi_sols_data = FALSE) {
  
  dir.create(path = paste0("data/", datasrc, "_output"), recursive = TRUE, showWarnings = FALSE)
  dir.create(path = paste0("plots/", datasrc, "_plots/", signatures_type), recursive = TRUE, showWarnings = FALSE)
  # CN_features <- extractRelativeCopynumberFeatures(copy_numbers_input)
  CN_features <- extractCopynumberFeatures(copy_numbers_input)
  # Notes:
  # For the customSignatures option
  # An explicit decision has been made not to re-generate the modelling of mixture components as the britroc sample group 
  # is so much larger and doing so again with just 22 samples will never work as well.
  # It also would not be a fair comparison with the Britroc samples modelling.
  # This option will likely need to be run twice. First to generate the num_sigs plots.
  # The second time to run it with your chosen number of signatures.
  # 
  # compareSignatures option:
  # The only difference with this option is that we do not choose or try to determine the 
  # optimal number of signatures. We just use the same number as was determined from the Britroc group.
  # 
  # for the fitMixtureModels() function
  # The segsize component wasn't structured properly (our listed segment sizes). Now it is. (We use the sample.seg files from ichorCNA)
  # These need to be variable. As in the continuous segments w/out CN change over the genome 
  #
  # This doesn't fix the errors associated with the changepoint modelling as well...
  # When I try modelling using a poisson distribution rather than a norm model it doesn't throw an error...
  # Why should this be a normal distribution?
  # I now have it as a poisson.... hopefully ok?
  #
  # cophenetic: after the trend begins increasing (after 2) choose the max value before it decreases again.
  # dispersion: higher == better, measures the reproducibility
  # silhouette: higher == better, measures consistency between clusters
  
  if (signatures_type == 'customSignatures') {
    all_components <- fitMixtureModels(CN_features)
    sample_by_component <- generateSampleByComponentMatrix(CN_features, all_components, cores = 1)
    component_by_signature <- generateSignatures(sample_by_component, nsigs)
    browser()
    num_sigs <- chooseNumberSignatures(sample_by_component)
    pdf(file = paste0('plots/', datasrc, "_plots/", signatures_type, "/", cn_caller, "_", window_size, "_custSigs_number_of_signatures.pdf"), 7, 7, onefile = FALSE)
    num_sigs
    dev.off()
    browser()
    # component_by_signature <- generateSignatures(sample_by_component, nsigs)
    saveRDS(component_by_signature, file = paste0("data/", datasrc, "_output/", cn_caller, "_", window_size, "_custSigs_component_by_signature.rds"))
    signatures <- quantifySignatures(sample_by_component, basis(component_by_signature))
    signatures <- make_diagnostics_cnsignatures_plots(signatures, metadata, component_by_signature, cn_caller, window_size, signatures_type, datasrc)
    
  } else if (signatures_type == 'compareSignatures') {
    sample_by_component <- generateSampleByComponentMatrix(CN_features, all_components)
    component_by_signature <- generateSignatures(sample_by_component, nsigs)
    saveRDS(component_by_signature, file = paste0("data/", datasrc, "_output/", cn_caller, "_", window_size, "_component_by_signature.rds"))
    signatures <- quantifySignatures(sample_by_component, basis(component_by_signature))
    signatures <- make_diagnostics_cnsignatures_plots(signatures, metadata, component_by_signature, cn_caller, window_size, signatures_type, datasrc)
    
  } else {
    sample_by_component <- generateSampleByComponentMatrix(CN_features, all_components)
    signatures <- quantifySignatures(sample_by_component, component_by_signature)
    signatures <- as.data.frame(signatures)
  }
  browser()
  pdf(file = paste0('plots/', datasrc, "_plots/", signatures_type, "/PC_heatmap_", cn_caller, "_", window_size, "_", datasrc, ".pdf"), 7, 7, onefile = FALSE)
  NMF::aheatmap(sample_by_component,Rowv=NULL, main="Component x Sample matrix")
  dev.off()
  
  return(signatures)
}

# Grab Britroc rCNs
copy_numbers_input <- readRDS(file = file.path("data", "input_copy_numbers", "britroc_relative_copynumber_30kb.rds"))
# Grab Britroc ACNs
copy_numbers_input <- readRDS(file = file.path("data", "input_copy_numbers", "britroc_absolute_copynumber.rds"))
# Grab QDNAseq output
copy_numbers_input <- readRDS(file = file.path("~/Documents", "projects", "cn_signatures_shallowWGS","qdnaseq_copy_number", "batch_1-13", "autosomes_only", "30kb_copyNumbersSegmented.rds"))
# Grab ACE output
copy_numbers_input <- readRDS(file = file.path("..", "ace_cn_scaler", "sWGS_pilot_output","chosenSegmentTablesList.rds"))
# Grab rascal output
copy_numbers_input <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_rascal_CN_Collapsed_segments_optimalMAD.rds')
copy_numbers_input <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1-13/autosomes_only/30kb_comCNVfilt_rascal_CN_Collapsed_segments_optimalVAF.rds')
# Grab multi-solution rascal output 
multi_sols_data <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/qdnaseq_copy_number/batch_1234/Xchr_included/100kb_rascal_multisolutions_CN_Collapsed_segments.rds')
  
# Read in Metadata
metadata <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv', sep = '\t', header = TRUE)
metrics <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/metrics_15kb.tsv', sep = '\t', header = TRUE)
outcomes <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/metadata/nsmp_outcomes.tsv', sep = '\t', header = TRUE)
biologic <- data.table::fread(file = '~/Downloads/swgs_sample_biologic_metadata.csv', sep = ',', header = TRUE)
# sort and filter metadata for just endometrium tissue from 
metadata <- metadata %>% left_join(metrics, by = c('sample_id'))
metadata <- metadata[order(metadata$sample_id), ]
metadata <- metadata %>% filter(rascal_solutions_num > 0, rascal_solutions_num < 7, num_segments < 225, tissue != 'ovary', tissue != 'cervical_stroma')
# metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
# metadata <- metadata[metadata$tissue == 'endometrium']
# metadata <- metadata %>% filter(cancer_type %in% c('endometrioid endometrial carcinoma'))
# metadata <- metadata %>% dplyr::filter(status %in% c('NSMP'))
# metadata <- metadata %>% left_join(outcomes, by = c('sample_id' = 'patient_id'))
metadata <- metadata %>% dplyr::select(sample_id, batch, path_AI_pred, tissue, status, hist, cancer_type, grade)
biologic <- biologic %>% dplyr::select(id, other_aliases, notes, quality_maxwell)
metadata <- metadata %>% left_join(biologic, by = c('sample_id' = 'other_aliases'))

# Restrict copy_number sample inputs to just the samples in the filtered metadata
copy_numbers_input <- copy_numbers_input[,sampleNames(copy_numbers_input) %in% metadata$sample_id]
multi_sols_data <- multi_sols_data[names(multi_sols_data) %in% metadata$sample_id]
multi_sols_data <- FALSE

cn_caller <- 'absolute'
window_size <- '30kb'
signatures_type <- 'vanSignatures'
datasrc <- 'pancan_UCEC'
nsigs <- 6
all_components <- 'vancouver_HQendo_VAFrascal'                                  # britroc_absCNs, vancouver_HQendo_VAFrascal
component_by_signature <- 'VAFrascalHQendoVan6sigs'                             # absbritroc, VAFrascalHQendoVan6sigs
output <- gen_signatures_for_custom_and_builtin(copy_numbers_input, metadata, cn_caller, window_size, signatures_type, nsigs, datasrc, all_components, component_by_signature)
# output <- gen_signatures_for_custom_and_builtin(copy_numbers_input, metadata, cn_caller, window_size, signatures_type, nsigs, datasrc, multi_sols_data)
# output <- gen_signatures_for_custom_and_builtin(sample_list, metadata, cn_caller='ichorCNA', window_size='1Mb', signatures_type='compareSignatures', nsigs, datasrc)



#################
### PLOTTING GGPLOT
#################

# Less optimal => ggplot as it doesn't include clustering and annotation bars 
save_file <- '~/repositories/cnsignatures/plots/pancan_UCEC_plots/britrocSignatures/autosomesOnly_pancanUCEC_7sigs_called_britroc_signatures_absolute_heatmap'

output1 <- output[c(2,3,4,5,6,7,1),]    # taken from hand-annotation of rBritroc sigs match to abs Britroc sigs assignment from paper
output1 <- output[c(2,6,5,4,7,3,1),]    # re-order Absolute Sigs based on Brenton paper order
# output <- output[arrange(metadata, cancer_type)$sample_id]
long_data <- gather(output1)
long_data$max_sig <- rep(apply(output1, 2, function(x) which.max(x)), times = 1, each = nsigs)
long_data$sigs <- rep(1:nsigs,dim(output1)[2])
colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
## snames <- map(snames, as.data.table)
# snames <- as.data.frame(long_data$key)
# plotting_data <- data.frame(X = snames, 
#                             Y = rownames(output), 
#                             Z = gather(output)$value)
long_data <- long_data %>% arrange(max_sig)
long_data$X <- factor(long_data$X, levels = unique(long_data$X))

g <- ggplot(long_data, aes(X, Y, fill=Z)) +  
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme(plot.title = element_text(size = 20),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(size = 15, angle = 75, vjust = 0.5, hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size = 18)) + 
  labs(fill = 'Signature \nExposure', x = "Samples", y = " ") + 
  scale_y_discrete(limits = c('S1', 'S2', 'S3', 'S4', 'S5', 'S6')) +
  ggtitle("Signature exposures called per sample")
g
ggsave(paste0(save_file,".png"), plot = g, width = 15, height = 10)

rownames(output1) <- c(1:7)
write.csv(output1, file = "~/repositories/cnsignatures/data/pancan_UCEC_output/signature_exposures_absolute_autosomesOnly_pancanUCEC_britrocSignatures.csv")

# Grab the sample names from a specific signature
uu <- unique((long_data %>% filter(max_sig == 7))$X)
uu <- as.character(uu)
vv <- metadata[metadata$sample_id %in% uu,]

#################
### PLOTTING PHEATMAP
#################
# Using pheatmap
# Specify colors
# palette <- c('#E5F5E0', '#C7E9C0', '#A1D99B', '#74C476', '#41AB5D', '#238B45', '#006D2C',
#              '#00441B', '#000000', '#67000D', '#A50F15', '#CB181D', '#EF3B2C', '#FB6A4A', '#FC9272', '#FCBBA1', '#FEE0D2')
palette <- viridis_pal()(10)
sample_type_pal <- c("#E6AB02", "#7570B3", "#A1D99B", "#CB181D", "#A6761D")
sampletype <- c(sample_type_pal[1:length(unique(annotation$status))])
names(sampletype) <- as.factor(unique(annotation$status))
ann_colors = list(
  SampleType = sampletype
)

# misc. parameters
clusterR <- FALSE
clusterC <- TRUE
fontsize <- 6
cell_size <- 10
height <- 8
width <- 6
plotTitle <- "Patient x Signature matrix"
plotName <- "pilotdata_plots/customSignatures/sigs_clustered_pheatmap_ACE_pilotdata.pdf"

pdf(file = plotName, width, height, onefile = FALSE)
pheatmap(signatures, cluster_cols = clusterC, 
         annotation_col = annotation.col,
         annotation_colors = ann_colors, color = palette, fontsize_col = fontsize,
         fontsize_row = fontsize, scale = "row", border_color = NA,
         clustering_method = "ward.D2", cluster_rows = clusterR, cellwidth = cell_size,
         cellheight = cell_size*3, breaks = seq(0,1,0.1), main = plotTitle)
dev.off()


temp <- output
temp <- as.data.frame(temp)
temp <- gather(temp, sample, sig_proportions, `VS07-16094-D14`:VOA626C, factor_key = TRUE)
temp <- mutate(temp, signatures = as.factor(rep(1:7,22)))
ggplot(temp, aes(x = sample, y = sig_proportions, fill = signatures)) + 
  geom_bar(stat = 'identity') +
  ggtitle("Signatures called per sample") +
  scale_fill_viridis(discrete=TRUE) +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5, hjust=0.5))



#################################
### Create ACE samples format input
#################################
sample <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/se_ichorCNA_outs_lowNormFrac_noSubCloneFracMax_agctNormPanel/VOA427C.cna.seg', 
                            sep = '\t', header = TRUE)
sample <- as.data.frame(sample)
ace_sample <- sample[,c('chr', 'start', 'end', 'VOA427C.Corrected_Copy_Number')]
colnames(ace_sample) <- c('chromosome', 'start', 'end', 'segmented')
ace_sample <- rbind(c(1,1,1000000,2), ace_sample)
ace_sample[is.na(ace_sample)] <- 2
write.table(ace_sample, file = "~/Documents/projects/cn_signatures_shallowWGS/ace_sample_format_CN_calls/VOA427C_cna_segs.tsv", row.names=FALSE, sep="\t")

sample <- data.table::fread(file = '/Users/maxwell/Documents/projects/cn_signatures_shallowWGS/se_ichorCNA_outs_lowNormFrac_noSubCloneFracMax_agctNormPanel/VOA427C.seg.txt', 
                            sep = '\t', header = TRUE)
ace_sample <- as.data.frame(sample)
ace_sample <- ace_sample[,c('chrom', 'start', 'end', 'Corrected_Copy_Number')]
colnames(ace_sample) <- c('chromosome', 'start', 'end', 'segmented')
ace_sample <- rbind(c(1,1,1000000,2), ace_sample)
ace_sample[is.na(ace_sample)] <- 2
write.table(ace_sample, file = "~/Documents/projects/cn_signatures_shallowWGS/ace_sample_format_CN_calls/VOA427C_segs.tsv", row.names=FALSE, sep="\t")
