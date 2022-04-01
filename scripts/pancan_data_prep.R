suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(purrr)
  library(pheatmap)
  library(viridis)
  library(cowplot)
})

##### Data Prep from https://www.synapse.org/#!Synapse:syn1703335
samples <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/TCGA_pancan/pancan12.sample_info.txt', sep = '\t')
samples <- samples %>% dplyr::filter(disease == 'UCEC')
test <- as.data.frame(str_split(samples$tcga_id, pattern = '-',simplify = TRUE))
samples$bcr_id <- paste0(test[,1], '-', test[,2], '-', test[,3])

# pancan clinical indicators
clinical <- data.table::fread(file = '~/Downloads/tcga_UCEC_clinical.patient', sep = '\t')
clinical <- clinical %>% dplyr::select(bcr_patient_barcode, tumor_tissue_site, 
                                       tumor_grade, weight, height, age_at_initial_pathologic_diagnosis, 
                                       histological_type) %>% 
                        dplyr::filter(tumor_tissue_site == "ENDOMETRIAL")
samples <- samples %>% dplyr::left_join(clinical, by = c('bcr_id' = 'bcr_patient_barcode'))

copy_numbers <- data.table::fread(file = '~/Documents/projects/cn_signatures_shallowWGS/TCGA_pancan/data/pancan12_absolute.segtab.txt', sep = '\t')
test <- as.data.frame(str_split(copy_numbers$Sample, pattern = '-',simplify = TRUE))
copy_numbers$bcr_id <- paste0(test[,1], '-', test[,2], '-', test[,3])

##### TCGA specific data
tcga_endo_samples <- data.table::fread(file = '~/Documents/projects/cn_sigs_swgs/tcga_nature_paper_endometrial_data/datafile.S1.1.KeyClinicalData.csv', sep = ',')
samplelist <- samples %>% dplyr::filter(segment_count < 1000, abs_call %in% c('called', 'non-aneuploid', 'low purity', 'high non-clonal')) %>% 
                          dplyr::filter(bcr_id %in% copy_numbers$bcr_id)
segs_out <- vector(mode = "list")
for (i in samplelist$bcr_id) {
  s <- copy_numbers %>% dplyr::filter(bcr_id == i) %>% 
                        dplyr::mutate(segVal = Modal_HSCN_1 + Modal_HSCN_2) %>% 
                        dplyr::select(Chromosome, Start, End, segVal) %>% 
                        dplyr::rename(chromosome = Chromosome, start = Start, end = End)
  segs_out[[i]] <- as.data.frame(s)
}
saveRDS(segs_out, file = "~/Documents/projects/cn_signatures_shallowWGS/TCGA_pancan/data/absolute_called_segments_tcga_pancan_374.rds")


##### Pancer MAFs
# file.maf = system.file('../../Downloads/UCEC_all_248_cases.whitelist.maf', package = 'maftools')
maf = maftools::read.maf(maf = '~/Downloads/UCEC_all_248_cases.whitelist.maf')
maftab <- data.table::fread(file = '~/Downloads/UCEC_all_248_cases.whitelist.maf', sep = '\t')
tp53 <- maftab %>% dplyr::filter(Hugo_Symbol == 'TP53')


# Create clustered signatures heatmap
p1 <- plot_signature_exposures(signatures, save_file = FALSE, order = FALSE, transpose = TRUE)
p2 <- swgs_cnv_heatmaps(copy_numbers, save_path, obj_name)
p3 <- # vertical columns for metadata

##### Run CN-Signatures
copy_numbers_input <- readRDS(file = '~/Documents/projects/cn_signatures_shallowWGS/TCGA_pancan/data/absolute_called_segments_tcga_pancan_374.rds')
metadata <- samplelist
cn_caller <- 'absolute'
window_size <- ''
signatures_type <- 'britrocSignatures'
datasrc <- 'pancan_UCEC'
nsigs <- 7
all_components <- 'britroc_absCNs'                                              # britroc_absCNs, 
component_by_signature <- 'absbritroc'                                          # absbritroc, 
output <- gen_signatures_for_custom_and_builtin(copy_numbers_input, metadata, 
                                                cn_caller, window_size, signatures_type, 
                                                nsigs, datasrc, all_components, 
                                                component_by_signature)



