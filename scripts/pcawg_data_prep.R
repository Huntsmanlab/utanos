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

# pcawg data transformations/injestion/prep
# segement tables in folder (https://dcc.icgc.org/releases/PCAWG/clinical_and_histology):
# consensus.20170119.somatic.cna.annotated
#
# uuid named segment files need to be cross-referenced with the file (https://dcc.icgc.org/releases/PCAWG/consensus_cnv):
# pcawg_specimen_histology_August2016_v9.xlsx
# tumour_subtype_consolidation_map.xlsx perhaps for restricting to just Uterine