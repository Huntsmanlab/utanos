suppressPackageStartupMessages({
	library(stringr)
	library(dplyr)
	library(ggplot2)
	library(data.table)
	library(purrr)
	library(pheatmap)
	library(viridis)
	library(cowplot)
  library(NMF)
  library(ggpubr)
  library(maftools)
  library(vcfR)
})

setwd("C:/Users/admin/Desktop/CNsigs")
source('scripts/main_functions.R')
source('scripts/plotting.R')
source('scripts/utils.R')

# pcawg data transformations/injestion/prep
# segment tables in folder (https://dcc.icgc.org/releases/PCAWG/clinical_and_histology):
# consensus.20170119.somatic.cna.annotated
#
# uuid named segment files need to be cross-referenced with the file (https://dcc.icgc.org/releases/PCAWG/consensus_cnv):
# pcawg_specimen_histology_August2016_v9.xlsx
# tumour_subtype_consolidation_map.xlsx perhaps for restricting to just Uterine


cnDataDir = "~/Documents/projects/cn_sigs_swgs/PCAWG_endometrial_corpus/data/consensus.20170119.somatic.cna.annotated/"
metadataPath = "~/Documents/projects/cn_sigs_swgs/PCAWG_endometrial_corpus/metadata/pcawg_specimen_histology_August2016_v9.csv"
metadataPath2 = "~/Documents/projects/cn_sigs_swgs/PCAWG_endometrial_corpus/metadata/tcga_UCEC_clinical.patient.tsv"
metadataPath3 = "data/pancan12.sample_info.txt"
metadataPath4 = "data/icgc_sample_annotations_summary_table.txt"

tp53path1 = "tp53Data/consensus_snv_indel.maf"
tp53path2 = "tp53Data/driver_mutations.tsv"
savePath = "~/Documents/projects/cn_sigs_swgs/PCAWG_endometrial_corpus/processed/PCAWG_tables.rds"

pcawgMetadata = as.data.frame(data.table::fread(file = metadataPath, sep = ',', header = TRUE))
pcawgMetadata2 = as.data.frame(data.table::fread(file = metadataPath2, sep = '\t', header = TRUE))
pcawgMetadata3 = as.data.frame(data.table::fread(file = metadataPath3, sep = '\t', header = TRUE))
pcawgMetadata4 = as.data.frame(data.table::fread(file = metadataPath4, sep = '\t', header = TRUE))

tp53data1 = read.maf(tp53path1)
tp53data2 = as.data.frame(data.table::fread(file = tp53path2, sep = '\t', header = TRUE))



endoUUID = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA",]$tcga_sample_uuid
pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA",]$histology_tier4
subID = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA",]$submitted_sample_id
iter = 1
for(id in subID) {
  subs = strsplit(id, "-")[[1]][1:3]
  subID[iter] = paste(subs[1], subs[2], subs[3], sep = "-")
  iter = iter + 1
}

x = unique(pcawgMetadata$tcga_specimen_uuid)
y = pcawgMetadata2$bcr_aliquot_uuid
pcawgMetadata2[y %in% x,]$tumor_grade
pcawgMetadata[x %in% y,]$histology_abbreviation
pcawgMetadata2[pcawgMetadata2$bcr_aliquot_uuid %in% subID]
pcawgMetadata2[pcawgMetadata2$bcr_patient_uuid %in% subID]$tumour_grade
pcawgMetadata3[pcawgMetadata3$tcga_id %in% pcawgMetadata$tcga_sample_uuid]
pcawgMetadata3[pcawgMetadata3$tcga_id %in% pcawgMetadata$submitted_sample_uuid]
pcawgMetadata3[pcawgMetadata3$tcga_id %in% pcawgMetadata$submitted_specimen_uuid]


subID = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA",]["# icgc_specimen_id"]
pcawgMetadata4[pcawgMetadata4$icgc_sample_id %in% subID,]

# Apr 20
tp53DataFiles = list.files("data/consensus/indel")
tp53DataSamples = vector()
for(i in 1:length(tp53DataFiles)) {
  tp53DataSamples[i] = strsplit(tp53DataFiles[i], "[.]")[[1]][1]
}

# Apr 21
#mutVCF = read.vcfR("data/new/muse.vcf")
mut1 = as.data.frame(data.table::fread(file = "data/new/joint_fpkm.tsv", sep = '\t', header = TRUE))
test = mut1[mut1$feature %like% "ENSG00000141510",]
mutSamples = colnames(mut1)
names(tables) %in% mutSamples

mut2 = as.data.frame(data.table::fread(file = "data/new/mutOccurrences.tsv", sep = '\t', header = TRUE))
endoSA = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA",]
allSA = pcawgMetadata$tcga_sample_uuid
mut2SA = mut2$`Donor ID`

tp53UUID = pcawgMetadata[pcawgMetadata$icgc_donor_id %in% mut2SA ,]$tcga_sample_uuid
tp53NegUUID = endoSA[!endoSA$icgc_donor_id %in% mut2SA ,]$tcga_sample_uuid

for(i in tables) {
  print(length(i[["segVal"]]))
}

sig5samples = list(tablesMut[["00db4dc2-3ec7-4ff9-9233-d69c8c8a607f"]][["segVal"]], tablesMut[["b38d0777-4901-48b8-9cdc-33b7f13a424f"]][["segVal"]], tablesMut[["8dd14f0e-8601-4aa1-864c-3c49e768cdd1"]][["segVal"]], tablesMut[["d12cfd8b-682d-41df-acf8-ee7f68a6241c"]][["segVal"]])

sample_id = endoSA$tcga_sample_uuid
histotype = endoSA$histology_tier4
PCAWGData = as.data.frame(cbind(sample_id, histotype))
PCAWGData[, "tp53_mutant"] <- NA

for(i in 1:nrow(PCAWGData)) {
  if(PCAWGData[i,]$sample_id %in% tp53UUID) {
    PCAWGData[i,]$tp53_mutant = TRUE
  } else {
    PCAWGData[i,]$tp53_mutant = FALSE
  }
}

saveRDS(PCAWGData, "data/PCAWG_metadata.rds")

serousUUID = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA" & pcawgMetadata$histology_tier4 == "Serous cystadenocarcinoma",]$tcga_sample_uuid
endometrioidUUID = pcawgMetadata[pcawgMetadata$histology_abbreviation == "Uterus-AdenoCA" & pcawgMetadata$histology_tier4 == "Adenocarcinoma, endometrioid",]$tcga_sample_uuid

pcawgMetadata2[pcawgMetadata2$bcr_patient_uuid %in% subID]

getSegTablePCAWG = function(UUIDs, cnDataPath = "data.tcga") {
  files = list.files(cnDataPath)
  iter = 1
  tables = list()
  fnames = c()
  for(file in files) {
  	fname = strsplit(file, "[.]")[[1]][1]
  	if(fname %in% UUIDs) {
  		fnames[iter] = fname
    	tables[[iter]] = as.data.frame(data.table::fread(file = paste(cnDataPath, file, sep = "/"), sep = '\t', header = TRUE))
    	tables[[iter]] = tables[[iter]][, c("chromosome", "start", "end", "total_cn")]
    	colnames(tables[[iter]]) = c("chromosome", "start", "end", "segVal")
    	tables[[iter]] = tables[[iter]][!is.na(tables[[iter]]$segVal),]
    	iter = iter + 1
  	}
  }
  names(tables) = fnames
  return(tables)
  #saveRDS(tables, savePath)
}

tables = getSegTablePCAWG(endoUUID, cnDataDir)
tablesMut = getSegTablePCAWG(tp53UUID, cnDataDir)
endoTables = getSegTablePCAWG(endometrioidUUID, cnDataDir)
serTables = getSegTablePCAWG(serousUUID, cnDataDir)
test = extractCopynumberFeatures(tables)
byComponent = generateSampleByComponentMatrix(test, "britroc_absCNs")
sigs = as.data.frame(quantifySignatures(byComponent, component_by_signature = "absbritroc"))
saveRDS(tables, "data/PCAWGSegTable.rds")

g = plot_signature_exposures(sigs)
g
