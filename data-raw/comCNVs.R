# Grab common structural variants data from http://dgv.tcag.ca/dgv/app/downloads?ref= using curl or wget or the like

hg19_comCNV <- data.table::fread(file = "~/Desktop/GRCh37_hg19_variants_2020-02-25.txt",
                                 sep = "\t", header = TRUE)
hg19.comCNVs.P01S10KB <- hg19_comCNV %>%
  dplyr::select(chr, start, end, samplesize, observedlosses, observedgains) %>%
  dplyr::mutate(size = end-start) %>%
  dplyr::filter(((observedlosses + observedgains)/samplesize > 0.01) &
                  (observedgains+observedlosses > 9) & (size > 10000))
usethis::use_data(hg19.comCNVs.P01S10KB, overwrite = TRUE)

hg38_comCNV <- data.table::fread(file = "~/Desktop/GRCh38_hg38_variants_2020-02-25.txt",
                                 sep = "\t")
hg38.comCNVs.P01S10KB <- hg38_comCNV %>%
  dplyr::select(chr, start, end, samplesize, observedlosses, observedgains) %>%
  dplyr::mutate(size = end-start) %>%
  dplyr::filter(((observedlosses + observedgains)/samplesize > 0.01) &
                  (observedgains+observedlosses > 9) & (size > 10000))
usethis::use_data(hg38.comCNVs.P01S10KB, overwrite = TRUE)
