## code to prepare `gaps.genome` datasets goes here

# These tables began by downloading the centromere and telomere data from the UCSC table browser.
# Web location: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2244581546_h4zkv9l9WA00Td2ONQ0EgsAkR6C8
# The process differed slightly for each genome. For hg38 they began apart.
# These tables were then read into R. They were very briefly processed as outlined below.
#

# To prep hg19 table
gaps.hg19 <- gaps.hg19 %>% dplyr::select(V2, V3, V4, V7, V8)
colnames(gaps.hg19) <- c("chromosome", "start", "end", "size", "type")
gaps.hg19 <- gaps.hg19 %>% dplyr::filter(type %in% c("centromere", "telomere"))
# usethis::use_data(gaps.hg19, overwrite = TRUE)

# To prep hg38 gaps table
centromeres.hg38 <- centromeres.hg38 %>% dplyr::filter(type == "centromere") %>%
  group_by(chrom) %>% dplyr::summarise(start = min(chromStart),
                                       end = max(chromEnd),
                                       type = dplyr::first(type)) %>%
  dplyr::mutate(size = end - start)
telomeres.hg38 <- telomeres.hg38 %>% dplyr::select(chrom, chromStart,
                                              chromEnd, type, size) %>%
  dplyr::rename(start = chromStart,
                end = chromEnd)
gaps.hg38 <- rbind(centromeres.hg38, telomeres.hg38)
gaps.hg38 <- gaps.hg38[,c(1,2,3,5,4)]
colnames(gaps.hg38) <- c("chromosome", "start", "end", "size", "type")
# usethis::use_data(gaps.hg38, overwrite = TRUE)
