# survival.R
#
# Survival analysis of sWGS CN-Sigs samples grouped by predicted signature groups.
# 

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
  library(purrr)
  library(readr)
  library(survival)
  library(ranger)
  library(ggplot2)
  library(ggfortify)
})


###########################
### Functions
###########################

survival_data_age <- fread(file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/toIntegrateTransform/survival_data.tsv', sep = '\t')
exposures <- fread(file = '~/repositories/cnsignatures/data/batch1-13_output/signature_exposures_QDNA_30kb_autosomesOnly_p53_abn_VAFabsbritrocSignatures.csv', sep = ',')
exposures <- fread(file = '~/repositories/cnsignatures/data/batch1-13_output/signature_exposures_QDNA_30kb_autosomesOnly_p53_abn_absbritrocSignatures.csv', sep = ',')
metadata <- data.table::fread(file = '~/Documents/projects/cn_signatures_shallowWGS/metadata/metadata.tsv', sep = '\t', header = TRUE)
exposures[,1] <- NULL

metadata$sample_id <- str_replace_all(metadata$sample_id, "-", ".")
survival_data$study_id <- str_replace_all(survival_data$study_id, "-", ".")
survival_data$study_id[185] <- 'CC.VGH.0041'
survival_data$study_id[186] <- 'CC.VGH.0036'
survival_data$study_id[197] <- 'CC.VGH.0033a'

metadata <- metadata %>% dplyr::select(sample_id, batch, status, hist, cancer_type, grade)
survival_data_1 <- survival_data_1 %>% dplyr::left_join(metadata, by = c('study_id' = 'sample_id'))
survival_data_1$agebins <- '80+'
survival_data_1$agebins[(survival_data_1$age_surg >= 70) & (survival_data_1$age_surg < 80)] <- '70-80'
survival_data_1$agebins[(survival_data_1$age_surg >= 60) & (survival_data_1$age_surg < 70)] <- '60-70'
survival_data_1$agebins[survival_data_1$age_surg < 60] <- 'under 60'
survival_data_1$binaryAge <- ifelse(survival_data_1$age_surg >= 65, "65+", "under 65")
survival_data_1$os_bin <- ifelse(survival_data_1$os_sts == "os.event", 1, 0)
survival_data_1$dss_bin <- ifelse(survival_data_1$dss_sts == "os.event", 1, 0)
survival_data_1 <- survival_data_1 %>% dplyr::filter(status.y == 'p53_abn')

survival_data_1 <- survival_data_1 %>% dplyr::filter(survival_data_1$study_id %in% colnames(exposures))
# exposures <- exposures %>% dplyr::select( c(survival_data$study_id) )
sig_groups <- data.frame(study_id = colnames(exposures), max_sigs = as.character(apply(exposures, 2, function(x) which.max(x))), )
survival_data_1 <- survival_data_1 %>% dplyr::left_join(sig_groups, by = c('study_id'))

# Removing the signature classes that don't have enough observations for the km plot to make sense
survival_data_in <- survival_data_1 %>% dplyr::filter(max_sigs %in% c('1','3','4','5','7'))

km_trt_fit <- survfit(Surv(os_yrs, os_bin) ~ max_sigs, data=survival_data_in)
p <- autoplot(km_trt_fit, conf.int = FALSE, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients grouped by max. signature exposure") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))
p
p <- autoplot(km_trt_fit, conf.int = TRUE, xlab = 'time (years)') +
  labs(title = "Kaplan Meier survival curves for patients grouped by age") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12))
p
ggsave(filename = '~/Documents/projects/cn_signatures_shallowWGS/plotting/survival_plots/p53_abn/kapMeier_QDNA_30kb_autosomesOnly_p53_abn_MADabsbritrocSignature_exposures_3.png', 
       plot = p, width = 8, height = 8, units = 'in')


# Cox Proportional hazards
cox <- coxph(Surv(os_yrs, os_bin) ~ max_sigs + hist + binaryAge + grade, data = survival_data_in)
summary(cox)
cox_fit <- survfit(cox)
#plot(cox_fit, main = "cph model", xlab="Days")
autoplot(cox_fit)






metadata <- metadata %>% dplyr::filter(str_sub(sample_id,-1,-1) != 'N')
metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-2,-1) == '-T']
metadata$sample_id[str_sub(metadata$sample_id,-2,-1) == '-T'] <- str_sub(metadata$sample_id,1,8)[str_sub(metadata$sample_id,-2,-1) == '-T']
metadata$sample_id[str_sub(metadata$sample_id,-1,-1) == 'T'] <- str_sub(metadata$sample_id,1,11)[str_sub(metadata$sample_id,-1,-1) == 'T']
metadata$sample_id[25] <- 'VOA2050'

long_data <- gather(exposures)
long_data$max_sig <- rep(apply(exposures, 2, function(x) which.max(x)), times = 1, each = nsigs)
long_data$sigs <- rep(1:nsigs,dim(exposures)[2])
colnames(long_data) <- c('X', 'Z', 'max_sig', 'Y')
long_data <- long_data %>% arrange(max_sig)
long_data$X <- factor(long_data$X, levels = unique(long_data$X))

