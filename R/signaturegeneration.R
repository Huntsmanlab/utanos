# sWGS Utils package - utanos
# User facing functions used in the generation of copy-number signatures

###########################
### Functions
###########################
# CallSignatures()
###


#' Call/Determine Signature Exposures for a set of Samples
#'
#' @description
#'
#' CallSignatures() calculates copy-number signature exposures for a list of samples given their segmented copy-number profiles.
#' This function makes use of code found in the CN-Signatures bitbucket repo.
#' and described more robustly here https://www.nature.com/articles/s41588-018-0179-8.
#'
#'
#' @param copy_numbers_input A list of dataframes. Each dataframe should be the segmented copy-numbers of a sample.
#' @param component_models String. The mixture models to use in determining the components present in this dataset. \cr
#' Options: \cr
#' '30kb_ovarian/component_models_britroc_aCNs.rds' - Components determined by modelling the Britroc Ovarian sWGS samples. (Macintyre 2018) \cr
#' 'vancouver_HQendo_VAFrascal' - Components determined by modelling high-quality p53abn sWGS endometrial samples. (2022) \cr
#' @param signatures String. The signatures themselves. Please name those for which to calculate exposures. \cr
#' Options: \cr
#' '30kb_ovarian/component_by_signature_britroc_aCNs.rds' - The signatures created using the Britroc Ovarian sWGS samples. (Macintyre 2018) \cr
#' 'VAFrascalHQendoVan6sigs' - The signatures created using the high-quality p53abn sWGS endometrial samples. \cr
#' @param data_path String. The path where to find signatures objects and component models.
#' @param plot_savepath (optional) String. The path where to save a sample-by-component heatmap. Please provide a directory.
#' @param sigs_savepath (optional) String. The path where to save a csv of the calculated signature exposures. Please provide a directory.
#' @param relativeCN_data (optional) Logical. If using relative Copy-Number data as input, set to true. Otherwise, ignore.
#' @returns A dataframe of the exposures for each sample to each signature.
#' @details
#' ```
#' todo - make examples

#' ```
#'
#' @export
CallSignatures <- function (copy_numbers_input,
                            component_models = NULL,
                            signatures = NULL,
                            data_path = NULL,
                            plot_savepath = NULL,
                            sigs_savepath = NULL,
                            relativeCN_data = NULL) {

  stopifnot(!is.null(component_models))                                         # We want the user to explicitly declare the components they intend to use
  stopifnot(!is.null(signatures))                                               # We want the user to explicitly declare the signatures they intend to use
  stopifnot(!is.null(data_path))                                                # We want the user to explicitly declare the signatures they intend to use

  dir.create(path = sigs_savepath,
             recursive = TRUE, showWarnings = FALSE)
  dir.create(path = plot_savepath,
             recursive = TRUE, showWarnings = FALSE)

  if (!is.null(relativeCN_data)) {
    CN_features <- ExtractRelativeCopyNumberFeatures(copy_numbers_input)
  } else {
    CN_features <- ExtractCopyNumberFeatures(copy_numbers_input, 'hg19')
  }

  component_models <- paste0(data_path, '/', component_models)
  sample_by_component <- GenerateSampleByComponentMatrix(CN_features, component_models)
  signatures <- paste0(data_path, '/', signatures)
  sigex <- QuantifySignatures(sample_by_component, signatures)
  sigex <- as.data.frame(sigex)

  if (!is.null(plot_savepath)) {
    datetoday <- Sys.Date()
    pdf(file = paste0(plot_savepath, "/componentbysample_heatmap_", component_models,
                      "_", signatures, "_", datetoday, ".pdf"),
        7, 7, onefile = FALSE)
    NMF::aheatmap(sample_by_component,Rowv=NULL, main="Component x Sample matrix")
    dev.off()
  }

  if (!is.null(sigs_savepath)) {
    rownames(sigex) <- c(1:dim(sigex)[1])
    write.csv(sigex, file = paste0(sigs_savepath, "/signature_exposures_",
                                        component_models, "_", signatures, "_",
                                        datetoday, ".csv"))
  }

  return(sigex)
}

result <- CallSignatures(copy_numbers_input = CalculateACNs_output,
                         component_models = ComponentModelsBritrocACNs,
                         signatures = ComponentBySignatureBritrocACNs,
                         data_path = path,
                         plot_savepath = "~/projects/signature_and_component_models",
                         sigs_savepath = "~/projects/signature_and_component_models")


#CalculateACNs_output <- readRDS(file = "~/projects/replication_files/final_output_aCN.rds")
#ComponentModelsBritrocACNs <- "component_models_britroc_aCNs.rds"
#ComponentBySignatureBritrocACNs <- "component_by_signature_britroc_aCNs.rds"
#path <- "~/projects/signature_and_component_models"

