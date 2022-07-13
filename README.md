# cn-sigs-utils

### File list:
* __main_functions.R__ - Central logic functions used in the calling of copy-number signatures. Needs to be re-named.
* __helper_functions.R__ - Functions assisting in the calling of copy-number signatures.  Needs to be re-named.
* __utils.R__ - General purpose utility functions for anything related to the package.
* __plotting.R__ - General purpose functions for anything related to plotting data analyzed using the package.
* __rascal_batch_and_plots.R__ - Script that needs to be 'pipeline-i-fied' (2nd half). A plotting function needs to be extracted. The workflow described needs to be copied into a vignette and then described.
* __coverage_plots.R__ - First half of this script should probably be converted to a proper function and added to the 'plotting.R' file. Second half relates to targeted panel seq.
* __compare_signature_by_component.R__ - File contains a somewhat tricky function that tries to make a plot comparing component-by-signature matrices (described elsewhere). Likely worth fixing-up and adding to the 'plotting.R' file.
* __gainLoss_qdnaseq.R__ - First half of file creates summary plots for rCN data out of QDNAseq. Probably should be extracted and placed in a vignette. 2nd half of file is a modified version of a QDNAseq function (adds the ability to filter CN change calls that fall below a given threshold).
* __shallowHRD_hg19_1.13_QDNAseq_chrX.R__ - Script from the shallowHRD repo that was more useful to just copy over to this directory. (https://github.com/aeeckhou/shallowHRD). Might be worth just contacting the authors to see if they would be willing to allow us to integrate it into the package.
* __britrocSampleProcessing.R__ - Some scripting very early on in the sWGS CN-Sigs project, might not be terribly useful atm.
* __filter_segmented_CN_calls.R__ - The end needs to be converted into a vignette (demonstrate use). THe functions should already have been added to the 'utils.R' file.
* __getHumanReadableCNprofile.R__ - Functions to generate output tables useful to wetlab researchers.
* __makeFilteredCNSummaryTable.R__ - Functions to generate output tables useful to wetlab researchers.
* __make_qualityMetricsFile.R__ - Random scripting not ready for package integration. Just being kept here for now.
* __get_acns_from_vafs_rascal.R__ - Somewhat redundant set of functions - they should already exist in the 'utils.R' file.



_______________________________________________________________________________________________________________________
README for the CN-Signatures portion of the package

This repository contains all code and documentation necessary to reproduce the analysis undertaken in the 
manuscript [Copy-number signatures and mutational processes in ovarian carcinoma](https://www.biorxiv.org/content/early/2017/09/04/174201)


### Source code
The source code for copy-number signature analysis is split across two files:   
* main_functions.R contains a set of wrapper functions which allow identification and quantification of copy-number signatures.   
* helper_functions.R contains the underlying functions that carry out the copy-number signature analysis.  

These functions require the following packages to be installed:  
[NMF](https://cran.r-project.org/web/packages/NMF/index.html)  
[flexmix](https://cran.r-project.org/web/packages/flexmix/index.html)  
[QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html)  
[YAPSA](https://bioconductor.org/packages/devel/bioc/html/YAPSA.html)

A brief description of the wrapper functions:  
* __extractCopynumberFeatures__ This function takes as input a collection of absolute copy-number profiles and returns a list of copy-number features extracted from these samples. Copy-number profiles can be input as either a QDNAseq object, or as a list of segment tables. The segment tables (one for each sample) should have the following column headers: "chromosome", "start", "end", "segVal".  
* __fitMixtureModels__ This function takes as input a list of copy nnumber features and outputs a list of mixture model components fit to these features.  
* __generateSampleByComponentMatrix__ Given a set of extracted copy number features profiles this function returns a sum-of-posterior sample-by-component matrix. If the all_components argument is specified, then the sum-of-posteriors is calculated using these components, otherwise the component definitions from the manuscript are used.   
* __chooseNumberSignatures__ This function provides a wrapper to nmf functions required for determining the number of copy-number signatures in a dataset. It takes as input sample-by-component matrix and outputs a diagnostic plot for determining the optimal number of copy-number signatures.  
* __generateSignatures__ This function takes a sample-by-component matrix and a specified number of signatures and provides a wrapper to the nmf function for identifying copy-number signatures which returns a NMFfit object.  
* __quantifySignatures__ Given a sample-by-component matrix this function quantifies signature exposures using the LCD function from the YAPSA package, returning a normalised signature-by-sample matrix. If the component_by_signature matrix is specified then this matrix is used to define the signatures otherwise the signature definitions from the manuscript are used.

### Data files
The following data files, located in data/ are required for copy-number signature analysis:  
* __component_parameters.rds__ A RDS object containing the mixture model compenent definitions from the manuscript.  
* __feat_sig_mat.rds__ A RDS object containing the signature definition matrix from the manuscript.  
* __gap_hg19.txt__ An [annotation](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz) file containing the centromere locations for human hg19 genome build.  
* __hg19.chrom.sizes.txt__ An [annotation](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes) file for the hg19 build chromosome lengths.

### How to quantify existing signatures in your own samples
Once you have formatted your copy number data into a list of segment tables, run the following functions:

* CN_features<-extractCopynumberFeatures(<segment tables\>)
* sample_by_component<-generateSampleByComponentMatrix(CN_features)
* quantifySignatures(sample_by_component)

### How to generate your own signatures (advanced user)
Once you have formatted your copy number data into a list of segment tables, run the following functions:

* CN_features<-extractCopynumberFeatures(<segment tables\>)
* all_components<-fitMixtureModels(CN_features) (this step is optional)
* sample_by_component<-generateSampleByComponentMatrix(CN_features,all_components)
* chooseNumberSignatures(sample_by_component)
* component_by_signature<-generateSignatures(sample_by_component,<num signatures\>)
* quantifySignatures(sample_by_component,component_by_signature)

### Manuscript analysis Rmarkdown
The manuscript_Rmarkdown/ directory contains a collection of markdown documents which when compiled, reproduce the entire analysis carried out in the manuscript. To compile this document run knitr to html on the CN_signature_main.R file.

Packages required to compile the Rmarkdown document include:
knitr
RColorBrewer
ggplot2
plyr
dplyr
QDNAseq
Biobase
survival
survcomp
flexmix
DOSE
ppcor
VariantAnnotation
org.Hs.eg.db
TxDb.Hsapiens.UCSC.hg19.knownGene
rtracklayer
ReactomePA
Hmisc
corrplot
ggpubr
NbClust
kableExtra

