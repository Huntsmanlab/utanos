Create New Copy-Number Signatures
================

This vignette serves to demonstrate how to use copy-number (CN)
abnormalities, which could represent the imprints of distinct mutational
processes, and create signatures from these CN-features. In order to
generate new copy-number signatures, a set of samples are needed for
which you have segmented copy-number profiles. These profiles should be
formatted as a list of named dataframes. One for each sample.

### I - The input data…

``` r
library(utanos)
library(QDNAseq)
library(data.table)
library(dplyr)

# Reading in test data
segmented_aCNs <- readRDS("~/projects/cn-signatures-project/segmented_aCNs.rds")
# QDNAseq object
dplyr::glimpse(head(segmented_aCNs, 3))

List of 3
 $ CC.CHM.1341:'data.frame':    210 obs. of  4 variables:
  ..$ chromosome: chr [1:210] "1" "1" "1" "1" ...
  ..$ start     : int [1:210] 2850001 32100001 44670001 116790001 117960001 118740001 119010001 145530001 171090001 171660001 ...
  ..$ end       : int [1:210] 32100000 44490000 116790000 117960000 118740000 119010000 120270000 171090000 171660000 172530000 ...
  ..$ segVal    : num [1:210] 1.12 4.68 1.41 5.12 4.15 ...
 $ CC.CHM.1347:'data.frame':    683 obs. of  4 variables:
  ..$ chromosome: chr [1:683] "1" "1" "1" "1" ...
  ..$ start     : int [1:683] 2850001 3450001 4020001 5970001 15630001 15960001 18030001 19200001 22020001 23340001 ...
  ..$ end       : int [1:683] 3120000 3810000 5970000 15630000 15960000 18030000 19200000 21840000 23340000 23910000 ...
  ..$ segVal    : num [1:683] 3.01 5.19 3.27 5.48 1.39 ...
 $ CC.CHM.1355:'data.frame':    193 obs. of  4 variables:
  ..$ chromosome: chr [1:193] "1" "1" "1" "10" ...
  ..$ start     : int [1:193] 2850001 145530001 201600001 150001 30090001 31500001 33690001 34230001 42960001 106230001 ...
  ..$ end       : int [1:193] 120270000 201600000 249210000 30090000 31500000 33690000 34230000 38430000 106230000 111540000 ...
  ..$ segVal    : num [1:193] 1.94 2.01 3.04 1.98 3.04 ...
```

### II - Do some modelling…

To begin…  
1. First, extract CN-features using the ExtractCopynumberFeatures()
function.  
2. Fit mixture components.  
3. Calculate a patient-by-component sum-of-posterior probabilities
matrix.  
4. Search for the optimal number of signatures (interval of 3-12 is used
by default), running NMF 1000 times with different random seeds for each
signature number.

``` r
CN_features <- ExtractCopynumberFeatures(copy_numbers_input)
component_models <- FitMixtureModels(CN_features)
sample_by_component <- GenerateSampleByComponentMatrix(CN_features, component_models)
ChooseNumberSignatures(sample_by_component)
```

\#3 Expanded: For each copy-number event in each sample, compute the
posterior probability of belonging to a component. For each sample, sum
these posterior event vectors to a sum-of-posterior probabilities
vector. Combine sum-of-posterior vectors into a patient-by-component
sum-of-posterior probabilities matrix.

The output from the ChooseNumberSignatures() function is a figure
plotting several metrics that aid in the selection of an optimal number
of signatures. These metrics include cophenetic, dispersion, silhouette,
and sparseness.

What to look for in the values of each metric:  
cophenetic: after the trend begins increasing (after 2) choose the max
value before it decreases again.  
dispersion: higher == better, measures the reproducibility  
silhouette: higher == better, measures consistency between clusters  
sparseness: higher == better, try to choose the maximum sparsity which
can be achieved without exceeding that which was observed in the
randomly permuted matrices

Here is an example of said figure.  
<img src="../images/choose_optimal_Nsignatures_metrics.png" width="943" />

### III - Stop and evaluate

At this point it is a good idea to stop and reflect about the quality of
the data being used to create your new signatures. If the metrics
provided above for the input matrices do not appear to be outperforming
the randomized matrices, you may wish to set more restrictive demands on
what samples can be used as input data. It also could be helpful to
examine the sample-by-component matrix as a heatmap to see if there are
major contributions from many different samples to each component.  
If a clear Factorization rank (number of signatures) emerges, continue
on…

### IV - Create the signatures

``` r
signatures <- GenerateSignatures(sample_by_component,<num signatures>)
saveRDS(signatures, file = paste0("path_to_signatures/signatures_object.rds"))

# Also save the component models created earlier for re-use
saveRDS(component_models, file = paste0("path_to_signatures/component_models_object.rds"))
```

Create the signatures using the GenerateSignatures() function, and then
save them. Evaluate your handi-work by visualizing both the:  
-\> component-by-signature matrix as a heatmap  
-\> sample-by-signature matrix as a heatmap (contribution of each sample
to each signature)

``` r
# Using the aheatmap functions built into NMF package...
library(NMF)

pdf(file = '~/projects/cn-signatures-project/plots/sample-by-signature_heatmap.pdf', 7, 7, onefile = FALSE)
coefmap(signatures, Colv="consensus", tracks=c("basis:"), main="Sample x Signature matrix")
dev.off()

pdf(file = file = '~/projects/cn-signatures-project/plots/component-by-signature_heatmap.pdf', 7, 7, onefile = FALSE)
basismap(signatures, Rowv = NA, main = "Signature x Component matrix")
dev.off()
```

### V - Example figures

Signature by Component matrix  
<img src="../images/signature-by-component_heatmap.png" width="672" />

Sample by Signature matrix  
<img src="../images/sample-by-signature_heatmap.png" width="671" />
