# utanos

### Package Goal: 
To be a swiss army knife in analyzing shallow/low-pass WGS data

### Package Description/Details: 
Includes functions for plotting, quality evaluation, data re-structuring, and copy-number signature inference. 
Additionally there are tools that export analyses/datsets in formats useful for human examination. (i.e. for non-computational purposes).  

### Installation Instructions:
Install R (version > 4.2)  
[Cran would be the place to find the right R](https://cran.r-project.org/index.html)  

Ensure the following basic packages are installed.
```R
install.packages("librarian")
BiocManager::install(version = "3.16")
BiocManager::install("Biobase")
```

Install the bulk of this package's dependencies making use of librarian:
```R
librarian::shelf(caret, CGHcall, data.table, DBI, DescTools, devtools, doMC, dplyr, EnsDb.Hsapiens.v75, flexmix, GenomicRanges, ggalt, ggplot2, ggpubr, ggrepel, gridExtra, hrbrthemes, ks, magrittr, NMF, pheatmap, plyr, purrr, QDNAseq, readr, RMySQL, stringr, tidyr, viridis, YAPSA)
```

Install annotables from github:
```R
devtools::install_github("stephenturner/annotables")
```

Finally, install utanos:
```R
install_github("Huntsmanlab/utanos")
```


_____________________________________________________________________________________________________________________________________


### File list:
* __utils.R__ - General purpose utility functions for anything related to the package.
* __plotting.R__ - General purpose functions for anything related to plotting data analyzed using the package.
* __samplequality.R__ - Functions related to evaluating the quality of each sample being analyzed.
* __main_functions.R__ - Central logic functions used in the calling of copy-number signatures. Needs to be re-named.
* __helper_functions.R__ - Functions assisting in the calling of copy-number signatures.  Needs to be re-named.
* __makeFilteredCNSummaryTable.R__ - Contains functions for taking a Relative Copy Number (CN) QDNAseq object and writing out a TSV table summarizing the observed CN changes.
* __signaturegeneration.R__ - Code to calculate copy-number signature exposures for a list of samples given their segmented copy-number profiles.
* __wisecondorXtoQDNAseq.R__ - Code to build a QDNAseq object from rCN WisecondorX output files.



