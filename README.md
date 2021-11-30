# shinyMicrobiota

### Summary

shinyMicrobiota is a data visualisation app created in R using the `shiny` library framework. This app has been designed for interactive visualisation of 16S rRNA gene sequencing data. Visualisations include:
- Interactive Data Tables for metadata and alpha diversity values
- Taxonomic composition barplots at each taxonomic level
- Alpha diversity metrics as boxplots or jitter plots
- Beta diversity Non-Metric Multidimensional Scaling (NMDS) plots
- *Bacteroidetes* to *Firmicutes* ratios as boxplots

Each of these plots can be customised using Shiny's reactive programming to display data based on variables present in the metadata file. Data can also be subset to select for a particular variable, for example if you wish to look at the alpha diversity for an individual over various time points, or to select only patients and exclude control samples. Statistical testing is performed for alpha diversity and BF ratio plots by comparing means between groups using the `stat_compare_means()` function from the `ggpubr` package in R. This can be customised to compare means overall between all variables, or by setting a particular group as the reference group (using the radio buttons and drop-down menu), against which all other groups will be compared using a Wilcoxon rank sum test to produce p-values between groups. 

### Quick Start

To launch shinyMicrobiota, download the `app.R` script and run this from RStudio; this will launch the app in your browser. See the 'test_data/plug_and_play' directory to download some test data that can be immediately uploaded to the app. 


### Data processing
If you wish to use your own amplicon data, ata must be pre-processed in R before uploading to shinyMicrobiota. This pre-processing allows for the reactivity of the app to function successfully and ensures that plots are made correctly. The data to be uploaded is: metadata, phyloseq, alpha diversity, two beta diversity files, and the bacteroidetes to firmicutes ratio. These should all be imported into R (for example, using `qiime2R::read_qza()` for files outputted from the QIIME2 workflow, or `dplyr::read_tsv()` for metadata files). The QIIME2 files (metadata, features, taxonomy and tree) are combined to form a `phyloseq` object. Metadata and diversity files should be pre-processed to alter any variables as needed, ensure that the required column names (sampleID, individualID, Timepoint, Group, and Condition) exist in both the metadata and `phyloseq` objects, remove `q2-types` row, and perform any other data wrangling as required for the dataset. After sufficient pre-processing has been performed, all R objects should be saved in RDS format using `base::saveRDS()` which can then be uploaded to shinyMicrobiota. The app's "Data" page provides brief instructions on data requirements and a detailed example of the pre-processing workflow is provided in the `exampleBuild.R` script in the GitHub repository.

### Requirements

shinyMicrobiota requires the following R packages to run: `ggplot2`, `dplyr`, `microbiome`, `ggpubr`, `phyloseq`, `RColorBrewer`, `shinyBS` and `shiny`, which can be installed with:

```
packages <- c("shiny", "shinyBS", "RColorBrewer", "ggpubr","microbiome", "phyloseq", "dplyr", "ggplot2")
if (!require(packages)) install.packages(packages)
```

shinyMicrobiota built with the following software versions:

```
R version 4.0.2 (2020-06-22)
shinyBS_0.61      RColorBrewer_1.1-2  ggpubr_0.4.0
microbiome_1.12.0 phyloseq_1.34.0     dplyr_1.0.7
ggplot2_3.3.5     shiny_1.6.0
```
