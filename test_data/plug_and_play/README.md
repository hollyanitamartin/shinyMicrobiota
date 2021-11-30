# shinyMicrobiota 'Plug and Play' dataset (incomplete)

Here we have included a 'Plug and Play' dataset that can be uploaded to shinyMicrobiota to allow you to have a play around with with the Shiny app without performing your own data processing.

If you prefer to try your own R processing of the QIIME2 output files, see the 'shinyMicrobiota/test_data/unprocessed' directory.

### Data
This dataset consists of 6 files generated from 16S V4 rRNA amplicon sequencing data from [Jangi *et al.*, 2016](https://doi.org/10.1038/ncomms12015). The raw sequencing data is available on NCBI's Sequencing Read Archive (SRA) with the accession number [PRJNA321051](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA321051/).The samples were collected from healthy controls and patients with Multiple Sclerosis. 

### Instructions

1. Download the following files to your local computer:
   - `1_metadata.rds`: Processed metadata file       
   - `2_physeq.rds`: `phyloseq` object created from QIIME2 output files    
   - `3_shannon.rds`: Alpha diversity values (Shannon's Index)
   - `4_betaCounts.rds`: Beta diversity counts
   - `5_betaOrd.rds`: Beta diversity ordination
   - `6_bfratio.rds`: *Bacteroidetes* to *Firmicutes* Ratio
   - `app.R`: shinyMicrobiota script (located in the root directory of this repo)
2. Open `app.R` in RStudio and click the `Run App` button in the top right of the `app.R` window. This will launch shinyMicrobiota in your browser. 
3. Upload the above files to shinyMicrobiota using the instructions provided on the 'Data' landing page.


### Citations

Original Paper:
```
@article{Jangi2016,
abstract = {The gut microbiome plays an important role in immune function and has been implicated in several autoimmune disorders. Here we use 16S rRNA sequencing to investigate the gut microbiome in subjects with multiple sclerosis (MS, n=60) and healthy controls (n=43). Microbiome alterations in MS include increases in Methanobrevibacter and Akkermansia and decreases in Butyricimonas, and correlate with variations in the expression of genes involved in dendritic cell maturation, interferon signalling and NF-kB signalling pathways in circulating T cells and monocytes. Patients on disease-modifying treatment show increased abundances of Prevotella and Sutterella, and decreased Sarcina, compared with untreated patients. MS patients of a second cohort show elevated breath methane compared with controls, consistent with our observation of increased gut Methanobrevibacter in MS in the first cohort. Further study is required to assess whether the observed alterations in the gut microbiome play a role in, or are a consequence of, MS pathogenesis.},
author = {Jangi, Sushrut and Gandhi, Roopali and Cox, Laura M. and Li, Ning and {Von Glehn}, Felipe and Yan, Raymond and Patel, Bonny and Mazzola, Maria Antonietta and Liu, Shirong and Glanz, Bonnie L. and Cook, Sandra and Tankou, Stephanie and Stuart, Fiona and Melo, Kirsy and Nejad, Parham and Smith, Kathleen and Top{\c{c}}uolu, Beg{\"{u}}m D. and Holden, James and Kivis{\"{a}}kk, Pia and Chitnis, Tanuja and {De Jager}, Philip L. and Quintana, Francisco J. and Gerber, Georg K. and Bry, Lynn and Weiner, Howard L.},
doi = {10.1038/ncomms12015},
file = {:home/hollymartin/Documents/Mendeley/Alterations of the human gut microbiome in multiple sclerosis.pdf:pdf},
issn = {20411723},
journal = {Nature Communications},
month = {nov},
number = {1},
pages = {12015},
pmid = {27352007},
title = {{Alterations of the human gut microbiome in multiple sclerosis}},
url = {http://www.nature.com/articles/ncomms12015},
volume = {7},
year = {2016}
}
```

shinyMicrobiota software:
```
@software{Martin_shinyMicrobiota,
author = {Martin, Holly},
title = {{shinyMicrobiota}}
}
```
