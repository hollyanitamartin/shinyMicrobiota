# Load Libraries
library(phyloseq)
library(qiime2R)
library(tidyverse)

##### Metadata #####

# Import Metadata. If comma-separated, replace 'tsv' with 'csv'.
metadata <- readr::read_tsv("/path/to/metadata.tsv")
# Select relevant metadata columns (if required).
metadata <- dplyr::select(metadata, c("Colname1", "Colname2", "Colname3", "Colname4"))
# If applicable, create case/control or other metadata column
metadata <- mutate(metadata, Condition = if_else(metadata$Group == "Control", "control", "case"))
# Rename columns to ensure necessary columns labelled correctly (must be these titles). Ensure labels match columns.
colnames(metadata) <- c("sampleID", "individualID", "Timepoint", "Group")
# Convert all metadata columns to factors for plotting. 
# Set labels and levels manually if relevant (i.e. to ensure variables are plotted in a logical order).
metadata$sampleID <- factor(metadata$sampleID)
metadata$individualID <- factor(metadata$individualID)
metadata$Timepoint <- factor(metadata$Timepoint)
metadata$Group <- factor(metadata$Group)
# Write metadata object as a tsv file for use in creating the phyloseq object.
write_tsv(metadata, "/path/to/metadata_R.tsv")

##### Phyloseq #####

# Build phyloseq object from QIIME2 output and metadata table as created above.
physeq <- qiime2R::qza_to_phyloseq(features = "/path/to/feature-table.qza",
                                   tree = "/path/to/rooted-tree.qza",
                                   taxonomy = "/path/to/taxonomy.qza",
                                   metadata = "/path/to/metadata_R.tsv")

##### Alpha Diversity #####

# Import Shannon diversity metrics and add metadata (example for .qza QIIME2 format).
shannon <- qiime2R::read_qza("/path/to/shannon.qza")
shannon <- shannon$data %>% tibble::rownames_to_column("sampleID")
shannon <- left_join(metadata, shannon)

# Import Pielou's evenness diversity metrics and add metadata (example for .tsv format).
evenness <- readr::read_tsv("/path/to/pielou_evenness.tsv")
colnames(evenness) <- c("sampleID", "pielou_evenness") 
evenness <- left_join(metadata, evenness)

# Import Chao1 Index diversity metrics.
chao1 <- qiime2R::read_qza("/path/to/chao1.qza")
chao1 <- chao1$data %>% tibble::rownames_to_column("sampleID")
chao1 <- left_join(metadata, chao1)

##### Beta Diversity #####

# Build beta diversity objects for plotting.
# Transform phyloseq object to account for the compositional nature of the data.
pseq <- microbiome::transform(physeq, "compositional")
# Transform sample counts.
betaCounts <- phyloseq::transform_sample_counts(pseq, function(otu) otu/sum(otu))
# Create ordination from sample counts using bray-curtis distance metric.
betaOrd <- phyloseq::ordinate(betaCounts, method="NMDS", distance="bray")

##### Bacteroidetes to Firmicutes Ratio #####
bfratio <- microbiome::bfratio(physeq) %>% as.data.frame() %>% rownames_to_column()
colnames(bfratio) <- c("sampleID", "bfratio")
bfratio <- left_join(metadata, bfratio)

##### Save all R object as .rds files for upload to shinyMicrobiota.
saveRDS(metadata, file = "/path/to/save/location/1_metadata.rds")
saveRDS(physeq, file = "/path/to/save/location/2_physeq.rds")
saveRDS(shannon, file = "/path/to/save/location/3_shannon.rds")
saveRDS(betaCounts, file = "/path/to/save/location/4_betaCounts.rds")
saveRDS(betaOrd, file = "/path/to/save/location/5_betaOrd.rds")
saveRDS(bfratio, file = "/path/to/save/location/6_bfratio.rds")
saveRDS(evenness, file = "/path/to/save/location/7_evenness.rds")
saveRDS(chao1, file = "/path/to/save/location/8_chao1.rds")