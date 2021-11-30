# Load Libraries
library(phyloseq)
library(qiime2R)
library(tidyverse)

##### Metadata #####

# Import Metadata. If comma-separated, replace 'tsv' with 'csv'.
metadata <- readr::read_tsv("/media/Datas/Bioinformatics/MS_testdata/subset_metadata_MS.tsv")
# Remove #q2:types row (must only perform once)
# metadata <- metadata[-1,]
# Reformat/improve poor data entries
metadata$host_phenotype <- str_replace_all(metadata$host_phenotype, 
                                           "Multiple_Sclerosis", "MS")
metadata$host_phenotype <- str_replace_all(metadata$host_phenotype, 
                                           "_-_", "_")
metadata$host_phenotype <- str_replace_all(metadata$host_phenotype, 
                                           "copaxone", "Copaxone")
metadata$host_phenotype <- str_replace_all(metadata$host_phenotype, 
                                           "untreated", "Untreated")
metadata$Host_disease <- str_replace_all(metadata$Host_disease, 
                                           "Multiple_sclerosis", "MS")
# If applicable, create case/control or other metadata column
metadata <- mutate(metadata, Condition = if_else(metadata$Host_disease == "Healthy_Control", "control", "case"))
metadata <- mutate(metadata, individualID = metadata$`sample-id`)
# Rename columns to ensure necessary columns labelled correctly (must be these titles). Ensure labels match columns.
colnames(metadata) <- c("sampleID", "Group", "Phenotype", "Condition", "individualID")
# Convert all metadata columns to factors for plotting. 
# Set labels and levels manually if relevant (i.e. to ensure variables are plotted in a logical order).
metadata$sampleID <- factor(metadata$sampleID)
metadata$individualID <- factor(metadata$individualID)
metadata$Phenotype <- factor(metadata$Phenotype, 
                             levels = c("Healthy_Control", "MS_Untreated",
                                        "MS_Copaxone", "MS_Interferon"),
                             labels = c("Healthy_Control", "MS_Untreated",
                                        "MS_Copaxone", "MS_Interferon"))
metadata$Group <- factor(metadata$Group,
                         levels = c("Healthy_Control", "MS_in_remission"),
                         labels = c("Healthy_Control", "MS_in_remission"))
metadata$Condition <- factor(metadata$Condition)
# Write metadata object as a tsv file for use in creating the phyloseq object.
write_tsv(metadata, "/media/Datas/Bioinformatics/MS_testdata/metadata_MS_R.tsv")

##### Phyloseq #####

# Build phyloseq object from QIIME2 output and metadata table as created above.
physeq <- qiime2R::qza_to_phyloseq(features = "/media/Datas/Bioinformatics/MS_testdata/qiime/filt-feature-table.qza",
                                   tree = "/media/Datas/Bioinformatics/MS_testdata/qiime/rooted-tree.qza",
                                   taxonomy = "/media/Datas/Bioinformatics/MS_testdata/qiime/hybrid-taxonomy.qza",
                                   metadata = "/media/Datas/Bioinformatics/MS_testdata/metadata_MS_R.tsv")

##### Alpha Diversity #####

# Import Shannon diversity metrics and add metadata (example for .qza QIIME2 format).
shannon <- qiime2R::read_qza("/media/Datas/Bioinformatics/MS_testdata/qiime/core-diversity/shannon_vector.qza")
shannon <- shannon$data %>% tibble::rownames_to_column("sampleID")
shannon <- left_join(metadata, shannon)

# Import Pielou's evenness diversity metrics and add metadata (example for .tsv format).
evenness <- qiime2R::read_qza("/media/Datas/Bioinformatics/MS_testdata/qiime/core-diversity/evenness_vector.qza")
evenness <- evenness$data %>% tibble::rownames_to_column("sampleID")
evenness <- left_join(metadata, evenness)

##### Beta Diversity #####

# Build beta diversity objects for plotting.
# Transform phyloseq object to account for the compositional nature of the data.
pseq <- microbiome::transform(physeq, "compositional")
# Transform sample counts.
betaCounts <- phyloseq::transform_sample_counts(pseq, function(otu) otu/sum(otu))
# Create ordination from sample counts using bray-curtis distance metric. Can use any method supported by ordinate().
betaOrd <- phyloseq::ordinate(betaCounts, method="NMDS", distance="bray")

##### Bacteroidetes to Firmicutes Ratio #####
bfratio <- microbiome::bfratio(physeq) %>% as.data.frame() %>% rownames_to_column()
colnames(bfratio) <- c("sampleID", "bfratio")
bfratio <- left_join(metadata, bfratio)

##### Save all R object as .rds files for upload to shinyMicrobiota.
saveRDS(metadata, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/1_metadata.rds")
saveRDS(physeq, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/2_physeq.rds")
saveRDS(shannon, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/3_shannon.rds")
saveRDS(betaCounts, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/4_betaCounts.rds")
saveRDS(betaOrd, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/5_betaOrd.rds")
saveRDS(bfratio, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/6_bfratio.rds")
saveRDS(evenness, file = "/media/Datas/Bioinformatics/shiny/data/MS_dataset/7_evenness.rds")