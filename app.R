library(tidyverse)
library(qiime2R)
library(microbiome)
library(phyloseq)
library(ggpubr)
library(RColorBrewer)
library(shinyBS)
library(shiny)

### Possible names for the app:
# shinyMicRobiota
# micRobiota
# shinyAmplicon
# microbialVisualisR
# microbiotR
# microbiotaViz
# ############################

ui <- fluidPage(titlePanel("shinyMicRobiota", windowTitle = "16S Microbial Analysis"),
                tabsetPanel(              
                  tabPanel(title = "Data",
                           fluidRow(column(2, h2("Data set up"),
                                           h4(tags$b("Data Upload")),
                                           p(fileInput("metadata", label = "Upload Metadata", 
                                                       multiple = FALSE, accept = ".rds")),
                                           p(fileInput("physeq", label = "Upload phyloseq",
                                                       multiple = FALSE, accept = ".rds")),
                                           p(fileInput("alphadiv", label = "Upload Alpha Diversity", 
                                                       multiple = FALSE, accept = ".rds")),
                                           p(fileInput("betaCounts", label = "Upload Beta Diversity (betaCounts)", 
                                                       multiple = FALSE, accept = ".rds")),
                                           p(fileInput("betaOrd", label = "Upload Beta Diversity (betaOrd)", 
                                                        multiple = FALSE, accept = ".rds")),
                                           p(fileInput("bfratio", label = "Upload Bacteroidetes to Firmicutes ratio", 
                                                       multiple = FALSE, accept = ".rds"))),
                                    column(5, h3("How to upload your data"),"Text Instructions",
                                           tags$ol(
                                             tags$li("Upload data from pre-saved RDS files using the browse buttons to the left.
                                             These files are created using the buildPhylo.R script from QIIME2 output. Metadata column names and phyloseq sample_data must include: 
                                             'sampleID', 'individualID', 'Timepoint', 'Group', and 'Condition', but may include additional column names."), 
                                             tags$li("Once uploaded, the metadata table will be shown on the right and can be interactively browsed using the 
                                                     headers and search functions. The metadata file must be uploaded to view any other file, the remaining files can be uploaded as required."), 
                                             tags$li("The phyloseq file will be summarised below, including information on min/max reads and metadata information."),
                                             tags$li("Alpha diversity can be in the form of Shannon diversity or Pielou's evenness."),
                                             tags$li("Beta diversity can be in the form of bray-curtis or jaccard metrics. This should be uploaded in two separate files; 
                                                     betaCounts is a count matrix created using transform_sample_counts(), betaOrd is an ordination
                                                     of this betaCounts object. Both files are required to create the beta diversity plot with plot_ordination()"),
                                             tags$li("The phyloseq summary section will show information about the phyloseq object once it has been created")
                                             ),
                                           h3("Phyloseq summary"), textOutput("physeq_sum")),
                                    column(4, h3("Metadata"),
                                           dataTableOutput("metadata")),
                           )),
                  tabPanel(title = "Taxa Barplots",
                           fluidRow(column(2, h3("Input Settings"),
                                           textInput("barplotTitle", label = "Plot title"),
                                           p(selectInput("taxLevel", 
                                                          label = "Taxonomic Level",
                                                          choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                                          multiple = FALSE, selected = "Phylum")),
                                           selectInput("taxaPalette",
                                                       label = "Colour palette",
                                                       choices = c("Set3", "Paired", "Spectral", "RdYlGn"),
                                                       selected = "Paired"),
                                           downloadButton("expBarplot", "Export Barplot")),
                                    column(10, h3("Barplot of taxonomic composition"),
                                           plotOutput("barplot"))
                           )),
                  tabPanel(title = "Alpha Diversity",
                           fluidRow(column(2, h3("Input Settings"),
                                           textInput("alphaTitle", label = "Plot title"),
                                           p(radioButtons("alphadivGeom", label = "Type of plot",
                                                          choiceNames = c("Jitter", "Boxplot"),
                                                          choiceValues = c("jitter", "boxplot"),
                                                          inline = TRUE)),
                                           uiOutput("ADgroup_ui"),
                                           uiOutput("ADcol_ui"),
                                           radioButtons("checkRef", label = "Compare means:",
                                                        choiceNames = c("Overall", "With reference group"),
                                                        choiceValues = c("overall", "refGroup"),
                                                        inline = TRUE),
                                           conditionalPanel("input.checkRef == 'refGroup'",
                                                            uiOutput("ADref_ui")),
                                           checkboxInput("ADsubsetTick", "Subset data?"),
                                           conditionalPanel(condition = "input.ADsubsetTick" == "TRUE",
                                                            uiOutput("ADsubsetCol_ui"),
                                                            uiOutput("ADsubsetVar_ui")),
                                           selectInput("divPalette",
                                                       label = "Colour palette",
                                                       choices = c("Dark2", "Set1", "Set2", "Set3", 
                                                                   "Pastel2", "Accent", "Paired", "Spectral"),
                                                       selected = "Dark2"),
                                           p(downloadButton("expAlpha", "Export Plot"))),
                                    column(width = 5, plotOutput("alphadivPlot")),
                                    column(width = 4, h3("Alpha Diversity table"), dataTableOutput("alphadiv"))),
                           ),
                  tabPanel(title = "Beta Diversity",
                           fluidRow(column(2, h3("Input Settings"),
                                           textInput("betaTitle", label = "Plot title"),
                                           uiOutput("BDcol_ui"),
                                           selectInput("betaPalette",
                                                       label = "Colour palette",
                                                       choices = c("Set1", "Dark2", "Paired"),
                                                       selected = "Set1"),
                                           p(downloadButton("expBeta", label = "Export plot"))),
                                    column(10, h3("Beta diversity plot"),
                                           plotOutput("betaPlot")))
                           ),
                  tabPanel(title = "Bacteroidetes to Firmicutes",
                           fluidRow(column(2, h3("Input Settings"),
                                           textInput("bfTitle", label = "Plot title"),
                                           uiOutput("BFgroup_ui"),
                                           uiOutput("BFcol_ui"),
                                           radioButtons("BFcheckRef", label = "Compare means:",
                                                        choiceNames = c("Overall", "With reference group"),
                                                        choiceValues = c("overall", "refGroup"),
                                                        inline = TRUE),
                                           conditionalPanel("input.BFcheckRef == 'refGroup'",
                                                            uiOutput("BFref_ui")),
                                           checkboxInput("BFsubsetTick", "Subset data?"),
                                           conditionalPanel(condition = "input.BFsubsetTick" == "TRUE",
                                             uiOutput("BFsubsetCol_ui"),
                                             uiOutput("BFsubsetVar_ui")),
                                           selectInput("bfPalette",
                                                       label = "Colour palette",
                                                       choices = c("Dark2", "Set1", "Set2", "Set3", 
                                                                   "Pastel2", "Accent", "Paired", "Spectral"),
                                                       selected = "Dark2"),
                                           p(downloadButton("expBF", "Export Plot"))),
                                    column(width = 5, plotOutput("bfPlot")))
                  )
              )
)

server <- function(input, output) {
  output$metadata <- renderDataTable({
    file <- input$metadata
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    
    readRDS(file$datapath)
  }, options = list(pageLength = 10))
    output$physeq_sum <- renderText({
      file <- input$physeq
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "rds", "Please upload an RDS file"))
      
      physeq <- readRDS(file$datapath)
      
      microbiome::summarize_phyloseq(physeq) %>% as.list() %>% as.character()
   })
    output$ADgroup_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      selectInput("group", "Group by", choices = ColNames,
                  selected = "Timepoint")
    })
    output$ADcol_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      selectInput("colour", "Colour by", choices = ColNames,
                  selected = "Timepoint")
    })
    output$ADref_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      Groups <- metadata$Group %>% unique() %>% as.character()
      Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
      Conditions <- metadata$Condition %>% unique() %>% as.character()
      Individuals <- metadata$individualID %>% unique() %>% as.character()
      

        tipify(selectInput("alphaRefGroup", "Reference group",
                           choices = if (input$group == "individualID") { Individuals }
                           else if (input$group == "Timepoint") { Timepoints }
                           else if (input$group == "Group") { Groups }
                           else if (input$group == "Condition") { Conditions }), 
               title = "Hint: To compare means ensure that the Group and Colour options are the same.", 
               placement = "top", trigger = "hover")
    })
    output$ADsubsetCol_ui <- renderUI({
      req(input$ADsubsetTick)
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      Groups <- metadata$Group %>% unique() %>% as.character()
      Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
      Conditions <- metadata$Condition %>% unique() %>% as.character()
      Individuals <- metadata$individualID %>% unique() %>% as.character()
      
      selectInput("ADsubsetCol", "Subset data by",
                  choices = ColNames, selected = "Group")
    })
    output$ADsubsetVar_ui <- renderUI({
      req(input$ADsubsetTick)
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      Groups <- metadata$Group %>% unique() %>% as.character()
      Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
      Conditions <- metadata$Condition %>% unique() %>% as.character()
      Individuals <- metadata$individualID %>% unique() %>% as.character()
      Samples <- metadata$sampleID %>% unique() %>% as.character()
      
      selectInput("ADsubset", "Select variable to filter",
                  choices = if (input$ADsubsetCol == "individualID") { Individuals }
                  else if (input$ADsubsetCol == "Timepoint") { Timepoints }
                  else if (input$ADsubsetCol == "Group") { Groups }
                  else if (input$ADsubsetCol == "Condition") { Conditions }
                  else if (input$ADsubsetCol == "sampleID") { Samples })
    })
    output$barplot <- renderPlot({
      file <- input$physeq
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "rds", "Please upload an RDS file"))
      physeq <- readRDS(file$datapath)
      
      if (input$taxLevel == "Kingdom") {
        kingdom_ps <- microbiome::transform(physeq, "compositional")
        kingdom_ps <- aggregate_rare(kingdom_ps, level = "Kingdom", detection = 1/100, prevalence = 50/100)
        kingdom_ps <- aggregate_taxa(kingdom_ps, level = "Kingdom")
        
        plot_composition(kingdom_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Bacteria") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = "Kingdom-level Relative abundance",
               x = "Sample", y = "Relative Abundance",
               fill = "Kingdom") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      } else if (input$taxLevel == "Phylum") {
          phylum_ps <- microbiome::transform(physeq, "compositional")
          phylum_ps <- aggregate_rare(phylum_ps, level = "Phylum", detection = 1/100, prevalence = 50/100)
          phylum_ps <- phylum_ps %>% aggregate_taxa(level = "Phylum") %>%  
            microbiome::transform(transform = "compositional")
          
          plot_composition(phylum_ps, x.label = "Sample",
                           plot.type = "barplot",
                           sample.sort = "Bacteroidetes") +
            theme(legend.position = "bottom", 
                  legend.title = element_text(size = 14), 
                  axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(title = input$barplotTitle,
                 x = "Sample", y = "Relative Abundance",
                 fill = "Phylum") +
            scale_fill_brewer(palette = input$taxaPalette) +
            theme(axis.text.x = element_text(angle = 90, size = 14), 
                  text = element_text(size = 14))
      } else if (input$taxLevel == "Class") {
        class_ps <- microbiome::transform(physeq, "compositional")
        class_ps <- aggregate_rare(class_ps, level = "Class", detection = 1/100, prevalence = 50/100)
        class_ps <- class_ps %>% aggregate_taxa(level = "Class") %>%  
          microbiome::transform(transform = "compositional")
        
        plot_composition(class_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Bacteroidia") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Class") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      } else if (input$taxLevel == "Order") {
        order_ps <- microbiome::transform(physeq, "compositional")
        order_ps <- aggregate_rare(order_ps, level = "Order", detection = 1/100, prevalence = 50/100)
        order_ps <- aggregate_taxa(order_ps, level = "Order") %>%  
          microbiome::transform(transform = "compositional")
        
        plot_composition(order_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Bacteroidales") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Order") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      } else if (input$taxLevel == "Family") {
        family_ps <- microbiome::transform(physeq, "compositional")
        family_ps <- aggregate_rare(family_ps, level = "Family", detection = 1/100, prevalence = 50/100)
        family_ps <- aggregate_taxa(family_ps, level = "Family") %>%  
          microbiome::transform(transform = "compositional")
        
        plot_composition(family_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Bacteroidaceae") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Family") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      } else if (input$taxLevel == "Genus") {
        genus_ps <- microbiome::transform(physeq, "compositional")
        genus_ps <- aggregate_rare(genus_ps, level = "Genus", detection = 1/100, prevalence = 50/100)
        genus_ps <- aggregate_taxa(genus_ps, level = "Genus") %>%  
          microbiome::transform(transform = "compositional")
        
        plot_composition(genus_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Bacteroides") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Genus") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      } else if (input$taxLevel == "Species") {
        species_ps <- microbiome::transform(physeq, "compositional")
        species_ps <- aggregate_rare(species_ps, level = "Species", detection = 1/100, prevalence = 50/100)
        species_ps <- aggregate_taxa(species_ps, level = "Species") %>%  
          microbiome::transform(transform = "compositional")
        
        plot_composition(species_ps, x.label = "Sample",
                         plot.type = "barplot",
                         sample.sort = "Unknown") +
          theme(legend.position = "bottom", 
                legend.title = element_text(size = 14), 
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Species") +
          scale_fill_brewer(palette = input$taxaPalette) +
          theme(axis.text.x = element_text(angle = 90, size = 14), 
                text = element_text(size = 14))
      }
    }, height = 600, width = 1000 )
    output$expBarplot <- downloadHandler(
      filename = function() {
        paste(input$barplotTitle)
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = last_plot(), device = "png")
      })
    output$alphadiv <- renderDataTable({
    file <- input$alphadiv
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    readRDS(file$datapath)
  }, options = list(pageLength = 10))
  output$alphadivPlot <- renderPlot({
    file <- input$alphadiv
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    shannon <- readRDS(file$datapath)
    
    healthy_avg <- subset(shannon, shannon$Condition == "control") %>% 
      summarize(avg = median(shannon_entropy)) %>%
      pull(avg)
    
    p <- theme_bw(base_size = 14, 
               base_family = "Open Sans", 
               base_line_size = 0.09, 
               base_rect_size = 0.1) +
      theme(axis.text.x = element_text(size = 14), 
            text = element_text(size = 14),
            legend.position = "bottom")
    
    if (input$alphadivGeom == "jitter" && input$checkRef == "overall" && input$ADsubsetTick == "FALSE") {
      shannon %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          color = input$colour)) + 
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_color_brewer(palette = input$divPalette) +
        geom_hline(aes(yintercept = healthy_avg), color = "gray70", size = 0.6) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.03) +
        stat_summary(fun = "median", geom = "point", size = 3) +
        stat_compare_means() +
        p 
      } else if (input$alphadivGeom == "boxplot" && input$checkRef == "overall" && input$ADsubsetTick == "FALSE") {
      shannon %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          fill = input$colour)) +
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_fill_brewer(palette = input$divPalette) +
        geom_boxplot() +
        stat_compare_means() +
        p
      } else if (input$alphadivGeom == "jitter" && input$checkRef == "refGroup" && input$ADsubsetTick == "FALSE") {
      shannon %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          color = input$colour)) + 
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_color_brewer(palette = input$divPalette) +
        stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_hline(aes(yintercept = healthy_avg), color = "gray70", size = 0.6) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.03) +
        stat_summary(fun = "median", geom = "point", size = 3) +
        p 
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "refGroup") {
      shannon %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          fill = input$colour)) +
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_fill_brewer(palette = input$divPalette) +
        stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_boxplot() +
        p
    } else  if (input$alphadivGeom == "jitter" && input$checkRef == "overall" && input$ADsubsetTick == "TRUE") {
      subsetAD <- shannon %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          color = input$colour)) + 
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_color_brewer(palette = input$divPalette) +
        geom_hline(aes(yintercept = healthy_avg), color = "gray70", size = 0.6) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.03) +
        stat_summary(fun = "median", geom = "point", size = 3) +
        stat_compare_means() +
        p 
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "overall" && input$ADsubsetTick == "TRUE") {
      subsetAD <- shannon %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          fill = input$colour)) +
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_fill_brewer(palette = input$divPalette) +
        geom_boxplot() +
        stat_compare_means() +
        p
    } else if (input$alphadivGeom == "jitter" && input$checkRef == "refGroup" && input$ADsubsetTick == "TRUE") {
      subsetAD <- shannon %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          color = input$colour)) + 
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_color_brewer(palette = input$divPalette) +
        stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_hline(aes(yintercept = healthy_avg), color = "gray70", size = 0.6) +
        geom_jitter(size = 2, alpha = 0.5, width = 0.03) +
        stat_summary(fun = "median", geom = "point", size = 3) +
        p 
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "refGroup" && input$ADsubsetTick == "TRUE") {
      subsetAD <- shannon %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = "shannon_entropy",
                          fill = input$colour)) +
        labs(y = "Shannon Diversity",
             title = input$alphaTitle) +
        scale_fill_brewer(palette = input$divPalette) +
        stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_boxplot() +
        p
    }
  }, height = 600, width = 700 )
  output$expAlpha <- downloadHandler(
    filename = function() {
      paste(input$alphaTitle)
    },
    content = function(file) {
    ggplot2::ggsave(file, plot = last_plot(), device = "png")
    }
  )
  output$BDcol_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BDcolour", "Colour by", choices = ColNames,
                selected = "Timepoint")
  })
  output$betaPlot <- renderPlot({
    file1 <- input$betaCounts
    ext1 <- tools::file_ext(file1$datapath)
    
    req(file1)
    validate(need(ext1 == "rds", "Please upload an RDS file"))
    betaCounts <- readRDS(file1$datapath)
    
    file2 <- input$betaOrd
    ext2 <- tools::file_ext(file2$datapath)
    
    req(file2)
    validate(need(ext2 == "rds", "Please upload an RDS file"))
    betaOrd <- readRDS(file2$datapath)
    
    plot_ordination(betaCounts, betaOrd, 
                    color = input$BDcolour) +
      geom_point(size = 3) +
      labs(title = input$betaTitle) +
      scale_color_brewer(palette = input$betaPalette) +
      theme_bw(base_size = 24, 
               base_family = "Open Sans", 
               base_line_size = 0.09, 
               base_rect_size = 0.1) +
      theme(axis.text.x = element_text(size = 14), 
            text = element_text(size = 14))
  }, height = 700, width = 900 )
  output$expBeta <- downloadHandler(
    filename = function() {
      paste(input$betaTitle)
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = last_plot(), device = "png")
    }
  )
  output$BFgroup_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BFgroup", "Group by", choices = ColNames,
                selected = "Timepoint")
  })
  output$BFcol_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BFcolour", "Colour by", choices = ColNames,
                selected = "Timepoint")
  })
  output$BFref_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    Groups <- metadata$Group %>% unique() %>% as.character()
    Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
    Conditions <- metadata$Condition %>% unique() %>% as.character()
    Individuals <- metadata$individualID %>% unique() %>% as.character()
    
    tipify(selectInput("BFrefGroup", "Reference group",
                       choices = if (input$BFgroup == "individualID") { Individuals }
                       else if (input$BFgroup == "Timepoint") { Timepoints }
                       else if (input$BFgroup == "Group") { Groups }
                       else if (input$BFgroup == "Condition") { Conditions }), 
           title = "Hint: To compare means ensure that the Group and Colour options are the same.", 
           placement = "top", trigger = "hover")
  })
  output$BFsubsetCol_ui <- renderUI({
    req(input$BFsubsetTick)
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    Groups <- metadata$Group %>% unique() %>% as.character()
    Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
    Conditions <- metadata$Condition %>% unique() %>% as.character()
    Individuals <- metadata$individualID %>% unique() %>% as.character()

    selectInput("BFsubsetCol", "Subset data by",
                choices = ColNames, selected = "Group")
  })
  output$BFsubsetVar_ui <- renderUI({
    req(input$BFsubsetTick)
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    Groups <- metadata$Group %>% unique() %>% as.character()
    Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
    Conditions <- metadata$Condition %>% unique() %>% as.character()
    Individuals <- metadata$individualID %>% unique() %>% as.character()
    Samples <- metadata$sampleID %>% unique() %>% as.character()
    
    selectInput("BFsubset", "Select variable to filter",
                choices = if (input$BFsubsetCol == "individualID") { Individuals }
                          else if (input$BFsubsetCol == "Timepoint") { Timepoints }
                          else if (input$BFsubsetCol == "Group") { Groups }
                          else if (input$BFsubsetCol == "Condition") { Conditions }
                          else if (input$BFsubsetCol == "sampleID") { Samples })
  })
    output$bfPlot <- renderPlot({
      file <- input$bfratio
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "rds", "Please upload an RDS file"))
      bfratio <- readRDS(file$datapath)
      
    if (input$BFsubsetTick == "TRUE" && input$BFcheckRef == "overall") {
      subsetBF <- bfratio %>% dplyr::filter(.data[[input$BFsubsetCol]] == as.character(input$BFsubset))
      subsetBF %>%
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                             fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        scale_fill_brewer(input$bfPalette) +
        stat_compare_means() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
              text = element_text(size = 14),
              line = element_line(size = 0.1), 
              rect = element_rect(size = 0.1))
    } else if (input$BFsubsetTick == "FALSE" && input$BFcheckRef == "overall") {
        bfratio %>% 
          ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                            fill = input$BFcolour)) + 
          geom_boxplot() +
          labs(title = input$bfTitle,
               y = "Bacteroidetes:Firmicutes Ratio") +
          scale_fill_brewer(input$bfPalette) +
          stat_compare_means() +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
                text = element_text(size = 14),
                line = element_line(size = 0.1), 
                rect = element_rect(size = 0.1))
    } else if (input$BFsubsetTick == "TRUE" && input$BFcheckRef == "refGroup") {
      subsetBF <- bfratio %>% dplyr::filter(.data[[input$BFsubsetCol]] == as.character(input$BFsubset))
      subsetBF %>%
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        scale_fill_brewer(input$bfPalette) +
        stat_compare_means(ref.group = input$BFrefGroup) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
              text = element_text(size = 14),
              line = element_line(size = 0.1), 
              rect = element_rect(size = 0.1))
    } else if (input$BFsubsetTick == "FALSE" && input$BFcheckRef == "refGroup") {
      bfratio %>% 
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        scale_fill_brewer(input$bfPalette) +
        stat_compare_means(ref.group = input$BFrefGroup) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
              text = element_text(size = 14),
              line = element_line(size = 0.1), 
              rect = element_rect(size = 0.1))
    }
    }) 
}

shinyApp(server = server, ui = ui)