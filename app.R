library(ggplot2)
library(dplyr)
library(microbiome)
library(ggpubr)
library(phyloseq)
library(RColorBrewer)
library(shinyBS)
library(shiny)

ui <- fluidPage(titlePanel("shinyMicrobiota: A Shiny app for 16S sequencing data visualisation", windowTitle = "shinyMicrobiota"),
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
                                    column(5, h3("How to upload your data"),
                                           tags$ol(
                                             tags$li("Upload data from pre-saved RDS files using the browse buttons to the left.
                                             These files are created using the exampleBuild.R script from QIIME2 output. Metadata column names and phyloseq sample_data must include: 
                                             'sampleID', 'individualID', 'Timepoint', 'Group', and 'Condition', but may include additional columns."), 
                                             tags$li("Once uploaded, the metadata table will be shown on the right and can be interactively browsed using the 
                                                     headers and search functions. The metadata file must be uploaded to view any other file, the remaining files can be uploaded as required."), 
                                             tags$li("The phyloseq file will be summarised below, including information on min/max reads and metadata information."),
                                             tags$li("Alpha diversity can be any numeric measure of alpha diversity (e.g. Shannon index, Chao1 index or Pielou's evenness). The column containing alpha diversity values 
                                                     must be the final column in the alpha diversity table."),
                                             tags$li("Beta diversity should be uploaded in two separate files; 
                                                     betaCounts is a count matrix created using transform_sample_counts(), betaOrd is an ordination
                                                     of this betaCounts object. Both files are required to create the beta diversity plot with plot_ordination()"),
                                             tags$li("The phyloseq summary section will show information about the phyloseq object once it has been created"),
                                             tags$li("For colour-blind friendly plots, select one of the following colour palettes from the drop-down menus: 
                                                     'Paired', 'Dark2', 'Set2', 'RdYlBu'.")
                                    ),
                                           h3("Phyloseq summary"), textOutput("physeq_sum")),
                                    column(4, h3("Metadata"),
                                           dataTableOutput("metadata")),
                           )),
                  tabPanel(title = "Taxonomic Composition",
                           fluidRow(column(2, h3("Input Settings"),
                                           shinyBS::tipify(textInput("barplotTitle", label = "Plot title"),
                                                  title = "Hint: Add your own title before exporting your plot.", 
                                                  placement = "top", trigger = "hover"),
                                           p(selectInput("taxLevel", 
                                                          label = "Taxonomic Level",
                                                          choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                                          multiple = FALSE, selected = "Phylum")),
                                           selectInput("taxaPalette",
                                                       label = "Colour palette",
                                                       choices = c("Set3", "Paired", "Spectral", "RdYlGn", "RdYlBu"),
                                                       selected = "Spectral"),
                                           checkboxInput("subsetSubject_tick", "Subset by individual? [Under Construction]"),
                                           conditionalPanel(condition = "input.subsetSubject_tick" == "TRUE",
                                             uiOutput("subsetSubject_ui"))),
                                           #radioButtons("taxExpType", label = "Export plot as",
                                            #            choiceNames = c("PDF", "PNG"),
                                             #           choiceValues = c("pdf", "png"),
                                              #          inline = TRUE)),
                                           #downloadButton("expBarplot", "Export Barplot")),
                                    column(10, plotOutput("barplot"))
                           )),
                  tabPanel(title = "Alpha Diversity",
                           fluidRow(column(2, h3("Input Settings"),
                                           shinyBS::tipify(textInput("alphaTitle", label = "Plot title"), 
                                                  title = "Hint: Add your own title before exporting your plot.", 
                                                  placement = "top", trigger = "hover"),
                                           p(radioButtons("alphadivGeom", label = "Type of plot",
                                                          choiceNames = c("Jitter", "Boxplot"),
                                                          choiceValues = c("jitter", "boxplot"),
                                                          inline = TRUE)),
                                           uiOutput("ADgroup_ui"),
                                           uiOutput("ADcol_ui"),
                                           radioButtons("checkRef", label = "Compare means:",
                                                        choiceNames = c("Overall", "Against reference group"),
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
                                                                   "Pastel2", "Accent", "Paired", "Spectral", "RdYlBu"),
                                                       selected = "Spectral")),
                                          # radioButtons("alphaExpType", label = "Export plot as",
                                           #             choiceNames = c("PDF", "PNG"),
                                            #            choiceValues = c("pdf", "png"),
                                             #           inline = TRUE)),
                                           #p(downloadButton("expAlpha", "Export Plot"))),
                                    column(width = 5, plotOutput("alphadivPlot")),
                                    column(width = 4, h3("Alpha Diversity table"), dataTableOutput("alphadiv"))),
                           ),
                  tabPanel(title = "Beta Diversity",
                           fluidRow(column(2, h3("Input Settings"),
                                           shinyBS::tipify(textInput("betaTitle", label = "Plot title"),
                                                  title = "Hint: Add your own title before exporting your plot.", 
                                                  placement = "top", trigger = "hover"),
                                           uiOutput("BDcol_ui"),
                                           selectInput("betaPalette",
                                                       label = "Colour palette",
                                                       choices = c("Set1", "Dark2", "Paired", "Spectral", "RdYlBu"),
                                                       selected = "Dark2"),
                                           selectInput("BD_point",
                                                       label = "Point size",
                                                       choices = c(1,2,3,4,5,6),
                                                       selected = 5)),
                                           #radioButtons("betaExpType", label = "Export plot as",
                                          #              choiceNames = c("PDF", "PNG"),
                                           #             choiceValues = c("pdf", "png"),
                                            #            inline = TRUE)),
                                           #p(downloadButton("expBeta", label = "Export plot"))),
                                    column(10, plotOutput("betaPlot")))
                           ),
                  tabPanel(title = "Bacteroidetes to Firmicutes",
                           fluidRow(column(2, h3("Input Settings"),
                                           shinyBS::tipify(textInput("bfTitle", label = "Plot title"),
                                                  title = "Hint: Add your own title before exporting your plot.", 
                                                  placement = "top", trigger = "hover"),
                                           uiOutput("BFgroup_ui"),
                                           uiOutput("BFcol_ui"),
                                           radioButtons("BFcheckRef", label = "Compare means:",
                                                        choiceNames = c("Overall", "Against reference group"),
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
                                                                   "Pastel2", "Accent", "Paired", "Spectral", "RdYlBu"),
                                                       selected = "Spectral")),
                                           #radioButtons("bfExpType", label = "Export plot as",
                                           #             choiceNames = c("PDF", "PNG"),
                                          #              choiceValues = c("pdf", "png"),
                                           #             inline = TRUE)),
                                           #p(downloadButton("expBF", "Export Plot"))),
                                    column(width = 8, plotOutput("bfPlot")))
                  )
              )
)

server <- function(input, output) {
  output$metadata <- renderDataTable({
    file <- input$metadata
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    
    base::readRDS(file$datapath)
  }, options = list(pageLength = 10))
    output$physeq_sum <- renderText({
      file <- input$physeq
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "rds", "Please upload an RDS file"))
      
      physeq <- base::readRDS(file$datapath)
      
      microbiome::summarize_phyloseq(physeq) %>% as.list() %>% as.character()
   })
  output$subsetSubject_ui <- renderUI({
    req(input$subsetSubject_tick)
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- base::readRDS(fileM$datapath)
    Individuals <- metadata$individualID %>% unique() %>% as.character()
    
    selectInput("subsetSubject", "Select Individual ID",
                choices = Individuals)
  })
   plot_barplot <- reactive({
      file <- input$physeq
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "rds", "Please upload an RDS file"))
      physeq <- base::readRDS(file$datapath)
      p <- theme(legend.position = "bottom",
              title = element_text(size = 20),
              legend.title = element_text(size = 14),
              axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95),
              text = element_text(size = 14),panel.background = element_rect(fill="gray95"))

      taxaPlot <- if (input$taxLevel == "Kingdom") {
        kingdom_ps <- microbiome::transform(physeq, "compositional")
        kingdom_ps <- microbiome::aggregate_rare(kingdom_ps, level = "Kingdom", detection = 1/100, prevalence = 50/100)
        kingdom_ps <- microbiome::aggregate_taxa(kingdom_ps, level = "Kingdom")
        
        taxaPlot <- microbiome::plot_composition(kingdom_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "Bacteria") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Kingdom") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Phylum") {
        phylum_ps <- microbiome::transform(physeq, "compositional")
        phylum_ps <- microbiome::aggregate_rare(phylum_ps, level = "Phylum", detection = 1/100, prevalence = 50/100)
        phylum_ps <- phylum_ps %>% microbiome::aggregate_taxa(level = "Phylum") %>%  
          microbiome::transform(transform = "compositional")
        if (input$subsetSubject_tick == TRUE) {
          ## Needs work. If "AYA_0740" instead  of .data[[input$subsetSubject]] or "input.subsetSubject" it works.
          # even with ID written non-dynamically, order of samples not logical (e.g D33, Dx, D79).
          phylum_ps <- phyloseq::subset_samples(phylum_ps, individualID == "input.subsetSubject")
        }
        taxaPlot <- microbiome::plot_composition(phylum_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "Bacteroidetes") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Phylum") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Class") {
        class_ps <- microbiome::transform(physeq, "compositional")
        class_ps <- microbiome::aggregate_rare(class_ps, level = "Class", detection = 1/100, prevalence = 50/100)
        class_ps <- class_ps %>% microbiome::aggregate_taxa(level = "Class") %>%  
          microbiome::transform(transform = "compositional")
        
        taxaPlot <- microbiome::plot_composition(class_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "Bacteroidia") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Class") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Order") {
        order_ps <- microbiome::transform(physeq, "compositional")
        order_ps <- microbiome::aggregate_rare(order_ps, level = "Order", detection = 1/100, prevalence = 50/100)
        order_ps <- microbiome::aggregate_taxa(order_ps, level = "Order") %>%  
          microbiome::transform(transform = "compositional")
        
        taxaPlot <- microbiome::plot_composition(order_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "Bacteroidales") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Order") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Family") {
        family_ps <- microbiome::transform(physeq, "compositional")
        family_ps <- microbiome::aggregate_rare(family_ps, level = "Family", detection = 1/100, prevalence = 50/100)
        family_ps <- microbiome::aggregate_taxa(family_ps, level = "Family") %>%  
          microbiome::transform(transform = "compositional")
        
        taxaPlot <- microbiome::plot_composition(family_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "neatmap") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Family") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Genus") {
        genus_ps <- microbiome::transform(physeq, "compositional")
        genus_ps <- microbiome::aggregate_rare(genus_ps, level = "Genus", detection = 1/100, prevalence = 50/100)
        genus_ps <- microbiome::aggregate_taxa(genus_ps, level = "Genus") %>%  
          microbiome::transform(transform = "compositional")
        
        taxaPlot <- microbiome::plot_composition(genus_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "neatmap") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Genus") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      } else if (input$taxLevel == "Species") {
        species_ps <- microbiome::transform(physeq, "compositional")
        species_ps <- microbiome::aggregate_rare(species_ps, level = "Species", detection = 1/100, prevalence = 50/100)
        species_ps <- microbiome::aggregate_taxa(species_ps, level = "Species") %>%  
          microbiome::transform(transform = "compositional")
        
        taxaPlot <- microbiome::plot_composition(species_ps, x.label = "Sample",
                                                 plot.type = "barplot",
                                                 sample.sort = "neatmap") +
          labs(title = input$barplotTitle,
               x = "Sample", y = "Relative Abundance",
               fill = "Species") +
          scale_fill_brewer(palette = input$taxaPalette) +
          p
        taxaPlot
      }
      taxaPlot
    })
    output$barplot <- renderPlot({
      plot_barplot()
    }, height = 800, width = 1200 )
    output$expBarplot <- downloadHandler(
      filename = function() {
        paste(input$barplotTitle)
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = plot_barplot(), device = input$taxExpType, scale = 2, limitsize = FALSE)
      })
    output$ADgroup_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- base::readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      selectInput("group", "Group by", choices = ColNames,
                  selected = "Group")
    })
    output$ADcol_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- base::readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      selectInput("colour", "Colour by", choices = ColNames,
                  selected = "Group")
    })
    output$ADref_ui <- renderUI({
      fileM <- input$metadata
      extM <- tools::file_ext(fileM$datapath)
      req(fileM)
      validate(need(extM == "rds", "Please upload an RDS file"))
      
      metadata <- base::readRDS(fileM$datapath)
      ColNames <- colnames(metadata)
      Groups <- metadata$Group %>% unique() %>% as.character()
      Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
      Conditions <- metadata$Condition %>% unique() %>% as.character()
      Individuals <- metadata$individualID %>% unique() %>% as.character()
      
      shinyBS::tipify(selectInput("alphaRefGroup", "Reference group",
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
      
      metadata <- base::readRDS(fileM$datapath)
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
      
      metadata <- base::readRDS(fileM$datapath)
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
    output$alphadiv <- renderDataTable({
    file <- input$alphadiv
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    base::readRDS(file$datapath)
  }, options = list(pageLength = 10))
  plot_div <- reactive({
    file <- input$alphadiv
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    alpha <- base::readRDS(file$datapath)
    
    p <- theme_classic(base_size = 14, base_line_size = 0.09, base_rect_size = 0.1) +
      theme(title = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.02, hjust=0.95, size = 14),
            text = element_text(size = 14), line = element_line(size = 0.1), 
            rect = element_rect(size = 0.1), axis.title.y = element_text(vjust = 1), 
            plot.margin = margin(2,2,2,0.5, "cm"), legend.position = "bottom")
    
    if (input$alphadivGeom == "jitter" && input$checkRef == "overall" && input$ADsubsetTick == "FALSE") {
      alphaPlot <- alpha %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          color = input$colour)) + 
        labs(title = input$alphaTitle) +
        geom_jitter(size = 4, alpha = 0.5, width = 0.03) +
        ggplot2::stat_summary(fun = "median", geom = "point", size = 6) +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means() +
        p
      alphaPlot
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "overall" && input$ADsubsetTick == "FALSE") {
      alphaPlot <- alpha %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          fill = input$colour)) +
        labs(title = input$alphaTitle) +
        geom_boxplot() +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means() +
        p
      alphaPlot
    } else if (input$alphadivGeom == "jitter" && input$checkRef == "refGroup" && input$ADsubsetTick == "FALSE") {
      alphaPlot <- alpha %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          color = input$colour)) + 
        labs(title = input$alphaTitle) +
        ggpubr::stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_jitter(size = 4, alpha = 0.5, width = 0.03) +
        ggplot2::stat_summary(fun = "median", geom = "point", size = 6) +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        p
      alphaPlot
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "refGroup" && input$ADsubsetTick == "FALSE") {
      alphaPlot <- alpha %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          fill = input$colour)) +
        labs(title = input$alphaTitle) +
        ggpubr::stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_boxplot() +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        p
      alphaPlot
    } else if (input$alphadivGeom == "jitter" && input$checkRef == "overall" && input$ADsubsetTick == "TRUE") {
      subsetAD <- alpha %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      alphaPlot <- subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          color = input$colour)) + 
        labs(title = input$alphaTitle) +
        geom_jitter(size = 4, alpha = 0.5, width = 0.03) +
        ggplot2::stat_summary(fun = "median", geom = "point", size = 6) +
        ggpubr::stat_compare_means() +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        p
      alphaPlot
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "overall" && input$ADsubsetTick == "TRUE") {
      subsetAD <- alpha %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      alphaPlot <- subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          fill = input$colour)) +
        labs(title = input$alphaTitle) +
        geom_boxplot() +
        ggpubr::stat_compare_means() +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        p
      alphaPlot
    } else if (input$alphadivGeom == "jitter" && input$checkRef == "refGroup" && input$ADsubsetTick == "TRUE") {
      subsetAD <- alpha %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      alphaPlot <- subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          color = input$colour)) + 
        labs(title = input$alphaTitle) +
        ggpubr::stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_jitter(size = 4, alpha = 0.5, width = 0.03) +
        ggplot2::stat_summary(fun = "median", geom = "point", size = 6) +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        p
      alphaPlot
    } else if (input$alphadivGeom == "boxplot" && input$checkRef == "refGroup" && input$ADsubsetTick == "TRUE") {
      subsetAD <- alpha %>% dplyr::filter(.data[[input$ADsubsetCol]] == as.character(input$ADsubset))
      alphaPlot <- subsetAD %>%
        ggplot(aes_string(x = input$group,
                          y = colnames(alpha[,ncol(alpha)]),
                          fill = input$colour)) +
        labs(title = input$alphaTitle) +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$divPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means(ref.group = input$alphaRefGroup) +
        geom_boxplot() +
        p
      alphaPlot
    }
  })
  output$alphadivPlot <- renderPlot({
    plot_div()
  }, height = 800, width = 800 )
  output$expAlpha <- downloadHandler(
    filename = function() {
      paste(input$alphaTitle)
    },
    content = function(file) {
    ggplot2::ggsave(file, plot = plot_div(), device = input$alphaExpType, scale = 2, limitsize = FALSE)
    }
  )
  output$BDcol_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- base::readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BDcolour", "Colour by", choices = ColNames,
                selected = "Group")
  })
  plot_beta <- reactive({
    file1 <- input$betaCounts
    ext1 <- tools::file_ext(file1$datapath)
    
    req(file1)
    validate(need(ext1 == "rds", "Please upload an RDS file"))
    betaCounts <- base::readRDS(file1$datapath)
    
    file2 <- input$betaOrd
    ext2 <- tools::file_ext(file2$datapath)
    
    req(file2)
    validate(need(ext2 == "rds", "Please upload an RDS file"))
    betaOrd <- base::readRDS(file2$datapath)
    
    betaPlot <- phyloseq::plot_ordination(betaCounts, betaOrd, 
                                color = input$BDcolour) +
      geom_point(size = as.numeric(input$BD_point)) +
      labs(title = input$betaTitle) +
      suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$betaPalette),times=6),
                                           aesthetics = c("colour", "fill"))) +
      theme_bw(base_size = 24, 
               base_line_size = 0.09, 
               base_rect_size = 0.1) +
      theme(title = element_text(size = 20),
            axis.text.x = element_text(size = 14, vjust = 0.2, hjust=0.95), 
            text = element_text(size = 14))
    betaPlot
  })
  output$betaPlot <- renderPlot({
    plot_beta()
  }, height = 700, width = 900 )
  output$expBeta <- downloadHandler(
    filename = function() {
      paste(input$betaTitle)
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_beta(), device = input$betaExpType, scale = 2, limitsize = FALSE)
    }
  )
  output$BFgroup_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- base::readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BFgroup", "Group by", choices = ColNames,
                selected = "Group")
  })
  output$BFcol_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- base::readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    selectInput("BFcolour", "Colour by", choices = ColNames,
                selected = "Group")
  })
  output$BFref_ui <- renderUI({
    fileM <- input$metadata
    extM <- tools::file_ext(fileM$datapath)
    req(fileM)
    validate(need(extM == "rds", "Please upload an RDS file"))
    
    metadata <- base::readRDS(fileM$datapath)
    ColNames <- colnames(metadata)
    Groups <- metadata$Group %>% unique() %>% as.character()
    Timepoints <- metadata$Timepoint %>% unique() %>% as.character()
    Conditions <- metadata$Condition %>% unique() %>% as.character()
    Individuals <- metadata$individualID %>% unique() %>% as.character()
    
    shinyBS::tipify(selectInput("BFrefGroup", "Reference group",
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
    
    metadata <- base::readRDS(fileM$datapath)
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
    
    metadata <- base::readRDS(fileM$datapath)
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
  plot_bf <- reactive({
    file <- input$bfratio
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "rds", "Please upload an RDS file"))
    bfratio <- base::readRDS(file$datapath)
    bfTheme <- theme_classic() +
      theme(title = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95, size = 14),
            text = element_text(size = 14),
            line = element_line(size = 0.1), 
            rect = element_rect(size = 0.1), 
            plot.margin = margin(2,2,2,1, "cm"))
    
    if (input$BFsubsetTick == "TRUE" && input$BFcheckRef == "overall") {
      subsetBF <- bfratio %>% dplyr::filter(.data[[input$BFsubsetCol]] == as.character(input$BFsubset))
      BFplot <- subsetBF %>%
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$bfPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means() +
        bfTheme
      BFplot
    } else if (input$BFsubsetTick == "FALSE" && input$BFcheckRef == "overall") {
      BFplot <- bfratio %>% 
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$bfPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means() +
        bfTheme
      BFplot
    } else if (input$BFsubsetTick == "TRUE" && input$BFcheckRef == "refGroup") {
      subsetBF <- bfratio %>% dplyr::filter(.data[[input$BFsubsetCol]] == as.character(input$BFsubset))
      BFplot <- subsetBF %>%
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$bfPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means(ref.group = input$BFrefGroup) +
        bfTheme
      BFplot
    } else if (input$BFsubsetTick == "FALSE" && input$BFcheckRef == "refGroup") {
      BFplot <- bfratio %>% 
        ggplot(aes_string(x = input$BFgroup, y = "bfratio",
                          fill = input$BFcolour)) + 
        geom_boxplot() +
        labs(title = input$bfTitle,
             y = "Bacteroidetes:Firmicutes Ratio") +
        suppressWarnings(scale_colour_manual(values=rep(brewer.pal(12,input$bfPalette),times=6),
                                             aesthetics = c("colour", "fill"))) +
        ggpubr::stat_compare_means(ref.group = input$BFrefGroup) +
        bfTheme
      BFplot
    }
  })
  output$bfPlot <- renderPlot({
      plot_bf()
    }, height = 768, width = 1024, res = 92)
  output$expBF <- downloadHandler(
    filename = function() {
      paste(input$bfTitle)
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_bf(), device = input$bfExpType, scale = 2, limitsize = FALSE)
    }
  )
}

shinyApp(server = server, ui = ui)