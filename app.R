library(tidyverse)
library(qiime2R)
library(microbiome)
library(shiny)

ui <- fluidPage(titlePanel("16S Microbial Analysis", windowTitle = "16S Microbial Analysis"),
                tabsetPanel(              
                  tabPanel(title = "Data",
                           wellPanel(fluidRow(h3("Information about the page"), 
                                              "Brief Description of data expected/phyloseq object/qiime2 data")),
                           fluidRow(column(5, h3("Data set up"), 
                                           br("File Upload or dataset selection from pre-loaded data"),
                                           p(fileInput("metadata", label = "Upload Metadata", 
                                                       multiple = FALSE, accept = ".tsv")),
                                           p(fileInput("featureTable", label = "Upload Feature Table", 
                                                 multiple = FALSE, accept = ".qza")),
                                           p(fileInput("tree", label = "Upload Rooted Tree", 
                                                         multiple = FALSE, accept = ".qza")),
                                           p(fileInput("taxonomy", label = "Upload Taxonomy", 
                                                         multiple = FALSE, accept = ".qza")),
                                           p(checkboxInput("compTransform", label = "Transform data (compositional)")),
                                           p(actionButton(inputId = "buildps", label = "Build Phyloseq"))),
                                    column(7, h3("How to upload your data"),"Text Instructions")),
                           fluidRow(column(12, h3("Phyloseq summary"), p("Phyloseq Object Summary with microbiome::summarize_phyloseq()"),
                                           verbatimTextOutput("physeq_sum"))),
                           fluidRow(column(12, h3("Metadata"),
                                    dataTableOutput("metadata")))
                           ),
                  tabPanel(title = "Taxa Barplots",
                           wellPanel(fluidRow(h3("Information about the page"), "Brief Description of taxonomic composition, 
                                              relative abundances etc.")),
                           fluidRow(column(5, h3("Input Settings"), "[tax level, colour palette, plot title,
                                           order by, export button, update button]"), 
                                    column(7, h3("Barplot of taxonomic composition")))),
                  tabPanel(title = "Alpha Diversity",
                           wellPanel(fluidRow(h3("Information about the page"), 
                                              "Brief Description of alpha diversity, definitions, metrics")),
                           fluidRow(column(width = 5, h3("Input Settings"), "[plot title, metric (shannon, chao1), 
                                           geom (boxplot, jitter), metadata to plot (x axis, facet_grid),
                                           colour by, stat_compare_means tick box (ref group?),
                                           filter by metadata value (e.g. select only controls), 
                                           export plot button (PNG, pdf), update button]"), 
                                    column(width = 7, h3("Alpha diversity plot"), "boxplot or geom_jitter with median plotted"))),
                  tabPanel(title = "Beta Diversity",
                           wellPanel(fluidRow(h3("Information about the page"), 
                                              "Brief Description of beta diversity, definitions, metrics")),
                           fluidRow(column(5, h3("Input Settings"), "type (NMDS, PCA), colour by (metadata), plot title, 
                                           export plot, update button"), 
                                    column(7, h3("Beta diversity plot with legend"))))
              )
)

server <- function(input, output) {
  output$metadata <- renderDataTable({
    file <- input$metadata
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "tsv", "Please upload a tsv file"))
    
    read_tsv(file$datapath)
  })
  physeq <- observeEvent(input$buildps, {
    qiime2R::qza_to_phyloseq(features = input$featureTable,
                             tree = input$tree,
                             taxonomy = input$taxonomy,
                             metadata = input$metadata)
  })
  output$physeq_sum <- renderText({
    microbiome::summarize_phyloseq(physeq)
  })
}

shinyApp(server = server, ui = ui)