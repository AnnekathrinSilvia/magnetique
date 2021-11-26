# loading libraries -------------------------------------------------------
library("GeneTonic")
library("DESeq2")
library("shiny")
library("shinydashboard")
library("shinycssloaders")
library("visNetwork")
library("bs4Dash")
library("ggplot2")
library("ggrepel")
library("igraph")
library("plotly")

options(spinner.type = 6)

# sourcing external files -------------------------------------------------
source("volcano_plot.R")

source("DTU/plots_with_se_obj.R") # this one takes long! promiseable?
source("utils.R")

# this one is for loading data, candidate to be `promise`d
source("load_data_magnetique.R")

# ui definition -----------------------------------------------------------

magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  
  header = shinydashboard::dashboardHeader(title = "magnetique"),
  # header = bs4Dash::bs4DashNavbar(
    # controlbarIcon = icon("cogs")
  # ),

  # sidebar definition ------------------------------------------------------
  sidebar = shinydashboard::dashboardSidebar(
    title = "Options",
    
    selectInput("selected_contrast",
                label = "Contrast id",
                choices = c("DCMvsHCM", 
                            "DCMvsNFD", 
                            "HCMvsNFD"),
                selected = "DCMvsHCM"),
    selectInput("selected_ontology",
                label = "Ontology",
                choices = c("BP", "MF", "CC"),
                selected = "BP"),
    numericInput("number_genesets", 
                 "Number of genesets",
                 value = 15,
                 min = 0),
    selectInput("color_by",
                "Color by",
                choices = c("z_score",
                            "gs_pvalue"), 
                selected = "z_score")
    
  ),
    

  # body definition ---------------------------------------------------------
  body = shinydashboard::dashboardBody(
    tabBox(
      width = 12,
      shiny::tabPanel(
        title = "Welcome!", icon = icon("magnet"), value = "tab-welcome",
        fluidRow(
          column(
            width = 12,
            includeMarkdown("data/overview.md")
          )
        )
      ),
      shiny::tabPanel(
        title = "DE!", icon = icon("heartbeat"), value = "tab-de",
        fluidRow(
          column(
            width = 5,
            DT::dataTableOutput("de_table")
          ),
          column(
            width = 4,
            withSpinner(
              plotlyOutput("de_volcano")
            )
          )
        )
      ),
      shiny::tabPanel(
        title = "Enrichment map!", icon = icon("project-diagram"), value = "tab-emap",
        fluidRow(
          column(
            width = 8,
            withSpinner(
              visNetworkOutput("visnet_em")
            ) # ,
            # withSpinner(
            #   visNetworkOutput("visnet_ggs")
            # )
          ),
          column(
            width = 4,
            plotOutput("emap_signature")
          )
        )
      ),
      shiny::tabPanel(
        title = "DTU!", icon = icon("flask"), value = "tab-dtu",
        fluidRow(
          column(
            width = 12,
            selectizeInput(
                "gene_name", "Choose one gene:", choices = NULL
              ),
            withSpinner(
              plotOutput("dtu_plot")
              ),
            tableOutput("dtu_table")
          )
        )
      ),
      shiny::tabPanel(
        title = "Carnival!", icon = icon("sitemap"), value = "tab-carnival",
        fluidRow(
          column(
            width = 8,
            withSpinner(
              visNetworkOutput("visnet_igraph")
            )
          ),
          column(
            width = 4,
            plotOutput("carnival_counts")
          )
        )
      ),
      shiny::tabPanel(
        title = "About us", icon = icon("users"), value = "tab-aboutus",
        fluidRow(
          column(
            width = 12,
            h2('Project members (alphabetical order)'),
            tableOutput("team_list")
          )
        )
      )
    )
  )
)


# server definition -------------------------------------------------------

magnetique_server <- function(input, output, session) {
  
  rvalues <- reactiveValues()
  rvalues$mygtl <- NULL
  rvalues$myigraph <- NULL
  

  # selector of gtl object --------------------------------------------------
  rvalues$mygtl <- reactive({
    
    message(input$selected_contrast)
    message(input$selected_ontology)
    
    all_gtls[[input$selected_contrast]][[input$selected_ontology]]
    # all_gtls[["DCMvsHCM"]][["BP"]]
  })
  
  rvalues$myigraph <- reactive({
    all_igraph[[input$selected_contrast]]
    # all_igraph[["DCMvsHCM"]]
  })
  
  # DE related content ---------------------------------------------------------
  output$de_table <- DT::renderDataTable({
    mygtl <- rvalues$mygtl()
    myde <- mygtl$res_de
    myde <- GeneTonic::deseqresult2df(myde)
    ensembl_url <- "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
    rownames(myde) <- lapply(rownames(myde), function(x) format_url(ensembl_url, x))
    myde <- myde[c("SYMBOL", "log2FoldChange", "padj")]
    DT::datatable(myde, escape = FALSE, options = list(scrollX = TRUE))  %>% 
      DT::formatRound(columns=c('log2FoldChange', 'padj'), digits=3)
  })
  
  output$de_volcano <- renderPlotly({
    mygtl <- rvalues$mygtl()
    myde <- mygtl$res_de
    p <- volcano_plot(myde, mygtl$annotation_obj, volcano_labels = 0)
    plotly::ggplotly(p, tooltip = "text")
  })
  
  # enrichment map related content ---------------------------------------------
  emap_graph <- reactive({
    emg <- enrichment_map(
      gtl = rvalues$mygtl(),
      n_gs = input$number_genesets,
      overlap_threshold = 0.1,
      scale_edges_width = 200,
      color_by = input$color_by
      # color_by = "z_score"
    )
    # rank_gs <- rank(V(emg)$name)
    # emg <- permute.vertices(emg, rank_gs)
    return(emg)
  })
  

  output$visnet_em <- renderVisNetwork({
    visNetwork::visIgraph(emap_graph()) %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        nodesIdSelection = TRUE
      ) %>%
      visExport(
        name = "emap_network",
        type = "png",
        label = "Save enrichment map"
      )
  })
  
  output$emap_signature <- renderPlot({
    cur_gsid <- rvalues$mygtl()$res_enrich$gs_id[match(input$visnet_em_selected, rvalues$mygtl()$res_enrich$gs_description)]
    validate(need(!is.na(cur_gsid),
                  message = "Please select a gene set from the Enrichment Map."
    ))
    
    # if (!is.null(input$exp_condition)) {
    # message(cur_gsid)
    gs_heatmap(
      se = vst(rvalues$mygtl()$dds) ,
      gtl = rvalues$mygtl(),
      geneset_id = cur_gsid,
      FDR = 0.05,
      de_only = FALSE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      center_mean = TRUE,
      scale_row = TRUE,
      anno_col_info = "Etiology"
    )
    # } else {
    #   gs_heatmap(
    #     myvst,
    #     res_de,
    #     res_enrich,
    #     annotation_obj = annotation_obj,
    #     geneset_id = cur_gsid,
    #     FDR = input$de_fdr,
    #     de_only = FALSE,
    #     cluster_rows = TRUE,
    #     cluster_columns = TRUE,
    #     center_mean = TRUE,
    #     scale_row = TRUE
    #   )
    # }
  })
  
  
  # DTU related content --------------------------------------------------------
  genes_dtu <- unique(rowData(se)$gene_name)
  updateSelectizeInput(
    session, "gene_name",
    choices = genes_dtu, server = TRUE, selected=NULL
  )
  
  output$dtu_plot <- renderPlot({
    req(input$gene_name)
    gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
    plot_dtu(mcols(gtf_gene)[["gene_id"]], se, gtf)
    
  })
  
  output$dtu_table <- renderTable({
    req(input$gene_name)
    gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
    results_table(mcols(gtf_gene)[["gene_id"]], se) 
    
  })
  
  # carnival-related content ---------------------------------------------------
  output$carnival_counts <- renderPlot({
    mygtl <- rvalues$mygtl()
    
    g <- rvalues$myigraph()
    cur_sel <- input$visnet_igraph_selected
    cur_node <- match(cur_sel, V(g)$name)
    cur_nodetype <- V(g)$nodetype[cur_node]
    # validate(need(cur_nodetype == "Feature",
    #               message = "" # "Please select a gene/feature."
    # ))
    
    validate(need(cur_sel != "",
                  message = "Please select a node from the graph to plot the expression values."
    ))
    
    cur_geneid <- mygtl$annotation_obj$gene_id[match(cur_sel, mygtl$annotation_obj$gene_name)]
    
    message(length(cur_sel))
    message(cur_geneid)
    
    genes_exp <- rownames(mygtl$dds)
    validate(need(cur_geneid %in% genes_exp,
                  message = "gene not found in expression matrix" 
    ))
    
    gene_plot(gtl = mygtl, gene = cur_geneid, 
              intgroup = "Etiology")
  })
  
  
  output$visnet_igraph <- renderVisNetwork({
    
    visNetwork::visIgraph(rvalues$myigraph()) %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        nodesIdSelection = TRUE
      ) %>%
      visHierarchicalLayout(levelSeparation = 100, nodeSpacing = 500,
        shakeTowards = "leaves") %>% # same as visLayout(hierarchical = TRUE) 
      visExport(
        name = "igraph",
        type = "png",
        label = "Save igraph graph"
      )
  })
  
  # Other content --------------------------------------------------------------
  output$de_volcano_signature <- renderPlot({
    signature_volcano(gtl = rvalues$mygtl(),
                      FDR = 0.05)
  })
  
  
  output$gtl_loaded <- renderText({
    describe_gtl(gtl = rvalues$mygtl())
  })
  
  output$team_list <- renderTable({
    team_df
  })
  
  
}

# Launching magnetique! --------------------------------------------------------
shinyApp(magnetique_ui, magnetique_server)

# same for the diff exp things (but they are anyway in the GTL)

# do compute the z score or so for the res_enrich

## and then have some emap interactive/ggs interactive on that? we need then a numericinput for the nr of genesets
