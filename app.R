library("GeneTonic")
library("DESeq2")

MAGNet_DCMvsHCM_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsHCM_GeneTonic.rds")
MAGNet_DCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsNFD_GeneTonic.rds")
MAGNet_HCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_HCMvsNFD_GeneTonic.rds")
MAGNet_DCMvsHCM_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_hcm.RData"))
MAGNet_DCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_nfd.RData"))
MAGNet_HCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_hcm_vs_nfd.RData"))


# DCM vs HCM --------------------------------------------------------------
gtl_DCMvsHCM_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de[order(MAGNet_DCMvsHCM_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de[order(MAGNet_DCMvsHCM_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de[order(MAGNet_DCMvsHCM_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)


# DCM vs NFD --------------------------------------------------------------
gtl_DCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de[order(MAGNet_DCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de[order(MAGNet_DCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de[order(MAGNet_DCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)


# HCM vs NFD --------------------------------------------------------------
gtl_HCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de[order(MAGNet_HCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de[order(MAGNet_HCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de[order(MAGNet_HCMvsNFD_GeneTonic$res_de$padj), ],
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

all_gtls <- list(
  DCMvsHCM = list(
    BP = gtl_DCMvsHCM_BP,
    MF = gtl_DCMvsHCM_MF,
    CC = gtl_DCMvsHCM_CC
  ),
  DCMvsNFD = list(
    BP = gtl_DCMvsNFD_BP,
    MF = gtl_DCMvsNFD_MF,
    CC = gtl_DCMvsNFD_CC
  ),
  HCMvsNFD = list(
    BP = gtl_HCMvsNFD_BP,
    MF = gtl_HCMvsNFD_MF,
    CC = gtl_HCMvsNFD_CC
  )
)

all_igraph <- list(DCMvsHCM = MAGNet_DCMvsHCM_igraph,
                   DCMvsNFD = MAGNet_DCMvsNFD_igraph,
                   HCMvsNFD = MAGNet_HCMvsNFD_igraph)

all_vst <- list(DCMvsHCM = vst(all_gtls[["DCMvsHCM"]][["BP"]]$dds),
                DCMvsNFD = vst(all_gtls[["DCMvsNFD"]][["BP"]]$dds),
                HCMvsNFD = vst(all_gtls[["HCMvsNFD"]][["BP"]]$dds))

rm(MAGNet_DCMvsHCM_GeneTonic,
   MAGNet_DCMvsNFD_GeneTonic,
   MAGNet_HCMvsNFD_GeneTonic)

rm(gtl_DCMvsHCM_BP,
   gtl_DCMvsHCM_MF,
   gtl_DCMvsHCM_CC,
   gtl_DCMvsNFD_BP,
   gtl_DCMvsNFD_MF,
   gtl_DCMvsNFD_CC,
   gtl_HCMvsNFD_BP,
   gtl_HCMvsNFD_MF,
   gtl_HCMvsNFD_CC)

rm(MAGNet_DCMvsHCM_igraph,
   MAGNet_DCMvsNFD_igraph,
   MAGNet_HCMvsNFD_igraph)


# loading libraries -------------------------------------------------------

library("shiny")
library("visNetwork")
library("bs4Dash")
library("shinydashboard")
library("ggplot2")
library("ggrepel")
library("igraph")
library("plotly")
library("shinycssloaders")
library("shinycustomloader")
options(spinner.type = 6)


# sourcing external files -------------------------------------------------

source("volcano_plot.R")
source("DTU/plots_with_se_obj.R")
source("utils.R")
se <- readRDS("MAGNetApp_cloud_data/data/summarized_experiment.RDS")


# ui definition -----------------------------------------------------------

magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  
  header = shinydashboard::dashboardHeader(title = "magnetique"),
  
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
      id = "tabs",
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
        title = "Gene View", icon = icon("heartbeat"), value = "tab-gene-view",
        fluidRow(
          column(
            width = 4,
            DT::dataTableOutput("de_table")
          ),
          column(
            width = 4,
            #withLoader(plotlyOutput("de_volcano"), type="image", loader="/heart.gif") # only works when starting app over RunApp Button (and doesn't look good)
            withSpinner(
              plotlyOutput("de_volcano")
            )
          ),
          column(
            width = 4,
            withSpinner(
              plotlyOutput("dtu_volcano")
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            withSpinner(
              plotOutput("dtu_plot")
            ),
            uiOutput("carnival_launch")
          )
        )
      ),
      shiny::tabPanel(
        id = "tab-geneset-view",
        title = "Geneset View", icon = icon("project-diagram"), value = "tab-geneset-view",
        fluidRow(
          column(
            width = 5,
            DT::dataTableOutput("enrich_table")
          ),
          column(
            width = 7,
            withSpinner(
              plotlyOutput("enriched_funcres")
            )
          )
        ),
        fluidRow(
          column(
            width = 5,
            withSpinner(
              visNetworkOutput("visnet_em")
            ) 
          ),
          column(
            width = 4,
            plotOutput("emap_signature")
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
  rvalues$myvst <- NULL
  
  
  # selector of gtl object --------------------------------------------------
  rvalues$mygtl <- reactive({
    
    message(input$selected_contrast)
    message(input$selected_ontology)
    
    all_gtls[[input$selected_contrast]][[input$selected_ontology]]
  })
  
  rvalues$myigraph <- reactive({
    all_igraph[[input$selected_contrast]]
  })
  
  rvalues$myvst <- reactive({
    all_vst[[input$selected_contrast]]
  })
  
  # DE related content ---------------------------------------------------------
  output$de_table <- DT::renderDataTable({
    mygtl <- rvalues$mygtl()
    myde <- mygtl$res_de
    myde <- GeneTonic::deseqresult2df(myde)
    ensembl_url <- "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g="
    rownames(myde) <- lapply(rownames(myde), function(x) format_url(ensembl_url, x))
    myde <- myde[c("SYMBOL", "log2FoldChange", "padj")]
    DT::datatable(myde, escape = FALSE, options = list(scrollX = TRUE), selection = 'single')  %>% 
      DT::formatRound(columns=c('log2FoldChange', 'padj'), digits=3)
  })
  
  
  output$de_volcano <- renderPlotly({
    mygtl <- rvalues$mygtl()
    myde <- mygtl$res_de
    p <- volcano_plot(myde,
                      mygtl$annotation_obj, 
                      volcano_labels = 0, 
                      plot_title = "Volcano Plot - Differentially expressed genes")
    plotly::ggplotly(p, tooltip = "text") %>% 
    toWebGL() 
  })
  
  output$dtu_volcano <- renderPlotly({
    genes_dtu <- unique(rowData(se)$gene_id)
    mygtl <- rvalues$mygtl()
    my_de <- mygtl$res_de
    filter <- rownames(my_de) %in% genes_dtu
    my_de <- my_de[filter, ]
    p <- volcano_plot(my_de, 
                      mygtl$annotation_obj, 
                      volcano_labels = 0,
                      plot_title = "Volcano Plot - Differential Transcript Usage")
    plotly::ggplotly(p, tooltip = "text") %>%
    toWebGL() 
  })
  
  # enrichment map related content ---------------------------------------------
  output$enrich_table <- DT::renderDataTable({
    mygtl <- rvalues$mygtl()
    myres_enrich <- mygtl$res_enrich
    df <- data.frame(description = myres_enrich$gs_description,
                     obs = myres_enrich$DE_count,
                     exp = myres_enrich$Expected,
                     padj = myres_enrich$gs_pvalue)
    df <- df[order(df$obs, decreasing = T), ]
    rownames(df) <- myres_enrich$gs_id
    colnames(df) <- c("Description", "Observed", "Expected", "padj")
    DT::datatable(df, escape = FALSE, options = list(scrollX = TRUE))  %>% 
      DT::formatRound(columns=c('padj'), digits=3)
  })
  
  emap_graph <- reactive({
    emg <- enrichment_map(
      gtl = rvalues$mygtl(),
      n_gs = input$number_genesets,
      overlap_threshold = 0.1,
      scale_edges_width = 200,
      color_by = input$color_by
    )
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
    
    gs_heatmap(
      se = rvalues$myvst() ,
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
  })
  
  output$enriched_funcres  <- renderPlotly({
    gtl <- rvalues$mygtl()
    ggplotly(enhance_table(gtl$res_enrich,
                           gtl$res_de,
                           annotation_obj = gtl$annotation_obj,
                           n_gs = input$number_genesets,
                           chars_limit = 50
    ))
  })
  
  # DTU related content --------------------------------------------------------
  genes_dtu <- unique(rowData(se)$gene_id)
  
  output$dtu_plot <- renderPlot({
    row <- input$de_table_rows_selected
    validate(need(!is.na(row),
                  message = "Please select an entry from the table."
    ))
    mygtl <- rvalues$mygtl()
    res_de <- mygtl$res_de
    entry <- rownames(res_de)[[row]]
    
    validate(need(entry %in% genes_dtu,
                  message = paste("The gene you selected, ", entry, ", is not a DTU gene. Please select another gene from the table")
    ))
    #gtf_gene <- subset(gtf, type == "gene" & gene_id == entry)
    #plot_dtu(mcols(gtf_gene)[["gene_id"]], se, gtf)
    plot_dtu(entry, se, gtf)
  })
  
  output$carnival_launch <- renderUI({
    tagList(
      actionButton(
        inputId = "btn_show_carnival",
        icon = icon("flask"),
        label = "Show Carnival View of selected DTU gene", style = .actionbutton_biocstyle
      ),
      actionButton(
        inputId = "btn_switch_emap",
        icon = icon("project-diagram"),
        label = "Jump to Enrichtment Map",
        style = .actionbutton_biocstyle
      )
    )
  })
  
  
  # carnival-related content ---------------------------------------------------
  #output$carnival_counts <- renderPlot({
  #mygtl <- rvalues$mygtl()
  
  #g <- rvalues$myigraph()
  #cur_sel <- input$visnet_igraph_selected
  #cur_node <- match(cur_sel, V(g)$name)
  #cur_nodetype <- V(g)$nodetype[cur_node]
  #validate(need(cur_nodetype == "Feature",
  #              message = "" # "Please select a gene/feature."
  #))
  
  #validate(need(cur_sel != "",
  #             message = "Please select a node from the graph to plot the expression values."
  #))
  
  #cur_geneid <- mygtl$annotation_obj$gene_id[match(cur_sel, mygtl$annotation_obj$gene_name)]
  
  #message(length(cur_sel))
  #message(cur_geneid)
  
  #genes_exp <- rownames(mygtl$dds)
  #validate(need(cur_geneid %in% genes_exp,
  #              message = "gene not found in expression matrix" 
  #))
  
  #gene_plot(gtl = mygtl, gene = cur_geneid, 
  #           intgroup = "Etiology")
  #})
  
  
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
  
  observeEvent(input$btn_show_carnival, {
    showModal(
      modalDialog(
        title = "Carnival View", size = "l", fade = TRUE,
        footer = NULL, easyClose = TRUE,
        visNetworkOutput("visnet_igraph"),
        #plotOutput("carnival_counts")
      )
    )
  })
  
  observeEvent(input$btn_switch_emap,{
    updateTabsetPanel(session, "tabs",selected = "tab-geneset-view")
  })
  
  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
  
  
}

# Launching magnetique! --------------------------------------------------------
shinyApp(magnetique_ui, magnetique_server)

