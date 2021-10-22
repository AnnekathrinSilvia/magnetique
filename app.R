library("GeneTonic")
library("DESeq2")

MAGNet_DCMvsHCM_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsHCM_GeneTonic.rds")
MAGNet_DCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsNFD_GeneTonic.rds")
MAGNet_HCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_HCMvsNFD_GeneTonic.rds")
MAGNet_DCMvsHCM_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_hcm.RData"))
MAGNet_DCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_nfd.RData"))
MAGNet_HCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_hcm_vs_nfd.RData"))

library(shinycssloaders)
options(spinner.type = 6)


# DCM vs HCM --------------------------------------------------------------
gtl_DCMvsHCM_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)


# DCM vs NFD --------------------------------------------------------------
gtl_DCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)


# HCM vs NFD --------------------------------------------------------------
gtl_HCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
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

library("shiny")
library("visNetwork")

ui <- fluidPage(
  fluidRow(
    column(
      width = 8,
      selectInput("selected_contrast", label = "Contrast id", choices = c("DCMvsHCM", "DCMvsNFD", "HCMvsNFD")),
      selectInput("selected_ontology", label = "Ontology", choices = c("BP", "MF", "CC")),
      
      actionButton("button_loadgtl", "load gtl data"),
      
      plotOutput("de_volcano"),
      
      verbatimTextOutput("gtl_loaded")
    ),
    column(
      width = 4,
      numericInput("number_genesets", "Number of genesets", value = 15, min = 0),
      selectInput("color_by", "Color by", choices = c("z_score", "gs_pvalue"), 
                  selected = "z_score")
    ),
    
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
        
    ),
    fluidRow(
      column(
        width = 8,
        withSpinner(
          visNetworkOutput("visnet_igraph")
        )
      )
    )
    
  )
  
)

server <- function(input, output, session) {
  
  rvalues <- reactiveValues()
  rvalues$mygtl <- NULL
  rvalues$myigraph <- NULL
  
  # observeEvent(input$button_loadgtl, {
    
    rvalues$mygtl <- reactive({
      
      message(input$selected_contrast)
      message(input$selected_ontology)
      
      all_gtls[[input$selected_contrast]][[input$selected_ontology]]
      
      # if(input$selected_contrast == "DCMvsHCM") {
      #   if(input$selected_ontology == "BP") {
      #     gtl_DCMvsHCM_BP
      #   } else if(input$selected_ontology == "MF") {
      #     gtl_DCMvsHCM_MF
      #   } else if(input$selected_ontology == "CC") {
      #     gtl_DCMvsHCM_CC
      #   }
      # }
    })
  # })
    
    rvalues$myigraph <- reactive({
      all_igraph[[input$selected_contrast]]
    })
  
  # selected_gtl <- reactive({
  #   message(input$selected_contrast)
  #   message(input$selected_ontology)
  #   
  #   mygtl <- 
  #     {
  #       if(input$selected_contrast == "DCMvsHCM") {
  #         if(input$selected_ontology == "BP") {
  #           gtl_DCMvsHCM_BP
  #         } else if(input$selected_ontology == "MF") {
  #           gtl_DCMvsHCM_MF
  #         } else if(input$selected_ontology == "CC") {
  #           gtl_DCMvsHCM_CC
  #         }
  #       }
  #       
  #       if(input$selected_contrast == "DCMvsNFD") {
  #         if(input$selected_ontology == "BP") {
  #           gtl_DCMvsNFD_BP
  #         } else if(input$selected_ontology == "MF") {
  #           gtl_DCMvsNFD_MF
  #         } else if(input$selected_ontology == "CC") {
  #           gtl_DCMvsNFD_CC
  #         }
  #       }
  #       
  #       if(input$selected_contrast == "HCMvsNFD") {
  #         if(input$selected_ontology == "BP") {
  #           gtl_HCMvsNFD_BP
  #         } else if(input$selected_ontology == "MF") {
  #           gtl_HCMvsNFD_MF
  #         } else if(input$selected_ontology == "CC") {
  #           gtl_HCMvsNFD_CC
  #         }
  #       }
  #     }
  #   return(mygtl)
  # })
  
  output$de_volcano <- renderPlot({
    signature_volcano(gtl = rvalues$mygtl(),
                      FDR = 0.05)
  })
  
  
  output$gtl_loaded <- renderText({
    describe_gtl(gtl = rvalues$mygtl())
  })
  
  # emap section ------------------------------------------------------------
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
  
  
  # ggs graph section -------------------------------------------------------
  # myggs_graph <- reactive({
  #   g <- ggs_graph(
  #     gtl = rvalues$mygtl(),
  #     n_gs = input$number_genesets,
  #     prettify = TRUE,
  #     geneset_graph_color = "gold"
  #   )
  #   # rank_gs <- rank(V(g)$name[V(g)$nodetype == "GeneSet"])
  #   # rank_feats <- rank(V(g)$name[V(g)$nodetype == "Feature"]) +
  #   #   length(rank_gs) # to keep the GeneSets first
  #   # g <- permute.vertices(g, c(rank_gs, rank_feats))
  #   # return(g)
  # })
  
  output$visnet_ggs <- renderVisNetwork({
    
    visNetwork::visIgraph(myggs_graph()) %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        nodesIdSelection = TRUE
      ) %>%
      visExport(
        name = "ggs_network",
        type = "png",
        label = "Save ggs graph"
      )
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
      visExport(
        name = "igraph",
        type = "png",
        label = "Save igraph graph"
      )
  })
}

shinyApp(ui, server)

# same for the diff exp things (but they are anyway in the GTL)

# do compute the z score or so for the res_enrich

## and then have some emap interactive/ggs interactive on that? we need then a numericinput for the nr of genesets


















