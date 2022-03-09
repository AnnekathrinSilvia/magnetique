# loading libraries -------------------------------------------------------
library(shiny)
library(future)
library(promises)
plan(multisession)

library(shinycssloaders)
library(plotly, warn.conflicts = FALSE)
library(reactable, warn.conflicts = FALSE)
library(bs4Dash, warn.conflicts = FALSE)
library(shinydashboard, warn.conflicts = FALSE)
library(visNetwork, warn.conflicts = FALSE)
library(rintrojs, warn.conflicts = FALSE)
library(shinyBS, warn.conflicts = FALSE)

options(spinner.type = 6)

# sourcing external files -------------------------------------------------
source("utils.R")

# ui definition -----------------------------------------------------------
magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  header = shinydashboard::dashboardHeader(title = "magnetique"),

  # sidebar definition ------------------------------------------------------
  sidebar = shinydashboard::dashboardSidebar(
    title = "Options",
    uiOutput("ui_sidebar")
  ),

  # body definition ---------------------------------------------------------
  body = shinydashboard::dashboardBody(
    introjsUI(),
    shiny::tags$script(
      HTML(
        "$(function(){
        $(document).keyup(function(e) {
        if (e.which == 17) {
          $('#bookmarker').click()
        }
        });
        })"
      )
    ),
    tabBox(
      width = 12,
      id = "magnetique_tab",
      shiny::tabPanel(
        title = "Welcome!", icon = icon("magnet"), value = "tab-welcome",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_firststeps",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_firststeps",
                "Click me to start a tour of this section!",
                "bottom",
                options = list(container = "body")
              ),
              style = "float:right"
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            includeMarkdown("data/overview.md")
          )
        ),
      ),
      shiny::tabPanel(
        title = "Gene View", icon = icon("heartbeat"), value = "tab-gene-view",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_geneview",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_geneview",
                "Click me to start a tour of this section!",
                "bottom",
                options = list(container = "body")
              ),
              style = "float:right"
            )
          )
        ),
        fluidRow(
          id = "geneview_row1",
          column(
            width = 6,
            withSpinner(
              reactableOutput("de_table")
            )
          ),
          column(
            width = 3,
            withSpinner(
              plotlyOutput("de_volcano")
            )
          ),
          column(
            width = 3,
            withSpinner(
              plotlyOutput("dtu_volcano")
            )
          )
        ),
        fluidRow(
          id = "geneview_row2",
          column(
            width = 4,
            withSpinner(
              plotlyOutput("gene_counts")
            )
          ),
          column(
            width = 4,
            withSpinner(
              plotOutput("transcript_proportion")
            )
          ),
          column(
            width = 4,
            withSpinner(
              plotOutput("gene_structure")
            )
          )
        ),
        fluidRow(
          column(
            width = 1,
          ),
          column(
            width = 10,
            withSpinner(
              plotlyOutput("wgcn_heatmap")
            )
          ),
          column(
            width = 1,
          ),
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("carnival_launch")
          )
        )
      ),
      shiny::tabPanel(
        id = "tab-geneset-view",
        title = "Geneset View", icon = icon("project-diagram"), value = "tab-geneset-view",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_genesetview",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_genesetview",
                "Click me to start a tour of this section!",
                "bottom",
                options = list(container = "body")
              ),
              style = "float:right"
            )
          )
        ),
        fluidRow(
          column(
            width = 5,
            withSpinner(
              reactableOutput("enrich_table")
            )
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
            width = 6,
            withSpinner(
              visNetworkOutput("visnet_em")
            )
          ),
          column(
            width = 6,
            withSpinner(
              plotOutput("emap_signature")
            )
          )
        )
      ),
      shiny::tabPanel(
        title = "Bookmarks", icon = icon("bookmark"), value = "tab-bookmark",
        fluidRow(
          column(
            width = 12,
            div(

              actionButton(
                "tour_bookmarks",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_bookmarks",
                "Click me to start a tour of this section!",
                "bottom",
                options = list(container = "body")
              ),
              style = "float:right"
            )
          )
        ),
        uiOutput("ui_bookmarks")
      ),
      shiny::tabPanel(
        title = "About us", icon = icon("users"), value = "tab-aboutus",
        fluidRow(
          column(
            width = 12,
            h2("Project members (alphabetical order)"),
            withSpinner(
              tableOutput("team_list")
            )
          )
        )
      )
    )
  )
)

# server definition -------------------------------------------------------
magnetique_server <- function(input, output, session) {
  showNotification("Connecting to the db.", id = "db_connect", duration=NULL)  

  # con <- DBI::dbConnect(str_replace
  #   RPostgres::Postgres(),
  #   dbname = "magnetique",
  #   host = "10.250.140.12",
  #   port = 5432,
  #   password = "wGpVDExWK2NppuWENFcjc9v3VKgL4h86ZBHF78pEFdqJwEQwfG",
  #   user = "magnetique_reader"
  # )
  con <- DBI::dbConnect(RSQLite::SQLite(), "magnetique.sqlite")
  
  session$onSessionEnded(function() DBI::dbDisconnect(con))
  removeNotification(id = "db_connect")

  showNotification("Loading magnet dataset", id = "dds_loading")
  
  ## Info on how to create dds_magnet (to be adapted and updated in the final version):
  # DCMvsHCM <- readRDS("~/Development/magnetique/MAGNetApp/DCMvsHCM.RDS")
  # dds_magnet <- DCMvsHCM$dds
  # saveRDS(dds_magnet, "dds_magnet.RDS")
  dds_magnet <- readRDS("dds_magnet.RDS")
  showNotification("Loaded magnet dataset!", type = "message")
  removeNotification(id = "dds_loading")
  
  showNotification("Loading libraries.", id = "lib_load", duration=NULL)  
  suppressPackageStartupMessages({
    library(dplyr, warn.conflicts = FALSE)
    library(tidyr, warn.conflicts = FALSE)
    library(ggplot2, warn.conflicts = FALSE)
  })
  removeNotification(id = "lib_load")

  # reactive objects and setup commands -------------------------------------
  rvalues <- reactiveValues()
  rvalues$mygtl <- NULL
  rvalues$key <- NULL
  rvalues$myvst <- NULL

  rvalues$mygenes <- c()
  rvalues$mygenesets <- c()


  # sidebar server-side -----------------------------------------------------
  output$ui_sidebar <- renderUI({
    tagList(
      sidebarMenu(
        selectInput("selected_contrast",
          label = "Contrast id",
          choices = c(
            "DCMvsHCM",
            "DCMvsNFD",
            "HCMvsNFD"
          ),
          selected = "DCMvsHCM"
        ),
        selectInput("selected_ontology",
          label = "Ontology",
          choices = c("BP", "MF", "CC"),
          selected = "BP"
        ),
        numericInput("number_genesets",
          "Number of genesets",
          value = 15,
          min = 0
        ),
        selectInput("color_by",
          "Color by",
          choices = c(
            "z_score",
            "gs_pvalue"
          ),
          selected = "z_score"
        ),
        actionButton("bookmarker",
          label = "Bookmark", icon = icon("heart"),
          style = "color: #ffffff; background-color: #ac0000; border-color: #ffffff"
        )
      )
    )
  })
  # selector trigger data loading
  rvalues$mygtl <- reactive({
    message(input$selected_contrast)
    message(input$selected_ontology)

    showNotification("Assembling the gtl object for you!",
      id = "info_assembling",
      duration = NULL
    )
    mygtl <- buildup_gtl(
      con,
      dds = dds_magnet,                   
      contrast = input$selected_contrast,
      ontology = input$selected_ontology
    )
    removeNotification(id = "info_assembling")
    showNotification("Done! gtl object ready!", type = "message")
    return(mygtl)
  })

  rvalues$myvst <- reactive({
    DESeq2::vst(rvalues$mygtl()$dds)
  })

  rvalues$key <- reactive({
    tbl(con, paste0("res_", input$selected_contrast)) %>%
      select(
        c(
          "gene_id",
          "SYMBOL",
          "padj",
          "log2FoldChange",
          "dtu_pvadj",
          "dtu_dif",
          "module",
          "rank"
        )
      ) %>%
      mutate_at(vars(padj, log2FoldChange, dtu_pvadj, dtu_dif), ~round(., 6)) %>%
      collect() %>%
      highlight_key(.)
  })

  # DE related content ---------------------------------------------------------
  output$de_table <- renderReactable({
    rvalues$key() %>%
      reactable(
        .,
        searchable = TRUE,
        striped = TRUE,
        defaultPageSize = 5,
        highlight = TRUE,
        selection = "single",
        onClick = "select",
        rowStyle = list(cursor = "pointer"),
        theme = reactableTheme(
          stripedColor = "#f6f8fa",
          highlightColor = "#f0f5f9",
          cellPadding = "8px 12px",
        ),
        defaultColDef = colDef(sortNALast = TRUE),
        list(
          gene_id = colDef(
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 140,
            header = with_tooltip("gene_id", "Link to Ensembl gene page")
          ),
          log2FoldChange = colDef(
            header = with_tooltip("log2FoldChange", "DESeq2 log2FoldChange")
          ),
          padj = colDef(
            header = with_tooltip("padj", "DESeq2 padj")
          ),
          dtu_pvadj = colDef(
            header = with_tooltip("dtu_pvadj", "DRIMseq minimum p-value")
          ),
          dtu_dif = colDef(
            header = with_tooltip("dtu_dif", "Differential isoform usage")
          ),
          module = colDef(
            header = with_tooltip("module", "Co-expressed genes module")
          ),
          rank = colDef(
            header = with_tooltip("rank", "Rank in module")
          )
        )
      )
  })

  output$de_volcano <- renderPlotly({
    rvalues$key() %>%
      plot_ly(., color = I("black"), showlegend = FALSE) %>%
      add_markers(
        x = ~log2FoldChange,
        y = ~ -log10(padj),
        type = "scatter",
        text = ~ paste0(
          "<b>", SYMBOL, "</b>",
          "<br><i>GeneID</i>: ", gene_id,
          "<br><i>Log2FC</i> = ", format(round(log2FoldChange, 2)),
          "<br><i>p-value (adjusted)</i> = ", format(round(padj, 2))
        ),
        hoverinfo = "text"
      ) %>%
      config(displayModeBar = FALSE) %>%
      toWebGL %>%
      layout(title = "Differentially expressed genes") %>%
      highlight(
        on = "plotly_click",
        off = "plotly_doubleclick",
        color = "red"
      )
  })

  output$dtu_volcano <- renderPlotly({
    rvalues$key() %>%
      plot_ly(., color = I("black"), showlegend = FALSE) %>%
      add_markers(
        x = ~dtu_dif,
        y = ~ -log10(dtu_pvadj),
        type = "scatter",
        text = ~ paste0(
          "<b>", SYMBOL, "</b>",
          "<br><i>GeneID</i>: ", gene_id,
          "<br><i>dif </i> = ", format(round(dtu_dif, 2)),
          "<br><i>p-value (adjusted)</i> = ", format(round(dtu_pvadj, 2))
        ),
        hoverinfo = "text"
      ) %>%
      config(displayModeBar = FALSE) %>%
      toWebGL %>%
      layout(title = "Genes with differential transcript usage") %>%
      highlight(
        on = "plotly_click",
        off = "plotly_doubleclick",
        color = "red"
      )
  })

  output$gene_counts <- renderPlotly({
    i <- getReactableState("de_table", "selected")
    validate(need(!is.na(i),
      message = "Please select an entry from the table."
    ))
    i <- rvalues$key()$data()[[i, "gene_id"]]

    counts <- con %>%
      tbl(paste0("counts_", input$selected_contrast)) %>%
      filter(row_names == !!i) %>%
      select(-row_names) %>%
      mutate(n = n()) %>%
      collect() %>%
      pivot_longer(-n)

    metadata <- con %>%
      tbl("metadata") %>%
      select(row_names, Etiology) %>%
      collect()

    left_join(counts, metadata, by = c("name" = "row_names")) %>%
      plot_ly(
        .,
        type = "box",
        x = ~Etiology,
        y = ~ log10(value),
        color = ~factor(Etiology, levels = c("NFD",  "DCM", "HCM")),
        colors = c(I("steelblue"), I("gold"), I("forestgreen"))
      ) %>%
      config(displayModeBar = FALSE) %>%
      layout(title = "Gene counts", xaxis = list(title = ""))
  })

  output$gene_structure <- renderPlot({
    i <- getReactableState("de_table", "selected")
    validate(
      need(!is.na(i),
        message = "Please select an entry from the table."
      )
    )
    dtu_tested <- rvalues$key()$data()[[i, "dtu_pvadj"]]
    validate(
      need(!is.na(dtu_tested),
        message = "Entry not tested for DTU."
      )
    )
    i <- rvalues$key()$data()[[i, "gene_id"]]
    gtf <- con %>%
      tbl("gtf") %>%
      filter(gene_id == !!i) %>%
      collect() %>%
      GenomicRanges::GRanges(.)
    gtf <- GenomicRanges::split(gtf, gtf$transcript_id)
    plot_gene_structure(gtf) + labs(title='Gene structure') + theme(plot.title = element_text(size=18))
  })

  output$transcript_proportion <- renderPlot({
    i <- getReactableState("de_table", "selected")
    validate(
      need(!is.na(i),
        message = "Please select an entry from the table."
      )
    )
    dtu_tested <- rvalues$key()$data()[[i, "dtu_pvadj"]]
    validate(
      need(!is.na(dtu_tested),
        message = "Entry not tested for DTU."
      )
    )
    i <- rvalues$key()$data()[[i, "gene_id"]]
    x <- con %>%
      tbl("gene2tx") %>%
      filter(gene_id == !!i) %>%
      left_join(
        tbl(con, "dtu_fit_proportions"),
        by = (c("transcript_id" = "row_names"))
      ) %>%
      collect()
    x <- x %>% pivot_longer(
      -c(gene_id, transcript_id),
      names_to = "Run",
      values_to = "proportion"
    )
    metadata <- con %>%
      tbl("metadata") %>%
      select(Run, Etiology) %>%
      collect()

    left_join(x, metadata, by = c("Run" = "Run")) %>%
      ggplot() +
      geom_jitter(aes(x = transcript_id, y = proportion, color = Etiology),
        position = position_jitterdodge( jitter.width = 0.2 ),
        alpha = 0.9, size = 2, show.legend = T, na.rm = TRUE
      ) +
      geom_boxplot(aes(x = transcript_id, y = proportion, fill = Etiology),
        outlier.size = 0, alpha = 0.4, lwd = 0.5, show.legend = F
      ) +
      scale_fill_manual(name = "Etiology", values = group_colors) +
      scale_colour_manual(name = "Etiology", values = group_colors) +
      coord_flip() +
      labs(y = "Transcript proportion", title = 'Transcript usage') +
      theme_minimal(20) +
      theme(
        axis.title.y = element_blank(),
        plot.title = element_text(size=18)
      )
  })


  # wgcn related content ---------------------------------------------
  output$wgcn_heatmap <- renderPlotly({
    y <- con %>%
      tbl("wgcn_hp") %>%
      data.frame() %>%
      tibble::column_to_rownames("row_names") %>%
      as.matrix()

    textMatrix <- ifelse(y < 0.01, signif(y, 1), "")
    x <- con %>%
      tbl("wgcn_hcor") %>%
      data.frame() %>%
      tibble::column_to_rownames("row_names") %>%
      as.matrix()
    a <- rep(0:(ncol(x) - 1), each = nrow(x))
    b <- rep(c(0:(nrow(x) - 1)), ncol(x))
    plot_ly(
      x = colnames(x),
      y = rownames(x),
      z = x,
      colors = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(256),
      zmin = -1, zmax = 1,
      type = "heatmap"
    ) %>%
      add_annotations(x = a, y = b, text = textMatrix, xref = "x", yref = "y", showarrow = FALSE, font = list(color = "black")) %>%
      config(displayModeBar = FALSE) %>%
      layout(title = "Module-trait correlation")
  })

  res_enrich <- reactive({
    tbl(con, paste0("res_enrich_", local(input$selected_contrast))) %>%
        filter(ontology == local(input$selected_ontology)) %>%
        rename(
          c(
            "id" = "gs_id",
            "description" = "gs_description", 
            "expected" = "Expected", 
            "observed"= "gs_de_count", 
            "pval" = "gs_pvalue"))  %>%
        select(id, description, expected, observed, pval) %>%
        arrange(pval) %>%
        collect()
  })


  # enrichment map related content ---------------------------------------------
  output$enrich_table <- renderReactable({
      reactable(
        res_enrich(),
        searchable = TRUE,
        striped = TRUE,
        defaultPageSize = 5,
        highlight = TRUE,
        selection = "single",
        onClick = "select",
        rowStyle = list(cursor = "pointer"),
        theme = reactableTheme(
          stripedColor = "#f6f8fa",
          highlightColor = "#f0f5f9",
          cellPadding = "8px 12px",
        ))
  })

  emap_graph <- reactive({
    mygtl <- rvalues$mygtl()
    emg <- GeneTonic::enrichment_map(
      gtl = mygtl,
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
    mygtl <- rvalues$mygtl()
    myvst <- rvalues$myvst()

    cur_gsid_em <- mygtl$res_enrich$gs_id[
      match(input$visnet_em_selected, mygtl$res_enrich$gs_description)
    ]
    
    # cur_gsid_tbl <- 
      
      
    validate(
      need(!is.na(cur_gsid_em),
        message = "Please select a gene set from the Enrichment Map."
      )
    )

    colnames(myvst) <- NULL
    GeneTonic::gs_heatmap(
      se = myvst,
      gtl = mygtl,
      geneset_id = cur_gsid_em,
      FDR = 0.05,
      de_only = FALSE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      center_mean = TRUE,
      scale_row = TRUE,
      anno_col_info = "Etiology"
    )
  })

  output$enriched_funcres <- renderPlotly({

    res_enrich <- tbl(con, paste0("res_enrich_", local(input$selected_contrast))) %>%
      filter(ontology == local(input$selected_ontology)) %>%
      collect()
    res_enrich <- as.data.frame(res_enrich)
    rownames(res_enrich) <- res_enrich$gs_id

    res_de <- tbl(con, paste0("res_", local(input$selected_contrast))) %>%
      collect()
    res_de <- res_de %>% 
      filter(!is.na(SYMBOL)) %>% 
      arrange(-desc(padj)) %>% 
      select(gene_id, log2FoldChange, padj, SYMBOL)
      
    res_de <- DESeq2::DESeqResults(
      S4Vectors::DataFrame(res_de))

    rownames(res_de) <- res_de$gene_id
    res_de$pvalue <- res_de$padj
    res_de$description <- ''

    annotation_obj <- tbl(con, "annotation_obj") %>%
      collect() %>%
      as.data.frame()
    rownames(annotation_obj) <- annotation_obj$gene_id 

    ggplotly(
      GeneTonic::enhance_table(
        res_enrich,
        res_de,
        annotation_obj = annotation_obj,
        n_gs = input$number_genesets,
        chars_limit = 50
      ) 
    )
  })

  # DTU related content --------------------------------------------------------
  output$carnival_launch <- renderUI({
    tagList(
      actionButton(
        inputId = "btn_show_carnival",
        icon = icon("flask"),
        label = "Show Carnival View", style = .actionbutton_biocstyle
      )
    )
  })


  output$visnet_igraph <- renderVisNetwork({
    con %>%
      tbl("carnival") %>%
      filter(contrast == local(input$selected_contrast)) %>%
      pull(igraph) %>%
      jsonlite::unserializeJSON() %>%
      visNetwork::visIgraph() %>%
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

  # Bookmarker -----------------------------------------------------------------
  output$ui_bookmarks <- renderUI({
    tagList(
      fluidRow(
        column(
          width = 6,
          h5("Bookmarked genes"),
          reactableOutput("bookmarks_genes")
        ),
        column(
          width = 6,
          h5("Bookmarked genesets"),
          reactableOutput("bookmarks_genesets")
        )
      )
    )
  })

  output$bookmarks_genes <- renderReactable({
    validate(
      need(
        length(rvalues$mygenes) > 0,
        "Please select at least one gene with the Bookmark button"
      )
    )

    annotation_obj <- tbl(con, "annotation_obj") %>%
      filter(gene_id %in% local(rvalues$mygenes)) %>%
      select(gene_id, gene_name) %>%
      collect() 

    reactable(annotation_obj)
  })

  output$bookmarks_genesets <- renderReactable({
    validate(
      need(
        length(rvalues$mygenesets) > 0,
        "Please select at least one geneset with the Bookmark button"
      )
    )
    res_enrich() %>%
      filter(id %in% rvalues$mygenesets) %>%
      select(id, description) %>%
      reactable(rownames = FALSE)
  })

  # observeEvent(input$selected_contrast, {
  #   message(input$selected_contrast)
  #   if (!is.null(rvalues$key)) {
  #     message(dim(rvalues$key()$data))
  #     # rvalues$key()$data <- NULL       
  #   }
  # })

  observeEvent(input$bookmarker, {
    if (input$magnetique_tab == "tab-welcome") {
      showNotification("Welcome to magnetique! Navigate to the main tabs of the application to use the Bookmarks functionality.")
    } else if (input$magnetique_tab == "tab-gene-view") {
      i <- getReactableState("de_table", "selected")
      if (!is.null(i)) {
        sel_gene <- rvalues$key()$data()[[i, "gene_id"]]
        if (!sel_gene %in% rvalues$mygenes) {
          rvalues$mygenes <- c(rvalues$mygenes, sel_gene)
          showNotification(
            sprintf(
              "The selected gene %s was added to the bookmarked genes.",
              sel_gene
            ),
            type = "default"
          )
        }
      }
    } else if (input$magnetique_tab == "tab-geneset-view") {
      i <- getReactableState("enrich_table", "selected")    
      if (!is.null(i)) {
        sel_gs <- res_enrich()[[i, "id"]]
        if (!sel_gs %in% rvalues$mygenesets) {
          rvalues$mygenesets <- c(rvalues$mygenesets, sel_gs)
          showNotification(
            sprintf(
              "The selected geneset %s was added to the bookmarked genesets.",
              sel_gs
            ),
            type = "default"
          )
        }
      }
    }
  })

  # Other content --------------------------------------------------------------
  output$team_list <- renderTable({
    make_team_df
  })

  # Tours observers
  observeEvent(input$tour_firststeps, {
    tour <- read.delim("tours/intro_firststeps.txt",
      sep = ";", stringsAsFactors = FALSE, row.names = NULL, quote = ""
    )
    introjs(session, options = list(steps = tour))
  })

  observeEvent(input$tour_geneview, {
    tour <- read.delim("tours/intro_geneview.txt",
      sep = ";", stringsAsFactors = FALSE, row.names = NULL, quote = ""
    )
    introjs(session, options = list(steps = tour))
  })

  observeEvent(input$tour_genesetview, {
    tour <- read.delim("tours/intro_genesetview.txt",
      sep = ";", stringsAsFactors = FALSE, row.names = NULL, quote = ""
    )
    introjs(session, options = list(steps = tour))
  })

  observeEvent(input$tour_bookmarks, {
    tour <- read.delim("tours/intro_bookmarks.txt",
      sep = ";", stringsAsFactors = FALSE, row.names = NULL, quote = ""
    )
    introjs(session, options = list(steps = tour))
  })

  observeEvent(input$btn_show_carnival, {
    showModal(
      modalDialog(
        title = "Carnival View", size = "l", fade = TRUE,
        footer = NULL, easyClose = TRUE,
        visNetworkOutput("visnet_igraph")
      )
    )
  })

  observeEvent(input$btn_switch_emap, {
    updateTabsetPanel(session, "tabs", selected = "tab-geneset-view")
  })

  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
}
shinyApp(magnetique_ui, magnetique_server)
