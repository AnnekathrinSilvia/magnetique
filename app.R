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
library(igraph)

options(spinner.type = 6)

# sourcing external files -------------------------------------------------
source("utils.R")
logo_url <- "https://gist.githubusercontent.com/tbrittoborges/3c86ffbaa62e671771f443c65cb04fdc/raw/7ae0ea4a76e8f5464139ef34164c67de7a297ce8/baltica_logo.png"
# ui definition -----------------------------------------------------------
magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  header = shinydashboard::dashboardHeader(disable = TRUE),

  # sidebar definition ------------------------------------------------------
  sidebar = dashboardSidebar(
    img(src = logo_url, class="img-responsive"),
    h1("Magnetique", align='center'),
    h2("Options:", style = 'margin: 15px'),
    uiOutput("ui_sidebar")
  ),

  # body definition ---------------------------------------------------------
  body = dashboardBody(
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
      tabPanel(
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
      tabPanel(
        title = "Gene View", icon = icon("dna"), value = "tab-gene-view",
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
      ),
      tabPanel(
        id = "tab-geneset-view",
        title = "Gene set View", icon = icon("project-diagram"), value = "tab-geneset-view",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_view",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_view",
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
      tabPanel(
        id = "tab-carnival",
        title = "Carnival View", icon = icon("viruses"), value = "tab-carnival",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_carnivalview",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_carnivalview",
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
                visNetworkOutput("visnet_igraph", height = "500px")
            )
        )
      ),
      tabPanel(
        id = "tab-wcgn",
        title = "Network correlation View", icon = icon("network-wired"), value = "tab-wcgn",
        fluidRow(
          column(
            width = 12,
            div(
              actionButton(
                "tour_wcgnview",
                label = "", icon = icon("question-circle"),
                style = .helpbutton_biocstyle
              ),
              shinyBS::bsTooltip(
                "tour_wcgnview",
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
                plotlyOutput("wgcn_heatmap")
            )
        )
      ),
      tabPanel(
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
      tabPanel(
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

  con <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = "magnetique",
    host = Sys.getenv("PGHOST"),
    port = Sys.getenv("PGPORT"),
    password = Sys.getenv("PGPASSWORD"),
    user = Sys.getenv("PGUSER")
  ) 
  session$onSessionEnded(function() DBI::dbDisconnect(con))
  removeNotification(id = "db_connect")

  showNotification("Loading libraries.", id = "lib_load", duration=NULL)  
  suppressPackageStartupMessages({
    library(dplyr, warn.conflicts = FALSE)
    library(tidyr, warn.conflicts = FALSE)
    library(ggplot2, warn.conflicts = FALSE)
  })
  removeNotification(id = "lib_load")

  # reactive objects and setup commands -------------------------------------
  rvalues <- reactiveValues()
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
          "Number of gene sets",
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
      mutate(dtu_dif = case_when(
        !is.na(dtu_pvadj) & is.na(dtu_dif) ~ 0,
        TRUE ~ dtu_dif)) %>%
      collect()
  })
  
  rvalues$metadata <- reactive({
    con %>%
      tbl("metadata") %>%
      collect()
  })

  rvalues$counts <- reactive({
    con %>%
      tbl(paste0("counts_", local(input$selected_contrast))) %>% 
      collect() 
  })

  rvalues$res_enrich <- reactive({
    tbl(con, paste0("res_enrich_", local(input$selected_contrast))) %>%
      filter(ontology == local(input$selected_ontology)) %>%
      collect() %>% 
      as.data.frame(.) %>%
      `rownames<-`(.$gs_id)
    })

  rvalues$annotation_obj <- reactive({
     tbl(con, "annotation_obj") %>%
      collect() %>%
      as.data.frame() %>% 
      `rownames<-`(.$gene_id)
  })

  # DE related content ---------------------------------------------------------
  output$de_table <- renderReactable({
    rvalues$key() %>%
      mutate_at(
        vars(padj, log2FoldChange, dtu_pvadj, dtu_dif), ~round(., 2)) %>%
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
          rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #FF0000")
        ),
        columns = list(
          gene_id = colDef(
            name = "gene_id",
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 140,
            header = with_tooltip("gene_id", "Link to Ensembl gene page")
          ),
          SYMBOL = colDef(
            name = "Gene name",
            header = with_tooltip("gene_name", "Gene symbol")
          ),
          log2FoldChange = colDef(
            name = "dge_log2fc",
            header = with_tooltip("dge_log2fc", "Fold change for the DESeq2 analysis") 
          ),
          padj = colDef(
            name = "dge_padj",
            header = with_tooltip("dge_padj", "Adjusted p-value for the DESeq2 analysis")
          ),
          dtu_pvadj = colDef(
            name = "dtu_padj",
            header = with_tooltip("dtu_padj", "Adjusted p-value for the DRIMseq analysis")
          ),
          dtu_dif = colDef(
            name = "dtu_dif",
            header = with_tooltip("dtu_dif", "Difference in isoform usage for the DRIMseq analysis")
          ),
          module = colDef(
            name = "Module",
            header = with_tooltip("net_module", "Co-expressed genes module for the network analysis")
          ),
          rank = colDef(
            name = "Rank",
            header = with_tooltip("net_rank", "Rank in module for the network analysis")
          )
        ),
        defaultColDef = colDef(sortNALast = TRUE)
      )
  })

  prev_selected <- reactive(getReactableState("de_table", "selected"))
  colors <- reactive(highlight_selected(prev_selected(), nrow(rvalues$key())))
  colors_dtu <- reactive({
    df <- rvalues$key()
    colors()[!is.na(df$dtu_pvadj)]
  })

  output$de_volcano <- renderPlotly({
    rvalues$key() %>%
      plot_ly(., color=I('black'), showlegend = FALSE) %>%
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
      layout(
        title = "Differentially expressed genes",
        yaxis = list(title = '-log10(dge_padj)'), 
        xaxis = list(title = 'dge_log2fc')
      ) 
  })


  output$dtu_volcano <- renderPlotly({  
    rvalues$key() %>%
      plot_ly(., color=I('black'), showlegend = FALSE) %>%
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
      layout(
        title = "Genes with differential transcript usage",
        yaxis = list(title = '-log10(dtu_padj)'), 
        xaxis = list(title = 'dtu_dif')
      )
  })

  observeEvent(colors(), {

    plotlyProxy("de_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list(colors()), 0)
    
    plotlyProxy("dtu_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list(colors()), 0)        
        
  })

  output$gene_counts <- renderPlotly({
    i <- getReactableState("de_table", "selected")
    validate(need(!is.na(i),
      message = "Please select an entry from the table."
    ))
    i <- rvalues$key()[[i, "gene_id"]]

    counts <- rvalues$counts() %>%
      filter(row_names == !!i) %>%
      select(-row_names) %>%
      mutate(n = n()) %>%
      collect() %>%
      pivot_longer(-n)

    metadata <- rvalues$metadata() %>%
      select(row_names, Etiology)

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
    dtu_tested <- rvalues$key()[[i, "dtu_pvadj"]]
    validate(
      need(!is.na(dtu_tested),
        message = "Entry not tested for DTU."
      )
    )
    i <- rvalues$key()[[i, "gene_id"]]
    gtf <- con %>%
      tbl("gtf") %>%
      filter(gene_id == !!i) %>%
      collect() %>%
      GenomicRanges::GRanges(.)
    gtf <- GenomicRanges::split(gtf, gtf$transcript_id)
    suppressMessages({
      plot_gene_structure(gtf) + labs(title='Gene structure') + theme(plot.title = element_text(size=18))
    })
    
  })

  output$transcript_proportion <- renderPlot({
    i <- getReactableState("de_table", "selected")
    validate(
      need(!is.na(i),
        message = "Please select an entry from the table."
      )
    )
    dtu_tested <- rvalues$key()[[i, "dtu_pvadj"]]
    validate(
      need(!is.na(dtu_tested),
        message = "Entry not tested for DTU."
      )
    )
    i <- rvalues$key()[[i, "gene_id"]]
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
    metadata <- rvalues$metadata() %>%
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
      type = "heatmap",
    ) %>%
      add_annotations(x = a, y = b, text = textMatrix, xref = "x", yref = "y", showarrow = FALSE, font = list(color = "black")) %>%
      config(displayModeBar = FALSE) %>%
      layout(title = "Module-trait correlation")
  })


  # enrichment map related content ---------------------------------------------
  output$enrich_table <- renderReactable({
    rvalues$res_enrich() %>%
        rename(
          c(
            "id" = "gs_id",
            "description" = "gs_description", 
            "expected" = "Expected", 
            "observed"= "gs_de_count", 
            "pval" = "gs_pvalue"))  %>%
        select(id, description, pval, expected, observed) %>%
        arrange(pval) %>%
      reactable(
        .,
        searchable = TRUE,
        striped = TRUE,
        rownames = FALSE,
        defaultPageSize = 5,
        highlight = TRUE,
        selection = "single",
        onClick = "select",
        rowStyle = list(cursor = "pointer"),
        theme = reactableTheme(
          stripedColor = "#f6f8fa",
          highlightColor = "#f0f5f9",
          cellPadding = "8px 12px",
          rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #FF0000")
        ),
        columns = list(
          id = colDef(
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'http://amigo.geneontology.org/amigo/term/' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 100),
          description = colDef(
            width = 250)
      )
    )
  })

  emap_graph <- reactive({

    res_enrich <- rvalues$res_enrich() %>% as.data.frame(.)
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

    annotation_obj <- rvalues$annotation_obj()

    emg <- GeneTonic::enrichment_map(
      res_enrich,
      res_de,
      annotation_obj,
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
    i <- getReactableState("enrich_table", "selected")    
    validate(
      need(!is.na(i),
        message = "Please select a gene set from the table."
      )
    )

    res_enrich <- rvalues$res_enrich() %>% as.data.frame(.)
    sel_gs <- res_enrich[[i, "gs_id"]]

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

    annotation_obj <- rvalues$annotation_obj()
    genes <- res_enrich[i, 'gs_genes']
    genes <- unlist(strsplit(genes, ","))
    genes <- annotation_obj[
      match(genes, annotation_obj$gene_name), ]$gene_id
    counts <- rvalues$counts() %>%
      filter(row_names %in% local(genes))

    metadata <- rvalues$metadata()

    common <- intersect(
      colnames(counts)[2: ncol(counts)],
      metadata$row_names)

    counts_rnames <- counts$row_names
    counts <- as.data.frame(counts[, common])
    rownames(counts) <- counts_rnames

    metadata <- metadata %>% 
      filter(row_names %in% common)

    dds <- DESeq2::DESeqDataSetFromMatrix(
      counts,
      colData=metadata,
      design = ~Etiology + Race + Sex + Age + SV1 + SV2)

    colnames(dds) <- NULL
    show_row_names <- TRUE
    if(length(rownames(dds)) > 30){
      show_row_names <- FALSE
    }
  
    GeneTonic::gs_heatmap(
      dds,
      res_de,
      res_enrich,
      annotation_obj = annotation_obj,
      genelist = genes,
      FDR = 0.05,
      de_only = FALSE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      center_mean = TRUE,
      scale_row = TRUE,
      anno_col_info = c("Etiology", "Race", "Sex", "Age", "SV1", "SV2"),
      show_row_names = show_row_names
    )
  })

  output$enriched_funcres <- renderPlotly({

    res_enrich <- rvalues$res_enrich()
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

    annotation_obj <- rvalues$annotation_obj()

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


  # Carnival related content ---------------------------------------------------
  output$visnet_igraph <- renderVisNetwork({
    con %>%
      tbl("carnival") %>%
      filter(contrast == local(input$selected_contrast)) %>%
      pull(igraph) %>%
      jsonlite::unserializeJSON() %>%
      igraph::upgrade_graph(.) %>%
      permute.vertices(., Matrix::invPerm(order(V(.)$name))) %>%
      visNetwork::visIgraph() %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        nodesIdSelection = TRUE
      ) %>%
      visLegend(addEdges = ledges, addNodes = lnodes, useGroups = FALSE) %>% 
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
          h5("Bookmarked gene sets"),
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
    rvalues$res_enrich() %>%
      filter(gs_id %in% rvalues$mygenesets) %>%
      select(gs_id, gs_description) %>%
      reactable(rownames = FALSE)
  })

  observeEvent(input$bookmarker, {
    if (input$magnetique_tab == "tab-welcome") {
      showNotification("Welcome to magnetique! Navigate to the main tabs of the application to use the Bookmarks functionality.")
    } else if (input$magnetique_tab == "tab-gene-view") {
      i <- getReactableState("de_table", "selected")
      if (!is.null(i)) {
        sel_gene <- rvalues$key()[[i, "gene_id"]]
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
        sel_gs <- rvalues$res_enrich()[i, ]$gs_id
        if (!sel_gs %in% rvalues$mygenesets) {
          rvalues$mygenesets <- c(rvalues$mygenesets, sel_gs)
          showNotification(
            sprintf(
              "The selected geneset %s was added to the bookmarked gene sets.",
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
  
  observeEvent(input$tour_carnivalview, {
    tour <- read.delim("tours/intro_carnivalview.txt",
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

  .actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
}
shinyApp(magnetique_ui, magnetique_server)
