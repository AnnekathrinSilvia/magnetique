# loading libraries -------------------------------------------------------
library(shiny)
library(future)
library(promises)
plan(multisession)

library(shinycssloaders)
library(shinyjs)
library(plotly, warn.conflicts = FALSE)
library(reactable, warn.conflicts = FALSE)
library(bs4Dash, warn.conflicts = FALSE)
library(shinydashboard, warn.conflicts = FALSE)
library(visNetwork, warn.conflicts = FALSE)
library(rintrojs, warn.conflicts = FALSE)
library(shinyBS, warn.conflicts = FALSE)
library(igraph)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)

options(spinner.type = 6)

# sourcing external files -------------------------------------------------
source("utils.R")

# ui definition -----------------------------------------------------------
magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  header = shinydashboard::dashboardHeader(disable = TRUE),

  # sidebar definition ------------------------------------------------------
  sidebar = dashboardSidebar(
    img(src = "magnetique_logo.png", class="img-responsive"),
    h1("Magnetique", align='center'),
    hidden(h2("Options:", style = 'margin: 15px', id='options')),
    uiOutput("ui_sidebar")
  ),

  # body definition ---------------------------------------------------------
  body = dashboardBody(
    introjsUI(),
    ## handling the overflow in vertical direction for the tabBox
    shiny::tags$head(
      tags$link(rel="shortcut icon", href="favicon.ico"),
      shiny::tags$style(
        HTML("#myScrollBox{
                overflow-y: scroll;
              }")
      )
    ),
    ## using shinyjs inside the app
    shinyjs::useShinyjs(),
    ## left ctrl to trigger the bookmarking functionality (via button)
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
    
    div(
      id = "myScrollBox",
      ## will adjust indentation later ##
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
            width = 6,
            includeMarkdown("www/overview.md"),
            actionButton(inputId = "popup_about_us",
                         label = "More about the development team",
                         icon = icon("users"), style = .actionbutton_biocstyle
            )
          ),
          column(
            width = 6,
            h2("Patient characteristics and technical covariates"),
            includeHTML("www/metadata_table.html")
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
        title = "Carnival View", icon = icon("network-wired"), value = "tab-carnival",
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
            width = 3,
            img(src = "https://user-images.githubusercontent.com/4517913/166660064-8b171940-b2ce-46e5-bfc6-c2c6685a5bd0.jpeg", class="img-responsive")
            ),
          column(
            width = 9,
            visNetworkOutput("visnet_carnival", height = "500px")
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
    ) ## end of the scrollbox
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


  # reactive objects and setup commands -------------------------------------
  rvalues <- reactiveValues()
  rvalues$mygenes <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("gene_id", "gene_name"))))
  rvalues$mygenesets <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("gs_id", "gs_description"))))


  # sidebar server-side -----------------------------------------------------
  output$ui_sidebar <- renderUI({
    tagList(
      sidebarMenu(
        hidden(
          selectInput("selected_contrast",
            label = "Contrast id",
            choices = c(
              "DCMvsHCM",
              "DCMvsNFD",
              "HCMvsNFD"
            ),
            selected = "DCMvsHCM"
          )
        ),
        shinyjs::hidden(
          tagList(
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
            )
          )
        ),
        hidden(
          actionButton("bookmarker",
          label = "Bookmark", icon = icon("heart"),
          style = "color: #ffffff; background-color: #ac0000; border-color: #ffffff"
          )
        )
      )
    )
  })
  
  observeEvent(input$magnetique_tab, {
    if (input$magnetique_tab == "tab-geneset-view") {
      shinyjs::show("selected_contrast")
      shinyjs::show("selected_ontology")
      shinyjs::show("number_genesets")
      shinyjs::show("color_by")
      shinyjs::show("bookmarker")
      shinyjs::show("options")
    } else if (
      input$magnetique_tab == "tab-gene-view" 
      || input$magnetique_tab == "tab-carnival") {
      shinyjs::show("selected_contrast")
      shinyjs::hide("selected_ontology")
      shinyjs::hide("number_genesets")
      shinyjs::hide("color_by")
      shinyjs::show("bookmarker")
      shinyjs::show("options")
    } else {
      shinyjs::hide("selected_contrast")
      shinyjs::hide("selected_ontology")
      shinyjs::hide("number_genesets")
      shinyjs::hide("color_by")
      shinyjs::hide("bookmarker")
      shinyjs::hide("options")
    }
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
          "dtu_dif"
        )
      ) %>%
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

  rvalues$vst <- reactive({
    con %>% tbl("vst")
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
        wrap = FALSE,
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
          )
        ),
        defaultColDef = colDef(sortNALast = TRUE)
      )
  })

  prev_selected <- reactive(getReactableState("de_table", "selected"))

  output$de_volcano <- renderPlotly({
    p <- rvalues$key() %>%
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
    
    if(length(prev_selected()) == 1) {
      df <- rvalues$key() %>%
        slice(prev_selected())

      p <- add_trace(
        p, 
        x = df[[1, "log2FoldChange"]],
        y = ~ -log10(df[[1, "padj"]]),
        type = "scatter",
        mode = "markers", 
        color = I("red"), 
        inherit = FALSE)
    }
    p
    
  })


  output$dtu_volcano <- renderPlotly({  
    p <- rvalues$key() %>%
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


  if(length(prev_selected()) == 1) {
    df <- rvalues$key() %>%
      slice(prev_selected())

    p <- add_trace(
      p, 
      x = df[[1, "dtu_dif"]],
      y = -log10(df[[1, "dtu_pvadj"]]),
      type = "scatter",
      mode = "markers", 
      color = I("red"), 
      inherit = FALSE)
  }
  p
 
  })

  observeEvent(prev_selected(), {

    plotlyProxy("de_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list('rgba(0,0,0,.10)'), 0)
    
    plotlyProxy("dtu_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list('rgba(0,0,0,.10)'), 0)        
        
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
        alpha = 0.9, size = 2, show.legend = TRUE, na.rm = TRUE
      ) +
      geom_boxplot(aes(x = transcript_id, y = proportion, fill = Etiology),
        outlier.size = 0, alpha = 0.4, lwd = 0.5, show.legend = FALSE
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
        mutate_at(vars(pval), ~round(., 2)) %>%
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
        wrap = FALSE,
        rowStyle = list(cursor = "pointer"),
        theme = reactableTheme(
          stripedColor = "#f6f8fa",
          highlightColor = "#f0f5f9",
          cellPadding = "8px 12px",
          rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #FF0000")
        ),
        defaultColDef = colDef(width = 70),
        columns = list(
          id = colDef(
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'http://amigo.geneontology.org/amigo/term/' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 100,
            width = 120,
            header = with_tooltip("gs_id", "Gene set ID and link to AMIGO db")
            ),
          description = colDef(
            width = 250,
            header = with_tooltip("gs_description", "Description for the gene set")),
          pval = colDef(
            header = with_tooltip("pvalue", "p-value for the TopGO enrichment test")),
          expected = colDef(
            header = with_tooltip("expec.", "Expected number of genes in the gene set")),
          observed = colDef(
            header = with_tooltip("obser.", "Number of genes observed in the gene set"))
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
    validate(
      need({ecount(emap_graph()) > 0}, 
           message = paste0(
             "No edges detected in the enrichment map. ",
             "Please select a larger number of genesets to generate a full graph")
      )
    )
    visNetwork::visIgraph(emap_graph()) %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1
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
    gs_description <- res_enrich[[i, "gs_description"]]

    res_de <- con %>%
      tbl(paste0("res_", local(input$selected_contrast))) %>%
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

    metadata <- rvalues$metadata()

    counts <- rvalues$vst() %>% 
      filter(row_names %in% local(genes)) %>% 
      collect() %>% 
      as.data.frame(.)

    rownames(counts) <- counts[, 1]
    counts[, 1] <- NULL

    se <- SummarizedExperiment::SummarizedExperiment(
      counts,
      colData = metadata,
      rowData = rownames(counts))

    show_row_names <- ifelse(nrow(se) > 30, FALSE, TRUE)  
    GeneTonic::gs_heatmap(
      se,
      res_de,
      res_enrich,
      annotation_obj = annotation_obj,
      genelist = genes,
      FDR = 0.05,
      de_only = FALSE,
      plot_title = stringr::str_glue(
        "Gene signature for {gs_description} geneset"), 
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      center_mean = TRUE,
      scale_row = TRUE,
      anno_col_info = c("Etiology", "Race", "Sex", "Age", "SV1", "SV2"),
      show_row_names = show_row_names,
      show_column_names = FALSE
    )
  })

  output$enriched_funcres <- renderPlotly({
    
    i <- getReactableState("enrich_table", "selected")    
    if(!is.null(i)) {
      res_enrich <- rvalues$res_enrich() %>% 
        slice(i)
      plot_title <- "Enrichment for genes comprising  
        {res_enrich[1, 'gs_id']} geneset \n ({input$selected_contrast})"

    } else {
      res_enrich <- rvalues$res_enrich()

      plot_title <- "Enrichment overview for top genesets 
          ({input$selected_ontology} and {input$selected_contrast})"
    }
        
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
        chars_limit = 50,
        plot_title = stringr::str_glue(plot_title)
      ) 
    )
  })


  # Carnival related content ---------------------------------------------------
  output$visnet_carnival <- renderVisNetwork({
    con %>%
      tbl("carnival") %>%
      filter(contrast == local(input$selected_contrast)) %>%
      pull(igraph) %>%
      jsonlite::unserializeJSON() %>%
      igraph::upgrade_graph(.) %>%
      permute.vertices(., Matrix::invPerm(order(V(.)$name))) %>%
      visNetwork::visIgraph() %>%
      visNodes(font = list(background = "white"))  %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1
        ),
        nodesIdSelection = TRUE
      ) %>%
      visInteraction(navigationButtons = TRUE) %>%
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
          uiOutput("ui_bookmarks_genes")
        ),
        column(
          width = 6,
          h5("Bookmarked gene sets"),
          uiOutput("ui_bookmarks_genesets")
        )
      )
    )
  })
  
  output$ui_bookmarks_genes <- renderUI({
    validate(
      need(
        nrow(rvalues$mygenes) > 0,
        "Please select at least one gene with the Bookmark button"
      )
    )
    tagList(
      reactableOutput("bookmarks_genes"),
      downloadButton("download_bookmarks_genes", "Download as csv")
    )
  })

  output$ui_bookmarks_genesets <- renderUI({
    validate(
      need(
        nrow(rvalues$mygenesets) > 0,
        "Please select at least one geneset with the Bookmark button"
      )
    )
    tagList(
      reactableOutput("bookmarks_genesets"),
      downloadButton("download_bookmarks_genesets", "Download as csv")
    )
  })
  
  output$bookmarks_genes <- renderReactable({
    reactable(rvalues$mygenes)
  })
  
  output$bookmarks_genesets <- renderReactable({
    reactable(rvalues$mygenesets, rownames = FALSE)
  })

  output$download_bookmarks_genes <- downloadHandler(
    filename = "magnetique_genes_bookmark.csv",
    content = function(file){
      write.csv(rvalues$mygenes, file=file)
      }
  )

  output$download_bookmarks_genesets <- downloadHandler(
    filename = "magnetique_genesets_bookmark.csv",
    content = function(file){
      write.csv(rvalues$mygenesets, file=file)
      }
  )

  observeEvent(input$bookmarker, {
    if (input$magnetique_tab == "tab-welcome") {
      showNotification("Welcome to magnetique! Navigate to the main tabs of the application to use the Bookmarks functionality.")
    } else if (input$magnetique_tab == "tab-gene-view") {
      i <- getReactableState("de_table", "selected")
      
      if (is.null(i)) {
        showNotification("Select a row in the main table to bookmark it", type = "warning")
      } else {
        key <- rvalues$key()
        df <- key[i, c("gene_id", "SYMBOL")] %>%
          rename(gene_name=SYMBOL)

        sel_gene_id <- df[[1, "gene_id"]]
        sel_gene <- df[[1, "gene_name"]]
        if (sel_gene_id %in% rvalues$mygenes$gene_id) {
          showNotification(sprintf("The selected gene %s (%s) is already in the set of the bookmarked genes.", sel_gene, sel_gene_id), type = "default")
        } else {
          rvalues$mygenes <- rbind(rvalues$mygenes, df)
          showNotification(sprintf("Added %s (%s) to the bookmarked genes. The list contains now %d elements", sel_gene, sel_gene_id, nrow(rvalues$mygenes)), type = "message")
        }
      }
    } else if (input$magnetique_tab == "tab-geneset-view") {
      # handling bookmarks from the table
      i <- getReactableState("enrich_table", "selected")
      # as well as from the interactive graph
      cur_em <- input$visnet_em_selected
      
      if (is.null(i) & (cur_em == "")) {
        showNotification("Select a row in the enrichment table or a node from the enrichment map to bookmark it", type = "warning")
      } else {
        # handling bookmarks from the table
        if (!is.null(i)) {
          df <- rvalues$res_enrich() %>%
            slice(i) %>%
            select(gs_id, gs_description)
          
          sel_gs_id <- df[[1, "gs_id"]]
          sel_gs <- df[[1, "gs_description"]]
          
          if (sel_gs_id %in% rvalues$mygenesets$gs_id) {
            showNotification(sprintf("The selected gene set, %s (%s), is already in the set of the bookmarked genesets.", sel_gs, sel_gs_id), type = "default")
          } else {
            rvalues$mygenesets <- rbind(rvalues$mygenesets, df)
            showNotification(sprintf("Added %s (%s) to the bookmarked genesets. The list contains now %d elements", sel_gs, sel_gs_id, nrow(rvalues$mygenesets)), type = "message")
          }
        }
        
        # as well as from the interactive graph
        if (cur_em != "") {
          re <- rvalues$res_enrich()
          cur_em_id <- re$gs_id[match(cur_em, re$gs_description)]
          
          df <- data.frame(
            gs_id = cur_em_id,
            gs_description = cur_em
          )

          if (cur_em_id %in% rvalues$mygenesets$gs_id) {
            showNotification(sprintf("The selected gene set, %s (%s), is already in the set of the bookmarked genesets.", cur_em, cur_em_id), type = "default")
          } else {
            rvalues$mygenesets <- rbind(rvalues$mygenesets, df)
            showNotification(sprintf("Added %s (%s) to the bookmarked genesets. The list contains now %d elements", cur_em, cur_em_id, nrow(rvalues$mygenesets)), type = "message")
          }
        }
      }
      
    } else if(input$magnetique_tab == "tab-carnival") {
      
      sel_gene <- input$visnet_carnival_selected
      annotation_obj <- rvalues$annotation_obj()
      
      if (sel_gene == "") {
        showNotification("Select a node in the carnival graph to bookmark it", type = "warning")
      } else {
        sel_gene_id <- annotation_obj$gene_id[match(sel_gene, annotation_obj$gene_name)]
        df <- data.frame(
          gene_id = sel_gene_id,
          gene_name = sel_gene
        )
        
        if (sel_gene_id %in% rvalues$mygenes$gene_id) {
          showNotification(sprintf("The selected gene %s (%s) is already in the set of the bookmarked genes.", sel_gene, sel_gene_id), type = "default")
        } else {
          rvalues$mygenes <- rbind(rvalues$mygenes, df)
          showNotification(sprintf("Added %s (%s) to the bookmarked genes. The list contains now %d elements", sel_gene, sel_gene_id, nrow(rvalues$mygenes)), type = "message")
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

  observeEvent(input$popup_about_us, {
    showModal(
      modalDialog(
        title = "Project members (alphabetical order)",
        size = "l",
        tableOutput("team_list"), 
        easyClose = TRUE,
        footer = ""
      )
    )
  })
}
shinyApp(magnetique_ui, magnetique_server)
