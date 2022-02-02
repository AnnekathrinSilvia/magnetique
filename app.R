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
library(visNetwork)

options(spinner.type = 6)

# sourcing external files -------------------------------------------------
source("utils.R")
source("data_preparation.R")
source("heatmap.R")


# ui definition -----------------------------------------------------------
magnetique_ui <- shinydashboard::dashboardPage(
  title = "magnetique",
  header = shinydashboard::dashboardHeader(title = "magnetique"),

  # sidebar definition ------------------------------------------------------
  sidebar = shinydashboard::dashboardSidebar(
    title = "Options",
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
  ),

  # body definition ---------------------------------------------------------
  body = shinydashboard::dashboardBody(
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
            includeMarkdown("data/overview.md")
          )
        )
      ),
      shiny::tabPanel(
        title = "Gene View", icon = icon("heartbeat"), value = "tab-gene-view",
        fluidRow(
          column(
            width = 6,
            withSpinner(
              reactableOutput("de_table")
            )
          ),
          column(
            width = 3,
            # withLoader(plotlyOutput("de_volcano"), type="image", loader="/heart.gif") # only works when starting app over RunApp Button (and doesn't look good)
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
            reactableOutput("enrich_table")
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
            plotOutput("emap_signature")
          )
        )
      ),
      shiny::tabPanel(
        title = "Bookmarks", icon = icon("bookmark"), value = "tab-bookmark",
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
  progress <- Progress$new(session)
  progress$set(value = 0.0, message = "Loading.")
  library(dplyr, warn.conflicts = FALSE)
  library(magrittr, warn.conflicts = FALSE)
  library("GeneTonic")
  library("ggplot2")
  library("ggrepel")
  library(igraph, warn.conflicts = FALSE)
  library("pheatmap")
  suppressPackageStartupMessages({
    library("ComplexHeatmap")
    library("DESeq2")
    library("RColorBrewer")
  })
  progress$set(value = 0.1, message = "Loading.")
  res_dtu <- readRDS("MAGNetApp/data/DTU/summarized_experiment.RDS")
  progress$set(value = .33, message = "Loading.")
  gtf <- readRDS("MAGNetApp/data/DTU/gtf.RDS")

  # reactive objects and setup commands -------------------------------------
  rvalues <- reactiveValues()
  rvalues$mygtl <- NULL
  rvalues$key <- NULL
  rvalues$myvst <- NULL

  # selector trigger data loading
  rvalues$data <- eventReactive(
    c(input$selected_contrast, input$selected_ontology),
    {
      future({
        prepare_data(input$selected_contrast, input$selected_ontology)
      })
    }
  )

  rvalues$mygtl <- reactive({
    rvalues$data() %...>% {
      data <- .
      extract2(data, "genetonic")
    }
  })

  rvalues$myvst <- reactive({
    rvalues$data() %...>% {
      data <- .
      data %>%
        extract2("genetonic") %>%
        extract2("dds") %>%
        vst(.)
    }
  })

  rvalues$key <- reactive({
    rvalues$data() %...>% {
      data <- .
      data %>%
        extract2("res") %>%
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
        highlight_key()
    }
  })

  # DE related content ---------------------------------------------------------
  output$de_table <- renderReactable(
    {
      rvalues$key() %...>% {
        data <- .
        reactable(
          data,
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
          list(
            log2FoldChange = colDef(
              cell = function(value) format(round(value, 2))
              ),
            padj = colDef(
              cell = function(value) format(round(value, 2))
            ),
            dtu_pvadj = colDef(
              cell = function(value) format(round(value, 2))
            ),
            dtu_dif = colDef(
              cell = function(value) format(round(value, 2))
            )
          )
        )
      }
    }
  )

  output$de_volcano <- renderPlotly({
    rvalues$key() %...>% {
      data <- .
      data %>%
        plot_ly(., color = I("black"), showlegend = FALSE) %>%
        add_markers(
          x = ~log2FoldChange,
          y = ~ -log10(padj),
          type = "scatter",
          text = ~ paste0(
            "<b>", SYMBOL, "</b>",
            "<br><i>GeneID</i>: ", gene_id,
            "<br><i>Log2FC</i> = ", format(round(log2FoldChange, 2), nsmall = 2),
            "<br><i>p-value (adjusted)</i> = ", format(padj),
            nsmall = 2
          ),
          hoverinfo = "text"
        ) %>%
        config(displayModeBar = FALSE) %>%
        layout(title = "Differentially expressed genes") %>%
        toWebGL() %>%
        highlight(
          on = "plotly_click",
          off = "plotly_doubleclick",
          color = "red"
        )
    }
  })

  output$dtu_volcano <- renderPlotly({
    rvalues$key() %...>% {
      data <- .
      data %>%
        plot_ly(., color = I("black"), showlegend = FALSE) %>%
        add_markers(
          x = ~dtu_dif,
          y = ~ -log10(dtu_pvadj),
          # size = ~n_transcript,
          type = "scatter",
          text = ~ paste0(
            "<b>", SYMBOL, "</b>",
            "<br><i>GeneID</i>: ", gene_id,
            "<br><i>max(abs(dif))</i> = ", format(round(log2FoldChange, 2), nsmall = 2),
            "<br><i>p-value (adjusted)</i> = ", format(padj),
            nsmall = 2
          ),
          hoverinfo = "text"
        ) %>%
        config(displayModeBar = FALSE) %>%
        layout(title = "Genes with differential transcript usage") %>%
        toWebGL() %>%
        highlight(
          on = "plotly_click",
          off = "plotly_doubleclick",
          color = "red"
        )
    }
  })

  # enrichment map related content ---------------------------------------------
  output$enrich_table <- renderReactable({
    rvalues$mygtl() %...>% {
      mygtl <- .
      myres_enrich <- mygtl$res_enrich
      df <- data.frame(
        description = myres_enrich$gs_description,
        obs = myres_enrich$DE_count,
        exp = myres_enrich$Expected,
        padj = myres_enrich$gs_pvalue
      )
      df <- df[order(df$obs, decreasing = T), ]
      rownames(df) <- myres_enrich$gs_id
      colnames(df) <- c("Description", "Observed", "Expected", "padj")
      reactable(
        df,
        columns = list(
          padj = colDef(
            cell = function(value) format(round(value, 2))
          )
        )
      )
    }
  })

  output$visnet_em <- renderVisNetwork({
    rvalues$mygtl() %...>% {

      mygtl <- .
      emg <- enrichment_map(
        gtl = mygtl,
        n_gs = input$number_genesets,
        overlap_threshold = 0.1,
        scale_edges_width = 200,
        color_by = input$color_by
      )

    visNetwork::visIgraph(emg) %>%
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
    }
  })

  output$emap_signature <- renderPlot({
    rvalues$data() %...>% {
      data <- .
      mygtl <- extract2(data, "genetonic")
      myvst <- extract2(mygtl, "dds")
      
      cur_gsid <- mygtl$res_enrich$gs_id[
        match(input$visnet_em_selected, mygtl$res_enrich$gs_description)]
      validate(
        need(!is.na(cur_gsid),
        message = "Please select a gene set from the Enrichment Map."
        ))
    heatmap(
      se = myvst,
      gtl = mygtl,
      geneset_id = cur_gsid,
      FDR = 0.05,
      de_only = FALSE,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      center_mean = TRUE,
      scale_row = TRUE,
      anno_col_info = "Etiology"
    )}
  })

  output$enriched_funcres <- renderPlotly({
    rvalues$mygtl() %...>% {
      gtl <- .
      ggplotly(
        enhance_table(
          gtl$res_enrich,
          gtl$res_de,
          annotation_obj = gtl$annotation_obj,
          n_gs = input$number_genesets,
          chars_limit = 50
        )
      )
    }
  })

  # DTU related content --------------------------------------------------------
  output$dtu_plot <- renderPlot({
    row <- req(getReactableState("de_table", "selected"))
    validate(
      need(!is.na(row),
        message = "Please select an entry from the table."
      )
    )

    rvalues$data() %...>% {
      data <- .
      res <- data %>%
        extract2("res") 

      genes_dtu <- unique(rowData(res_dtu)$gene_id)
      row <- res[row, "gene_id"]
      validate(need(
        row %in% genes_dtu,
          message = paste("The gene you selected, ", row, ", is not a DTU gene. Please select another gene from the table")
        ))
      plot_dtu(row, res_dtu, gtf)
    }
  })

  output$carnival_launch <- renderUI({
    tagList(
      actionButton(
        inputId = "btn_show_carnival",
        icon = icon("flask"),
        label = "Show Carnival View", style = .actionbutton_biocstyle
      ),
      actionButton(
        inputId = "btn_switch_emap",
        icon = icon("project-diagram"),
        label = "Jump to Enrichtment Map",
        style = .actionbutton_biocstyle
      )
    )
  })


  output$visnet_igraph <- renderVisNetwork({
    rvalues$data() %...>% {
      data <- .
      data %>%
        extract2("igraph") %>%
        visNetwork::visIgraph(.) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = TRUE
        ) %>%
        visHierarchicalLayout(
          levelSeparation = 100,
          nodeSpacing = 500,
          shakeTowards = "leaves"
        ) %>% # same as visLayout(hierarchical = TRUE)
        visExport(
          name = "igraph",
          type = "png",
          label = "Save igraph graph"
        )
    }
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
    book_df_genes <- data.frame(
      gene_id = rvalues$mygenes,
      gene_symbols = rvalues$mygenes
    )
    reactable(book_df_genes, rownames = FALSE)
  })

  output$bookmarks_genesets <- renderReactable({
    book_df_genesets <- data.frame(
      geneset_id = rvalues$mygenesets,
      geneset_description = rvalues$mygenesets
    )
    reactable(book_df_genesets, rownames = FALSE)
  })


  # Other content --------------------------------------------------------------
  output$de_volcano_signature <- renderPlot({
    signature_volcano(
      gtl = rvalues$mygtl(),
      FDR = 0.05
    )
  })


  output$gtl_loaded <- renderText({
    describe_gtl(gtl = rvalues$mygtl())
  })

  output$team_list <- renderTable({
    make_team_df
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
# Launching magnetique! --------------------------------------------------------
shinyApp(magnetique_ui, magnetique_server)
