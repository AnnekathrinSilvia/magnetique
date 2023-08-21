# loading libraries -------------------------------------------------------
library(shiny)
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
    img(src = "magnetique_logo.png", class = "img-responsive"),
    hidden(h2("Options:", style = "margin: 15px", id = "options")),
    uiOutput("ui_sidebar")
  ),
  # body definition ---------------------------------------------------------
  body = dashboardBody(
    introjsUI(),
    ## handling the overflow in vertical direction for the tabBox
    shiny::tags$head(
      shiny::tags$style(
        HTML("
              /* fix overflow */
              #myScrollBox{ overflow-y: scroll; }
              /* side bar color */
              .main-sidebar { background-color: white !important; }
              .skin-blue .main-sidebar .sidebar{ color: black; }
              /* body */
              .content-wrapper, .right-side {background-color: #d3d3d3 ;}")
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
              actionButton(
                inputId = "popup_about_us",
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
                tagList(
                  reactableOutput("de_table"),
                  tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('de_table', 'de_table.csv')")
                )
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
              width = 2,
              withSpinner(
                plotlyOutput("gene_counts")
              )
            ),
            column(
              width = 4,
              withSpinner(
                plotlyOutput("transcript_proportion")
              )
            ),
            column(
              width = 6,
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
                tagList(
                  reactableOutput("enrich_table"),
                  tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('enrich_table', 'enrich_table.csv')")
                )
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
                tagList(
                  visNetworkOutput("visnet_em"),
                  shinydashboard::box(
                    title = "About this enrichment map",
                    width = 12,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    htmlOutput("visnet_explanation")
                  )
                )
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
              id = "carnival_legend",
              width = 3,
              img(src = "https://user-images.githubusercontent.com/4517913/166660064-8b171940-b2ce-46e5-bfc6-c2c6685a5bd0.jpeg", class = "img-responsive")
            ),
            column(
              id = "carnival_network",
              width = 9,
              visNetworkOutput("visnet_carnival", height = "750px")
            )
          )
        ),
        tabPanel(
          title = "RBP:RNA View", icon = icon("arrows-alt-h"), value = "tab-rbp-view",
          fluidRow(
            column(
              width = 12,
              div(
                actionButton(
                  "tour_rbpview",
                  label = "", icon = icon("question-circle"),
                  style = .helpbutton_biocstyle
                ),
                shinyBS::bsTooltip(
                  "tour_rbpview",
                  "Click me to start a tour of this section!",
                  "bottom",
                  options = list(container = "body")
                ),
                style = "float:right"
              )
            )
          ),
          fluidRow(
            id = "rbpview_row1",
            column(
              width = 8,
              withSpinner(
                visNetworkOutput("rbp_network")
              )
            ),
          ),
          fluidRow(
            id = "rbpview_row2",
            column(
              width = 8,
              withSpinner(
                tagList(
                  reactableOutput("rbp_table"),
                  tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('rbp_table', 'rbp_table.csv')")
                )
              )
            ),
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
        )
      )
    ) ## end of the scrollbox
  )
)

# server definition -------------------------------------------------------
magnetique_server <- function(input, output, session) {
  # reactive objects and setup commands -------------------------------------
  rvalues <- reactiveValues()
  rvalues$mygenes <- data.frame(matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("gene_id", "gene_name"))))
  rvalues$mygenesets <- data.frame(matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("gs_id", "gs_description"))))


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
              label = "Gene Ontology (GO)",
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
      input$magnetique_tab == "tab-gene-view" ||
        input$magnetique_tab == "tab-carnival") {
      shinyjs::show("selected_contrast")
      shinyjs::hide("selected_ontology")
      shinyjs::hide("number_genesets")
      shinyjs::hide("color_by")
      shinyjs::show("bookmarker")
      shinyjs::show("options")
    } else if (input$magnetique_tab == "tab-rbp-view") {
      shinyjs::show("selected_contrast")
      shinyjs::hide("bookmarker")
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
    res[[input$selected_contrast]] %>%
      select(
        c(
          "gene_id",
          "SYMBOL",
          "padj",
          "log2FoldChange",
          "dtu_pvadj",
          "dtu_dif"
        )
      )
  })

  rvalues$metadata <- reactive({
    metadata
  })

  rvalues$counts <- reactive({
    counts
  })

  rvalues$vst <- reactive({
    vst
  })

  rvalues$res_enrich <- reactive({
    res_enrich[[input$selected_contrast]] %>%
      filter(ontology == input$selected_ontology) %>%
      as.data.frame(.) %>%
      `rownames<-`(.$gs_id)
  })

  rvalues$annotation_obj <- reactive({
    annotation_obj %>%
      as.data.frame() %>%
      `rownames<-`(.$gene_id)
  })

  # DE related content ---------------------------------------------------------
  output$de_table <- renderReactable({
    rvalues$key() %>%
      mutate_at(vars(log2FoldChange, dtu_dif), ~ round(., 2)) %>%
      mutate_at(vars(padj, dtu_pvadj), ~ round(-log10(.), 2)) %>%
      reactable(
        .,
        language = reactableLang(
          filterPlaceholder = "Filter"
        ),
        searchable = FALSE,
        striped = TRUE,
        showPageSizeOptions = TRUE,
        defaultSorted = c("padj"),
        defaultPageSize = 5,
        pageSizeOptions = c(5, 10, 25, 50),
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
            filterable = TRUE,
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
            filterable = TRUE,
            header = with_tooltip("gene_name", "Gene symbol")
          ),
          log2FoldChange = colDef(
            name = "dge_log2fc",
            header = with_tooltip("dge_log2fc", "Fold change for the DESeq2 analysis")
          ),
          padj = colDef(
            name = "-log10_dge_padj",
            defaultSortOrder = "desc",
            minWidth = 120,
            filterable = TRUE,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] >= filterValue
                })
                }"),
            header = with_tooltip("-log10_dge_padj", "Adjusted p-value for the DESeq2 analysis (-log10 transformed)")
          ),
          dtu_pvadj = colDef(
            name = "-log10_dtu_padj",
            filterable = TRUE,
            width = 130,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] >= filterValue
                })
                }"),
            header = with_tooltip("-log10_dtu_padj", "Adjusted p-value for the DRIMseq analysis (-log10 transformed)")
          ),
          dtu_dif = colDef(
            name = "dtu_dif",
            header = with_tooltip("dtu_dif", "Difference in isoform usage for the DRIMseq analysis")
          )
        ),
        defaultColDef = colDef(sortNALast = TRUE),
        elementId = "de_table"
      )
  })

  prev_selected <- reactive(getReactableState("de_table", "selected"))

  output$de_volcano <- renderPlotly({
    p <- rvalues$key() %>%
      plot_ly(., color = I("black"), showlegend = FALSE) %>%
      add_markers(
        x = ~log2FoldChange,
        y = ~ -log10(padj),
        type = "scatter",
        marker = list(opacity = 0.2),
        text = ~ paste0(
          "<b>", SYMBOL, "</b>",
          "<br><i>GeneID</i>: ", gene_id,
          "<br><i>Log2FC</i> = ", format(round(log2FoldChange, 2)),
          "<br><i>p-value (adjusted)</i> = ", format(round(padj, 2))
        ),
        hoverinfo = "text"
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c(
          "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d",
          "hoverCompareCartesian", "hoverClosestCartesian",
          "select2d", "lasso2d", "zoom2d"
        ),
        toImageButtonOptions = list(
          format = "svg",
          width = 700,
          height = 500,
          filename = stringr::str_glue("magnetique_dge_volcano_{input$selected_contrast}")
        )
      ) %>%
      toWebGL() %>%
      layout(
        title = "Differentially expressed genes",
        yaxis = list(title = "-log10(dge_padj)"),
        xaxis = list(title = "dge_log2fc")
      )

    if (length(prev_selected()) == 1) {
      df <- rvalues$key() %>%
        slice(prev_selected())

      p <- add_trace(
        p,
        x = df[[1, "log2FoldChange"]],
        y = ~ -log10(df[[1, "padj"]]),
        type = "scatter",
        mode = "markers",
        color = I("red"),
        inherit = FALSE
      )
    }
    p
  })


  output$dtu_volcano <- renderPlotly({
    p <- rvalues$key() %>%
      plot_ly(., color = I("black"), showlegend = FALSE) %>%
      add_markers(
        x = ~dtu_dif,
        y = ~ -log10(dtu_pvadj),
        type = "scatter",
        marker = list(opacity = 0.2),
        text = ~ paste0(
          "<b>", SYMBOL, "</b>",
          "<br><i>GeneID</i>: ", gene_id,
          "<br><i>dif </i> = ", format(round(dtu_dif, 2)),
          "<br><i>p-value (adjusted)</i> = ", format(round(dtu_pvadj, 2))
        ),
        hoverinfo = "text"
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c(
          "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d",
          "hoverCompareCartesian", "hoverClosestCartesian",
          "select2d", "lasso2d", "zoom2d"
        ),
        toImageButtonOptions = list(
          format = "svg",
          width = 700,
          height = 500,
          filename = stringr::str_glue("magnetique_dge_volcano_{input$selected_contrast}")
        )
      ) %>%
      toWebGL() %>%
      layout(
        title = "Genes with differential transcript usage",
        yaxis = list(title = "-log10(dtu_padj)"),
        xaxis = list(title = "dtu_dif")
      )


    if (length(prev_selected()) == 1) {
      df <- rvalues$key() %>%
        slice(prev_selected())

      p <- add_trace(
        p,
        x = df[[1, "dtu_dif"]],
        y = -log10(df[[1, "dtu_pvadj"]]),
        type = "scatter",
        mode = "markers",
        color = I("red"),
        inherit = FALSE
      )
    }
    p
  })
  observeEvent(prev_selected(), {
    plotlyProxy("de_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list("rgba(0,0,0,.10)"), 0)

    plotlyProxy("dtu_volcano", session) %>%
      plotlyProxyInvoke("restyle", "marker.color", list("rgba(0,0,0,.10)"), 0)
  })

  output$gene_counts <- renderPlotly({
    i <- getReactableState("de_table", "selected")
    validate(need(!is.na(i),
      message = "Please select an entry from the table."
    ))
    i <- rvalues$key()[[i, "gene_id"]]

    counts <- rvalues$counts() %>%
      select(-contrast) %>%
      filter(row_names == !!i) %>%
      dplyr::select(-row_names) %>%
      mutate(n = n()) %>%
      pivot_longer(-n)

    metadata <- rvalues$metadata() %>%
      select(row_names, Etiology, Race, Age, Sex, SV1, SV2) %>%
      tibble::rownames_to_column("name")

    left_join(counts, metadata, by = "name") %>%
      plot_ly(
        .,
        type = "box",
        x = ~Etiology,
        y = ~ log10(value),
        boxpoints = "all",
        jitter = 0.8,
        pointpos = 0,
        hoverinfo = "text",
        marker = list(opacity = 0.3),
        text = ~ paste0(
          "<br><i>Sex</i>: ", Sex,
          "<br><i>Race</i> = ", Race,
          "<br><i>Age</i> = ", Age,
          "<br><i>SV1</i> = ", round(SV1, 2),
          "<br><i>SV2</i> = ", round(SV2, 2)
        ),
        color = ~ factor(Etiology, levels = c("NFD", "DCM", "HCM")),
        colors = c(I("steelblue"), I("gold"), I("forestgreen"))
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c(
          "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d",
          "hoverCompareCartesian", "hoverClosestCartesian",
          "select2d", "lasso2d", "zoom2d"
        ),
        toImageButtonOptions = list(
          format = "svg",
          width = 700,
          height = 500,
          filename = stringr::str_glue("magnetique_gene_counts_{i}")
        )
      ) %>%
      layout(
        title = "Gene counts",
        xaxis = list(title = "", showticklabels = FALSE),
        yaxis = list(title = "log10(normalized counts)"),
        legend = list(orientation = "h")
      )
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
    this_gtf <- gtf %>% subset(gene_id == i)

    this_gtf <- GenomicRanges::split(this_gtf, gtf$transcript_id)
    suppressMessages({
      plot_gene_structure(this_gtf) + labs(title = "Gene structure") + theme(plot.title = element_text(size = 18))
    })
  })

  output$transcript_proportion <- renderPlotly({
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
    x <- gene2tx %>%
      filter(gene_id == !!i) %>%
      left_join(
        dtu_fit_proportions,
        by = (c("transcript_id" = "row_names"))
      )

    x <- x %>% pivot_longer(
      -c(gene_id, transcript_id),
      names_to = "Run",
      values_to = "proportion"
    )
    metadata <- rvalues$metadata()

    left_join(x, metadata, by = c("Run" = "Run")) %>%
      plot_ly(
        type = "box",
        boxpoints = "all",
        jitter = 1,
        pointpos = 0,
        y = ~transcript_id,
        x = ~proportion,
        hoverinfo = "text",
        text = ~ paste0(
          "<br><i>Sex</i>: ", Sex,
          "<br><i>Race</i> = ", Race,
          "<br><i>Age</i> = ", Age,
          "<br><i>SV1</i> = ", round(SV1, 2),
          "<br><i>SV2</i> = ", round(SV2, 2)
        ),
        color = ~ factor(Etiology, levels = c("NFD", "DCM", "HCM")),
        colors = c(I("steelblue"), I("gold"), I("forestgreen")),
        orientation = "h",
        marker = list(opacity = 0.3)
      ) %>%
      layout(
        boxmode = "group",
        title = "Transcript proportion",
        yaxis = list(title = ""),
        showlegend = FALSE
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c(
          "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d",
          "hoverCompareCartesian", "hoverClosestCartesian",
          "select2d", "lasso2d", "zoom2d"
        ),
        toImageButtonOptions = list(
          format = "svg",
          width = 700,
          height = 500,
          filename = stringr::str_glue("magnetique_transcript_proportion_{i}")
        )
      )
  })

  # enrichment map related content ---------------------------------------------
  output$enrich_table <- renderReactable({
    rvalues$res_enrich() %>%
      rename(
        c(
          "id" = "gs_id",
          "description" = "gs_description",
          "expected" = "Expected",
          "observed" = "gs_de_count",
          "pval" = "gs_pvalue"
        )
      ) %>%
      select(id, description, pval, expected, observed) %>%
      arrange(pval) %>%
      reactable(
        .,
        elementId = "enrich_table",
        showPageSizeOptions = TRUE,
        defaultPageSize = 5,
        pageSizeOptions = c(5, 10, 25, 50),
        searchable = FALSE,
        striped = TRUE,
        rownames = FALSE,
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
        language = reactableLang(
          filterPlaceholder = "Filter"
        ),
        columns = list(
          id = colDef(
            html = TRUE,
            filterable = TRUE,
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
            filterable = TRUE,
            header = with_tooltip("gs_description", "Description for the gene set")
          ),
          pval = colDef(
            filterable = TRUE,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] <= filterValue
                })
                }"),
            header = with_tooltip("pval", "p-value for the TopGO enrichment test")
          ),
          expected = colDef(
            header = with_tooltip("expec.", "Expected number of genes in the gene set")
          ),
          observed = colDef(
            header = with_tooltip("obser.", "Number of genes observed in the gene set")
          )
        )
      )
  })

  emap_graph <- reactive({
    res_enrich <- rvalues$res_enrich() %>% as.data.frame(.)
    rownames(res_enrich) <- res_enrich$gs_id

    res_de <- res[[input$selected_contrast]]
    res_de <- res_de %>%
      filter(!is.na(SYMBOL)) %>%
      arrange(-dplyr::desc(padj)) %>%
      select(gene_id, log2FoldChange, padj, SYMBOL)

    res_de <- DESeq2::DESeqResults(
      S4Vectors::DataFrame(res_de)
    )

    rownames(res_de) <- res_de$gene_id
    res_de$pvalue <- res_de$padj
    res_de$description <- ""

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
      need(
        {
          ecount(emap_graph()) > 0
        },
        message = paste0(
          "No edges detected in the enrichment map. ",
          "Please select a larger number of genesets to generate a full graph"
        )
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
        type = "pdf",
        label = "Save enrichment map"
      )
  })

  output$visnet_explanation <- renderUI({
    validate(
      need(
        {
          ecount(emap_graph()) > 0
        },
        message = ""
      )
    )

    HTML(
      paste0(
        "If selecting <b>color</b> by <em>z_score</em>, red nodes display gene sets with positive values, ",
        "while blue nodes represent gene sets with negative values. ",
        "If showing the gene set <em>p-value</em>, darker colors represent gene sets with smaller ",
        "p-values, indicating a more significant enrichment.",
        "<br>The <b>size</b> of the node is representing the size of the gene set (number of genes ",
        "assigned to it)."
      )
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

    res_de <- res[[input$selected_contrast]]

    res_de <- res_de %>%
      filter(!is.na(SYMBOL)) %>%
      arrange(-dplyr::desc(padj)) %>%
      select(gene_id, log2FoldChange, padj, SYMBOL)

    res_de <- DESeq2::DESeqResults(
      S4Vectors::DataFrame(res_de)
    )

    rownames(res_de) <- res_de$gene_id
    res_de$pvalue <- res_de$padj
    res_de$description <- ""

    annotation_obj <- rvalues$annotation_obj()
    genes <- res_enrich[i, "gs_genes"]
    genes <- unlist(strsplit(genes, ","))
    genes <- annotation_obj[
      match(genes, annotation_obj$gene_name),
    ]$gene_id

    metadata <- rvalues$metadata()

    counts <- rvalues$vst()[genes, ] %>%
      as.data.frame(.)

    se <- SummarizedExperiment::SummarizedExperiment(
      counts,
      colData = metadata,
      rowData = rownames(counts)
    )

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
        "Gene signature for {gs_description} gene set"
      ),
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
    if (!is.null(i)) {
      res_enrich <- rvalues$res_enrich() %>%
        slice(i)
      plot_title <- "Enrichment for genes comprising
      {res_enrich[1, 'gs_id']} gene set ({input$selected_contrast})"
    } else {
      res_enrich <- rvalues$res_enrich()

      plot_title <- "Enrichment overview
          ({input$selected_ontology} and {input$selected_contrast})"
    }

    res_enrich <- as.data.frame(res_enrich)
    rownames(res_enrich) <- res_enrich$gs_id

    res_de <- res[[input$selected_contrast]]

    res_de <- res_de %>%
      filter(!is.na(SYMBOL)) %>%
      arrange(-dplyr::desc(padj)) %>%
      select(gene_id, log2FoldChange, padj, SYMBOL)

    res_de <- DESeq2::DESeqResults(
      S4Vectors::DataFrame(res_de)
    )

    rownames(res_de) <- res_de$gene_id
    res_de$pvalue <- res_de$padj
    res_de$description <- ""

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
    ) %>% config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c(
        "zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d",
        "hoverCompareCartesian", "hoverClosestCartesian",
        "select2d", "lasso2d", "zoom2d"
      ),
      toImageButtonOptions = list(
        format = "svg",
        width = 700,
        height = 500,
        filename = gsub(x = stringr::str_glue(plot_title), " ", "_")
      )
    )
  })


  # Carnival related content ---------------------------------------------------
  output$visnet_carnival <- renderVisNetwork({
    g <- carnival %>%
      filter(contrast == input$selected_contrast) %>%
      pull(igraph) %>%
      jsonlite::unserializeJSON() %>%
      igraph::upgrade_graph(.) %>%
      permute.vertices(., Matrix::invPerm(order(V(.)$name)))
    V(g)$title <- V(g)$name


    visNetwork::visIgraph(g) %>%
      visNodes(font = list(background = "white")) %>%
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
        type = "pdf",
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
          uiOutput("ui_bookmarks_genes"),
          id = "bookmarks_genes"
        ),
        column(
          width = 6,
          h5("Bookmarked gene sets"),
          uiOutput("ui_bookmarks_genesets"),
          id = "bookmarks_genesets"
        )
      )
    )
  })

  # RBP related content --------------------------------------------------------
  rvalues$rbp_table <- reactive({
    dtu_sig <- rvalues$key() %>%
      filter(dtu_pvadj < 0.05) %>%
      pull(SYMBOL)


    rbp %>%
      group_by(gene_id_regulator) %>%
      mutate(FDR = p.adjust(Pvalue, method = "fdr")) %>%
      filter(FDR <= 0.05 & gene_name %in% dtu_sig) %>%
      ungroup()
  })

  output$rbp_network <- renderVisNetwork({
    i <- getReactableState("rbp_table", "selected")
    validate(
      need(!is.na(i),
        message = "Please select an entry from the table to display the regulator subgraph."
      )
    )

    df <- rvalues$rbp_table()
    g <- create_graph_rbp(df)

    g %>%
      make_ego_graph(
        .,
        order = 2,
        nodes = which(V(g)$name == df[[i, "gene_name_regulator"]])
      ) %>%
      magrittr::extract2(1) %>%
      visIgraph() %>%
      visIgraphLayout() %>%
      visGroups(
        groupname = "regulator",
        color = "#EBECF0",
        shape = "square"
      ) %>%
      visGroups(
        groupname = "target",
        color = "#E5C494",
        shape = "ellipse"
      ) %>%
      visLegend(
        width = 0.1,
        position = "right",
        main = ""
      ) %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1
        ),
        nodesIdSelection = TRUE
      ) %>%
      visInteraction(navigationButtons = TRUE) %>%
      visExport(
        name = "rbp_network",
        type = "pdf",
        label = "Save graph"
      )
  })

  output$rbp_gene <- renderPlot({
    validate(
      need(input$rbp_network_selected != "", message = "Select a node in the RBP graph")
    )

    g <- rvalues$rbp_graph()
    cur_sel <- input$rbp_network_selected
    cur_node <- match(cur_sel, V(g)$name)
    cur_nodetype <- V(g)$nodetype[cur_node]


    annotation_obj <- rvalues$annotation_obj()


    plot(1, 1, main = paste(cur_sel, cur_node, cur_nodetype))
  })

  output$rbp_dtu <- renderPrint({
    validate(
      need(input$rbp_network_selected != "", message = "Select a node in the RBP graph")
    )

    g <- rvalues$rbp_graph()
    cur_sel <- input$rbp_network_selected
    cur_node <- match(cur_sel, V(g)$name)
    cur_nodetype <- V(g)$nodetype[cur_node]


    annotation_obj <- rvalues$annotation_obj()


    paste("Select some DTU-centered content for", cur_sel, cur_node, cur_nodetype)
  })

  output$rbp_table <- renderReactable({
    rvalues$rbp_table() %>%
      mutate(
        FDR = round(-log10(FDR), 2),
        Association = ifelse(Association == 1, "positive", "negative")
      ) %>%
      select(
        gene_name_regulator,
        gene_id_regulator,
        transcript_name,
        transcript_id,
        transcript_biotype,
        FDR,
        Association
      ) %>%
      reactable(
        .,
        language = reactableLang(
          filterPlaceholder = "Filter"
        ),
        filterable = TRUE,
        striped = TRUE,
        defaultSorted = c("FDR"),
        showPageSizeOptions = TRUE,
        defaultPageSize = 5,
        pageSizeOptions = c(5, 10, 25, 50),
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
        columnGroups = list(
          colGroup(name = "Regulator", columns = c("gene_name_regulator", "gene_id_regulator")),
          colGroup(name = "Target", columns = c("transcript_name", "transcript_biotype"))
        ),
        columns = list(
          gene_name_regulator = colDef(
            name = "gene_name",
            html = TRUE,
            header = with_tooltip("gene_name", "Regulator gene name")
          ),
          gene_id_regulator = colDef(
            name = "gene_id",
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 140,
            header = with_tooltip("gene_id", "Link to Ensembl gene page")
          ),
          transcript_id = colDef(
            name = "transcript_id",
            html = TRUE,
            cell = JS("function(cellInfo) {
              const url = 'https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=' + cellInfo.value
              return '<a href=\"' + url + '\" target=\"_blank\">' + cellInfo.value + '</a>'
            }"),
            minWidth = 140,
            header = with_tooltip("transcript_id", "Link to Ensembl transcript page")
          ),
          transcript_name = colDef(
            name = "transcript_name",
            header = with_tooltip("transcript_name", "Target transcript name")
          ),
          transcript_biotype = colDef(
            name = "transcript_biotype",
            header = with_tooltip("transcript_biotype", "Biotype of the target transcript")
          ),
          FDR = colDef(
            name = "-log10_FDR",
            defaultSortOrder = "desc",
            header = with_tooltip("-log10_FDR", "FDR adjusted p-value for the reverse global test"),
            filterable = TRUE,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] >= filterValue
                })
                }"),
          ),
          Association = colDef(
            name = "association",
            header = with_tooltip("association", "Direction of change in transcript expression.")
          )
        ),
        defaultColDef = colDef(sortNALast = TRUE)
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
      tags$button(
        "Download as CSV",
        onclick = "Reactable.downloadDataCSV('bookmarks_genes', 'bookmarks_genes.csv')"
      )
    )
  })

  output$ui_bookmarks_genesets <- renderUI({
    validate(
      need(
        nrow(rvalues$mygenesets) > 0,
        "Please select at least one gene set with the Bookmark button"
      )
    )
    tagList(
      reactableOutput("bookmarks_genesets"),
      tags$button(
        "Download as CSV",
        onclick = "Reactable.downloadDataCSV('bookmarks_genesets', 'bookmarks_genesets.csv')"
      )
    )
  })

  output$bookmarks_genes <- renderReactable({
    genes_tbl <- rvalues$key() %>%
      mutate_at(vars(log2FoldChange, dtu_dif), ~ round(., 2)) %>%
      mutate_at(vars(padj, dtu_pvadj), ~ round(-log10(.), 2))

    genes_tbl[match(rvalues$mygenes$gene_id, genes_tbl$gene_id), ] %>%
      reactable(
        .,
        columns = list(
          gene_id = colDef(
            name = "gene_id",
          ),
          SYMBOL = colDef(
            name = "gene_name",
          ),
          log2FoldChange = colDef(
            name = "dge_log2fc",
          ),
          padj = colDef(
            name = "-log10_dge_padj",
          ),
          dtu_pvadj = colDef(
            name = "-log10_dtu_padj",
          ),
          dtu_dif = colDef(
            name = "dtu_dif",
          )
        ),
        defaultColDef = colDef(sortNALast = TRUE)
      )
  })

  output$bookmarks_genesets <- renderReactable({
    rvalues$res_enrich() %>%
      filter(gs_id %in% rvalues$mygenesets$gs_id) %>%
      select(-row_names) %>%
      reactable(rownames = FALSE)
  })

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
          rename(gene_name = SYMBOL)

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
            showNotification(sprintf("The selected gene set, %s (%s), is already in the set of the bookmarked genes ets.", sel_gs, sel_gs_id), type = "default")
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
    } else if (input$magnetique_tab == "tab-carnival") {
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

  observeEvent(input$tour_rbpview, {
    tour <- read.delim("tours/intro_rbpview.txt",
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
        title = "Project members",
        size = "l",
        renderTable(make_team_df(),
          sanitize.text.function = function(x) x
        ),
        easyClose = TRUE,
        footer = ""
      )
    )
  })
}
shinyApp(magnetique_ui, magnetique_server)
