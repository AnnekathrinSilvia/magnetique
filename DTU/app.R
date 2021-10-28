library(shiny)


#' Summarized experiment for MAGnet dataset
#' Results from DRIMSeq
#' Each row is a transcript
#' There are three assays:
#' - counts: read counts per transcript
#' - proportions: proportions of reads per gene
#' - fit_full: fitted proportions
#' colData has the sample information
#' rowData has the transcript annotation
#' Modelling results are accesible via rowData:
#' rowData "DRIMSeq_HCM_vs_NFD"
#' rowData "DRIMSeq_DCM_vs_NFD"
#' rowData "DRIMSeq_DCM_vs_HCM"
se <- readRDS("/beegfs/prj/MAGE/analysis/DTU/summarized_experiment.RDS")

gene_names <- unique(rowData(se)[["gene_name"]])

source("../plots_with_se_obj.R")

ui <- fluidPage(
  titlePanel("MAGnet DTU prototype"),
  ui <- fluidPage(
    selectizeInput(
      "gene_name", "Choose one gene:", choices = NULL
    ),
    textOutput("result"),
    plotOutput("plot"),
    tableOutput("table")
  )
)

server <- function(input, output, session) {
  
  
  updateSelectizeInput(
    session, "gene_name",
    choices = gene_names, server = TRUE, selected='OGT'
  )

  output$plot <- renderPlot({
    if (!is.null(input$gene_name)){
        gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
        plot_dtu(mcols(gtf_gene)[["gene_id"]], se, gtf)
    }
  })
  
  output$table <- renderTable({
    if (!is.null(input$gene_name)){
      gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
     results_table(mcols(gtf_gene)[["gene_id"]], se) 
    }
  })
}
shinyApp(ui = ui, server = server)
