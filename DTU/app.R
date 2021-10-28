library(shiny)

source("plots_with_se_obj.R")
se <- readRDS("summarized_experiment.RDS")

gene_names <- unique(rowData(se)[["gene_name"]])


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
    req(input$gene_name)
    gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
    plot_dtu(mcols(gtf_gene)[["gene_id"]], se, gtf)
    
  })
  
  output$table <- renderTable({
    req(input$gene_name)
    gtf_gene <- subset(gtf, type == "gene" & gene_name == input$gene_name)
    results_table(mcols(gtf_gene)[["gene_id"]], se) 
    
  })
}
shinyApp(ui = ui, server = server)
