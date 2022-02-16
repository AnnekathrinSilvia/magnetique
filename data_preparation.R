prepare_data <- function(contrast, ontology = 'BP') {
  progress$set(value = .66, message = 'Loading.')
  data <- readRDS(paste0("MAGNetApp/", contrast, ".RDS"))
  progress$set(value = 1, message = 'Loading.')
  progress$close()
  return(
    list(
        igraph = data$igraph,
        res = data$res,
        gtf = gtf,
        res_dtu = res_dtu,
        genetonic = GeneTonic_list(
          dds = data$dds,
          res_de = data$res_de[order(data$res_de$padj), ],
          res_enrich = data$res_enrich[[ontology]],
          annotation_obj =data$annotation_obj
          )
        )
      )
    }