system(
  "wget --no-clobber --output-document  data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive' ")
system("unzip -o data.zip")

library(dplyr)
dirloc <- file.path("MAGNetApp", "data", fsep=.Platform$file.sep)
# load DGE results
cdir <- file.path(dirloc, "DGE", fsep=.Platform$file.sep)
dge <- purrr::map(dir(cdir, "*.rds"), ~readRDS(file.path(cdir, .x)))
names(dge) <- lapply(dge, \(x) x$project_id) %>% unlist()

# DTU results
dtu <- readRDS(file.path(dirloc, "DTU", "summarized_experiment.RDS", fsep=.Platform$file.sep))

# annotation
gtf <- readRDS(file.path(dirloc, "DTU", "gtf.RDS", fsep=.Platform$file.sep))

# carnival
cdir <- file.path(dirloc, "networks", fsep=.Platform$file.sep)
carnival <- purrr::map(dir(cdir, "*.RData"), ~get(load(file.path(cdir, .x))))

get_names <- function(files) {
  l <- strsplit(files, "_", fixed = T)
  lapply(l, function(.x) {
    paste(toupper(unlist(.x)[c(2,4)]), collapse='vs')
  })
}
names(carnival) <- get_names(dir(cdir, "*.RData"))

# metadata
metadata <- read.csv(file.path(dirloc, "colData.txt", fsep=.Platform$file.sep))

## Preparing tables ####
counts <- dge$DCMvsHCM$dds %>%
  DESeq2::counts(.) %>%
  as_tibble(., rownames = "row_names")

message('Loading VST-transformed')
vst <- dge$DCMvsHCM$dds %>% 
  DESeq2::vst(., blind=FALSE) %>% 
  SummarizedExperiment::assay(.)

rowdata <- SummarizedExperiment::rowData(dtu)
contrasts <- names(rowdata@listData)[grep("^DRIMSeq", names(rowdata@listData))]
stopifnot(all(rownames(SummarizedExperiment::colData(dtu)) == metadata$Run))

compute_usage_dif <- function(transcript, dataset, c, .type = "fit_full") {
  ind <- which(SummarizedExperiment::rowData(dataset)[, "transcript_id"] == transcript)
  x <- SummarizedExperiment::assays(dataset)[[.type]][ind, ]
  stopifnot(all(metadata$Run == names(x)))
  x <- x %>% as.numeric()
  mean(x[metadata$Etiology==c[1]])-mean(x[metadata$Etiology==c[2]])
}

res <- purrr::map(contrasts, function(.x) {
  contrast <- unlist(sapply(strsplit(.x, "_", fixed = T), "[", c(2, 4)))
  res_dtu <- rowdata[[.x]] %>% 
    arrange(adj_pvalue) %>%
    group_by(gene_id) %>% 
    mutate(dtu_dif = compute_usage_dif(feature_id, dtu, contrast)) %>%
    summarise(
      n_transcript = n(),
      adj_pvalue = dplyr::first(adj_pvalue),
      dtu_dif = dplyr::first(dtu_dif)) %>%
    dplyr::rename(dtu_pvadj = adj_pvalue) %>%
    dplyr::select(gene_id, n_transcript, dtu_pvadj, dtu_dif)
  
  res_de <- dge[[paste(contrast, collapse='vs')]]$res_de %>% as.data.frame() %>%
    dplyr::select(log2FoldChange, padj, SYMBOL) %>% tibble::rownames_to_column("gene_id")
  res_de <- merge(res_de, res_dtu, all.x=T, by='gene_id')
  res_de
})
names(res) <- get_names(contrasts)

# DRIMSeq transcript proportions
dtu_fit_proportions <- dtu@assays@data$fit_full %>%
  data.frame() 

# Loads the gene ontology for Human
library("org.Hs.eg.db")

hs_go <-  toTable(org.Hs.egGO)
hs_go <- hs_go %>%
  group_by(go_id) %>%
  summarise(n_gene_p_term = n_distinct(gene_id)) %>%
  filter(n_gene_p_term < 300)

# enrichment results per contrast...
res_enrich <- purrr::map2(dge, names(dge), function(.x, .y) {
  do.call("rbind", purrr::map2(.x$res_enrich, names(.x$res_enrich), function(.xi, .yi) {
    .xi %>% data.frame() %>% 
      mutate(ontology=.yi) %>% 
      filter(gs_id %in% hs_go$go_id)
  }))
})

# gtf gene to transcript
gene2tx <- gtf[gtf$type == "transcript"] %>% 
  S4Vectors::mcols(.) %>% 
  as.data.frame() %>% 
  dplyr::select(gene_id, transcript_id)

# annotation 
annotation_obj <- dge$DCMvsNFD$annotation_obj

# RBPs
rbp <- read.csv(file.path(dirloc, "RBP", "mirna_mrna_revgt_interactions.csv", fsep=.Platform$file.sep))

# Carnival data
carnival <- list()

load(file.path(dirloc, "networks", "igraph_dcm_vs_hcm_hierarchic.RData"))
carnival[["DCMvsHCM"]] <-  jsonlite::serializeJSON(gg)

load(file.path(dirloc, "networks", "igraph_dcm_vs_nfd_hierarchic.RData"))
carnival[["DCMvsNFD"]] <-  jsonlite::serializeJSON(gg)

load(file.path(dirloc, "networks", "igraph_hcm_vs_nfd_hierarchic.RData"))
carnival[["HCMvsNFD"]] <-  jsonlite::serializeJSON(gg)

carnival <- data.frame(
  contrast = names(carnival),
  igraph = as.character(carnival)
)

system("rm data.zip")
system("rm -rf MAGNetApp/")