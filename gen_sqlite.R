#! /usr/bin/env Rscript

# Generate sqlite DB to match app.R 
# See CHANGELOG.

# NOTE: The nomenclature was changed in some branches. I follow recent changes to avoid breaking the app, but 
# we should avoid changing the names/variables! e.g. in main we use the GeneTonic nomenclature (dds, res_de, etc.),
# but here this is changed to counts, and res.

# NOTE: data downloaded and unpacked under local
# wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
# unzip -o data.zip
# rm data.zip
setwd('~/repos/magnetique/')
library(dplyr)

db_name <- "magnetique.sqlite"

# I/O ------------------------------------------------------

dirloc <- file.path("MAGNetApp", "data", fsep=.Platform$file.sep)

# load DGE results
cdir <- file.path(dirloc, "DGE", fsep=.Platform$file.sep)
dge <- purrr::map(dir(cdir, "*.rds"), ~readRDS(file.path(cdir, .x)))
names(dge) <-  sapply(strsplit(dir(cdir, "*.rds"), "_", fixed = T), "[[", 2)

# DTU results
dtu <- readRDS(file.path(dirloc, "DTU", "summarized_experiment.RDS", fsep=.Platform$file.sep))

# annotation
gtf <- readRDS(file.path(dirloc, "DTU", "gtf.RDS", fsep=.Platform$file.sep))

# TODO: carnival
cdir <- file.path(dirloc, "networks", fsep=.Platform$file.sep)
carnival <- purrr::map(dir(cdir, "*.RData"), ~get(load(file.path(cdir, .x))))

get_names <- function(files) {
  l <- strsplit(files, "_", fixed = T)
  lapply(l, function(.x) {
    paste(toupper(unlist(.x)[c(2,4)]), collapse='vs')
  })
}
names(carnival) <- get_names(dir(cdir, "*.RData"))


# colData.txt
coldata <- read.csv(file.path(dirloc, "colData.txt", fsep=.Platform$file.sep))


# Create DB ------------------------------------------------

# Generate contrast-specific tables + long tables (contrast as column) i
# case we want to use one table per analysis, and select the contrast on the fly.

library(DBI)
db <- dbConnect(RSQLite::SQLite(), db_name)

# write separate dds counts...
wrt <- purrr::map2(dge, names(dge), function(.x, .y) {
  DESeq2::counts(.x$dds) %>% 
    data.frame() %>% 
    dbWriteTable(db, paste("counts", .y, sep="_"), ., overwrite=TRUE, row.names=TRUE)
})

# ... and a long table
wrt <- purrr::map2(dge, names(dge), function(.x, .y) {
  DESeq2::counts(.x$dds) %>% 
    data.frame() %>% 
    mutate(contrast=.y)
})
do.call("rbind", wrt) %>% dbWriteTable(db, "counts", ., overwrite=TRUE, row.names=TRUE)

message('Loading VST-transformed counts into DB')
dds_object <- dge$DCMvsNFD$dds
transformed <- DESeq2::vst(dds_object, blind=FALSE)
vst <- SummarizedExperiment::assay(transformed)
vst <- as_tibble(vst, rownames='row_names')
dbWriteTable(db, "vst",vst , overwrite=TRUE)

# write DGE results, but first add DTU

# DRIMSeq can generate a single p-value per gene, which tests whether there is any differential transcript usage within the gene, but
# I don't know how this differs from the stageR procedure. Anyway, Thiago advised to use the min. transcript p-value per gene, and the 
# corresponding usage for this transcript.

rowdata <- SummarizedExperiment::rowData(dtu)
contrasts <- names(rowdata@listData)[grep("^DRIMSeq", names(rowdata@listData))]
stopifnot(all(rownames(SummarizedExperiment::colData(dtu)) == coldata$Run))

compute_usage_dif <- function(transcript, dataset, c, .type = "fit_full") {
  ind <- which(SummarizedExperiment::rowData(dataset)[, "transcript_id"] == transcript)
  x <- SummarizedExperiment::assays(dataset)[[.type]][ind, ]
  stopifnot(all(coldata$Run == names(x)))
  x <- x %>% as.numeric()
  mean(x[coldata$Etiology==c[1]])-mean(x[coldata$Etiology==c[2]])
}


wrt <- purrr::map(contrasts, function(.x) {
  contrast <- unlist(sapply(strsplit(.x, "_", fixed = T), "[", c(2, 4)))
  res_dtu <- rowdata[[.x]] %>% 
    arrange(adj_pvalue) %>%
    group_by(gene_id) %>% 
    mutate(dtu_dif=compute_usage_dif(feature_id, dtu, contrast)) %>%
    summarise(
      n_transcript = n(), 
      adj_pvalue = first(adj_pvalue),
      dtu_dif = first(dtu_dif)) %>%
    dplyr::rename(dtu_pvadj=adj_pvalue) %>%
    dplyr::select(gene_id, n_transcript, dtu_pvadj, dtu_dif)
  
  res_de <- dge[[paste(contrast, collapse='vs')]]$res_de %>% as.data.frame() %>%
    dplyr::select(log2FoldChange, padj, SYMBOL) %>% tibble::rownames_to_column("gene_id")
  res_de <- merge(res_de, res_dtu, all.x=T, by='gene_id')
  res_de
})
names(wrt) <- get_names(contrasts)
tmp <- purrr::map2(wrt, names(wrt), function(.x, .y) {
  .x %>% dbWriteTable(db, paste("res", .y, sep="_"), ., overwrite=TRUE, row.names=TRUE)
})

# ... also add one long table
wrt <- purrr::map2(wrt, names(wrt), function(.x, .y) {
  .x %>% mutate(contrast=.y)
})
do.call("rbind", wrt) %>% dbWriteTable(db, "res", ., overwrite=TRUE, row.names=TRUE)


# DRIMSeq transcript proportions
dtu@assays@data$fit_full %>%
  data.frame() %>% 
  dbWriteTable(db, "dtu_fit_proportions", .,  overwrite=TRUE, row.names=TRUE)

# Loads the gene ontology for Human
library(org.Hs.eg.db)

hs_go <-  toTable(org.Hs.egGO)
# filter terms with more than 300 genes
hs_go <- hs_go %>% 
  group_by(go_id) %>% 
  summarise(n_gene_p_term = n_distinct(gene_id)) %>%
  filter(n_gene_p_term < 300)

# enrichment results per contrast...
wrt <- purrr::map2(dge, names(dge), function(.x, .y) {
  do.call("rbind", purrr::map2(.x$res_enrich, names(.x$res_enrich), function(.xi, .yi) {
    .xi %>% data.frame() %>% 
      mutate(ontology=.yi) %>% 
      filter(gs_id %in% hs_go$go_id)
  })) %>%
    dbWriteTable(db, paste("res_enrich", .y, sep="_"), ., overwrite=TRUE, row.names=TRUE)
})

# ... and in one long table  
wrt <- purrr::map2(dge, names(dge), function(.x, .y) {
  do.call("rbind", purrr::map2(.x$res_enrich, names(.x$res_enrich), function(.xi, .yi) {
    .xi %>% data.frame() %>% 
      mutate(ontology=.yi)
      
  })) %>% 
    mutate(contrast=.y) %>% 
    filter(gs_id %in% hs_go$go_id)
})
do.call("rbind", wrt) %>% dbWriteTable(db, "res_enrich", ., overwrite=TRUE, row.names=TRUE)

# metadata
dbWriteTable(db, 'metadata', coldata, ., overwrite=TRUE, row.names=TRUE)

# gtf gene to transcript
tx <- gtf[gtf$type == "transcript"]
gene2tx <- as.data.frame(
  S4Vectors::mcols(tx)[
    , c("gene_id", "transcript_id")
  ]
)
dbWriteTable(db, 'gene2tx', gene2tx, overwrite=TRUE, row.names=FALSE)
gtf %>% data.frame() %>% dbWriteTable(db, 'gtf', ., overwrite=TRUE, row.names=FALSE)

# annotation_obj - pick one
dge$DCMvsNFD$annotation_obj %>%  dbWriteTable(db, "annotation_obj", ., overwrite=TRUE, row.names=TRUE)

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

dbWriteTable(db, "carnival", carnival, overwrite=TRUE)

dbWriteTable(db, "rbp", rbp, overwrite=TRUE)

dbDisconnect(db)

base_url <- "https://data.dieterichlab.org/public.php/webdav/"
library(httr)
PUT(
  paste0(base_url, basename(db_name)), 
  authenticate('3gCTLGT4DaAqaqb', ''), 
  body = upload_file(db_name),
  add_headers('X-Requested-With' = 'XMLHttpRequest'))
