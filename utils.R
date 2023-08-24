get_group_colors <- function(){
  group_colors <- c("steelblue", "gold", "forestgreen")
  names(group_colors) <- c("NFD",  "DCM", "HCM")
  group_colors  
}

group_colors <- get_group_colors()

# from https://glin.github.io/reactable/articles/cookbook/cookbook.html#tooltips
with_tooltip <- function(value, tooltip) {
  tags$abbr(style = "text-decoration: underline; text-decoration-style: dotted; cursor: help",
            title = tooltip, value)
}

#' Compute the mean proportion difference between groups
get_gid2name <- function(gtf) {
  setNames(
    nm = gtf$gene_id,
    gtf$gene_name)
}


#' Produce table with team info
make_team_df <- function() {
  df <- data.frame(
    name = c(
      "Annekathrin Ludt",
      "Christoph Dieterich",
      "Enio Gjerga",
      "Etienne Boileau",
      "Federico Marini",
      "Thiago Britto-Borges"),
    twitter = c(
      "https://twitter.com/AnnekathrinLudt",
      "",
      "https://twitter.com/e_ni_o",
      "",
      "https://twitter.com/FedeBioinfo",
      "https://twitter.com/tbrittoborges"),
    orcid = c(
      "https://orcid.org/0000-0002-2475-4945",
      "https://orcid.org/0000-0001-9468-6311",
      "https://orcid.org/0000-0001-8042-0395",
      "https://orcid.org/0000-0001-9355-0973",
      "https://orcid.org/0000-0003-3252-7758",
      "https://orcid.org/0000-0002-6218-4429")
  )
  df$twitter <- paste0("<a href='", df$twitter, "' target='_blank'>", df$twitter, "</a>")
  df$orcid <- paste0("<a href='", df$orcid, "' target='_blank'>", df$orcid, "</a>")
  return(df)
  
  }



#' Compute the mean proportion difference between groups
#' @param gene selection
#' @param dataset SummarizedExperiment object with the dataset
#' @param .type either `proportions` or  `fit_full` for fitted proportions
#' @export
compute_usage_dif <- function(gene, dataset, .type = "fit_full") {
  stop("Not implemented", call. = FALSE)
}

#' Get gene expression for `gene`
#' @param gene selection
#' @param dataset SummarizedExperiment object with the dataset
#' @export
gene_proportions <- function(gene, dataset, .type = "fit_full") {
  ind <- which(rowData(dataset)[, "gene_id"] == gene)
  x <- SummarizedExperiment::assays(dataset)[[.type]][ind, ]
  stopifnot(all(rownames(
    SummarizedExperiment::colData(dataset)
  ) == colnames(x)))

  col_names <- rownames(x)
  x <- suppressMessages(reshape2::melt(x))
  colnames(x) <- c("transcript_id", "sample_id", "value")
  x["group"] <- SummarizedExperiment::colData(dataset)[
    x[["sample_id"]], "Etiology"
  ]
  return(x)
}

plot_gene_structure <- function(gtf) {
    gtf$type <- plyr::revalue(
      gtf$type, c("five_prime_utr" = "utr", "three_prime_utr" = "utr")
    )
    p <- ggplot() +
      ggbio::geom_alignment(gtf,
        fill = "black", cds.rect.h = .3, utr.rect.h = .2,
        exon.rect.h = .2, label = T
      ) +
      theme_minimal(20) +
      theme(plot.margin = margin())
    return(p)
}

create_graph_rbp <- function(tbl_rbp) {
    
  df <- tbl_rbp %>%
   select(gene_name_regulator, transcript_name, everything())

  g <- graph.data.frame(df, directed = TRUE) 
  g <- permute.vertices(g, Matrix::invPerm(order(V(g)$name)))
  rbp_graph_color <- "gold"
  V(g)$group <- ifelse(grepl("-", names(V(g))), 'target', 'regulator')
  V(g)$title <- V(g)$name
    
  return(g)
}

createLinkGO <- function(val) {
  sprintf(
    "<a href=\"http://amigo.geneontology.org/amigo/term/%s\" target=\"_blank\" class=\"btn btn-primary\">%s</a>",
    val, val
  )
}

.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"

# Launching magnetique! --------------------------------------------------------
#' Build up a GeneTonicList, from the magnetique DB 
#'
#' @param con The DB connection (an SQLiteConnection object)
#' @param contrast The contrast, as specified e.g. in the app
#' @param ontology The ontology to focus upon (BP, MF, CC)
#' @param verbose Logical, whether to display messages while constructing
#'
#' @return A GeneTonicList object, to be used in concert with GeneTonic's function
#' @export
#'
#' @examples
#' mygtl <- buildup_gtl(con, "DCMvsHCM", "BP")
buildup_gtl <- function(con,
                        contrast,
                        ontology,
                        verbose = TRUE) {
  
  coldata <- metadata

  if(verbose) message("... building annotation...")
  annotation <- annotation_obj %>% 
    select(c("gene_id", "gene_name")) %>% as.data.frame()
  rownames(annotation) <- annotation$gene_id 

  if(verbose) message("Done!")
  
  if(verbose) message("... building counts...")
  
  counts <- as.data.frame(counts)
  counts_rownames <- counts$row_names
  counts$row_names <- NULL
  counts$contrast <- NULL
  rownames(counts) <- counts_rownames

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~Etiology + Race + Sex + Age + SV1 + SV2)
  
  if(verbose) message("Done!")
  
  if(verbose) message("... building DE table...")
  tbl_de <- res[[contrast]]
  tbl_de <- S4Vectors::DataFrame(tbl_de)
  rownames(tbl_de) <- rownames(annotation)
  
  res_de <- DESeq2::DESeqResults(tbl_de)
  if(verbose) message("Done!")
  
  if(verbose) message("... building enrichment table...")
  tbl_enrich <- res_enrich[[contrast]] %>%
    filter(ontology == ontology) %>%
    select(-ontology) %>%
    as.data.frame()
  rownames(tbl_enrich) <- tbl_enrich$gs_id
  
  res_enrich <- tbl_enrich
  if(verbose) message("Done!")
  
  gtl <- GeneTonic::GeneTonic_list(
    dds = dds,
    res_de = res_de,
    res_enrich = res_enrich,
    annotation_obj = annotation
  )
  return(gtl)
} 

highlight_selected <- function(selected, nelements) {     
    if(!is.null(selected)){
      colors <- rep('rgba(0,0,0,.10)', nelements)
      colors[selected] <- 'rgba(255,0,0,1)'
    } else {
      colors <- 'rgb(0,0,0)'
    }
    return(colors)
}