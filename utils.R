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
make_team_df <- data.frame(
    name = c(
      "Annekathrin Ludt",
      "Christoph Dieterich",
      "Enio Gjerga",
      "Etienne Boileau",
      "Federico Marini",
      "Thiago Britto-Borges"),
    twitter = c(
      "https://twitter.com/AnnekathrinLudt",
      NA,
      "https://twitter.com/e_ni_o",
      NA,
      "https://twitter.com/FedeBioinfo",
      "https://twitter.com/tbrittoborges"),
    orcid = c(
      "https://orcid.org/0000-0002-2475-4945",
      "https://orcid.org/0000-0001-9468-6311",
      "https://orcid.org/0000-0001-8042-0395",
      "https://orcid.org/0000-0001-9355-0973",
      "https://orcid.org/0000-0003-3252-7758",
      "https://orcid.org/0000-0002-8984-9084")
  )


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

rbp_results <- "MAGNetApp/data/RBP/mirna_mrna_revgt_interactions.csv"
create_graph_rbp <- function(rbp_results, filter_adj_pval = 0.03) {
  tbl_rbp <- read.csv(rbp_results)
  
  table(tbl_rbp$gene_name_regulator)
  
  # find out what is connected with what
  all_rbps <- sort(unique(tbl_rbp$gene_name_regulator))
  all_mrnas <- sort(unique(tbl_rbp$gene_name))
  
  table(tbl_rbp$gene_name_regulator == tbl_rbp$gene_name, useNA = "always")
  sum(is.na(tbl_rbp$gene_name_regulator == tbl_rbp$gene_name))
  
  tbl_rbp$adj_pval <- NA
  
  # adjusting the p-value, "by RBPs"
  for (i in all_rbps) {
    interactors_with_rbp <- (tbl_rbp$gene_name_regulator == i)
    this_rbp <- tbl_rbp[interactors_with_rbp, ]
    # message(nrow(this_rbp))
    tbl_rbp$adj_pval[interactors_with_rbp] <- p.adjust(this_rbp$Pvalue, method = "BH")
  }
  
  # plot(tbl_rbp$Pvalue, tbl_rbp$adj_pval)
  # might need all p-values to correct adequately? i.e. before filtering?
  
  # there are some NAs, removing them
  tbl_rbp <- na.omit(tbl_rbp)
  
  # all rbps are also mrnas
  intersect(all_rbps, all_mrnas)
  # some rbps act also on themselves?
  sum(tbl_rbp$gene_name_regulator == tbl_rbp$gene_name)
  
  # enhancing the info encoding in the table
  tbl_rbp$Association[tbl_rbp$Association == 1] <- "pos"
  tbl_rbp$Association[tbl_rbp$Association == 0] <- "neg"
  
  
  # filtering the interactions 
  tbl_rbp <- tbl_rbp[tbl_rbp$adj_pval <= filter_adj_pval, ]
  
  tbl_rbp <- dplyr::select(tbl_rbp, gene_name_regulator, gene_name, everything())
  
  g <- graph.data.frame(tbl_rbp, directed = FALSE)
  
  nodeIDs_rbps <- which(names(V(g)) %in% all_rbps)
  nodeIDs_genes <- which((names(V(g)) %in% setdiff(all_mrnas, all_rbps)))
  
  V(g)$nodetype <- NA
  V(g)$nodetype[nodeIDs_rbps] <- "RBP"
  V(g)$nodetype[nodeIDs_genes] <- "mRNA"
  
  rbp_graph_color = "gold"
  
  V(g)$shape <- c("box", "ellipse")[factor(V(g)$nodetype, levels = c("RBP", "mRNA"))]
  
    
  # to be reused if we want tooltips
  # V(g)$title <- NA
  # V(g)$title[nodeIDs_rbps] <- paste0(
  #   "<h4>",
  #   sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>', enriched_gsids[nodeIDs_rbps], enriched_gsids[nodeIDs_rbps]), "</h4><br>",
  #   V(g)$name[nodeIDs_rbps], "<br><br>",
  #   paste0(strwrap(enriched_gsdescs[nodeIDs_rbps], 50), collapse = "<br>")
  # )
  # V(g)$title[nodeIDs_genes] <- paste0(
  #   "<h4>", V(g)$name[nodeIDs_genes], "</h4><br>",
  #   "logFC = ", format(round(fcs_genes, 2), nsmall = 2)
  # )
  
  
  V(g)$color[nodeIDs_genes] <- "#EBECF0"
  V(g)$color[nodeIDs_rbps] <- "#E5C494"
  
  # re-sorting the vertices alphabetically
  rank_rbps <- rank(V(g)$name[V(g)$nodetype == "RBP"])
  rank_mrnas <- rank(V(g)$name[V(g)$nodetype == "mRNA"]) +
    length(rank_rbps) # to keep the RBPs first
  g <- permute.vertices(g, c(rank_rbps, rank_mrnas))
  
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
  
  coldata <- tbl(con, "metadata") %>% collect()

  if(verbose) message("... building annotation...")
  annotation <- tbl(con, "annotation_obj") %>% 
    select(c("gene_id", "gene_name")) %>% collect() %>% as.data.frame()
  rownames(annotation) <- annotation$gene_id 

  if(verbose) message("Done!")
  
  if(verbose) message("... building counts...")
  
counts <- tbl(con, "counts") %>%
  filter(contrast == contrast) %>% 
  mutate(row_names = str_replace_all(row_names, 'DCMvsHCM.', "")) %>% 
  collect()

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
  tbl_de <- tbl(con, paste0("res_", local(contrast))) %>% 
    collect()
  tbl_de <- S4Vectors::DataFrame(tbl_de)
  rownames(tbl_de) <- rownames(annotation)
  
  res_de <- DESeq2::DESeqResults(tbl_de)
  if(verbose) message("Done!")
  
  if(verbose) message("... building enrichment table...")
  tbl_enrich <- tbl(con, paste0("res_enrich_", local(contrast))) %>%
    filter(ontology == ontology) %>%
    select(-ontology) %>%
    collect() %>%
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

ledges <- data.frame(
  color = c("black", "red"),
  label = c("\nActivatory Interaction", "\nInhibitory Interaction"), 
  arrows =c("to", "to"))

lnodes <- data.frame(
  label = c("Down-regulation", "Up-regulation", "Pertubation", "Inferred Protein", "Transcription Factor"),
  shape = c("ellipse", "ellipse", "ellipse", "circle", 'square'), 
  color = c("#0000FF", "#FF0000",  "#00FF00", "lightgrey", "lightgrey"),
  id = 1:5)

highlight_selected <- function(selected, nelements) {     
    if(!is.null(selected)){
      colors <- rep('rgba(0,0,0,.10)', nelements)
      colors[selected] <- 'rgba(255,0,0,1)'
    } else {
      colors <- 'rgb(0,0,0)'
    }
    return(colors)
}