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
      "https://twitter.com/fedebioinfo",
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
gene_counts <- function(gene, dataset = se) {
  ind <- which(rowData(dataset)[, "gene_id"] == gene)

  x <- SummarizedExperiment::assays(dataset)[["counts"]][ind, ]
  x["transcript_id"] <- rownames(x)

  x <- suppressMessages(reshape2::melt(x))

  colnames(x) <- c("transcript_id", "sample_id", "value")

  x["group"] <- SummarizedExperiment::colData(dataset)[
    x[["sample_id"]], "Etiology"
  ]

  return(x)
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

#' Get gene expression for `gene`
#' @param gene selection
#' @param dataset SummarizedExperiment object with the dataset
#' @param gtf GRanges with gene structure
#' @export
gene_structure <- function(gene, gtf, dataset) {
  x <- BiocGenerics::subset(gtf, gene_id %in% c(gene))
  x <- BiocGenerics::subset(x, transcript_id %in% rownames(dataset))
  x$type <- plyr::revalue(
    x$type, c("five_prime_utr" = "utr", "three_prime_utr" = "utr")
  )
  x <- split(x, x$transcript_id)
  return(x)
}

#' Get gene expression for `gene`
#' @param gene selection
#' @param dataset SummarizedExperiment object with the dataset
#' @param gtf GRanges with gene structure
#' @export
plot_dtu <- function(gene, dataset, gtf) {
  library(patchwork)
  geneid2name <- setNames(
    nm = gtf$gene_id,
    gtf$gene_name
  )

  x_tx <- gene_proportions(gene, dataset)
  x_structure <- gene_structure(gene, gtf, dataset)

  p2 <- x_tx %>%
    ggplot() +
    geom_jitter(aes(x = transcript_id, y = value, color = group),
      position = position_jitterdodge(),
      alpha = 0.9, size = 2, show.legend = T, na.rm = TRUE
    ) +
    geom_boxplot(aes(x = transcript_id, y = value, fill = group),
      outlier.size = 0, alpha = 0.4, lwd = 0.5, show.legend = F
    ) +
    scale_fill_manual(name = "Etiology", values = group_colors) +
    scale_colour_manual(name = "Etiology", values = group_colors) +
    coord_flip() +
    labs(y = "transcript proportion") + 
    theme_minimal(20) +
    theme(
      axis.title.y = element_blank()
    )
  

  p3 <- suppressMessages({
    ggplot() +
    ggbio::geom_alignment(x_structure,
      fill = "black", cds.rect.h = .3, utr.rect.h = .2,
      exon.rect.h = .2, label = T
    ) +
    theme_minimal(20) +
    theme(plot.margin = margin())
    })

  p_final <- (p3 | p2) +
    plot_annotation(
      theme = theme(plot.margin = margin()),
      title = stringr::str_glue("{geneid2name[gene]} - {gene}"),
      caption = "Source: MAGNet / GRCh38.96"
    )
  return(p_final)
}
