#' Plot a volcano plot of a geneset
#'
#' Plot a volcano plot for the geneset of the provided data, with the remaining
#' genes as shaded dots in the background of the plot.
#'
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to 0.05.
#' @param color Character string to specify color of filtered points in the plot.
#' Defaults to #1a81c2 (shade of blue).
#' @param volcano_labels Integer, maximum number of labels for the gene sets to be
#' plotted as labels on the volcano scatter plot. Defaults to 25.
#' @param plot_title Character string, to specify the title of the plot,
#' displayed over the volcano plot. If left to `NULL` as by default, it tries to use
#' the information on the geneset identifier provided.
#'
#' @return A plot returned by the [ggplot()] function
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#' library("apeglm")
#' library("ggplot2")
#' library("ggrepel")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'     keys = rownames(dds_macrophage),
#'     column = "SYMBOL",
#'     keytype = "ENSEMBL"
#'   ),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' # res object
#' dds_macrophage <- DESeq(dds_macrophage)
#' res_de <- results(dds_macrophage, name = "condition_IFNg_vs_naive", alpha = 0.05)
#' res_de <- lfcShrink(dds_macrophage, coef = "condition_IFNg_vs_naive", type = "apeglm", res = res_de) 
#'
#' volcano_plot(res_de,
#'   anno_df,
#' )
volcano_plot <- function(res_de,
                         annotation_obj,
                         FDR = 0.05,
                         color = "#1a81c2",
                         alpha = 0.10,
                         volcano_labels = 25,
                         plot_title = NULL) {
  
  # Prepare the data
  gene_ids <- rownames(res_de)
  gene_names <-
    annotation_obj$gene_name[match(gene_ids, annotation_obj$gene_id)]
  
  padj_complete <- res_de[gene_ids, "padj"]
  filter <- sapply(padj_complete, function(x) x <= FDR)
  padj_complete <- sapply(padj_complete, function(x) -log10(x))
  
  log2FoldChange_complete <- res_de[gene_ids, "log2FoldChange"]
  
  complete_data <- data.frame(
    gene_ids,
    padj_complete,
    log2FoldChange_complete,
    filter
  )
  
  colnames(complete_data) <- c(
    "genes",
    "logTransformedpvalue",
    "log2FoldChange",
    "significant"
  )
  
  
  # Prepare plotting
  volcano_df <- complete_data
  volcano_df$gene_names <- gene_names
  max_x <- max(abs(range(complete_data["log2FoldChange"])))
  limit_x <- max_x * c(-1, 1)
  
  
  # Prepare plot title
  if (is.null(plot_title)) {
    title <- paste0("Volcano Plot")
  } else {
    title <- plot_title
  }
  
  # handling the tooltips (works if plotlyfied)
  volcano_df$gene_info <- paste0(
    "<b>",volcano_df$gene_names, "</b>",
    "<br><i>GeneID</i>: ", volcano_df$genes,
    "<br><i>Log2FC</i> = ", format(round(volcano_df$log2FoldChange, 2), nsmall = 2),
    "<br><i>p-value (adjusted)</i> = ", format(res_de$padj))
  
  # Plot data
  p <- ggplot(
    volcano_df,
    aes_string(x = "log2FoldChange", 
               y = "logTransformedpvalue",
               text = "gene_info")
  ) +
    geom_point(aes_string(
      color = "significant",
    ), alpha=alpha) +
    labs(
      x = "log2 Fold Change",
      y = "-log10 p-value",
      color = paste0("pvalue <= ", FDR)
    ) +
    scale_x_continuous(limits = limit_x) +
    scale_color_manual(
      labels = c("significant", "not significant"),
      breaks = c("TRUE", "FALSE"),
      values = c(color, "grey25")
    ) + 
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold")
    )
  
  
  
  # adding labels to the significant points
  if (volcano_labels > 0) {
    p <- p + geom_text_repel(
      data = subset(volcano_df, filter),
      aes_string(label = "gene_names"),
      size = 4,
      max.overlaps = volcano_labels
    )
  }
  
  return(p)
}
