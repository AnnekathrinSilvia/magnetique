#' Plot a heatmap of the gene signature on the data
#'
#' Plot a heatmap for the selected gene signature on the provided data, with the possibility to compactly display also DE only genes
#'
#' @param se A `SummarizedExperiment` object, or an object derived from this class,
#' such as a `DESeqTransform` object (variance stabilized transformed data, or
#' regularized logarithm transformed), in where the transformation has been applied
#' to make the data more homoscedastic and thus a better fit for visualization.
#' @param res_de A `DESeqResults` object.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param geneset_id Character specifying the gene set identifier to be plotted
#' @param genelist A vector of character strings, specifying the identifiers
#' contained in the row names of the `se` input object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to 0.05.
#' @param de_only Logical, whether to include only differentially expressed genes
#' in the plot
#' @param cluster_rows Logical, determining if rows should be clustered, as
#' specified by [pheatmap::pheatmap()]
#' @param cluster_columns Logical, determining if columns should be clustered, as
#' specified by [pheatmap::pheatmap()]
#' @param center_mean Logical, whether to perform mean centering on the row-wise
#' @param scale_row Logical, whether to standardize by row the expression values
#' @param anno_col_info A character vector of names in `colData(dds)` to use for
#' decorating the heatmap as annotation.
#' @param plot_title Character string, to specify the title of the plot,
#' displayed over the heatmap. If left to `NULL` as by default, it tries to use
#' the information on the geneset identifier provided
#' 
#' @return A plot returned by the [pheatmap::pheatmap()] function
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' vst_macrophage <- vst(dds_macrophage)
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
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#'
#' heatmap(vst_macrophage,
#'   res_de,
#'   res_enrich,
#'   anno_df,
#'   geneset_id = res_enrich$gs_id[1],
#'   cluster_columns = TRUE,
#'   anno_col_info = "condition"
#' )
heatmap <- function(se,
                    res_de,
                    res_enrich,
                    annotation_obj = NULL,
                    gtl = NULL,
                    geneset_id = NULL,
                    genelist = NULL,
                    FDR = 0.05,
                    de_only = FALSE,
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    center_mean = TRUE,
                    scale_row = FALSE,
                    anno_col_info = NULL,
                    plot_title = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  # check that the data would ideally be a DST, so that it is not the counts/normalized?
  mydata <- assay(se)

  if (!is.null(geneset_id)) {
    if (geneset_id %in% res_enrich[["gs_id"]]) {
      thisset_name <- res_enrich[geneset_id, "gs_description"]
      thisset_members <-
        unlist(strsplit(res_enrich[geneset_id, "gs_genes"], ","))
      thisset_members_ids <-
        annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]
    }
  } else {
    # overridable via a list
    if (!all(genelist %in% rownames(se))) {
      not_there <- genelist[!(genelist %in% rownames(se))]
      warning(
        "Some of the provided gene ids were not found in the SummarizedExperiment",
        "\nNot found: ",
        not_there
      )
    }
    thisset_members_ids <- intersect(rownames(se), genelist)
    thisset_name <- "Custom list"
  }
  
  sig_to_keep <- (thisset_members_ids %in% rownames(se)) #
  thisset_members_ids_available <- thisset_members_ids[sig_to_keep]
  
  mydata_sig <- mydata[thisset_members_ids_available, ]
  
  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove, ]
  
  hm_name <- "Expression \nvalues"
  
  if (center_mean) {
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)
    hm_name <- "Expression \nvalues"
  }
  
  if (scale_row) {
    mydata_sig <- t(scale(t(mydata_sig)))
    hm_name <- "Z-scores \nExpression \nvalues"
  }
  
  if (de_only) {
    de_res <- deseqresult2df(res_de, FDR)
    de_genes <- de_res$id
    de_to_keep <- rownames(mydata_sig) %in% de_genes
    mydata_sig <- mydata_sig[de_to_keep, , drop = FALSE]
  }
  
  if (is.null(plot_title)) {
    title <-
      paste0("Signature heatmap - ", thisset_name, " - ", geneset_id)
  } else {
    title <- plot_title
  }
  
  anno_col_info <- anno_col_info[anno_col_info %in% colnames(colData(se))]
  sample_decoration <- as.data.frame(colData(se))[, anno_col_info, drop = FALSE]
  anns_colors <- list(Etiology = c(DCM = "gold", HCM = "forestgreen", NFD = "steelblue"))

  
  p <- ComplexHeatmap::pheatmap(mydata_sig,
            cluster_rows = cluster_rows,
            cluster_cols = cluster_columns,
            scale = ifelse(scale_row, "row", "none"),
            main = title,
            labels_row = annotation_obj[rownames(mydata_sig), ]$gene_name,
            annotation_col = sample_decoration,
            annotation_colors = anns_colors)
}
