
library(SummarizedExperiment)

d <- readRDS("drimseq_age_sex_race_sva.RDS")

# DCMvsNFD
if (file.exists("DTU_mage_result_DCMvsNFD.RDS")) {
  d1 <- readRDS("DTU_mage_result_DCMvsNFD.RDS")
} else {
  d1 <- dmTest(d, coef = 2, verbose = 1, BPPARAM = BPPARAM)
  saveRDS(d1, "DTU_mage_result_DCMvsNFD.RDS")
}
# HCMvsNFD
if (file.exists("DTU_mage_result_HCMvsNFD.RDS")) {
  d2 <- readRDS("DTU_mage_result_HCMvsNFD.RDS")
} else {
  d2 <- dmTest(d, coef = 3, verbose = 1, BPPARAM = BPPARAM)
  saveRDS(d2, "DTU_mage_result_HCMvsNFD.RDS")
}
# DCMvsHCM
if (file.exists("DTU_mage_result_DCMvsHCM.RDS")) {
  d3 <- readRDS("DTU_mage_result_DCMvsHCM.RDS")
} else {
  d3 <- dmTest(d, contrast = c(0, 1, -1, 0, 0, 0, 0, 0), verbose = 1, BPPARAM = BPPARAM)
  saveRDS(d3, "DTU_mage_result_DCMvsHCM.RDS")
}

gtf <- rtracklayer::import.gff2("/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38.102.gtf")
tx <- gtf[gtf$type == "transcript"]
tx2gene <- as.data.frame(
  mcols(tx)[
    , c("gene_name", "gene_id", "gene_biotype", "transcript_name", "transcript_id", "transcript_biotype")
  ]
)


row_data <- tx2gene[ match(proportions(d)[, 2], tx2gene$transcript_id), ]
rownames(row_data) <- proportions(d)[, 2]

se <- SummarizedExperiment(
  assays = list(
    counts = counts(d)[3: length(counts(d))],
    proportions = proportions(d)[3: length(counts(d))],
    fit_full = d@fit_full@unlistData
  ),
  rowData = row_data,
  colData = d@samples
)

rowData(se)[["DRIMSeq_HCM_vs_NFD"]] <- DRIMSeq::results(d1, level = "feature")
rowData(se)[["DRIMSeq_DCM_vs_NFD"]] <- DRIMSeq::results(d2, level = "feature")
rowData(se)[["DRIMSeq_DCM_vs_HCM"]] <- DRIMSeq::results(d3, level = "feature")

saveRDS(se, "summarized_experiment.RDS")
