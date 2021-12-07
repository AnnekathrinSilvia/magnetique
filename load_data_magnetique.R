# load in main objects ----------------------------------------------------

MAGNet_DCMvsHCM_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsHCM_GeneTonic.rds")
MAGNet_DCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsNFD_GeneTonic.rds")
MAGNet_HCMvsNFD_GeneTonic <- readRDS("MAGNetApp_cloud_data/data/DGE/MAGNet_HCMvsNFD_GeneTonic.rds")
MAGNet_DCMvsHCM_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_hcm.RData"))
MAGNet_DCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_nfd.RData"))
MAGNet_HCMvsNFD_igraph <- get(load("MAGNetApp_cloud_data/data/networks/igraph_hcm_vs_nfd.RData"))


# DCM vs HCM --------------------------------------------------------------
gtl_DCMvsHCM_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)

gtl_DCMvsHCM_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsHCM_GeneTonic$dds,
  res_de = MAGNet_DCMvsHCM_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsHCM_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsHCM_GeneTonic$annotation_obj
)


# DCM vs NFD --------------------------------------------------------------
gtl_DCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)

gtl_DCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_DCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_DCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_DCMvsNFD_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_DCMvsNFD_GeneTonic$annotation_obj
)


# HCM vs NFD --------------------------------------------------------------
gtl_HCMvsNFD_BP <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$BP,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_MF <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$MF,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

gtl_HCMvsNFD_CC <- GeneTonic_list(
  dds = MAGNet_HCMvsNFD_GeneTonic$dds,
  res_de = MAGNet_HCMvsNFD_GeneTonic$res_de,
  res_enrich = MAGNet_HCMvsNFD_GeneTonic$res_enrich$CC,
  annotation_obj = MAGNet_HCMvsNFD_GeneTonic$annotation_obj
)

all_gtls <- list(
  DCMvsHCM = list(
    BP = gtl_DCMvsHCM_BP,
    MF = gtl_DCMvsHCM_MF,
    CC = gtl_DCMvsHCM_CC
  ),
  DCMvsNFD = list(
    BP = gtl_DCMvsNFD_BP,
    MF = gtl_DCMvsNFD_MF,
    CC = gtl_DCMvsNFD_CC
  ),
  HCMvsNFD = list(
    BP = gtl_HCMvsNFD_BP,
    MF = gtl_HCMvsNFD_MF,
    CC = gtl_HCMvsNFD_CC
  )
)

all_igraph <- list(DCMvsHCM = MAGNet_DCMvsHCM_igraph,
                   DCMvsNFD = MAGNet_DCMvsNFD_igraph,
                   HCMvsNFD = MAGNet_HCMvsNFD_igraph)

rm(MAGNet_DCMvsHCM_GeneTonic,
   MAGNet_DCMvsNFD_GeneTonic,
   MAGNet_HCMvsNFD_GeneTonic)

rm(gtl_DCMvsHCM_BP,
   gtl_DCMvsHCM_MF,
   gtl_DCMvsHCM_CC,
   gtl_DCMvsNFD_BP,
   gtl_DCMvsNFD_MF,
   gtl_DCMvsNFD_CC,
   gtl_HCMvsNFD_BP,
   gtl_HCMvsNFD_MF,
   gtl_HCMvsNFD_CC)

rm(MAGNet_DCMvsHCM_igraph,
   MAGNet_DCMvsNFD_igraph,
   MAGNet_HCMvsNFD_igraph)


# DTU data loading --------------------------------------------------------

se_dtu <- readRDS("MAGNetApp_cloud_data/data/summarized_experiment.RDS")
gtf <- rtracklayer::import.gff2("MAGNetApp_cloud_data/data/GRCh38.96.gtf.gz")


# WGCNA data loading ------------------------------------------------------


