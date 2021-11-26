#!/bin/sh
mkdir -p MAGNetApp_cloud_data/data/DGE

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2FDGE&files=MAGNet_DCMvsHCM_GeneTonic.rds&downloadStartSecret=ddvmnsspk96" -O MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsHCM_GeneTonic.rds 

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2FDGE&files=MAGNet_DCMvsNFD_GeneTonic.rds&downloadStartSecret=dhjqyerzsi" -O MAGNetApp_cloud_data/data/DGE/MAGNet_DCMvsNFD_GeneTonic.rds

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2FDGE&files=MAGNet_HCMvsNFD_GeneTonic.rds&downloadStartSecret=jwfm2l43h1k" -O MAGNetApp_cloud_data/data/DGE/MAGNet_HCMvsNFD_GeneTonic.rds

mkdir -p MAGNetApp_cloud_data/data/networks

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2Fnetworks&files=igraph_dcm_vs_hcm.RData&downloadStartSecret=z1kpacl8cbd" -O MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_hcm.RData

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2Fnetworks&files=igraph_dcm_vs_nfd.RData&downloadStartSecret=c6cpfceqazg" -O MAGNetApp_cloud_data/data/networks/igraph_dcm_vs_nfd.RData

wget "https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download?path=%2Fdata%2Fnetworks&files=igraph_hcm_vs_nfd.RData" -O MAGNetApp_cloud_data/data/networks/igraph_hcm_vs_nfd.RData

wget "https://data.dieterichlab.org/s/dMzSgF98DRQ7yYc/download/summarized_experiment.RDS" -O MAGNetApp_cloud_data/data/summarized_experiment.RDS

wget "http://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz" -O data/GRCh38.96.gtf.gz

