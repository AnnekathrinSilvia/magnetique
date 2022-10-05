# `magnetique`

This repository contains the source for the `magnetique` application, presented in the work "Magnetique: An interactive web application to explore transcriptome signatures of heart failure".

The data required for running this app is located at: https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb.

A separate repository (https://github.com/dieterich-lab/magnetiqueCode2022) contains the source codes for the analyses of the MAGNet RNA-seq dataset, as they are made available via the `magnetique` application.  

## Detailing files
- `Dockerfile`: contains the instructions to create the docker image
- `Rprofile.site`: Shiny options
- `renv.loc`: R package enviroment with provenance 

## Magnetique DTU view

Here I describe the process to plot the DTU results.

The annotation we use is here: http://ftp.ensembl.org/pub/release-102/gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz

The summarized experiment object for DRIMseq results with the MAGnet dataset goes here: https://data.dieterichlab.org/s/dMzSgF98DRQ7yYc. 
These files  should be in the data directory `DTU/data/`.

Each row is in the object is a transcript, and it has three assays:
- counts: read counts per transcript
- proportions: proportions of reads per gene
- fit_full: fitted proportions
I used the fitter proportion for the plots. 

In addition, the colData and the rowData have the sample information and the transcript annotation, respectively.

Modelling results are accesible via rowData and named as "DRIMSeq_HCM_vs_NFD", "DRIMSeq_DCM_vs_NFD" and "DRIMSeq_DCM_vs_HCM"

## Files:
- `app.R`: quick and dirty prototype to check how the plots look into the browser
- `plots_with_se_obj.R`: defines the functions used app
- `write_summarized_experiment.R`: defines how I create a SE object from the DRIMSeq result
- `summarized_experiment.RDS`: the SE object define above

## DTU TODO:
- This is missing a table with the main functional differences from the transcript isoforms
- Also, it would be interesting to make the proportion plot interactive, showing the information patient information + phenotype (hypertension, diabetes etc.) because some of group show proportions that are bimodal
