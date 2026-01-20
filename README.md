# CAF Neighborhood Analysis in Pancreatic Ductal Adenocarcinoma (PDAC)

## Overview

This repository contains R scripts for custom analyses of cancer-associated fibroblast (CAF) neighborhoods and their impact on pancreatic ductal adenocarcinoma (PDAC). 

The code covers cell-type phenotype visualization, CAF neighborhood analyses, and validation of PhenoGraph IMC clusters using Pixie, as presented in the associated manuscript.

## Associated paper

This repository is associated with the following manuscript:

> **Cancer associated fibroblasts promote epithelial to mesenchymal transition and classical to basal change in pancreatic cancer cells in association with loss of IL-8 expression**  

## Input data

All raw MCD data files, and the `backup_output.rds` files containing the fully annotated data frame required to reproduce the manuscript figures, are available at 10.5281/zenodo.14871170 

`Raw TIFF images` used for cell-type validation with Pixie are also available at doi.org/10.5281/zenodo.17887333

## Scripts

Scripts for cell-type phenotype visualization and neighborhood analyses are organized as separate modules. To run either module independently, first execute main.R to load the required data and libraries.

For example, you may skip 01_celltypePhenotype.R and proceed directly to the neighborhood analysis scripts if cell-type visualization is not required.

To run validation of PhenoGraph IMC clusters using Pixie, please ensure that all required TIFF files from Zenodo have been downloaded and run the scripts sequentially.


