#load libraries 
library(readxl)
library(dplyr)
library(pals)
library(ComplexHeatmap)
library(matrixStats)
library(Rphenograph)
library(Hmisc)
library(reshape2)
library(dittoSeq)
library(cowplot)
library(ggpubr)
library(tibble)
library(data.table)
library(spatstat)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(readr)
library(tidyverse)
library(limma)
library(ggrepel)
library(stringr)  
library(viridis)

#Set up levels ======
samplevels=c("317_S_hi_i_hi","317_S_hi_i_lo","317_T_hi_i_lo","317_T_hi_i_hi","348_T_hi_i_hi","348_S_hi_i_hi","348_T_hi_i_lo",
             "348_S_hi_i_lo","352_S_hi_i_hi","352_T_hi_i_hi","352_S_hi_i_lo","352_T_hi_i_lo","357_T_hi_i_lo","357_S_hi_i_hi",
             "357_S_hi_i_lo","361_T_hi_i_hi","361_S_hi_i_hi","361_T_hi_i_lo","361_S_hi_i_lo","362_T_hi_i_lo","362_T_hi_i_hi",
             "362_S_hi_i_lo","362_S_hi_i_hi","368_S_hi_i_hi","368_T_hi_i_hi","368_T_hi_i_lo","368_S_hi_i_lo","369_S_hi_i_lo",
             "369_S_hi_i_hi","369_T_hi_i_lo","369_T_hi_i_hi","372_S_hi_i_hi","372_T_hi_i_hi","372_S_hi_i_lo","380_S_hi_i_hi",
             "380_T_hi_i_hi","380_S_hi_i_lo","380_T_hi_i_lo","383_S_hi_i_lo","383_S_hi_i_hi","383_T_hi_i_hi","387_S_hi_i_lo",
             "387_T_hi_i_hi","387_S_hi_i_hi","387_T_hi_i_lo","388_T_hi_i_lo","388_S_hi_i_lo","388_T_hi_i_hi","388_S_hi_i_hi",
             "390_T_hi_i_hi","390_T_hi_i_lo","390_S_hi_i_lo","417_T_hi_i_lo","417_T_hi_i_hi","417_S_hi_i_lo","417_S_hi_i_hi")


tumorlevels=c("Stroma high immune high",
              "Stroma high immune low",
              "Tumor high immune low",
              "Tumor high immune high")

caselevels=as.character(c(1:15))
patientlevels=c("JHH317","JHH348","JHH352","JHH357","JHH361","JHH362","JHH368","JHH369","JHH372","JHH380","JHH383","JHH387","JHH388","JHH390","JHH417")

## load backup output RDS =======
output <- readRDS("./backup/global_data_01_095_fixed.RDS")
data_full <- data.frame(output[1])
data <- data.matrix(output[2])
data01 <- data.frame(output[3])
csv_full <- output[4]

clusterMergeFile = "./Config/merge_fulldata.xlsx" 
cluster_merging <- read_excel(clusterMergeFile)
mm1 <- match(data_full$cluster, cluster_merging$original_cluster)
data_full$cluster1m <- cluster_merging$new_cluster[mm1]

#clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx('./Config/cleanpanel.xlsx')
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]

## Read-in metadata and clean =======
metaDataFile = "./Config/metadata.xlsx"
ifelse(grepl(metaDataFile,pattern='.xlsx'),
       md <- read_excel(metaDataFile),
       md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv

md$file_name <- factor(md$file_name)
md$File_order <- factor(md$File_order)
md$Tumor <- factor(md$Tumor, levels = tumorlevels) #### check high/low -> then low/high level order 
colnames(md)[7]<-"Patient"


## Read markers =========
#sort panels into different categories
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
Tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
Epi <- cleanpanel$clean_names[cleanpanel$epi == 1]
Stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]
CD45markers <- cleanpanel$clean_names[cleanpanel$cd45 == 1]



####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

###CAF =========================================================================
cafsubset_data<- readRDS('./backup/CAFsubset_phenographoutput_01_095_fixed.RDS')
data_caf<- cafsubset_data[[1]]

## Load merge file
clusterMergeFile = "./Config/merge_CAF_broad.xlsx"
cluster_merging <- read_excel(clusterMergeFile)
mm1 <- match(data_caf$cluster_str, cluster_merging$original_cluster)
data_caf$cluster1m_refined <- cluster_merging$new_cluster[mm1]

clusterMergeFile = "./Config/merge_CAF.xlsx"
cluster_merging <- read_excel(clusterMergeFile)
mm1 <- match(data_caf$cluster_str, cluster_merging$original_cluster)
data_caf$cluster2m_refined <- cluster_merging$new_cluster[mm1]


#### from data_tumor ===========================================================
tumorsubset_data<- readRDS('./backup/tumorsubset_phenographoutput_01_095_fixed.RDS')
data_tumor<-tumorsubset_data[[1]]

### from data_CD45 =============================================================
cd45subset_data<- readRDS('./backup/CD45subset_phenographoutput_01_095_fixed.RDS')
data_cd45<- cd45subset_data[[1]]

### from data_tcell =============================================================
tcellsubset_data<- readRDS('./backup/tcellsubset_phenographoutput_01_095_fixed.RDS')
data_tcell<- tcellsubset_data[[1]]

####new ANNOTATIONS AND VISUALIZATION OF CLUSTERS####
## Load merge file
clusterMergeFile = "./Config/tcell_annotation.xlsx"
cluster_merging <- read_excel(clusterMergeFile)

mm1 <- match(data_tcell$cluster_tcell, cluster_merging$original_cluster)
data_tcell$cluster1m_refined <- cluster_merging$new_cluster_updated[mm1]


# add broad cell type cluster
clusterMergeFile = "./Config/tcell_broad.xlsx"
cluster_merging <- read_excel(clusterMergeFile)

mm1 <- match(data_tcell$cluster_tcell, cluster_merging$original_cluster)
data_tcell$cluster1m_broad <- cluster_merging$new_cluster[mm1]



