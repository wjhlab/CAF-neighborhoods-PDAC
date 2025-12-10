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

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()
metaDataFile = paste0(work,"/Config/metadata.xlsx")
dataDirectory = paste0(work,"/Data")

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



## load output RDS =======
output <- readRDS("./global_data_01_095_fixed.RDS")
data_full <- data.frame(output[1])
data <- data.matrix(output[2])
data01 <- data.frame(output[3])
csv_full <- output[4]




## Read-in metadata and clean =======
ifelse(grepl(metaDataFile,pattern='.xlsx'),
       md <- read_excel(metaDataFile),
       md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv

md$file_name <- factor(md$file_name)
md$File_order <- factor(md$File_order)
md$Tumor <- factor(md$Tumor, levels = tumorlevels) #### check high/low -> then low/high level order 
colnames(md)[7]<-"Patient"



##input image id into metadata
image_id<-c()
for (i in 1:length(md$file_name)){
  tempfile <- read.csv(paste0(dataDirectory,"/",md$file_name[i]))
  df<- as.data.frame(cbind(paste0(md$file_name[i]), unique(tempfile$ImageId)))
  image_id<-rbind(image_id,df)
}
md$ImageId <- image_id$V2[match(md$file_name, image_id$V1)]  # should set md as the reference (come first)

## Make sure all files in metadata present in datadirectory
if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                     pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}

## Read csv into csv_raw =========
csv_raw <- lapply(paste0(dataDirectory,"/",md$file_name),read.csv)
csv_raw_full <- plyr::ldply(csv_raw, rbind)
csv_raw_full$ImageId <- md$sample_id[match(csv_raw_full$ImageId,md$ImageId)]

#export raw channel names to clean/rename panel + rename ImageId as sample_id
rawcolnames <- c()
rawcolnames$name <- colnames(csv_raw_full)
#rawcolnames$sum <- colSums(csv_raw_full)
#write.csv(rawcolnames, 'rawpanel.csv')

#clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx('./Config/cleanpanel.xlsx')
colnames(csv_raw_full) <- cleanpanel$clean_names
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]

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
## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_fulldata.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

#clusterlevels = c(1:5) #before annotations
clusterlevels=c("CD45",
                "unassigned",
                "CAF",
                "Endothelial",
                "Tumor")

mm1 <- match(data_full$cluster, cluster_merging$original_cluster)
data_full$cluster1m <- cluster_merging$new_cluster[mm1]


#subsetting different clusters from broadly annotated data
#CAF
data_caf <- data_full[data_full$cluster1m=="CAF",]

data_tumor <- data_full[data_full$cluster1m=="Tumor",]

data_cd45 <- data_full[data_full$cluster1m=="CD45",]

###CAF =========================================================================
data_caf1 <- data.matrix(data_caf[,Stromamarkers])
data_caf2 <- asinh(data_caf1 / 0.8)


#phenograph clustering of data
rng <- colQuantiles(data_caf2, probs = c(0.01, 0.95))
data_caf01 <- t((t(data_caf2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_caf01[data_caf01 < 0] <- 0; data_caf01[data_caf01 > 1] <- 1
#set.seed(1234)
#phenographout<-Rphenograph(data_caf01)
#data_caf$cluster_str<-factor(membership(phenographout[[2]]))

#save phenographoutput for subset
#save as RDS file
#cafsubset_data <- list(data_caf, data, data01, csv_full)
#saveRDS(cafsubset_data, "CAFsubset_phenographoutput_01_095_fixed.RDS")

cafsubset_data<- readRDS('./CAFsubset_phenographoutput_01_095_fixed.RDS')
data_caf<- cafsubset_data[[1]]

## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_CAF_broad.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("CAF")

mm1 <- match(data_caf$cluster_str, cluster_merging$original_cluster)
data_caf$cluster1m_refined <- cluster_merging$new_cluster[mm1]


cluster_mean_merged <- data.frame(data_caf01, cluster = data_caf$cluster1m_refined, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.data.frame(cluster_mean_merged[,Stromamarkers])
rownames(cluster_mean_merged_mat)<- cluster_mean_merged$cluster


## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusterlevels))
names(colorassigned)<-clusterlevels
color_list = list(clusters=colorassigned)

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_caf$cluster1m_refined)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_cafmergedbroad.pdf",width=10,height=10)
Caf_Heatmap<- Heatmap(as.matrix(cluster_mean_merged_mat),
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_mean_merged_mat)), to = round(max(cluster_mean_merged_mat)))),
        width = ncol(cluster_mean_merged_mat)*unit(4, "mm"), 
        height = nrow(cluster_mean_merged_mat)*unit(4, "mm"))
draw(Caf_Heatmap)
dev.off()




#### from data_tumor ===========================================================
data_tumor1 <- data.matrix(data_tumor[,Epi])
data_tumor2 <- asinh(data_tumor1 / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data_tumor2, probs = c(0.01, 0.95))
data_tumor01 <- t((t(data_tumor2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_tumor01[data_tumor01 < 0] <- 0; data_tumor01[data_tumor01 > 1] <- 1
#set.seed(1234)
#phenographout<-Rphenograph(data_tumor01)
#data_tumor$cluster_epi<-factor(membership(phenographout[[2]]))


#save phenographoutput
#tumorsubset_data <- list(data_tumor, data, data01, csv_full)
#saveRDS(tumorsubset_data, "tumorsubset_phenographoutput_01_095_fixed.RDS")
tumorsubset_data<- readRDS('tumorsubset_phenographoutput_01_095_fixed.RDS')
data_tumor<-tumorsubset_data[[1]]

## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_tumor_basalclassicalbroad.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("Tumor")

mm1 <- match(data_tumor$cluster_epi, cluster_merging$original_cluster)
data_tumor$cluster1m_refined <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data_tumor01, cluster = data_tumor$cluster1m_refined, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_merged_mat<-as.data.frame(cluster_mean_merged[,Epi])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster


## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusterlevels))
names(colorassigned)<-clusterlevels
color_list = list(clusters=colorassigned)

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_tumor$cluster1m_refined)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_tumormerged_broad.pdf",width=10,height=10)
tum_heatmap<-Heatmap(as.matrix(cluster_mean_merged_mat),
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_mean_merged_mat)), to = round(max(cluster_mean_merged_mat)))),
        width = ncol(cluster_mean_merged_mat)*unit(4, "mm"), 
        height = nrow(cluster_mean_merged_mat)*unit(4, "mm"))
draw(tum_heatmap)
dev.off()




### from data_CD45 =============================================================
data_immune1 <- data.matrix(data_cd45[,CD45markers])
data_immune2 <- asinh(data_immune1 / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data_immune2, probs = c(0.01, 0.95))
data_immune01 <- t((t(data_immune2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune01[data_immune01 < 0] <- 0; data_immune01[data_immune01 > 1] <- 1
#set.seed(1234)
#phenographout<-Rphenograph(data_immune01)
#data_cd45$cluster_cd45<-factor(membership(phenographout[[2]]))


#save RDS for CD45 subset
#cd45subset_data <- list(data_cd45, data, data01, csv_full)
#saveRDS(cd45subset_data, "CD45subset_phenographoutput_01_095_fixed.RDS")

cd45subset_data<- readRDS('CD45subset_phenographoutput_01_095_fixed.RDS')
data_cd45<- cd45subset_data[[1]]

####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####
## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_immune.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("T cell",
                "Immune unassigned",
                "CD57",
                "B cell",
                "Myeloid")

mm1 <- match(data_cd45$cluster_cd45, cluster_merging$original_cluster)
data_cd45$cluster1m_refined <- cluster_merging$new_cluster[mm1]


## For subsetting T cells ####
data_tcell <- data_cd45[data_cd45$cluster1m_refined=="T cell",]
### from data_CD3
data_tcell1 <- data.matrix(data_tcell[,Tcellmarkers])
data_tcell2 <- asinh(data_tcell1 / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data_tcell2, probs = c(0.01, 0.95))
data_tcell01 <- t((t(data_tcell2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_tcell01[data_tcell01 < 0] <- 0; data_tcell01[data_tcell01 > 1] <- 1
set.seed(1234)
#phenographout<-Rphenograph(data_tcell01)
#data_tcell$cluster_tcell<-factor(membership(phenographout[[2]]))

tcellsubset_data<- readRDS('tcellsubset_phenographoutput_01_095_fixed.RDS')
data_tcell<- tcellsubset_data[[1]]


####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####
## Load merge file
clusterMergeFile = paste0(work,"/Config/KYLIEmerge_tcell.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("CD8",
                "CD4",
                "Treg")

mm1 <- match(data_tcell$cluster_tcell, cluster_merging$original_cluster)
data_tcell$cluster1m_refined <- cluster_merging$new_cluster[mm1]


######RBINDING THE DATA INTO ONE DATAFRAME##### ================================
#clean up col names
data_caf_col<- colnames(data_caf)[colnames(data_caf) %nin% c("cluster","cluster1m","cluster_str")]
data_tumor_col<- colnames(data_tumor)[colnames(data_tumor) %nin% c("cluster","cluster1m","cluster_epi")]
data_tcell_col<- colnames(data_tcell)[colnames(data_tcell) %nin% c("cluster","cluster1m","cluster_cd45","cluster_tcell")]
data_cd45_col<- colnames(data_cd45)[colnames(data_cd45) %nin% c("cluster","cluster1m","cluster_cd45")]
final_colnames <- data_caf_col

#data clean up
data_caf_clean <- data_caf[,data_caf_col]
data_tumor_clean <- data_tumor[,data_tumor_col]
data_tcell_clean <- data_tcell[,data_tcell_col]
data_cd45_clean <- data_cd45[,data_cd45_col]
colnames(data_tcell_clean)<- final_colnames

#data_tcell_clean$cluster1m_refined2<-data_tcell_clean[] HEREEEEEE
names(data_tcell_clean)[73] <- "cluster1m_refined2"

#merge t cell annotations into cd45 dataset
tcell_a<- paste(data_tcell_clean$ImageId, data_tcell_clean$CellId, sep="_")
cd45_a<- paste(data_cd45_clean$ImageId, data_cd45_clean$CellId, sep="_")
data_cd45_clean$cluster1m_refined[match(tcell_a,cd45_a)] <- data_tcell_clean$cluster1m_refined

#data rbind
data_full_anno_b <- rbind(data_caf_clean, data_tumor_clean)
data_full_anno_b <- rbind(data_full_anno_b, data_cd45_clean)

unique(data_full_anno_b$cluster1m_refined)

#data_full_anno_b[grepl("unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Immune undefined"
#data_full_anno_b[grepl("Unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Tumor undefined"

#####from data_tumor
data_full_anno_b1 <- data.matrix(data_full_anno_b[,union(subtype_markers,functional_markers)])
data_full_anno_b2 <- asinh(data_full_anno_b1 / 0.8)

##clustering of data
rng <- colQuantiles(data_full_anno_b2, probs = c(0.01, 0.95))
data_full_anno_b01 <- t((t(data_full_anno_b2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_full_anno_b01[data_full_anno_b01 < 0] <- 0; data_full_anno_b01[data_full_anno_b01 > 1] <- 1

##CLUSTERING ANNOTATED FOR FINAL HEATMAP
cluster_mean_merged <- data.frame(data_full_anno_b01, cluster = data_full_anno_b$cluster1m_refined, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_merged_mat)<-1:nrow(cluster_mean_merged_mat)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-dittoColors()[1:(length(unique(cluster_mean_merged$cluster)))]

# set color order 
# CAF, Tumor, CD45 
colorlevels <- c(unique(data_caf_clean$cluster1m_refined), 
                 unique(data_tumor_clean$cluster1m_refined),
                 unique(data_cd45$cluster1m_refined),
                 unique(data_tcell$cluster1m_refined))
colorlevels <- setdiff(colorlevels, "T cell")
#colorlevels[colorlevels=="Unassigned"]<- "Tumor undefined"
#colorlevels[colorlevels=="unassigned"]<- "Immune undefined"

names(colorassigned)<-colorlevels
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster
color_list = list(clusters=colorassigned)

# arrange the row of cluster_mean_merged_mat into colorlevels 
selected_markers <- setdiff(union(subtype_markers, functional_markers), "CD45") # remove CD45

cluster_mean_merged_mat<-as.data.frame(cluster_mean_merged[ ,selected_markers])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster

# order heatmap row 
row_order <- c("Immune unassigned", "Myeloid", "B cell", "CD57", "Tumor" , "CAF","CD4", "CD8", "Treg")


cluster_mean_merged_mat_level<- cluster_mean_merged_mat[row_order, ] # arrange the heatmap in groups of similar celltypes 

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full_anno_b$cluster1m_refined)[rownames(cluster_mean_merged_mat_level)]),
                    gp = gpar(fill=colorassigned[rownames(cluster_mean_merged_mat_level)]),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))


### Supplemental Figure 4D ####
pdf("SuppFig4D_clusterheatmap_broad.pdf",width=10,height=10)
fig1C<-Heatmap(as.matrix(cluster_mean_merged_mat_level),
        column_title="Heatmap Annotated Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F, 
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=11),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_mean_merged_mat)), to = round(max(cluster_mean_merged_mat)))),
        width = ncol(cluster_mean_merged_mat)*unit(4, "mm"), 
        height = nrow(cluster_mean_merged_mat)*unit(4, "mm"))
draw(fig1C)
dev.off()

