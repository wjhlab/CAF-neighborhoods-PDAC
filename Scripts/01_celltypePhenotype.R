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
output <- readRDS("./backup/global_data_01_095_fixed.RDS")
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

cafsubset_data<- readRDS('./backup/CAFsubset_phenographoutput_01_095_fixed.RDS')
data_caf<- cafsubset_data[[1]]

## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_CAF_broad.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("iCAF",
                "myCAF",
                "CD105+ CAF",
                "ApCAF",
                "CXCL12+ CAF",
                "FAP+ PDPN+ HLADR+ CAF",
                "CAF undefined")

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

pdf("./output/clusterheatmap_cafmergedbroad.pdf",width=10,height=10)
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
tumorsubset_data<- readRDS('./backup/tumorsubset_phenographoutput_01_095_fixed.RDS')
data_tumor<-tumorsubset_data[[1]]

## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_tumor_basalclassicalbroad.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("Unassigned",
                "Basal",
                "Classical",
                "Mixed")

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

pdf("./output/clusterheatmap_tumormerged_broad.pdf",width=10,height=10)
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

cd45subset_data<- readRDS('./backup/CD45subset_phenographoutput_01_095_fixed.RDS')
data_cd45<- cd45subset_data[[1]]

####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####
## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_immune.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("T cell",
                "unassigned",
                "CD57+",
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

tcellsubset_data<- readRDS('./backup/tcellsubset_phenographoutput_01_095_fixed.RDS')
data_tcell<- tcellsubset_data[[1]]

####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####
## Load merge file
clusterMergeFile = paste0(work,"/Config/KYLIEmerge_tcell.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("CD8+",
                "CD4+",
                "Cytotoxic CD8",
                "Activated CD4",
                "Memory CD4",
                "Naïve CD4",
                "Naïve CD8",
                "Treg")

mm1 <- match(data_tcell$cluster_tcell, cluster_merging$original_cluster)
data_tcell$cluster1m_refined <- cluster_merging$KYLIE_new_cluster[mm1]



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

#merge t cell annotations into cd45 dataset
tcell_a<- paste(data_tcell_clean$ImageId, data_tcell_clean$CellId, sep="_")
cd45_a<- paste(data_cd45_clean$ImageId, data_cd45_clean$CellId, sep="_")
data_cd45_clean$cluster1m_refined[match(tcell_a,cd45_a)] <- data_tcell_clean$cluster1m_refined

#data rbind
data_full_anno_b <- rbind(data_caf_clean, data_tumor_clean)
data_full_anno_b <- rbind(data_full_anno_b, data_cd45_clean)

unique(data_full_anno_b$cluster1m_refined)

data_full_anno_b[grepl("unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Immune undefined"
data_full_anno_b[grepl("Unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Tumor undefined"

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
colorlevels[colorlevels=="Unassigned"]<- "Tumor undefined"
colorlevels[colorlevels=="unassigned"]<- "Immune undefined"

names(colorassigned)<-colorlevels
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster
color_list = list(clusters=colorassigned)

# arrange the row of cluster_mean_merged_mat into colorlevels 
selected_markers <- setdiff(union(subtype_markers, functional_markers), "CD45") # remove CD45

cluster_mean_merged_mat<-as.data.frame(cluster_mean_merged[ ,selected_markers])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster

# order heatmap row 
row_order <- c("Basal", "Mixed", "Classical" , "Tumor undefined" , 
                "CD57+",  "B cell",  "Myeloid", "Immune undefined",
               "CD8+", "CD4+", "Cytotoxic CD8", "Activated CD4", 
               "Memory CD4", "Naïve CD4", "Naïve CD8", "Treg", 
               "iCAF" ,  "myCAF", "CD105+ CAF",  "ApCAF", "CXCL12+ CAF", 
               "FAP+ PDPN+ HLADR+ CAF" ,  "CAF undefined")


cluster_mean_merged_mat_level<- cluster_mean_merged_mat[row_order, ] # arrange the heatmap in groups of similar celltypes 

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full_anno_b$cluster1m_refined)[rownames(cluster_mean_merged_mat_level)]),
                    gp = gpar(fill=colorassigned[rownames(cluster_mean_merged_mat_level)]),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))


### Fig 1C ####
pdf("./output/Fig1C_clusterheatmap_full_anno.pdf",width=10,height=10)
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


### Figure 2 - CAF/tumor subset ####

caf_tumors <- c("Basal", "Mixed", "Classical" , "Tumor undefined" , 
               "iCAF" ,  "myCAF", "CD105+ CAF",  "ApCAF", "CXCL12+ CAF", 
               "FAP+ PDPN+ HLADR+ CAF" ,  "CAF undefined")

caf_tumor_markers <- union(Stromamarkers, Epi)

data_full_anno_b1 <- data.matrix(data_full_anno_b[data_full_anno_b$cluster1m_refined%in%caf_tumors,
                                                  union(subtype_markers,functional_markers)])
data_full_anno_b2 <- asinh(data_full_anno_b1 / 0.8)

##clustering of data
rng <- colQuantiles(data_full_anno_b2, probs = c(0.01, 0.95))
data_full_anno_b01 <- t((t(data_full_anno_b2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_full_anno_b01[data_full_anno_b01 < 0] <- 0; data_full_anno_b01[data_full_anno_b01 > 1] <- 1

##CLUSTERING ANNOTATED FOR FINAL HEATMAP
clusters<- data_full_anno_b[data_full_anno_b$cluster1m_refined%in%caf_tumors,]$cluster1m_refined
cluster_mean_merged <- data.frame(data_full_anno_b01, cluster = clusters, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster


cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full_anno_b$cluster1m_refined)[caf_tumors]),
                    gp = gpar(fill=colorassigned[caf_tumors]),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

### Fig 2A ####

pdf("./output/Fig2A.pdf",width=10,height=10)
fig2A<-Heatmap(cluster_mean_merged_mat[caf_tumors,caf_tumor_markers],
        column_title="Heatmap Annotated Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_mean_merged_mat)), to = round(max(cluster_mean_merged_mat)))),
        width = ncol(cluster_mean_merged_mat)*unit(2, "mm"), 
        height = nrow(cluster_mean_merged_mat)*unit(3.5, "mm"))
draw(fig2A)
dev.off()


###Abundance Plots
counts_table_data_anno_b <- table(data_full_anno_b$cluster1m_refined, data_full_anno_b$ImageId)
props_table_data_anno_b<- t(t(counts_table_data_anno_b) / colSums(counts_table_data_anno_b)) * 100
counts <- as.data.frame.matrix(counts_table_data_anno_b)
props <-as.data.frame.matrix(props_table_data_anno_b)

#created dummy column
areas <- tibble(sample_id=colnames(counts))

#Densities
raw_areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'./output/counts_data_anno.csv')
write.csv(props,'./output/props_data_anno.csv')
write.csv(densities, './output/densities_data_anno.csv')

ggdf <- reshape2::melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
             id.vars = "cluster", value.name = "counts", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$case <- factor(md$Case[match(ggdf$sample_id, md$sample_id)], levels = caselevels)
ggdf$tumor <- factor(md$Tumor[match(ggdf$sample_id, md$sample_id)], levels = tumorlevels)
ggdf$patient <- factor(md$Patient[match(ggdf$sample_id, md$sample_id)], levels = patientlevels)
ggdf$cluster<- factor(ggdf$cluster, levels=colorlevels)

p <- ggplot(ggdf, aes(x = patient, y = counts, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorassigned) +
  guides(fill = guide_legend(ncol = 1))

# Extract the legend
legend <- get_legend(p)

# Plot the legend only
### Fig 1D Legend ####
pdf('./output/Fig1Dlegend.pdf',width=4,height=6)
fig1D_legend<- plot_grid(legend, ncol = 1)
print(fig1D_legend)
dev.off()

#proportion of cell breakdown by patient filled stack
p1<- ggplot(ggdf, aes(x = patient, y = counts, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +  # Setting position to "fill"
  scale_fill_manual(values = color_list$clusters) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.text.y = element_text(color="black"), 
        legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "Patient", y = "Proportion", fill = "Cluster")  # Correcting labels

### Fig 2B - abundance plot for CAF and tumors only ####
pdf('./output/Fig2B.pdf',width=6,height=3)
fig2B<- ggplot(ggdf[ggdf$cluster%in%caf_tumors, ], aes(x = patient, y = counts, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_fill_manual(values = color_list$clusters) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.text.y = element_text(color="black"), 
        legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  guides(fill = guide_legend(title = NULL))+
  labs(x = "Patient", y = "Proportion", fill = "Cluster") 
print(fig2B)
dev.off()


#counts of cell breakdown by patient filled stack
ggdf <- reshape2::melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
             id.vars = "cluster", value.name = "counts", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$case <- factor(md$Case[match(ggdf$sample_id, md$sample_id)], levels = caselevels)
ggdf$tumor <- factor(md$Tumor[match(ggdf$sample_id, md$sample_id)], levels = tumorlevels)
ggdf$patient <- factor(md$Patient[match(ggdf$sample_id, md$sample_id)], levels = patientlevels)
ggdf$cluster<- factor(ggdf$cluster, levels=colorlevels)


p2<- ggplot(ggdf, aes(x = patient, y = counts, fill = cluster)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = color_list$clusters) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(color="black"), 
        legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "", y = "Count", fill = "Cluster")  

### Fig 1D ####

pdf('./output/Fig1D.pdf',width=4,height=6)
fig1D<- plot_grid(p2,p1,ncol=1, align = "v", axis = "lr", rel_heights = c(0.9, 1))
print(fig1D)
dev.off()

### Fig 1D extended - proportion of patients for each celltype (cluster1m)
counts_table_data_1m <- table(data_full$cluster1m, data_full$ImageId)
counts_1m <- as.data.frame.matrix(counts_table_data_1m)

ggdf1m <- reshape2::melt(data.frame(cluster = rownames(counts_1m), counts_1m, check.names = FALSE),
                       id.vars = "cluster", value.name = "counts_1m", 
                       variable.name = "sample_id")
ggdf1m$sample_id <- factor(ggdf1m$sample_id, levels=samplevels)
ggdf1m$case <- factor(md$Case[match(ggdf1m$sample_id, md$sample_id)], levels = caselevels)
ggdf1m$tumor <- factor(md$Tumor[match(ggdf1m$sample_id, md$sample_id)], levels = tumorlevels)
ggdf1m$Patient <- factor(md$Patient[match(ggdf1m$sample_id, md$sample_id)], levels = patientlevels)

ggdf1m[grepl("Unassigned", ggdf1m$cluster), "cluster"] <- "Tumor_Unassigned"
ggdf1m$cluster<- factor(ggdf1m$cluster, levels=c("CAF", "CD45", "Endothelial", "Tumor", "Tumor_Unassigned"))

color_celltype <- dittoColors()[1:length(unique(ggdf1m$cluster))]
names(color_celltype)<- levels(ggdf1m$cluster)

pdf('./output/Patient_proportion_perCelltype.pdf', width=6, height=4)
p3<- ggplot(ggdf1m, aes(x = Patient, y = counts_1m, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_fill_manual(values = color_celltype) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.text = element_text(color="black"), 
        legend.position = "right")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "", y = "Proportion", fill = "Cell Type")
print(p3)
dev.off()

### Average collagen counts per patient in each cell type 
data_full_collagen <- data_full[,c("ImageId","Collagen", "cluster1m")]

collagen_CAF_Tumor <- data_full_collagen[data_full_collagen$cluster1m %in% c("CAF", "Tumor"),]

ggdfNEW <- collagen_CAF_Tumor %>%
  group_by(ImageId, cluster1m) %>%
  summarise(mean_collagen= mean(Collagen))

pdf("./output/average_collagen_per_sample.pdf",width=4,height=4)
collagen_CAF<-ggplot(ggdfNEW, aes(x=cluster1m, y=mean_collagen, fill= cluster1m))+
  geom_boxplot(outlier.color = NA, width=0.75, linewidth=0.25)+
  geom_point(position = position_jitter(width = 0.2)) +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = color_celltype)+
  labs(x = "Cell Type", y = "Average Collagen Expression per sample", fill = "Cluster") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("CAF", "Tumor")),
    #label.x.npc =  c('center'), #for adjusting p value location 
    size = 3,
    color = "black",
    hide.ns = FALSE)
print(collagen_CAF)
dev.off()

densities_1m <- reshape2::melt(t(t(counts_1m)/areas$TotalArea))
names(densities_1m)<-c("cluster", "sample_id", "densities")

densities_1m$Treatment <- md$Treatment[match(densities_1m$sample_id, md$sample_id)]

color_treatment <- dittoColors()[1:length(unique(densities_1m$Treatment))]
names(color_treatment)<- unique(densities_1m$Treatment)

ggplot(densities_1m, aes(x=Treatment, y=densities, fill= Treatment))+
  geom_boxplot(outlier.color = NA, width=0.75, linewidth=0.25)+
  geom_point(position = position_jitter(width = 0.2)) +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank(),
        legend.position = "right")+
  scale_fill_manual(values = color_treatment)+
  labs(x = "Cell Type", y = "Cell type densities per sample", fill = "Cluster") +
  stat_compare_means(
    method = "wilcox.test",
    size = 3,
    color = "black",
    comparisons = list(c("Untreated", "FFX")),
    hide.ns = FALSE)+
  facet_wrap(~cluster)

