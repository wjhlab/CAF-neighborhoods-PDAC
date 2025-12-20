####################################################################
#   01. cell phenotype clustering heatmap and cell type abundance  #
####################################################################

######RBINDING THE DATA INTO ONE DATAFRAME##### ================================
final_colnames<- colnames(data_caf)[colnames(data_caf) %nin% c("cluster","cluster_str")]

#merge t cell annotations into cd45 dataset
tcell_a<- paste(data_tcell$ImageId, data_tcell$CellId, sep="_")
cd45_a<- paste(data_cd45$ImageId, data_cd45$CellId, sep="_")

data_cd45$cluster1m_refined[match(tcell_a,cd45_a)] <- data_tcell$cluster1m_refined

#data rbind
data_epi_caf <- rbind(data_caf[, final_colnames], data_tumor[,final_colnames])
data_full_anno_b <- rbind(data_epi_caf, data_cd45[, final_colnames])

data_full_anno_b[grepl("unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Immune undefined"
data_full_anno_b[grepl("Unassigned", data_full_anno_b$cluster1m_refined), "cluster1m_refined"] <- "Tumor undefined"

unique(data_full_anno_b$cluster1m_refined)

# asinh tansform
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

colorlevels<- c("iCAF","myCAF","CD105+ CAF","ApCAF","CXCL12+ CAF","FAP+ PDPN+ HLADR+ CAF","CAF undefined",
                "Tumor undefined","Mixed","Classical", "Basal",
                "Immune undefined","CD57+","B cell","Myeloid", "Treg",
                "CD8+","CD4+","Cytotoxic CD8","Memory CD4","Activated CD4",
                "Naïve CD4","Naïve CD8")

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


### Figure 1C - all cell type phenotype clustering heatmap ####
pdf("./output/Fig1C.pdf",width=10,height=10)
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

### Figure 2A ####
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
raw_areas <- read_xlsx('./Config/area.xlsx')
areas$TotalArea <- raw_areas$TotalArea[match(areas$sample_id, raw_areas$sample_id)]
densities <- t(t(counts)/areas$TotalArea)

#write.csv(counts,'./output/counts_data_anno.csv')
#write.csv(props,'./output/props_data_anno.csv')
#write.csv(densities, './output/densities_data_anno.csv')

ggdf <- reshape2::melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
             id.vars = "cluster", value.name = "counts", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$case <- factor(md$Case[match(ggdf$sample_id, md$sample_id)], levels = caselevels)
ggdf$tumor <- factor(md$Tumor[match(ggdf$sample_id, md$sample_id)], levels = tumorlevels)
ggdf$patient <- factor(md$Patient[match(ggdf$sample_id, md$sample_id)], levels = patientlevels)
ggdf$cluster<- factor(ggdf$cluster, levels=colorlevels)

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

### Figure 2B - abundance plot for CAF and tumors only ####
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
p <- ggplot(ggdf, aes(x = patient, y = counts, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorassigned) +
  guides(fill = guide_legend(ncol = 1))

# Extract the legend
legend <- get_legend(p)

# Plot the legend only
### Fig.S5A Legend ####
pdf('./output/FigS5A_legend.pdf',width=4,height=6)
figS5A_legend<- plot_grid(legend, ncol = 1)
print(figS5A_legend)
dev.off()

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

### Figure S5A ####
pdf('./output/FigS5A.pdf',width=4,height=6)
figS5A<- plot_grid(p2,p1,ncol=1, align = "v", axis = "lr", rel_heights = c(0.9, 1))
print(figS5A)
dev.off()

### Figure 1D - proportion of patients for each celltype (cluster1m)
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

### Figure 1D ####
pdf('./output/Fig1D.pdf', width=6, height=4)
fig1D<- ggplot(ggdf1m, aes(x = Patient, y = counts_1m, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +  
  scale_fill_manual(values = color_celltype) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.text = element_text(color="black"), 
        legend.position = "right")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "", y = "Proportion", fill = "Cell Type")
print(fig1D)
dev.off()


### broader subset clustering heatmap ####
row_order <- c("Immune unassigned", "Myeloid", "B cell", "CD57+", "Tumor" , "CAF","CD4", "CD8", "Treg")

#data rbind
data_epi_caf <- rbind(data_caf[, final_colnames], data_tumor[,final_colnames])

#broad cluster annotation
tcell_a<- paste(data_tcell$ImageId, data_tcell$CellId, sep="_")
cd45_a<- paste(data_cd45$ImageId, data_cd45$CellId, sep="_")

data_cd45$cluster1m_refined[match(tcell_a,cd45_a)] <- data_tcell$cluster1m_broad

# group into broader cluster 
data_cd45_broader <- data_cd45
data_cd45_broader[grepl("unassigned", data_cd45_broader$cluster1m_refined), "cluster1m_refined"] <- "Immune unassigned"
data_cd45_broader$cluster1m<- data_cd45_broader$cluster1m_refined

data_full_anno_b <- rbind(data_epi_caf, data_cd45_broader[, final_colnames])

unique(data_full_anno_b$cluster1m)

# asinh tansform
data_full_anno_b1 <- data.matrix(data_full_anno_b[,union(subtype_markers,functional_markers)])
data_full_anno_b2 <- asinh(data_full_anno_b1 / 0.8)

##clustering of data
rng <- colQuantiles(data_full_anno_b2, probs = c(0.01, 0.95))
data_full_anno_b01 <- t((t(data_full_anno_b2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_full_anno_b01[data_full_anno_b01 < 0] <- 0; data_full_anno_b01[data_full_anno_b01 > 1] <- 1

##CLUSTERING ANNOTATED FOR FINAL HEATMAP
cluster_mean_merged <- data.frame(data_full_anno_b01, cluster = data_full_anno_b$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-dittoColors()[1:(length(unique(cluster_mean_merged$cluster)))]
names(colorassigned)<-unique(data_full_anno_b$cluster1m)
color_list = list(clusters=colorassigned)

cluster_mean_merged_mat_level<- cluster_mean_merged_mat[row_order, selected_markers] # arrange the heatmap in groups of similar celltypes 

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full_anno_b$cluster1m)[rownames(cluster_mean_merged_mat_level)]),
                    gp = gpar(fill=colorassigned[rownames(cluster_mean_merged_mat_level)]),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

### Figure S4D ####
pdf("./output/FigS4D.pdf",width=15,height=10)
figS4D<-Heatmap(cluster_mean_merged_mat_level,
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
                width = ncol(cluster_mean_merged_mat)*unit(3.5, "mm"), 
                height = nrow(cluster_mean_merged_mat)*unit(3.5, "mm"))
draw(figS4D)
dev.off()

### Average collagen counts per patient in each cell type 
data_full_collagen <- data_full[,c("ImageId","Collagen", "cluster1m")]

collagen_CAF_Tumor <- data_full_collagen[data_full_collagen$cluster1m %in% c("CAF", "Tumor"),]

ggdfNEW <- collagen_CAF_Tumor %>%
  group_by(ImageId, cluster1m) %>%
  summarise(mean_collagen= mean(Collagen))

### Figure S8A ####
pdf("./output/FigS8A.pdf",width=4,height=4)
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
