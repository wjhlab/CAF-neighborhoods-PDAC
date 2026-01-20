#######################################################
#             02. CAF neighborhood analysis           #
#######################################################

# FFX treated & untreated patients 
treated_sampleId<- unique(md[md$Treatment=="FFX", ]$sample_id)
untreated_sampleId<- unique(md[md$Treatment=="Untreated", ]$sample_id)

## Prepare clean data (N1, N2, N3)
data_caf_clean_ext <- data_caf %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"))

data_tumor_clean_ext <- data_tumor %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"))

### Figure 2C - broad CAF subtypes & Tumor (basal, classical, Mixed) relationship 

data_t_basal <- data_tumor[data_tumor$cluster1m_refined == "Basal",]
data_t_classical <- data_tumor[data_tumor$cluster1m_refined == "Classical",]
data_t_mixed <- data_tumor[data_tumor$cluster1m_refined == "Mixed",]

basalcell_ind <- paste(data_t_basal$ImageId,data_t_basal$CellId,sep="_")
classicalcell_ind <- paste(data_t_classical$ImageId,data_t_classical$CellId,sep="_")
mixedcell_ind <- paste(data_t_mixed$ImageId,data_t_mixed$CellId,sep="_")

### treated and untreated patient 
data_caf_clean_treated <- data_caf_clean_ext[data_caf_clean_ext$ImageId%in%treated_sampleId, ]
data_caf_clean_untreated <- data_caf_clean_ext[data_caf_clean_ext$ImageId%in%untreated_sampleId, ]

get_near <- function(df, subtype_indices) {
  mask1 <- df$NN1_ext %in% subtype_indices
  mask2 <- df$NN2_ext %in% subtype_indices
  mask3 <- df$NN3_ext %in% subtype_indices
  c(df$cluster1m_refined[mask1],
    df$cluster1m_refined[mask2],
    df$cluster1m_refined[mask3])
}

# Apply the function for each T-cell type and store in a list
get_near_result <- function(data, index_list) {
  lapply(index_list, function(ind) get_near(data, ind))
}

# convert result list to summarized frequency data frame
summarize_near <- function(result_list) {
  # Count frequencies
  tables <- lapply(result_list, function(near) {
    near <- as.character(near)  
    as.data.frame(table(near), stringsAsFactors = FALSE)
  })
  all_types <- unique(unlist(lapply(tables, function(tbl) tbl$near)))
  # Create a data frame and fill with matched frequencies
  summary_df <- data.frame(row.names = all_types)
  for (name in names(tables)) {
    tbl <- tables[[name]]
    summary_df[[name]] <- tbl$Freq[match(rownames(summary_df), tbl$near)]
  }
  summary_df[is.na(summary_df)] <- 0
  summary_df <- summary_df[rownames(summary_df) != "Unassigned", ] # remove unassigned
  return(summary_df)
}


# Apply to each subset
cafs_near_basal_treated     <- get_near(data_caf_clean_treated, basalcell_ind)
cafs_near_classical_treated <- get_near(data_caf_clean_treated, classicalcell_ind)
cafs_near_mixed_treated     <- get_near(data_caf_clean_treated, mixedcell_ind)

cafs_near_basal_untreated     <- get_near(data_caf_clean_untreated, basalcell_ind)
cafs_near_classical_untreated <- get_near(data_caf_clean_untreated, classicalcell_ind)
cafs_near_mixed_untreated     <- get_near(data_caf_clean_untreated, mixedcell_ind)

## combined - all samples
cafs_near_basal_combined     <- get_near(data_caf_clean_ext, basalcell_ind)
cafs_near_classical_combined <- get_near(data_caf_clean_ext, classicalcell_ind)
cafs_near_mixed_combined     <- get_near(data_caf_clean_ext, mixedcell_ind)


make_caf_summary <- function(basal, classical, mixed) {
  all_names <- unique(c(basal, classical, mixed))
  df <- data.frame(row.names = all_names)
  df$basal     <- table(factor(basal, levels = all_names))
  df$classical <- table(factor(classical, levels = all_names))
  df$mixed     <- table(factor(mixed, levels = all_names))
  return(df)
}

# Generate summaries
caf_near_epi_sum_combined  <- make_caf_summary(cafs_near_basal_combined, cafs_near_classical_combined, cafs_near_mixed_combined)
caf_near_epi_sum_treated   <- make_caf_summary(cafs_near_basal_treated, cafs_near_classical_treated, cafs_near_mixed_treated)
caf_near_epi_sum_untreated <- make_caf_summary(cafs_near_basal_untreated, cafs_near_classical_untreated, cafs_near_mixed_untreated)

##Heatmap of CAF subtypes by tumor cell subtypes
### Figure 2C ####
fig2C_1<-Heatmap(t(scale(t(as.matrix(caf_near_epi_sum_combined)))),
        cluster_rows = T, cluster_columns = T, 
        name = "Scaled\nFreq",
        column_title = "Combined",
        width = ncol(caf_near_epi_sum_combined)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum_combined)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

fig2C_2<-Heatmap(t(scale(t(as.matrix(caf_near_epi_sum_untreated)))),
        cluster_rows = , cluster_columns = T, 
        name = "Scaled\nFreq",
        column_title = "Untreated",
        width = ncol(caf_near_epi_sum_untreated)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum_untreated)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

fig2C_3<-Heatmap(t(scale(t(as.matrix(caf_near_epi_sum_treated)))),
        cluster_rows = T, cluster_columns = T, 
        name = "Scaled\nFreq",
        column_title = "Treated",
        width = ncol(caf_near_epi_sum_treated)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum_treated)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

### Figure 2C and Figure S5B & C ####
pdf('./output/Fig2C_combined.pdf')
draw(fig2C_1)
draw(fig2C_2) # Fig S5B
draw(fig2C_3) # Fig S5C
dev.off()

# Convert data to format for ggplot
plot_data_long <- reshape2::melt(
  as.data.frame.matrix(caf_near_epi_sum_combined) %>% 
    tibble::rownames_to_column("Subtype"),
  id.vars = "Subtype")
                                    
names(plot_data_long)<-c("Subtype","TumorCellType", "Frequency")
plot_data_long$Frequency <- as.numeric(plot_data_long$Frequency)

colorassigned<- dittoColors()[1:length(unique(plot_data_long$Subtype))]
clusterlevels=c("iCAF",
                "myCAF",
                "CD105+ CAF",
                "ApCAF",
                "CXCL12+ CAF",
                "FAP+ PDPN+ HLADR+ CAF",
                "CAF undefined")
names(colorassigned)<- clusterlevels
plot_data_long$Subtype<- factor(plot_data_long$Subtype, levels=clusterlevels)

# Create a stacked bar chart with custom colors
### Figure 2D ####
pdf('./output/Fig2D.pdf', height=3.5, width=4)
fig2D<-ggplot(plot_data_long, aes(x = TumorCellType, y = Frequency, fill = Subtype)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "CAF Subtype")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned)
print(fig2D)
dev.off()

### Figure 3K, 3L ####
data_tumor_clean_ext2 <- data_tumor %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"),
         NN4_ext = paste(ImageId, NN4, sep = "_"),
         NN5_ext = paste(ImageId, NN5, sep = "_"),
         NN6_ext = paste(ImageId, NN6, sep = "_"),
         NN7_ext = paste(ImageId, NN7, sep = "_"),
         NN8_ext = paste(ImageId, NN8, sep = "_"),
         NN9_ext = paste(ImageId, NN9, sep = "_"),
         NN10_ext = paste(ImageId, NN10, sep = "_"),
         NN11_ext = paste(ImageId, NN11, sep = "_"),
         NN12_ext = paste(ImageId, NN12, sep = "_"))
keymarkerlist <- c("KI67","HLADR","VIM","NCAD","Ecad")

CAF_ind <- paste0(data_caf$ImageId, "_", data_caf$CellId)

##for CAF near/not near Basal
basalnearCAF <- c()

data_tumor_clean_ext2_subtype<- data_tumor_clean_ext2[data_tumor_clean_ext2$cluster1m_refined=="Basal", ]
#save the cell index for any basal cell that has a given all CAF as its nearest neighbor (any one of 12)
#if any of the 12 neighbors belong within the CAF indices and the number of TRUEs are >0, then return the cell index
for(i in 1:nrow(data_tumor_clean_ext2_subtype)){
  if(
    sum(data_tumor_clean_ext2_subtype[i,c('NN1_ext','NN2_ext','NN3_ext','NN4_ext','NN5_ext','NN6_ext',
                                          'NN7_ext','NN8_ext','NN9_ext','NN10_ext','NN11_ext','NN12_ext')] %in% CAF_ind) >0){
    basalnearCAF <- c(basalnearCAF,data_tumor_clean_ext2_subtype[i,]$Id_ext)
  }}

#create a ggdf based on the cell indices returned from the loop  
ggdfvio_caf<- rbind(data_tumor_clean_ext2[data_tumor_clean_ext2$Id_ext %in% basalnearCAF,c("Id_ext",keymarkerlist)],
                    data_tumor_clean_ext2[!data_tumor_clean_ext2$Id_ext %in% basalnearCAF,c("Id_ext",keymarkerlist)])
ggdfvio_caf$near <- 1
ggdfvio_caf[ggdfvio_caf$Id_ext %in% basalnearCAF,]$near<-"Near"
ggdfvio_caf[!ggdfvio_caf$Id_ext %in% basalnearCAF,]$near<-"Not_Near"
ggdfvio_caf[,keymarkerlist] <- asinh(ggdfvio_caf[,keymarkerlist] / 0.8)

### Figure 3K ####
pdf("./output/3M_basal_CAF_NN_Ecad.pdf",width=3,height=3)
fig3M<-ggplot(ggdfvio_caf, aes(x=near, y=Ecad, fill=near))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(outlier.color = NA, width=0.25, alpha=0.2, linewidth=0.25)+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank())+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(aes(label = paste0("p=", after_stat(p.format))))
print(fig3M)
dev.off()

### Figure 3L ####
pdf("./output/3N_basal_CAF_NN_Vim.pdf",width=3,height=3)
fig3N<-ggplot(ggdfvio_caf, aes(x=near, y=VIM, fill=near))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(outlier.color = NA, width=0.25, alpha=0.2, linewidth=0.25)+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank())+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(aes(label = paste0("p=", after_stat(p.format))))
print(fig3N)
dev.off()


## finer CAF subtypes & Tumor (basal, classical, Mixed) relationship
get_cafs_near_subtype_2m <- function(df, subtype_indices) {
  mask1 <- df$NN1_ext %in% subtype_indices
  mask2 <- df$NN2_ext %in% subtype_indices
  mask3 <- df$NN3_ext %in% subtype_indices
  c(df$cluster2m_refined[mask1],
    df$cluster2m_refined[mask2],
    df$cluster2m_refined[mask3])
}

cafs_near_basal_combined_2m     <- get_cafs_near_subtype_2m(data_caf_clean_ext, basalcell_ind)
cafs_near_classical_combined_2m <- get_cafs_near_subtype_2m(data_caf_clean_ext, classicalcell_ind)
cafs_near_mixed_combined_2m     <- get_cafs_near_subtype_2m(data_caf_clean_ext, mixedcell_ind)

# Generate summaries
caf_near_epi_sum_combined_2m  <- make_caf_summary(cafs_near_basal_combined_2m, cafs_near_classical_combined_2m, cafs_near_mixed_combined_2m)

## Figure 5H - Heatmap of CAF subtypes by tumor cell subtypes ####
pdf('./output/Fig5G.pdf')
fig5G<-Heatmap(t(scale(t(as.matrix(caf_near_epi_sum_combined_2m)))),
        name = "Scaled\nFreq",
        width = ncol(caf_near_epi_sum_combined_2m)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum_combined_2m)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
draw(fig5G)
dev.off()


## immune tumor nearest neighbor
data_cd45<- cd45subset_data[[1]]

broadImmune_tbl<- split(data_cd45, data_cd45$cluster1m_refined)
broadImmune_indices <- lapply(broadImmune_tbl, function(df) {
  paste(df$ImageId, df$CellId, sep = "_")
})

# Extend data_tumor with necessary columns for neighbor identification
tumor_near_results  <- get_near_result(data_tumor_clean_ext, broadImmune_indices)

# Convert each result to a data frame with frequency counts
tumor_near_tables <- summarize_near(tumor_near_results)

### Figure 6B  ####
pdf('./output/Fig6B.pdf')
fig6B<- Heatmap(t(scale(t(as.matrix(tumor_near_tables)))),
        name = "Scaled\nFreq",
        width = ncol(tumor_near_tables)*unit(7, "mm"), 
        height = nrow(tumor_near_tables)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
draw(fig6B)
dev.off()



## immune CAF nearest neighbor

# Generate CAFs near subtypes for each condition
cafs_near_by_subtype_combined  <- get_near_result(data_caf_clean_ext,        broadImmune_indices)
cafs_near_by_subtype_treated   <- get_near_result(data_caf_clean_treated,    broadImmune_indices)
cafs_near_by_subtype_untreated <- get_near_result(data_caf_clean_untreated,  broadImmune_indices)

# Construct immune summary matrices
caf_near_immune_sum_combined <- summarize_near(cafs_near_by_subtype_combined)
caf_near_immune_sum_treated  <- summarize_near(cafs_near_by_subtype_treated)
caf_near_immune_sum_untreated <- summarize_near(cafs_near_by_subtype_untreated)

## Heatmap of CAF subtypes by tumor cell subtypes
CAF_Imm_1<-Heatmap(t(scale(t(as.matrix(caf_near_immune_sum_combined)))),
        cluster_rows = T, cluster_columns = T,
        name = "Scaled\nFreq",
        column_title = "Combined",
        width = ncol(caf_near_immune_sum_combined)*unit(7, "mm"),
        height = nrow(caf_near_immune_sum_combined)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

CAF_Imm_2<-Heatmap(t(scale(t(as.matrix(caf_near_immune_sum_untreated)))),
        cluster_rows = , cluster_columns = T,
        name = "Scaled\nFreq",
        column_title = "Untreated",
        width = ncol(caf_near_immune_sum_untreated)*unit(7, "mm"),
        height = nrow(caf_near_immune_sum_untreated)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

CAF_Imm_3<-Heatmap(t(scale(t(as.matrix(caf_near_immune_sum_treated)))),
        cluster_rows = T, cluster_columns = T,
        name = "Scaled\nFreq",
        column_title = "Treated",
        width = ncol(caf_near_immune_sum_treated)*unit(7, "mm"),
        height = nrow(caf_near_immune_sum_treated)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))

### Figure 6C  ####
pdf('./output/Fig6C.pdf')
draw(CAF_Imm_1)
draw(CAF_Imm_2)
draw(CAF_Imm_3)
dev.off()

## Cluster heatmap annotation for T cell subtypes =======
data_tcell1 <- data.matrix(data_tcell[,Tcellmarkers])
data_tcell2 <- asinh(data_tcell1 / 0.8)

rng <- colQuantiles(data_tcell2, probs = c(0.01, 0.95))
data_tcell01 <- t((t(data_tcell2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_tcell01[data_tcell01 < 0] <- 0; data_tcell01[data_tcell01 > 1] <- 1

clusterlevels=c("CD8+",
                "CD4+",
                "Cytotoxic CD8",
                "Activated CD4",
                "Memory CD4",
                "Naïve CD4",
                "Naïve CD8",
                "Treg")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(data_tcell$cluster1m_refined)))
names(colorassigned)<-clusterlevels

cluster_mean_merged <- data.frame(data_tcell01, cluster = data_tcell$cluster1m_refined, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_merged_mat<-as.data.frame(cluster_mean_merged[,Tcellmarkers])
rownames(cluster_mean_merged_mat)<-cluster_mean_merged$cluster


## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusterlevels))
names(colorassigned)<-clusterlevels
color_list = list(clusters=colorassigned)

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_tcell$cluster1m_refined)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

### Figure 6D ####
pdf("./output/Fig6D.pdf",width=10,height=10)
fig6D <- Heatmap(as.matrix(cluster_mean_merged_mat),
                 column_title="Heatmap Annotated Clusters",
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
draw(fig6D)
dev.off()


## frequency bar plots =======

# Examining T cell subtypes by Tumor cell subtypes 
# Function to get indices based on cluster type

# T cell clusterlevels
clusterlevels=c("CD8+",
                "CD4+",
                "Cytotoxic CD8",
                "Activated CD4",
                "Memory CD4",
                "Naïve CD4",
                "Naïve CD8",
                "Treg")

get_Tcell_indices <- function(cluster_name) {
  data_subset <- data_tcell[data_tcell$cluster1m_refined== cluster_name, ]
  paste(data_subset$ImageId, data_subset$CellId, sep = "_")
}


# finer T cell indices 
Tcell_indices_list <- lapply(clusterlevels, get_Tcell_indices)
names(Tcell_indices_list) <- clusterlevels

data_tumor_clean_treated <- data_tumor_clean_ext[data_tumor_clean_ext$ImageId%in%treated_sampleId, ]
data_tumor_clean_untreated <- data_tumor_clean_ext[data_tumor_clean_ext$ImageId%in%untreated_sampleId, ]

tumor_near_results_combined  <- get_near_result(data_tumor_clean_ext,        Tcell_indices_list)
tumor_near_results_treated   <- get_near_result(data_tumor_clean_treated,    Tcell_indices_list)
tumor_near_results_untreated <- get_near_result(data_tumor_clean_untreated,  Tcell_indices_list)

tumor_near_Tcell_sum_combined  <- summarize_near(tumor_near_results_combined)
tumor_near_Tcell_sum_treated   <- summarize_near(tumor_near_results_treated)
tumor_near_Tcell_sum_untreated <- summarize_near(tumor_near_results_untreated)

to_long_df <- function(df, condition_label) {
  df %>%
    tibble::rownames_to_column(var = "cluster1") %>%
    reshape2::melt(id.vars = "cluster1", variable.name = "cluster2", value.name = "Frequency") %>%
    dplyr::mutate(Condition = condition_label)
}

plot_data_long <- bind_rows(
  to_long_df(tumor_near_Tcell_sum_combined,  "Combined"),
  to_long_df(tumor_near_Tcell_sum_treated,   "Treated"),
  to_long_df(tumor_near_Tcell_sum_untreated, "Untreated")
)

# Immune subtypes color:
colorassigned<- c("#FFBE2D",
                  "#80C7EF",
                  "#00F6B3",
                  "#06A5FF", 
                  "#F4EB71",
                  "#FF8320",
                  "#D99BBD", 
                  "#4D4D4D")
names(colorassigned)<- clusterlevels
plot_data_long$cluster2<- factor(plot_data_long$cluster2, levels=clusterlevels)

# Create a stacked bar chart with custom colors
fig6E_1<- ggplot(plot_data_long[plot_data_long$Condition=="Combined", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Combined") +
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

fig6E_2<-ggplot(plot_data_long[plot_data_long$Condition=="Treated", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Treated") +
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

fig6E_3<-ggplot(plot_data_long[plot_data_long$Condition=="Untreated", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Untreated") +
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

### Figure 6E & Figure S12A & B  ####
pdf('./output/Fig6E.pdf', height=3.5, width=3)
print(fig6E_1)
print(fig6E_3) # FigS12A
print(fig6E_2) # FigS12B
dev.off()


# Examining T cell subtypes by CAF subtypes
data_caf_clean_treated <- data_caf_clean_ext[data_caf_clean_ext$ImageId%in%treated_sampleId, ]
data_caf_clean_untreated <- data_caf_clean_ext[data_caf_clean_ext$ImageId%in%untreated_sampleId, ]

CAF_near_results_combined  <- get_near_result(data_caf_clean_ext,        Tcell_indices_list)
CAF_near_results_treated   <- get_near_result(data_caf_clean_treated,    Tcell_indices_list)
CAF_near_results_untreated <- get_near_result(data_caf_clean_untreated,  Tcell_indices_list)

CAF_near_Tcell_sum_combined  <- summarize_near(CAF_near_results_combined)
CAF_near_Tcell_sum_treated   <- summarize_near(CAF_near_results_treated)
CAF_near_Tcell_sum_untreated <- summarize_near(CAF_near_results_untreated)

plot_data_long <- bind_rows(
  to_long_df(CAF_near_Tcell_sum_combined,  "Combined"),
  to_long_df(CAF_near_Tcell_sum_treated,   "Treated"),
  to_long_df(CAF_near_Tcell_sum_untreated, "Untreated")
)

plot_data_long$cluster2<- factor(plot_data_long$cluster2, levels=clusterlevels)
plot_data_long<- plot_data_long[plot_data_long$cluster1!="CAF undefined", ] # remove CAF undefined
plot_data_long$cluster1 <- factor(plot_data_long$cluster1, levels=c("iCAF",
                                                                    "myCAF",
                                                                    "CD105+ CAF",
                                                                    "ApCAF",
                                                                    "CXCL12+ CAF",
                                                                    "FAP+ PDPN+ HLADR+ CAF")) # manual arrangement of CAF level

fig6F_1<- ggplot(plot_data_long[plot_data_long$Condition=="Combined", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Combined") +
  theme(axis.text.x = element_text(angle = 45,vjust =1, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

fig6F_2<-ggplot(plot_data_long[plot_data_long$Condition=="Treated", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Treated") +
  theme(axis.text.x = element_text(angle = 45,vjust =1, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

fig6F_3<-ggplot(plot_data_long[plot_data_long$Condition=="Untreated", ], aes(x = cluster1, y = Frequency, fill = cluster2)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  ggtitle("Untreated") +
  theme(axis.text.x = element_text(angle = 45,vjust =1, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 

### Figure 6F & Figure S12C & D  ####
pdf('./output/Fig6F.pdf', height=4, width=3)
print(fig6F_1)
print(fig6F_3) # FigS12C
print(fig6F_2) # FigS12D
dev.off()

## Assessing cytokine secretion from classical, basal, or changer CAFs ##

##### Looking at combined biological replicates ######

# Read data
profiler <- read_xlsx('./Config/revision_CAF_CM_Zscore_COMBO.xlsx')

# Set first column as rownames
profiler <- profiler %>% column_to_rownames(var = colnames(profiler)[1])

# Check if numeric
if(!all(apply(profiler, 2, is.numeric))) stop("Data contains non-numeric columns!")

# Scale rows
profiler_scaled <- as.data.frame(t(scale(t(profiler))))

#Take out NAs
profiler_clean <- profiler_scaled[rowSums(is.na(profiler_scaled)) == 0, ]

# order by changer 
profiler_scaled_ordered <- profiler_clean[order(profiler_clean$Changer, decreasing = TRUE), ]

mat_clean <- as.matrix(profiler_scaled_ordered)

### Figure S11A #### 
pS11a<- Heatmap(mat_clean,
        cluster_rows = FALSE, # do not cluster row if want to arrange by changer
        cluster_columns = TRUE,
        name = "Scaled",
        row_names_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 40),
        width = ncol(mat_clean) * unit(70, "mm"),
        height = nrow(mat_clean) * unit(7, "mm"),
        column_names_rot = 0, 
        column_names_centered = TRUE,
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 15), legend_width = unit(30, "cm")))

pdf('./output/ALL_ProteomeProf.pdf', height=35, width = 25)
draw(pS11a)                             
dev.off()
                                    
