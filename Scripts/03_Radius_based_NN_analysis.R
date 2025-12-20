#######################################################
#   03. Neighborhood analysis based on radius         #
#######################################################
## neighbors within radius 10um and 30um distance

## calculate nearest-neighbor distance 
# euclidean dist  sqrt((x1 - x2)^2 + (y1 - y2)^2))]
select_cols <- c("X_coord", "Y_coord", "cluster1m_refined", "ImageId", "CellId")

combined_data<- rbind(data_caf[, select_cols], data_tumor[, select_cols])

# add 2m refined cluster annotation for CAF
combined_data$cluster2m_refined<- combined_data$cluster1m_refined
combined_data[combined_data$cluster1m_refined%in%unique(data_caf$cluster1m_refined), "cluster2m_refined"]<- data_caf$cluster2m_refined

combined_data$Id_ext<- paste0(combined_data$ImageId, "_",combined_data$CellId)
unique(combined_data$cluster1m_refined)

# remove unassigned tumor 
combined_data_filtered <- combined_data%>% 
  filter(cluster1m_refined!="Unassigned")

# first create point pattern data form (separately for each image ID)
# finding dist to Tumor types from CAF point of view 
# CAF ids 
CAF_Id<- combined_data_filtered[combined_data_filtered$cluster1m_refined%in%unique(data_caf$cluster1m_refined), "Id_ext"]
# TUMOR ids 
TUM_Id<- combined_data_filtered[combined_data_filtered$cluster1m_refined%in%unique(data_tumor$cluster1m_refined), "Id_ext"]

# calculate euclidean distance
dist_df <- c()
for (i in samplevels){
  df_k <- combined_data_filtered[combined_data_filtered$ImageId==i, ] # per sample
  pp <- with(df_k,
             ppp(x = X_coord,
                 y = Y_coord,
                 window = owin(range(X_coord), range(Y_coord)),     # bounding box
                 marks  = Id_ext)                                 # keep IDs as marks
  )
  
  caf   <- subset(pp, marks(pp) %in% CAF_Id)
  tumor <- subset(pp, marks(pp) %in% TUM_Id)
  
  r <- 100 # set max distance between pairs of points 
  xcp <- crosspairs(caf, tumor, rmax = r, what = "all")   # return i,j indices, x,y-coordinates, distance 
  
  result_df <- data.table(from     = marks(caf)[xcp$i],
                          to       = marks(tumor)[xcp$j],
                          distance = xcp$d)
  dist_df <- rbind(dist_df, result_df)
}


dist_df <- merge(dist_df, combined_data_filtered[, c("Id_ext", "cluster1m_refined")], 
                 by.x = "from", by.y = "Id_ext", all.x = TRUE)
colnames(dist_df)[4]<- "from_celltype"


dist_df <- merge(dist_df, combined_data_filtered[, c("Id_ext", "cluster1m_refined")], 
                 by.x = "to", by.y = "Id_ext", all.x = TRUE)
colnames(dist_df)[5]<- "to_celltype"

# add 2m clustering annotation 
dist_df <- merge(dist_df, combined_data_filtered[, c("Id_ext", "cluster2m_refined")], 
                 by.x = "from", by.y = "Id_ext", all.x = TRUE)

# filter radius 10um (CAF as the reference)

dist_df_10um <- dist_df[dist_df$distance<=10, ]
interaction_count_table10<-table(dist_df_10um$from_celltype, dist_df_10um$to_celltype)

# 10um radius
### Figure S5D ####
pdf('./output/FigS5D.pdf')
FigS5D<- Heatmap(t(scale(t(as.matrix(interaction_count_table10)))),
        name = "Scaled\nFreq",
        width = ncol(interaction_count_table10)*unit(7, "mm"), 
        height = nrow(interaction_count_table10)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
draw(FigS5D)
dev.off()


# filter radius 30um (CAF as the reference)

dist_df_30um <- dist_df[dist_df$distance<=30, ]
interaction_count_table30<-table(dist_df_30um$from_celltype, dist_df_30um$to_celltype)

# 30um radius
### Figure S5E ####
pdf('./output/FigS5E.pdf')
FigS5E<- Heatmap(t(scale(t(as.matrix(interaction_count_table30)))),
        name = "Scaled\nFreq",
        width = ncol(interaction_count_table30)*unit(7, "mm"), 
        height = nrow(interaction_count_table30)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
draw(FigS5E)
dev.off()


### Fig S4C ####
### Visualize CAF undefined
custom_colors <- c("ApCAF" = "#E69F00", "myCAF" = "#D55E00", "iCAF" = "#AD7700", 
                   "CD105+ CAF"="#F0E442", "CXCL12+ CAF"="#CC79A7", "FAP+ CAF"="#B14380", "CAF undefined"="#666666",
                   "Basal"="#154d6e","Classical"="#0979ba", "Mixed"="#c4e1f2")  # modify as needed

colorlevel<- c("ApCAF", "myCAF", "iCAF", "CD105+ CAF", "CXCL12+ CAF", "FAP+ CAF",
               "CAF undefined","Basal","Classical", "Mixed")

# visualization of spatial organization
pdf('./output/CAF_Tumor_spatial organization.pdf', height=25, width = 28)
CAF_tumor<-ggplot(combined_data_filtered, aes(x = X_coord, y = Y_coord, color = factor(cluster1m_refined, levels=colorlevel))) +
  geom_point(size = 1, alpha=0.7) +
  labs(color = "Cluster") +
  scale_color_manual(values = custom_colors) +
  scale_y_reverse()+ # flip y-axis to match MCD image 
  theme_minimal()+
  facet_wrap(~ImageId)
print(CAF_tumor)
dev.off()


select_cols <- c("ImageId","X_coord", "Y_coord", "cluster1m_refined", "Collagen", "SMA")
vis_data<- data_caf[, select_cols]
unique(vis_data$cluster1m_refined)

pdf('./output/CAF_undefined_Collagen.pdf', height=25, width = 28)
CAF_undef<- ggplot(vis_data[vis_data$cluster1m_refined%in%c("myCAF","CAF undefined"), ], 
       aes(x = X_coord, y = Y_coord, color = Collagen)) +
  geom_point(size = 1) +
  scale_color_gradientn(
    colors = c("blue", "yellow", "red"),
    values = scales::rescale(c(min(vis_data$Collagen, na.rm = TRUE),
                               median(vis_data$Collagen, na.rm = TRUE),
                               max(vis_data$Collagen, na.rm = TRUE))),
    name = "Collagen"
  )+
  theme_minimal()+
  scale_y_reverse()+
  facet_wrap(~ImageId)
print(CAF_undef)
dev.off()


