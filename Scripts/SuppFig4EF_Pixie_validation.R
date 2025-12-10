##Validation of phenograph IMC clusters with pixie

#Nicole Gross 
#R version 4.4.3 (2025-02-28)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sequoia 15.5

rm(list = ls())

library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(readr)
library(readxl)
library(tidyverse)
library(limma)
library(ggrepel)
library(stringr)  
library(ggpubr)
library(reshape2)
library(tibble)
library(pals)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

#import counts docs 
pixie<- read.csv(paste0(work,'/config/pixel_meta_cluster_summary.csv'))
counts<- read.csv(paste0(work,'/config/counts_broad.csv'),check.names = FALSE, row.names = 1)

# rename pixie file names with sample_id using metadata

#set up metadata
pixie_md = paste0(work,"/config/pixie_metadata.xlsx") 
pixie_md <- read_excel(pixie_md)

md = paste0(work,"/config/metadata.xlsx") 
md <- read_excel(md)

#add sample_id to pixie df
pixie <- pixie %>%
  left_join(pixie_md %>% select(file_name, sample_id), by = "file_name") %>%
  select(-file_name) %>%  #remove old file_name column
  rename(sample_id = sample_id) #ensure column is named sample_id

pixie <- pixie %>%
  select(sample_id, everything()) #move sample_id to first column

write.csv(pixie, file='pixie.csv')

pixie_wide <- pixie %>%
  group_by(cluster, sample_id) %>%
  summarise(pixel_count = sum(pixel_count), .groups = "drop") %>%
  pivot_wider(
    names_from = sample_id,
    values_from = pixel_count,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "cluster")


write.csv(pixie_wide, file='pixie_wide_broad.csv')

#calculate proportions
props_tablep <- t(t(pixie_wide) / colSums(pixie_wide)) * 100
props_p <- as.data.frame.matrix(props_tablep)

props_table <- t(t(counts) / colSums(counts)) * 100
props <- as.data.frame.matrix(props_table)

#set up levels
clusterlevels=c("Immune unassigned",
                "Myeloid",
                "B cell",
                "CD57",
                "Tumor",
                "CAF",
                "CD4",
                "CD8",
                "Treg")

samplevels_orig=c("317_S_hi_i_hi",
              "317_S_hi_i_lo",
              "317_T_hi_i_hi",
              "317_T_hi_i_lo",
              "348_S_hi_i_hi",
              "348_S_hi_i_lo",
              "348_T_hi_i_hi",
              "348_T_hi_i_lo",
              "352_S_hi_i_hi",
              "352_S_hi_i_lo",
              "352_T_hi_i_hi",
              "352_T_hi_i_lo",
              "357_S_hi_i_hi",
              "357_S_hi_i_lo",
              "357_T_hi_i_lo",
              "361_S_hi_i_hi",
              "361_S_hi_i_lo",
              "361_T_hi_i_hi",
              "361_T_hi_i_lo",
              "362_S_hi_i_hi",
              "362_S_hi_i_lo",
              "362_T_hi_i_hi",
              "362_T_hi_i_lo",
              "368_S_hi_i_hi",
              "368_S_hi_i_lo",
              "368_T_hi_i_hi",
              "368_T_hi_i_lo",
              "369_S_hi_i_hi",
              "369_S_hi_i_lo",
              "369_T_hi_i_hi",
              "369_T_hi_i_lo",
              "372_S_hi_i_hi",
              "372_S_hi_i_lo",
              "372_T_hi_i_hi",
              "380_S_hi_i_hi",
              "380_S_hi_i_lo",
              "380_T_hi_i_hi",
              "380_T_hi_i_lo",
              "383_S_hi_i_hi",
              "383_S_hi_i_lo",
              "383_T_hi_i_hi",
              "387_S_hi_i_hi",
              "387_S_hi_i_lo",
              "387_T_hi_i_hi",
              "387_T_hi_i_lo",
              "388_S_hi_i_hi",
              "388_S_hi_i_lo",
              "388_T_hi_i_hi",
              "388_T_hi_i_lo",
              "390_S_hi_i_lo",
              "390_T_hi_i_hi",
              "390_T_hi_i_lo",
              "417_S_hi_i_hi",
              "417_S_hi_i_lo",
              "417_T_hi_i_hi",
              "417_T_hi_i_lo")

samplevels_pixie=c("317_S_hi_i_hi",
              "317_S_hi_i_lo",
              "317_T_hi_i_hi",
              "317_T_hi_i_lo",
              "348_S_hi_i_hi",
              "348_S_hi_i_lo",
              "348_T_hi_i_hi",
              "348_T_hi_i_lo",
              "352_S_hi_i_hi",
              "352_S_hi_i_lo",
              "352_T_hi_i_hi",
              "352_T_hi_i_lo",
              "357_S_hi_i_hi",
              "357_S_hi_i_lo",
              "357_T_hi_i_lo",
              "361_S_hi_i_hi",
              "361_S_hi_i_lo",
              "361_T_hi_i_hi",
              "361_T_hi_i_lo",
              "362_S_hi_i_hi",
              "362_S_hi_i_lo",
              "362_T_hi_i_hi",
              "362_T_hi_i_lo",
              "368_S_hi_i_hi",
              "368_S_hi_i_lo",
              "368_T_hi_i_hi",
              "368_T_hi_i_lo",
              "369_S_hi_i_hi",
              "369_S_hi_i_lo",
              "369_T_hi_i_hi",
              "369_T_hi_i_lo",
              "372_S_hi_i_hi",
              "372_S_hi_i_lo",
              "372_T_hi_i_hi",
              "380_S_hi_i_hi",
              "380_S_hi_i_lo",
              "380_T_hi_i_hi",
              "380_T_hi_i_lo",
              "383_S_hi_i_hi",
              "383_S_hi_i_lo",
              "383_T_hi_i_hi",
              "387_S_hi_i_hi",
              "387_S_hi_i_lo",
              "387_T_hi_i_hi",
              "387_T_hi_i_lo",
              "388_S_hi_i_hi",
              "388_S_hi_i_lo",
              "388_T_hi_i_hi",
              "388_T_hi_i_lo",
              "390_S_hi_i_lo",
              "390_T_hi_i_hi",
              "390_T_hi_i_lo",
              "417_S_hi_i_lo",
              "417_T_hi_i_hi",
              "417_T_hi_i_lo")

patientlevels=c("JHH317", "JHH348", "JHH352", "JHH357", "JHH361", "JHH362", "JHH368", "JHH369", "JHH372", "JHH380", "JHH383", "JHH387","JHH388", "JHH390", "JHH417")


## Set up the data frame for proportional plotting of ORIGINAL clusters
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels_orig)
ggdf$Name <- factor(md$Name[match(ggdf$sample_id,md$sample_id)], levels=patientlevels)

#set up data frame for proportional plotting of PIXIE clusters
ggdfp <- melt(data.frame(cluster = rownames(props_p),props_p, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdfp$sample_id <- factor(ggdfp$sample_id, levels=samplevels_pixie)
ggdfp$Name <- factor(md$Name[match(ggdfp$sample_id,md$sample_id)], levels=patientlevels)



# Load data
pixie_counts <- read.csv(paste0(work,'/pixie_wide_broad.csv'), row.names = 1, check.names = FALSE) 
pheno_counts <- read.csv(paste0(work,'/config/counts_broad.csv'), row.names = 1, check.names = FALSE)

#calculate proportions
props_tablep <- t(t(pixie_counts) / colSums(pixie_counts)) * 100
props_pixie <- as.data.frame.matrix(props_tablep)

props_table <- t(t(pheno_counts) / colSums(pheno_counts)) * 100
props_pheno <- as.data.frame.matrix(props_table)

# Find shared samples
shared_samples <- intersect(colnames(props_pixie), colnames(props_pheno))

# Find shared clusters
shared_clusters <- intersect(rownames(props_pixie), rownames(props_pheno))

# Subset both matrices to shared samples
props_pixie <- props_pixie[, shared_samples]
props_pheno <- props_pheno[, shared_samples]

# Subset both matrices to shared clusters
props_pixie <- props_pixie[shared_clusters,]
props_pheno <- props_pheno[shared_clusters,]

# Transpose (samples = rows, clusters = columns)
pixie_t <- t(props_pixie)
pheno_t <- t(props_pheno)

# Add method and sample_id columns
pixie_df <- as.data.frame(pixie_t)
pixie_df$method <- "Pixie"
pixie_df$sample_id <- rownames(pixie_df)
pixie_df$plot_id <- paste(pixie_df$sample_id, pixie_df$method, sep = "_")

pheno_df <- as.data.frame(pheno_t)
pheno_df$method <- "Phenograph"
pheno_df$sample_id <- rownames(pheno_df)
pheno_df$plot_id <- paste(pheno_df$sample_id, pheno_df$method, sep = "_")

# Combine into one data frame
combined <- rbind(pixie_df, pheno_df)

write_csv(combined,"combined_props_broad.csv")


##Correlation

#Get only the cluster columns (all numeric ones)
cluster_cols <- names(combined)[sapply(combined, is.numeric)]

#Split into Pixie and Phenograph data frames
pixie_df <- combined[combined$method == "Pixie", ]
pheno_df <- combined[combined$method == "Phenograph", ]

#Sort both by sample_id to match the order
pixie_df <- pixie_df[order(pixie_df$sample_id), ]
pheno_df <- pheno_df[order(pheno_df$sample_id), ]

#Extract cluster values and flatten to a vector
pixie_props <- as.vector(t(pixie_df[, cluster_cols]))
phenograph_props <- as.vector(t(pheno_df[, cluster_cols]))

# Create sample_id vector: repeat each sample_id for each cluster
sample_ids <- rep(pixie_df$sample_id, each = length(cluster_cols))

all.equal(length(pixie_props), length(phenograph_props))  #check vectors are equal in length

pearson_result<- cor.test(pixie_props, phenograph_props, method = "pearson")

r_p   <- round(pearson_result$estimate, 3)
p_p   <- signif(pearson_result$p.value, digits = 4)

#to annotate on plot
annotation_text <- paste0("Pearson r = ", r_p, ", p = ", p_p)

# Repeat cluster names for each sample (since you transposed then flattened)
cluster_names <- rep(cluster_cols, times = nrow(pixie_df))

comparison_df <- data.frame(
  Pixie = pixie_props,
  Phenograph = phenograph_props,
  cluster = cluster_names,
  sample_id = sample_ids)

comparison_df$patient_id <- substr(comparison_df$sample_id, 1, 3)

# Plot Pearson Corr

#color by cluster
pdf("SuppFig4E.pdf", width = 10, height = 10)
ggplot(comparison_df, aes(x = sqrt(Phenograph), y = sqrt(Pixie), color=cluster)) +
  geom_point(alpha = 0.5, size=3) +
  # Add smoothing lines with mapped aesthetics for legend
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black", linetype = "solid") +
  # Correlation annotation
  annotate("text", x = Inf, y = Inf, label = annotation_text,
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic") +
  labs(title = "Pixie vs Phenograph Cluster Proportions",
    x = "Phenograph Sqrt(% of cells)",
    y = "Pixie Sqrt(% of area)") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))
dev.off()

#color by patient
pdf("SuppFig4F.pdf", width = 10, height = 10)
ggplot(comparison_df, aes(x = sqrt(Phenograph), y = sqrt(Pixie), color = patient_id)) +
  geom_point(alpha = 0.5, size=3) +
  # Add smoothing lines with mapped aesthetics for legend
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black", linetype = "solid") +
  # Correlation annotation
  annotate("text", x = Inf, y = Inf, label = annotation_text,
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic") +
  labs(title = "Pixie vs Phenograph Cluster Proportions",
       x = "Phenograph Sqrt(% of cells)",
       y = "Pixie Sqrt(% of area)") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))
dev.off()
