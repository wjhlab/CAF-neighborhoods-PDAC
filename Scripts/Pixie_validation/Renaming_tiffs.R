#Nicole Gross
#R version 4.4.3 (2025-02-28)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sequoia 15.5

##Rename single page tiffs to simplified channel names
rm(list = ls())
library(readxl)
library(fs)       
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work <- getwd()

# Load channel names from Excel
channel_names <- read_excel("channels.xlsx", col_names = TRUE)
colnames(channel_names) <- c("Original", "New")

# Define input and output directories
input_dir <- "Original_TIFFs"
output_dir <- "Renamed_TIFFs"
dir_create(output_dir)

# List all subdirectories (one per sample)
sample_dirs <- dir_ls(input_dir, type = "directory")

# Loop through each sample folder
for (sample_path in sample_dirs) {
  sample_id <- path_file(sample_path)  # Just the folder name
  output_sample_path <- path(output_dir, sample_id)
  dir_create(output_sample_path)
  
  # List all TIFFs in the current sample folder
  tiff_files <- dir_ls(sample_path, regexp = "\\.ome\\.tiff$", type = "file")
  
  for (file_path in tiff_files) {
    original_name <- path_file(file_path)
    
    # Match with Excel sheet
    if (original_name %in% channel_names$Original) {
      new_name <- paste0(channel_names$New[channel_names$Original == original_name], ".tiff") #define what file extension to use
      new_path <- path(output_sample_path, new_name)
      
      file_copy(file_path, new_path, overwrite = TRUE)
    } else {
      warning(paste("No matching new name for", original_name))
    }
  }
}
