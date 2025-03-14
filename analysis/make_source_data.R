# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2025 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2025-03-10
#
# Script Name:    ~/R-projects/reheat2_pilot/make_source_data.R
#
# Script Description:
# helper functions to set up source data of figures 

# create folder in project repo
create_source_data_folder <- function(source_data_list, output_dir = "Source_Data") {
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  

}

create_source_data_folder()
save_source_data <- function(fig_type = TRUE, fig_number, panel_letter, data, 
                             output_dir = "Source_Data", 
                             description = NULL, 
                             bottom_description = NULL, 
                             row_names = FALSE) {
  # Convert boolean fig_type to appropriate label
  fig_label <- ifelse(fig_type, "fig", "suppfig")
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Construct the filename based on input parameters
  file_name <- paste0(fig_label, "_", fig_number, panel_letter, ".csv")
  file_path <- file.path(output_dir, file_name)
  
  # Add description as the first row if provided, followed by an empty line
  if (!is.null(description)) {
    writeLines(c(description, ""), file_path)
  }
  
  # Write the main data table
  write.table(data, file_path, sep = ",", row.names = row_names, col.names = TRUE, append = TRUE)
  
  # Add bottom description if provided, preceded by an empty line
  if (!is.null(bottom_description)) {
    writeLines(c("", bottom_description), file_path, append = TRUE)
  }
  
  message("Saved: ", file_path)
}

# Example usage:

#df <- data.frame(x = 1:10, y = rnorm(10))
#save_source_data(fig_type =T, fig_number =  1, "b", description = "test_title", df, row.names = T)


save_source_data <- function(fig_type = TRUE, fig_number, panel_letter, data, output_dir = "Source_Data", description = NULL, bottom_description = NULL, row_names = FALSE) {
  # Convert boolean fig_type to appropriate label
  fig_label <- ifelse(fig_type, "fig", "suppfig")
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Construct the filename based on input parameters
  file_name <- paste0(fig_label, "_", fig_number, panel_letter, ".csv")
  file_path <- file.path(output_dir, file_name)
  
  # Open a connection to the file
  con <- file(file_path, "wt")
  
  # Write description if provided
  if (!is.null(description)) {
    writeLines(description, con)
    writeLines("", con)  # Empty line
  }
  
  close(con)  # Close before appending data
  
  # Write the main data table
  write.table(data, file_path, sep = ",", row.names = row_names, col.names = TRUE, append = TRUE)
  
  # Write bottom description if provided
  if (!is.null(bottom_description)) {
    write.table(data.frame(""), file_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)  # Empty line
    write.table(data.frame(bottom_description), file_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  message("Saved: ", file_path)
}


# update the .xslx file to bundle the source data -------------------------


# Load required libraries
library(readr)
library(openxlsx)

# Define input folder and output file
input_folder <- "../reheat2_pilot/Source_Data/"  # Change this to your folder path
output_file <- "../reheat2_pilot/Source_Data/LanzerRamirez_HFpatientmap_sourcedata.xlsx"

# List all CSV files in the folder
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each CSV file and add it as a sheet
# Loop through each CSV file and add it as a sheet
for (file in csv_files) {
  base_name <- tools::file_path_sans_ext(basename(file))  # Get filename without extension
  sheet_name <- gsub("Figure", "fig_", base_name)  # Replace "Figure" with "fig_"
  sheet_name <- substr(sheet_name, 1, 31)  # Ensure sheet name does not exceed Excel’s 31-char limit
  
  data <- read_csv(file)  # Read CSV file
  addWorksheet(wb, sheet_name)  # Create sheet
  writeData(wb, sheet_name, data)  # Write data
}
# Save workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

