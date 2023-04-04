library(readxl)
library(tools)
library(dplyr)
library(tidyr)
library(data.table)
library(gdata)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(stringr)
library(rstudioapi)
library(R.utils)
library(tools)
library(readr)

# Functions ####################################################################

# This functions tries to find potential targets in the list of mutations for each sample
find_applicable_drugs <- function(data_drugs, samples_df){
  # Create an output data frame
  patients_df <- data.frame('Sample ID' = "X", 'Diagnosis' = 'Q', 'Target' = "Y", 'Drug' = 'Z')
  
  # For each SAMPLE
  for(i in 1:nrow(samples_df)){
    # Select the mutations
    mutations <- strsplit(samples_df[i,3], ', ')[[1]]
    found_drug = F
    # FOR each DRUG
    for(q in 1:nrow(data_drugs)){
      # Check if a drug is in the mutations list
      if(data_drugs[q, 2] %in% mutations) {
        # if it is there - add information about the finding
        patients_df <- rbind(patients_df, c(samples_df[i,1], samples_df[i,2], data_drugs[q, 2], data_drugs[q, 1]))
        found_drug = T
      }
    }
    if(found_drug == F){
      patients_df <- rbind(patients_df, c(samples_df[i,1], samples_df[i,2], "---", "---"))
    }
  }
  patients_df <- patients_df[-1,]
  return(patients_df)
}

# Main Body ####################################################################

# Set Dynamic Working Directory
directory <- getSourceEditorContext()$path
directory <- str_replace(directory, "/[^/]*$", "")
directory <- str_replace(directory, "/[^/]*$", "")
setwd(directory)

# Output and log files
outputfile_path <- paste(directory, "results/output.txt", sep = "/")
logfile_path <- paste(directory, "results/log.log", sep = "/")
file.create(outputfile_path)
file.create(logfile_path)
unlink(paste(directory, "data/temp", sep = "/"), recursive = TRUE)

# Drugs Data
data_drugs <- read.csv(paste(directory, "database/data_drugs.txt", sep = "/"), sep = ",")
data_drugs <- aggregate(target ~ drug, data_drugs, function(x) paste(x, collapse=", "))

# Get list of all files in data folder
data_files_list <- list.files(paste(directory, "data", sep = "/"))

# For each file in data folder
for (file_num in 1:length(data_files_list)){
  # load current data file
  data_file <- paste(directory, "data", data_files_list[file_num], sep = "/")
  # Get directory of the data file
  data_file_extension = file_ext(data_file)
  # If it is an archive - it should be unpacked
  if(data_file_extension == "gz")
  {
    tar_file <- gunzip(data_file)
  } else{
    tar_file <- data_file
  }
  # Get the name of tar file folder
  tar_file_folder <- basename(file_path_sans_ext(tar_file))
  # The directory where the tar filw will be extracted
  ex_dir = paste(directory, "data/temp", data_files_list[file_num], sep = "/")
  # Extract the tar file
  untar(tarfile = tar_file, exdir = ex_dir)
  # Directory of the folder containing mutations file
  tar_directory <- paste(directory, "data/temp", basename(tar_file), tar_file_folder, sep = "/")
  # Check if the Mutations file Exist
  if(file.exists(paste(tar_directory, "data_mutations.txt", sep = "/")) == T){
    write(paste("# OK: ", tar_file_folder, " - Data Mutations Exists", sep = ""), logfile_path, append=TRUE)
    # Read the data_muations file
    # data_df <- read.table(paste(tar_directory, "data_mutations.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    data_df <- read_tsv(paste(tar_directory, "data_mutations.txt", sep = "/"))
    if("Hugo_Symbol" %in% colnames(data_df) & "Tumor_Sample_Barcode" %in% colnames(data_df)){
      write("# OK: Hugo_Symbol and Tumor_Sample_Barcode Columns exist.", logfile_path, append=TRUE)
      # Select only the columns we need
      data_df <- data_df[,c("Hugo_Symbol", "Tumor_Sample_Barcode")] 
      # Check if the Clinical Sample file Exists
      if(file.exists(paste(tar_directory, "data_clinical_sample.txt", sep = "/")) == T){
        write(paste("# OK: ", tar_file_folder, " - Clinical Sample File Exists", sep = ""), logfile_path, append=TRUE)
        # Read Clinical Sample file
        data_clinical_df <- read.table(paste(tar_directory, "data_clinical_sample.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        if("SAMPLE_ID" %in% colnames(data_clinical_df) & "CANCER_TYPE_DETAILED" %in% colnames(data_clinical_df)){
          write("# OK: SAMPLE_ID and CANCER_TYPE_DETAILED Columns exist.", logfile_path, append=TRUE)
          # Select only the columns we need
          data_clinical_df <- data_clinical_df[,c("SAMPLE_ID", "CANCER_TYPE_DETAILED")]
          # Rename the Columns
          colnames(data_clinical_df) <- c("Tumor_Sample_Barcode", "Cancer_Type_Detailed")
          # Merge Data of Mutations with Clinical Samples
          data_merged_df <- merge(data_df, data_clinical_df, by="Tumor_Sample_Barcode")
          # Connect each sample with all mutations discovered
          samples_df <- aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode + Cancer_Type_Detailed, data_merged_df, function(x) paste(x, collapse=", "))
          # Rename the Columns
          colnames(samples_df) <- c("Sample_ID", "Cancer Type", "Mutated_Genes")
          ### Find Applicable Drugs
          output_df <- find_applicable_drugs(data_drugs, samples_df)
          print(output_df)
        }
        else{
          write("# ERROR: Some of the required Columns DO NOT Exist.", logfile_path, append=TRUE)
        }
      } else{
        write(paste("# ERROR", tar_file_folder, " - Clinical Sample File is MISSING", sep = ""), logfile_path, append=TRUE)
      }
    } else{
      write("# ERROR: Some of the required Columns DO NOT Exist.", logfile_path, append=TRUE)
    }
  } else {
    write(paste("# ERROR", tar_file_folder, " - Data Mutations File is MISSING", sep = ""), logfile_path, append=TRUE)
  }
  unlink(paste(directory, "data/temp", sep = "/"), recursive = TRUE)
}

################################################################################


# write(paste("##### Mutations of Diffuse Large B Cell Lymphoma ##### "), output_file_path, append=TRUE)
# write(paste("### Number of Samples:", length(unique(data_diffuse_large_b_cell_type$Tumor_Sample_Barcode)), sep = " "), output_file_path, append=TRUE)
# write.fwf(cancer_type_df, file=output_file_path, sep=" ", quote=F, rownames=F, colnames = T, na="NA", append=TRUE)



