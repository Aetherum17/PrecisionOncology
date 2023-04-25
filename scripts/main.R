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
library(htmltools)
library(htmlTable)
library(ggplot2)

# Functions ####################################################################

# This functions tries to find potential targets in the list of mutations for each sample
find_applicable_drugs <- function(data_drugs, samples_df){
  # Create an output data frame
  patients_df <- data.frame('Sample_ID' = "Z", 'Diagnosis' = 'X', 'Target' = "C", 'Drug' = 'V', 'Found_Drug' = 'B')
  
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
        patients_df <- rbind(patients_df, c(samples_df[i,1], samples_df[i,2], data_drugs[q, 2], data_drugs[q, 1], TRUE))
        found_drug = T
      }
    }
    if(found_drug == F){
      patients_df <- rbind(patients_df, c(samples_df[i,1], samples_df[i,2], "---", "---", FALSE))
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
summary_path <- paste(directory, "results/summary.html", sep = "/")
summary_txt_path <- paste(directory, "results/summary.txt", sep = "/")
file.create(outputfile_path)
write(c("Sample_ID;Diagnosis;Target;Drug;Found_Drug"), outputfile_path, append=TRUE)
file.create(logfile_path)
file.create(summary_path)
file.create(summary_txt_path)
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
    data_df <- read_tsv(paste(tar_directory, "data_mutations.txt", sep = "/"), show_col_types = FALSE)
    if("Hugo_Symbol" %in% colnames(data_df) & "Tumor_Sample_Barcode" %in% colnames(data_df)){
      write("# OK: Hugo_Symbol and Tumor_Sample_Barcode Columns exist.", logfile_path, append=TRUE)
      # Select only the columns we need
      data_df <- data_df[,c("Hugo_Symbol", "Tumor_Sample_Barcode")] 
      # Check if the Clinical Sample file Exists
      if(file.exists(paste(tar_directory, "data_clinical_sample.txt", sep = "/")) == T){
        write(paste("# OK: ", tar_file_folder, " - Clinical Sample File Exists", sep = ""), logfile_path, append=TRUE)
        # Read Clinical Sample file
        data_clinical_df <- read.table(paste(tar_directory, "data_clinical_sample.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        if("SAMPLE_ID" %in% colnames(data_clinical_df) & "CANCER_TYPE" %in% colnames(data_clinical_df) & "CANCER_TYPE_DETAILED" %in% colnames(data_clinical_df)){
          write("# OK: SAMPLE_ID and CANCER_TYPE_DETAILED Columns exist.", logfile_path, append=TRUE)
          # Select only the columns we need
          data_clinical_df <- data_clinical_df[,c("SAMPLE_ID", "CANCER_TYPE", "CANCER_TYPE_DETAILED")]
          # Rename the Columns
          colnames(data_clinical_df) <- c("Tumor_Sample_Barcode", "Cancer_Type", "Cancer_Type_Detailed")
          # Merge Cancer Type and Cancer Type Detailed
          data_clinical_df$Cancer_Type_Detailed <- paste(data_clinical_df$Cancer_Type, data_clinical_df$Cancer_Type_Detailed, sep = ": ")
          data_clinical_df <- subset(data_clinical_df, select = -c(Cancer_Type))
          # Merge Data of Mutations with Clinical Samples
          data_merged_df <- merge(data_df, data_clinical_df, by="Tumor_Sample_Barcode")
          # Connect each sample with all mutations discovered
          samples_df <- aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode + Cancer_Type_Detailed, data_merged_df, function(x) paste(x, collapse=", "))
          # Rename the Columns
          colnames(samples_df) <- c("Sample_ID", "Cancer Type", "Mutated_Genes")
          ### Find Applicable Drugs
          output_df <- find_applicable_drugs(data_drugs, samples_df)
          write.table(output_df, file=outputfile_path, sep=";", quote=F, row.names=F, col.names = F, na="NA", append=TRUE)
          write.table(samples_df, file=summary_txt_path, sep=";", quote=F, row.names=F, col.names = F, na="NA", append=TRUE)
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

# Load the Output data
output_df <- read.csv(outputfile_path, row.names=NULL, sep = ";")
head(output_df)

### How Many people for each Diagnosis can be assigned to Drug therapy?
output_df_1 <- output_df[,c("Sample_ID", "Diagnosis", "Found_Drug")]
# Some patients may be suited for multiple drug therapy, so we remove the duplicated lines
output_df_1 <- output_df_1[!duplicated(output_df_1),]
# Drop Sample ID column, as we do not need it anymore
output_df_1 <- subset(output_df_1, select = -c(Sample_ID))
# Create table
output_df_1_table <- table(output_df_1)
output_df_1_table_prop <- round(prop.table(table(output_df_1)), 2)*100
# Add some of rows and columns
output_df_1_table <- cbind(output_df_1_table, rowSums(output_df_1_table))
output_df_1_table <- rbind(output_df_1_table, c(colSums(output_df_1_table), sum(output_df_1_table)))
# Rename Rows and Columns
colnames(output_df_1_table) <- c("No Drug Found", "Drug Exists", "Sum")
rownames(output_df_1_table)[nrow(output_df_1_table)] <- "Sum"
colnames(output_df_1_table_prop) <- c("No Drug Found (%)", "Drug Exists (%)")
# Graph
output_df_1_table_graph <- as.data.frame(output_df_1_table)[-nrow(output_df_1_table),-ncol(output_df_1_table)]
output_df_1_table_graph$Diagnosis <- row.names(output_df_1_table_graph)
output_df_1_table_graph <- output_df_1_table_graph %>% pivot_longer(cols=c('No Drug Found', "Drug Exists"),
                                         names_to='Type',
                                         values_to='Patients')
output_df_1_table_graph$DiagnosisType <- paste(output_df_1_table_graph$Diagnosis, output_df_1_table_graph$Type, sep = ":\n")
# Make the ggplot
options(repr.plot.width =20, repr.plot.height =20)
output_df_1_table_graph$fraction = output_df_1_table_graph$Patients / sum(output_df_1_table_graph$Patients)
output_df_1_table_graph$ymax = cumsum(output_df_1_table_graph$fraction)
output_df_1_table_graph$ymin = c(0, head(output_df_1_table_graph$ymax, n=-1))
output_df_1_table_graph$labelPosition <- (output_df_1_table_graph$ymax + output_df_1_table_graph$ymin) / 2
ggplot(output_df_1_table_graph, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=DiagnosisType)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=Patients), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "bottom")+ scale_fill_brewer(palette="Dark2")+ggtitle("Patients with Mature B-Cell Neoplasms")

### What are the most common drugs used?
output_df_2 <- output_df[,c("Drug")]
output_df_2_table <- sort(table(output_df_2), decreasing=TRUE)
names(output_df_2_table)[1] <- "No Drug Found"
output_df_2 <- as.data.frame(output_df_2_table)
output_df_2$Percent <- round(output_df_2$Freq/sum(output_df_2$Freq)*100, 2)
colnames(output_df_2) <- c("Drug", "Number of Patients", "Number of Patients (%)")

### What are the most common targets used?
output_df_3 <- output_df[,c("Target")]
output_df_3_table <- sort(table(output_df_3), decreasing=TRUE)
names(output_df_3_table)[1] <- "No Target Found"
output_df_3 <- as.data.frame(output_df_3_table)
output_df_3$Percent <- round(output_df_3$Freq/sum(output_df_3$Freq)*100, 2)
colnames(output_df_3) <- c("Target", "Number of Patients", "Number of Patients (%)")

### What are the most frequent mutations in the samples?
output_df_4 <- data_df[,c("Hugo_Symbol")]
output_df_4_table <- sort(table(output_df_4), decreasing = TRUE)
output_df_4 <- as.data.frame(output_df_4_table)
output_df_4$Percent <- round(output_df_4$Freq/sum(output_df_4$Freq)*100, 2)
output_df_4 <- output_df_4[output_df_4$Freq>4,]
colnames(output_df_4) <- c("Mutated Gene", "Frequency", "In Number of Samples (%)")

################################################################################

summary_file <- file(summary_path, open = "w")

# Create HTML File
html <- tags$html(
  tags$body(
    tags$h1("Statistics of Targeted Drugs for Hematological Malignancies"),
    tags$h2("How Many Samples for each Diagnosis can be assigned to the Drug therapy?"),
    htmlTable(output_df_1_table, css.cell = "padding: 0 50px;"),
    htmlTable(output_df_1_table_prop, css.cell = "padding: 0 70px;"),
    tags$h2("What are the most common drugs used?"),
    htmlTable(output_df_2, css.cell = "padding: 0 50px;"),
    tags$h2("What are the most common targets used?"),
    htmlTable(output_df_3, css.cell = "padding: 0 50px;"),
    tags$h2("What are the most frequent mutations in the samples?"),
    htmlTable(output_df_4, css.cell = "padding: 0 50px;"),
  )
)

writeLines(as.character(html), summary_file)
close(summary_file)

found_drugs <- unique(output_df[,c("Drug")])
db_drugs <- unique(data_drugs$drug)
print(paste("Not applied drugs from DB: ", found_drugs[which(!(found_drugs %in% db_drugs))]))

################################################################################

print("DONE!")
