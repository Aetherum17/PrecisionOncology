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

# Main Body ####################################################################

# Set Dynamic Working Directory
directory <- getSourceEditorContext()$path
directory <- str_replace(directory, "/[^/]*$", "")
directory <- str_replace(directory, "/[^/]*$", "")
setwd(directory)

# BCGSC Database ###############################################################

data_mutations_bcgsc <- read.csv(paste(directory, "data/data_mutations_bcgsc.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_bcgsc <- data_mutations_bcgsc[,colSums(is.na(data_mutations_bcgsc))<nrow(data_mutations_bcgsc)]
data_mutations_bcgsc$Databse <- c("BCGSC")

data_clinical_sample_bcgsc <- read.csv(paste(directory, "data/data_clinical_sample_bcgsc.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_clinical_sample_bcgsc <- data_clinical_sample_bcgsc[,colSums(is.na(data_clinical_sample_bcgsc))<nrow(data_clinical_sample_bcgsc)]
# Remove first 4 rows
data_clinical_sample_bcgsc <- tail(data_clinical_sample_bcgsc, -4)
colnames(data_clinical_sample_bcgsc)[2] ="Tumor_Sample_Barcode"

df_bcgsc <- merge(data_mutations_bcgsc, data_clinical_sample_bcgsc,by="Tumor_Sample_Barcode")

# Broad Database ###############################################################

data_mutations_broad <- read.csv(paste(directory, "data/data_mutations_broad.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_broad <- data_mutations_broad[,colSums(is.na(data_mutations_broad))<nrow(data_mutations_broad)]
data_mutations_broad$Databse <- c("Broad")

data_clinical_sample_broad <- read.csv(paste(directory, "data/data_clinical_sample_broad.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_clinical_sample_broad <- data_clinical_sample_broad[,colSums(is.na(data_clinical_sample_broad))<nrow(data_clinical_sample_broad)]
# Remove first 4 rows
data_clinical_sample_broad <- tail(data_clinical_sample_broad, -4)
colnames(data_clinical_sample_broad)[2] ="Tumor_Sample_Barcode"

df_broad <- merge(data_mutations_broad, data_clinical_sample_broad,by="Tumor_Sample_Barcode")

# TCGA Database ################################################################

data_mutations_tcga <- read.csv(paste(directory, "data/data_mutations_tcga.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_tcga <- data_mutations_tcga[,colSums(is.na(data_mutations_tcga))<nrow(data_mutations_tcga)]
data_mutations_tcga$Databse <- c("TCGA")

data_clinical_sample_tcga <- read.csv(paste(directory, "data/data_clinical_sample_tcga.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_clinical_sample_tcga <- data_clinical_sample_tcga[,colSums(is.na(data_clinical_sample_tcga))<nrow(data_clinical_sample_tcga)]
# Remove first 4 rows
data_clinical_sample_tcga <- tail(data_clinical_sample_tcga, -4)
colnames(data_clinical_sample_tcga)[2] ="Tumor_Sample_Barcode"

df_tcga <- merge(data_mutations_tcga, data_clinical_sample_tcga,by="Tumor_Sample_Barcode")

################################################################################

df_broad_colnames <- colnames(df_broad)
df_bcgsc_colnames <- colnames(df_bcgsc)
df_tcga_colnames <- colnames(df_tcga)

common_column_names <- Reduce(intersect, list(df_broad_colnames, df_bcgsc_colnames, df_tcga_colnames))

df_bcgsc <- subset(df_bcgsc, select = common_column_names)
df_broad <- subset(df_broad, select = common_column_names)
df_tcga <- subset(df_tcga, select = common_column_names)

data <- rbind(df_bcgsc, df_broad, df_tcga)

write.table(data, file = paste(directory, "results/data.csv", sep = "/"), quote=F, sep = ";", row.names=FALSE)

################################################################################

person_mutation_df <- data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Cancer.Type.Detailed")]
person_mutation_df_wide <- aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode + Cancer.Type.Detailed, person_mutation_df, function(x) paste(x, collapse=", "))

data_drugs <- read.csv(paste(directory, "data/data_drugs.tsv", sep = "/"), sep = ",")
data_drugs <- aggregate(mutation ~ drug, data_drugs, function(x) paste(x, collapse=", "))

patients_df <- data.frame('Patient ID' = "X", 'Mutation' = "Y", 'Drug' = 'Z')

for(i in 1:nrow(person_mutation_df_wide)){
  mutations <- strsplit(person_mutation_df_wide[i,3], ', ')[[1]]
  found_drug = F
  for(q in 1:nrow(data_drugs)){
    if(data_drugs[q, 2] %in% mutations) {
      patients_df <- rbind(patients_df, c(person_mutation_df_wide[i,1], data_drugs[q, 2], data_drugs[q, 1]))
      found_drug = T
    }
  }
  if(found_drug == F){
    patients_df <- rbind(patients_df, c(person_mutation_df_wide[i,1], "---", "---"))
  }
}
patients_df <- patients_df[-1,]

################################################################################
################################################################################
################################################################################

# Find Unique Cancer subtypes
print(unique(data$Cancer.Type.Detailed))

# Germinal Center B-Cell Type
data_germinal_center_type <- data[data$Cancer.Type.Detailed == "Germinal Center B-Cell Type", ]
data_germinal_center_type_f_table <- table(data_germinal_center_type$Hugo_Symbol)
data_germinal_center_type_f_table <- data_germinal_center_type_f_table[order(data_germinal_center_type_f_table, decreasing = TRUE)]

# Activated B-cell Type
data_activated_b_cell_type <- data[data$Cancer.Type.Detailed == "Activated B-cell Type", ]
data_activated_b_cell_type_f_table <- table(data_activated_b_cell_type$Hugo_Symbol)
data_activated_b_cell_type_f_table <- data_activated_b_cell_type_f_table[order(data_activated_b_cell_type_f_table, decreasing = TRUE)]

# Diffuse Large B-Cell Lymphoma, NOS
data_diffuse_large_b_cell_type <- data[data$Cancer.Type.Detailed == "Diffuse Large B-Cell Lymphoma, NOS", ]
data_diffuse_large_b_cell_type_f_table <- table(data_diffuse_large_b_cell_type$Hugo_Symbol)
data_diffuse_large_b_cell_type_f_table <- data_diffuse_large_b_cell_type_f_table[order(data_diffuse_large_b_cell_type_f_table, decreasing = TRUE)]


# >10
data_germinal_center_type_f_table[data_germinal_center_type_f_table>10]
data_activated_b_cell_type_f_table[data_activated_b_cell_type_f_table>10]
data_diffuse_large_b_cell_type_f_table[data_diffuse_large_b_cell_type_f_table>10]

table(data_germinal_center_type$HGVSc)[table(data_germinal_center_type$HGVSc)>1]

### Determine what mutations have people ### ###################################
cancer_type_df <- data_germinal_center_type

output_file_path <- paste(directory, "Results", paste(unique(cancer_type_df$Cancer.Type.Detailed), "_summary.txt", sep = ""), sep = "/")
file.create(output_file_path)
write(paste("##### Summary of ", unique(cancer_type_df$Cancer.Type.Detailed), " #####", sep = ""), output_file_path, append=TRUE)
write(paste("### Number of Samples:", length(unique(data_germinal_center_type$Tumor_Sample_Barcode)), sep = " "), output_file_path, append=TRUE)

# Keep only related columns
mutations_samples_cancer_type_df <- cancer_type_df[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSc")]
# List of samples
saples_list <- unique(mutations_samples_cancer_type_df$Tumor_Sample_Barcode)
# For each sample
sample_df_genes <- list()

for(sample_num in 1:length(saples_list)){
  sample <- saples_list[sample_num]
  sample_df <- mutations_samples_cancer_type_df[mutations_samples_cancer_type_df$Tumor_Sample_Barcode == sample,]
  sample_df_hugo_symbol_f_table <- table(sample_df$Hugo_Symbol)
  sample_df_genes[[length(sample_df_genes)+1]] <- unique(sample_df$Hugo_Symbol)
}
print(sample_df_genes)

for(sample_df_genes_num in 1:length(sample_df_genes)){
  print((sample_df_genes))
}

print(Reduce(intersect, sample_df_genes))

### Find most common mutations #################################################

# 99 patients
length(unique(data_diffuse_large_b_cell_type$Tumor_Sample_Barcode))

# 17 patients
length(unique(data_activated_b_cell_type$Tumor_Sample_Barcode))

# 31 patients
length(unique(data_germinal_center_type$Tumor_Sample_Barcode))

# Work with data_diffuse_large_b_cell_type data
cancer_type_df <- data_diffuse_large_b_cell_type

# Keep only related columns
cancer_type_df <- cancer_type_df[, c("Hugo_Symbol", "Variant_Classification")]

cancer_type_df <- cancer_type_df %>%
  add_count(Hugo_Symbol, Variant_Classification) %>% 
  distinct(Hugo_Symbol, Variant_Classification, .keep_all = TRUE) %>%
  mutate(frequency = n/sum(n))

cancer_type_df <- cancer_type_df[order(cancer_type_df$n, decreasing = T),]

output_file_path <- paste(directory, "Results", "data_diffuse_large_b_cell_type_mutations_frequency.txt", sep = "/")
file.create(output_file_path)
write(paste("##### Mutations of Diffuse Large B Cell Lymphoma ##### "), output_file_path, append=TRUE)
write(paste("### Number of Samples:", length(unique(data_diffuse_large_b_cell_type$Tumor_Sample_Barcode)), sep = " "), output_file_path, append=TRUE)
write.fwf(cancer_type_df, file=output_file_path, sep=" ", quote=F, rownames=F, colnames = T, na="NA", append=TRUE)



