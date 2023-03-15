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

#write.fwf(data,file=paste(directory, "results/data.csv", sep = "/"),sep=";", quote=F, rownames=F, na="NA", append=F)