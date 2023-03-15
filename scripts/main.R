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

# BCGSC Database

data_mutations_bcgsc <- read.csv(paste(directory, "data/data_mutations_bcgsc.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_bcgsc <- data_mutations_bcgsc[,colSums(is.na(data_mutations_bcgsc))<nrow(data_mutations_bcgsc)]
data_mutations_bcgsc$Databse <- c("BCGSC")

# Broad Database

data_mutations_broad <- read.csv(paste(directory, "data/data_mutations_broad.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_broad <- data_mutations_broad[,colSums(is.na(data_mutations_broad))<nrow(data_mutations_broad)]
data_mutations_broad$Databse <- c("Broad")

# TCGA Database

data_mutations_tcga <- read.csv(paste(directory, "data/data_mutations_tcga.txt", sep = "/"), sep = "\t")
# Drop columns with all NA
data_mutations_tcga <- data_mutations_tcga[,colSums(is.na(data_mutations_tcga))<nrow(data_mutations_tcga)]
data_mutations_tcga$Databse <- c("TCGA")

data_mutations_broad_colnames <- colnames(data_mutations_broad)
data_mutations_bcgsc_colnames <- colnames(data_mutations_bcgsc)
data_mutations_tcga_colnames <- colnames(data_mutations_tcga)

common_column_names <- Reduce(intersect, list(data_mutations_broad_colnames, data_mutations_bcgsc_colnames, data_mutations_tcga_colnames))

data_mutations_bcgsc <- subset(data_mutations_bcgsc, select = common_column_names)
data_mutations_broad <- subset(data_mutations_broad, select = common_column_names)
data_mutations_tcga <- subset(data_mutations_tcga, select = common_column_names)



data <- rbind(data_mutations_bcgsc, data_mutations_broad, data_mutations_tcga)

head(data_mutations_tcga)