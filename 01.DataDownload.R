# find . -maxdepth 4 | grep CXR | grep dcm > fileList.txt

##################################
# Clinical Data Preprocessing - R#
##################################
rm(list = ls())
gc()

library(openxlsx)
library(pROC)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(xgboost)
library(plotROC)
library(gridExtra)

# 1_CRF_210416, 2_CRF_210416 files that can only be accessed from within the server
CRF1 <- openxlsx::read.xlsx("/home/Lunit/DATA/1_CRF_210416.xlsx", sheet=1)
CRF2 <- openxlsx::read.xlsx("/home/Lunit/DATA/2_CRF_210416.xlsx", sheet=1)

# Since column names are very complex strings, access only with numbers is inevitable.
CRF1 <- CRF1[,c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,20,21,22,23,73,74,75,76,77,78)]
CRF2 <- CRF2[,c(1,2,3,8,9,10,11,12,14,15,16,17,18,19,33,34,35,36,132,133,134,135,136,137)]

Columns <- c("ID", "Age.Group", "Sex", "Smoking", "HT", "DM", "CVD", "Cancer", "Fever", "Cough",
              "Sputum", "Dyspnea", "Myalgia","Sorethroat", "Lymphocyte", "Platelet", "CRP", "LDH",
              "O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality", "DischargeDate")

colnames(CRF1) <- Columns
colnames(CRF2) <- Columns

# Delete the first, second rows which is the detail description for the column
CRF1 <- CRF1[-c(1,2),]
CRF2 <- CRF2[-c(1,2),]

# Combine Whole data
CRF <- data.frame(rbind(CRF1, CRF2))

# Check the number of Not Available data according to the columns
data.frame(N_na=colSums(is.na(CRF)))

# Some NA variables is not actual NA, but 0 (False event)
CRF[, c(-1,-17)][is.na(CRF[, c(-1, -17)])] <- "0"

CRF <- data.frame(CRF)

# Remove the row in which the Age.group variable is a missing value.
CRF <- CRF[-which(is.na(CRF$Age.Group)),]

#which(is.na(CRF$Sex))

# Set the labels
Event = data.matrix(CRF[,c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality")])

# Data curators said that NA is not actual NA, but False event.
Event[is.na(Event)] <- 0
Event_vector <- !rowSums(Event) == 0

# Final data frame structure, and the type of variables
CRF <- data.frame(
    ID = as.character(CRF$ID),
    Age = as.numeric(CRF$Age.Group),
    Sex = CRF$Sex,
    HT = as.numeric(CRF$HT),
    DM = as.numeric(CRF$DM),
    CVD = as.numeric(CRF$CVD),
    Cancer = as.numeric(CRF$Cancer),
    Fever = as.numeric(CRF$Fever),
    Cough = as.numeric(CRF$Cough),
    Sputum = as.numeric(CRF$Sputum),
    Dyspnea = as.numeric(CRF$Dyspnea),
    Myalgia = as.numeric(CRF$Myalgia),
    Lymphocyte <- as.numeric(CRF$Lymphocyte),
    Platelet <- as.numeric(CRF$Platelet),
    CRP <- as.numeric(CRF$CRP),
    LDH <- as.numeric(CRF$LDH),
    Sorethroat = as.numeric(CRF$Sorethroat),
    Event,
    Event = Event_vector
)

# file List: DICOM X-ray file names
Files <- read.csv("/DATA/lunit/fileList.txt", sep="\t", header=F, stringsAsFactors=F)[[1]]
Files <- gsub("\\./", "/DATA/lunit/", Files)
Files_patient <- do.call(rbind, strsplit(Files, "/"))[,4]
Files <- data.frame(Files, ID = Files_patient)

Files$ID <- as.numeric(as.character(Files$ID))
CRF$ID <- as.numeric(as.character(CRF$ID))

# Merge Clinical information and X-ray file names by Patient ID
DF <- merge(Files, CRF, by="ID")
write.table(DF, file="~/CRF_Event.tsv", quote=F, row.names=F, sep="\t")


# The way to access server, patient data, to download the x-ray files

#sudo mount /dev/sdc /DATA/
#sudo azcopy copy 'https://smc01storage.blob.core.windows.net/lunit/?sv=2020-02-10&ss=b&srt=co&sp=rlx&se=2021-12-31T14:59:59Z&st=2021-03-18T06:58:36Z&spr=https&sig=k4Ibpbtyl9%2BAdzGehXknR93zGUqYoT3Vx3bUL02ccfw%3D' './' --recursive=true
#sudo docker run -it --rm --ipc=host -v /DATA:/DATA lunit/base_cxr:v13_base_pytorch_1.5_cuda10.1_cudnn7 /bin/bash