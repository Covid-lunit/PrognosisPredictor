##################################
# Clinical Data Preprocessing - R#
# 1. Azure 서버 계정 활용
# ID: lunit@smcdevcloud.onmicrosoft.com
# Password: sosal OTP
# smc storage id: https://smc01storage.dfs.core.windows.net/lunit

# ssh: 20.194.32.110
# id: Lunit
# pw: sosal OTP
##################################

# library load
library(reshape2)
library(ggplot2)
library(glmnet)
library(randomForestSRC)
library(ggpubr)
library(pROC)
library(Epi)
library(readr)

# set current path
setwd("C:\\Users\\sosal\\Documents\\Lunit\\COVID\\")
        
# Read clinical datas
CRF <- read.csv("DATA/CRF_Final.tsv", sep="\t", header=T, stringsAsFactors=F)

# Set patient ID to 5-digit number
IDX <- which(nchar(CRF$ID) == 5)
CRF$ID[IDX] <- paste0("0", CRF$ID[IDX])

# Delete duplicates
FileName <- do.call(rbind, strsplit(CRF$Files, "/"))
CRF$Patient <- FileName[,2]
CRF <- CRF[, c(-1, -2)]
CRF <- CRF[!duplicated(CRF),]

# read CRF file & Validator result
Validation <- read.csv("DATA/Validation_All.csv", sep=",", stringsAsFactors=F)
Validation <- Validation[, -which(data.frame(colSums(is.na(Validation))) > 0)]

# preprocessing CXR filenames
Validation$Patient <- substr(Validation$InputFileName, 1, 6)
Validation$Order   <- substr(Validation$InputFileName, 12, 13)
Validation_First <- do.call(rbind, by(Validation, INDICES=Validation$Patient, function(x) x[which.min(as.numeric(x$Order)),]))

DATA <- merge(Validation_First, CRF, by="Patient")
DATA$Cohort <- substr(DATA$Patient, 1, 2)
colnames(DATA)[1] <- "ID"

Event <- DATA[, c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality")]
DATA <- DATA[, c("ID", "Cohort", "Event", "Age", "Sex", "HT", "DM", "CVD", "Cancer", "Fever",
                 "Cough", "Sputum", "Dyspnea", "Myalgia", "Sorethroat", 
                 "Lymphocyte", "Platelet",  "CRP", "LDH",
                 "Nodule","Consolidation", "Pneumothorax", "Pleural_effusion","Cardiomegaly",
                 "Fibrosis", "Pneumoperitoneum", "Mediastinal_widening", "Calcification", "Atelectasis")]

do.call(rbind, by(Event, INDICES=DATA$Cohort, function(x) c(nrow(x), colSums(x) )))

# Basic statistics
# Check statistics of logistic regression
Multivariates <- coef(summary(glm(Event ~ ., data=DATA[,2:14], family="binomial")))[-1,c(1,4)]
Univariates <- list()
for(i in 3:14) Univariates[[i]] <- coef(summary(glm(Event ~ ., data=DATA[,c(2, i)], family="binomial")))[2,c(1,4)]

LogitRegression <- data.frame(cbind(Multivariates, do.call(rbind, Univariates)))

get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}

# show heatmap of inter-correlation
melted_cormat <- melt(get_upper_tri(cor(DATA[,3:14])), na.rm=T)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

# Feature selection by lasso regularization
lassoModel <- cv.glmnet(x=as.matrix(DATA[,3:14]), y=DATA[,2], family="binomial", alpha=1)
coef(lassoModel, "lambda.min")

# RFSRC: Random Forest Survival, Regression, Classification
# Default model of RFSRC & Importance

# Three type of multi-omics data: Clinical, Validator (CXR), AllData
Clinical  <- data.frame(Label = factor(DATA$Event), DATA[, c(3,4,9,12)])
Validator <- data.frame(Label = factor(DATA$Event), DATA[, c(15:24)])
AllData   <- data.frame(Label = factor(DATA$Event), DATA[, c(3,4,5,6,9,12,14, 15:24)])

# Train each single-modal models using randomForestSRC
RF_clinical  <- rfsrc(Label ~ ., data=Clinical,  ntree=400, importance = TRUE)
RF_validator <- rfsrc(Label ~ ., data=Validator, ntree=400, importance = TRUE)
RF_allData   <- rfsrc(Label ~ ., data=AllData,   ntree=400, importance = TRUE)

# Set 5-folds, for cross validation
set.seed(0)
CV_IDX <- sample(rep(1:5, 200)[1:length(unique(DATA$ID))])
CV_IDX_DF <- data.frame(ID=unique(DATA$ID), CV_IDX = CV_IDX)
CV_IDX <- CV_IDX_DF[match(DATA$ID, CV_IDX_DF$ID), "CV_IDX"]

# CXR & Event correlation map
DATA[,14:23]
DATA[,24:29]

corMat <- c()
for(i in 14:23){
    corVec <- c()
    for(j in 24:29) corVec[ length(corVec) + 1] <- cor(DATA[,i], DATA[,j])
    corMat <- rbind(corMat, corVec)
}

rownames(corMat) <- colnames(DATA[,14:23])
colnames(corMat) <- colnames(DATA[,24:29])


#############################
# 5 FOLD - CROSS VALIDATION #
#############################

#Training
RF_clinical_CV <- RF_validator_CV <- RF_allData_CV <- list()
for(CV in 1:5){
    RF_clinical_CV[[CV]]  <- rfsrc(Label ~ ., data=Clinical [CV_IDX != CV,], ntree=400, importance = TRUE)
    RF_validator_CV[[CV]] <- rfsrc(Label ~ ., data=Validator[CV_IDX != CV,], ntree=400, importance = TRUE)
    RF_allData_CV[[CV]]   <- rfsrc(Label ~ ., data=AllData  [CV_IDX != CV,], ntree=400, importance = TRUE)
}

#Validation
RF_res <- list()
for(CV in 1:5) RF_res[[CV]] <- list(
    predict(RF_clinical_CV [[CV]], newdata=Clinical [CV_IDX == CV,])$predicted,
    predict(RF_validator_CV[[CV]], newdata=Validator[CV_IDX == CV,])$predicted,
    predict(RF_allData_CV  [[CV]], newdata=AllData  [CV_IDX == CV,])$predicted
)

# Get AUROC of all folds according to the data type.
# (Clinical$Label variable does not change, because all of the single-modal labels are same.)
rocTable <- c()
for(CV in 1:5) rocTable <- rbind(rocTable,
    c(roc(Clinical$Label[CV_IDX == CV] ~ RF_res[[CV]][[1]][,2])$auc,
      roc(Clinical$Label[CV_IDX == CV] ~ RF_res[[CV]][[2]][,2])$auc,
      roc(Clinical$Label[CV_IDX == CV] ~ RF_res[[CV]][[3]][,2])$auc))

# Get Feature Importance using RFSRC
FeatureImportance <- list(
    apply(do.call(cbind,lapply(RF_clinical_CV , function(x) x$importance[,1])), 1, median),
    apply(do.call(cbind,lapply(RF_validator_CV, function(x) x$importance[,1])), 1, median),
    apply(do.call(cbind,lapply(RF_allData_CV  , function(x) x$importance[,1])), 1, median))
FeatureImportanceDF <- lapply(FeatureImportance, function(x) data.frame(Feature=as.character(names(x)), Importance=x))

# plot for importance graph
ggplot(FeatureImportanceDF[[1]], aes(x=reorder(Feature, Importance), y=Importance)) + geom_bar(stat='identity') + coord_flip() + theme_bw() + xlab("") + ylab("Relative Feature Importance")
ggplot(FeatureImportanceDF[[2]], aes(x=reorder(Feature, Importance), y=Importance)) + geom_bar(stat='identity') + coord_flip() + theme_bw() + xlab("") + ylab("Relative Feature Importance")
ggplot(FeatureImportanceDF[[3]], aes(x=reorder(Feature, Importance), y=Importance)) + geom_bar(stat='identity') + coord_flip() + theme_bw() + xlab("") + ylab("Relative Feature Importance")

# Show ROC Curves, and the Youden's index
Epi::ROC(RF_res[[CV]][[1]][,2], as.numeric(Clinical$Label[CV_IDX == CV])-1)
Epi::ROC(RF_res[[CV]][[2]][,2], as.numeric(Validator$Label[CV_IDX == CV])-1)
Epi::ROC(RF_res[[CV]][[3]][,2], as.numeric(AllData$Label[CV_IDX == CV])-1)


#############################
# USE First X-ray data #
#############################
setwd("C:\\Users\\sosal\\Documents\\Lunit\\COVID\\인수인계\\")

# data preprocessing for clinical data
CRF <- read.csv("DATA/CRF_Final.tsv", sep="\t", header=T, stringsAsFactors=F)
IDX <- which(nchar(CRF$ID) == 5)
CRF$ID[IDX] <- paste0("0", CRF$ID[IDX])

FileName <- do.call(rbind, strsplit(CRF$Files, "/"))
CRF$Patient <- FileName[,2]
CRF <- CRF[, c(-1, -2)]
CRF <- CRF[!duplicated(CRF),]

# read CRF file & Validator result
Validation <- read.csv("DATA/Validation_All.csv", sep=",", stringsAsFactors=F)
Validation <- Validation[, -which(data.frame(colSums(is.na(Validation))) > 0)]

Validation$Patient <- substr(Validation$InputFileName, 1, 6)
Validation$Order   <- substr(Validation$InputFileName, 12, 13)
Validation_First <- do.call(rbind, by(Validation, INDICES=Validation$Patient, function(x) x[which.min(as.numeric(x$Order)),]))

DATA <- merge(Validation_First, CRF, by="Patient")
DATA$Cohort <- substr(DATA$Patient, 1, 2)
colnames(DATA)[1] <- "ID"

# Set columns
Event <- DATA[, c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality")]
CRF <- DATA[, c("ID", "Cohort", "Event", "Age", "Sex", "HT", "DM", "CVD", "Cancer", "Fever",
                 "Cough", "Sputum", "Dyspnea", "Myalgia", "Sorethroat", 
                 "Lymphocyte", "Platelet",  "CRP", "LDH",
                 "Nodule","Consolidation", "Pneumothorax", "Pleural_effusion","Cardiomegaly",
                 "Fibrosis", "Pneumoperitoneum", "Mediastinal_widening", "Calcification", "Atelectasis")]

Event[is.na(Event)] <- 0
Event_vector <- !rowSums(Event) == 0

# Final Data frame
CRF <- data.frame(
    ID = as.character(CRF$ID),
    Cohort = substr(CRF$ID, 1, 2),
    Age = as.numeric(CRF$Age),
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
    Sorethroat = as.numeric(CRF$Sorethroat),
    Lymphocyte = as.numeric(CRF$Lymphocyte),
    Platelet = as.numeric(CRF$Platelet),
    CRP = as.numeric(CRF$CRP),
    LDH = as.numeric(CRF$LDH),
    Nodule=as.numeric(CRF$Nodule),
    Consolidation=as.numeric(CRF$Consolidation),
    Pneumothorax=as.numeric(CRF$Pneumothorax),
    Pleural_effusion=as.numeric(CRF$Pleural_effusion),
    Cardiomegaly=as.numeric(CRF$Cardiomegaly),
    Fibrosis=as.numeric(CRF$Fibrosis),
    Pneumoperitoneum=as.numeric(CRF$Pneumoperitoneum),
    Mediastinal_widening=as.numeric(CRF$Mediastinal_widening),
    Calcification=as.numeric(CRF$Calcification),
    Atelectasis=as.numeric(CRF$Atelectasis),
    Event,
    Event = Event_vector
)

# NA data for Some variables are not actual NA, but False event (No events occur)
CRF[, c("HT", "DM", "CVD", "Cancer", "Fever", "Cough", "Sputum", "Dyspnea", "Myalgia", "Sorethroat")][
    is.na(CRF[, c("HT", "DM", "CVD", "Cancer", "Fever", "Cough", "Sputum", "Dyspnea", "Myalgia", "Sorethroat")])] <- 0

# Set Event Vector, and check Distribution for Multi-center validation
EventVector <- c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality", "Event")
N_distribution = do.call(rbind, by(CRF, INDICES=DATA$Cohort, function(x) c(nrow(x), colSums(x[, EventVector]))))
N_distribution[order(N_distribution[,1], decreasing=T),]

# Set data according to the Events
DATAList <- list()
for(Event_idx in 1:6){
    DATAList[[ Event_idx ]] <- CRF[, c(colnames(CRF)[1:28], EventVector[Event_idx])]
    colnames(DATAList[[ Event_idx ]])[29] <- "Event"
}

# Basic statistics
LogitRegression <- list()
for(Event_idx in 1:6){
    Multivariates <- coef(summary(glm(Event ~ ., data=DATAList[[ Event_idx ]][,c(2:13, 24)], family="binomial")))[-1,c(1,4)]
    Univariates <- list()
    for(i in 2:13) Univariates[[i]] <- coef(summary(glm(Event ~ ., data=DATAList[[ Event_idx ]][,c(i, 24)], family="binomial")))[2,c(1,4)]
    LogitRegression[[Event_idx]] <- data.frame(cbind(Multivariates, do.call(rbind, Univariates)))
}

# Make 4-types of data -> Clinical only, Laboratory data only, Validator (CXR) data only, and fianlly multi-modal data
Clinical <- Labdata <- Validator <- AllData <- list()
for(Event_idx in 1:6){
    Clinical [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(3:14)])
    Labdata  [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(15:18)])
    Validator[[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(19:28)])
    AllData  [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(3:28)])
}

# Set NA to the False event (This is not actual NA)
for(Event_idx in 1:6){
    Clinical[[Event_idx]][is.na(Clinical[[Event_idx]])] <- 0
    AllData[[Event_idx]][,1:13][ is.na(AllData[[Event_idx]][,1:13]) ] <- 0
}

# Correlation analysis
cors <- list(
    cor(Clinical [[1]][which(rowSums(is.na(Labdata  [[1]])) == 0),-1], Labdata  [[1]][which(rowSums(is.na(Labdata  [[1]])) == 0),-1]),
    cor(Clinical [[1]][which(rowSums(is.na(Validator[[1]])) == 0),-1], Validator[[1]][which(rowSums(is.na(Validator[[1]])) == 0),-1]),
    cor(Validator[[1]][which(rowSums(is.na(Labdata[[1]])) == 0 & rowSums(is.na(Validator[[1]])) == 0),-1], Labdata[[1]][which(rowSums(is.na(Labdata[[1]])) == 0 & rowSums(is.na(Validator[[1]])) == 0),-1]))

# Visualization of correlation matrix according to data types (Clinical, Validatir (CXR), Laboratory data)
melted_cormat <- lapply(cors, function(x) melt(x, na.rm=T))
ggplot(data = melted_cormat[[1]], aes(Var2, Var1, fill = value)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed() + xlab("LabData") + ylab("Clinical Data")

ggplot(data = melted_cormat[[2]], aes(Var2, Var1, fill = value)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

ggplot(data = melted_cormat[[3]], aes(Var2, Var1, fill = value)) + geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

# Visualization and wilcoxon test according to laboratory data
tmp_ggplot = cbind(Event, AllData[[6]])
colnames(tmp_ggplot)[6] <- "Event"
tmp = melt(tmp_ggplot[, c(1:6, 21)],id.vars=c("CRP"))
ggplot(tmp, aes(x=value, y=CRP, fill=value)) + geom_boxplot() + stat_compare_means(comparison=list(c("TRUE", "FALSE")), method="wilcox.test") + facet_wrap(~variable, nrow=1) 
tmp = melt(tmp_ggplot[, c(1:6, 22)],id.vars=c("LDH"))
ggplot(tmp, aes(x=value, y=LDH, fill=value)) + geom_boxplot() + stat_compare_means(comparison=list(c("TRUE", "FALSE")), method="wilcox.test") + facet_wrap(~variable, nrow=1) 


#########################
# clinical data summary #
DATA_CRF <- melt(DATA[,c( 4:14, 25:30)], id.vars=c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality", "Event"))
DATA_CXR <- melt(DATA[,c(15:30)],        id.vars=c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality", "Event"))
for(i in 1:5) DATA_CRF[,i] <- factor(DATA_CRF[,i])
for(i in 1:5) DATA_CXR[,i] <- factor(DATA_CXR[,i])
ggplot(DATA_CRF, aes(x=variable, y=value, fill=O2supply)) + geom_boxplot()

# Visualization for the association between events and input variables
ggplot(DATA_CXR, aes(x=O2supply,              y=value, fill=O2supply             )) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")
ggplot(DATA_CXR, aes(x=Mechnical_Ventilation, y=value, fill=Mechnical_Ventilation)) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")
ggplot(DATA_CXR, aes(x=ECMO,                  y=value, fill=ECMO                 )) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")
ggplot(DATA_CXR, aes(x=ICU_admission,         y=value, fill=ICU_admission        )) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")
ggplot(DATA_CXR, aes(x=Mortality,             y=value, fill=Mortality            )) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")
ggplot(DATA_CXR, aes(x=Event,                 y=value, fill=Event                )) + geom_boxplot() + stat_compare_means(comparisons = list(c("0","1")), method = "wilcox.test") + facet_wrap(~variable, nrow=1) + theme(legend.position = "none")

DATA_CXR_melt <- melt(DATA_CXR, id.vars=c("variable", "value"))
colnames(DATA_CXR_melt) <- c("CXR", "featureValue", "variable", "Event")
DATA_CXR_melt$Event <- factor(DATA_CXR_melt$Event %in% c("1", "TRUE"))

DATA_CXR_melt_tmp <- DATA_CXR_melt[DATA_CXR_melt$CXR %in% c("Nodule", "Consolidation"),]
ggplot(DATA_CXR_melt_tmp, aes(x=Event, y=featureValue, fill=Event)) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=.2, size=.1) + facet_wrap(~CXR + variable, ncol=6) + 
    stat_compare_means(comparisons = list(c("FALSE","TRUE")), method = "wilcox.test", label.y=105) + ylim(c(0,110))


cors <- matrix(0, ncol=6, nrow=12)
for(eventIDX in 1:6) for(idx in 3:14) cors[idx-2, eventIDX] <- cor.test(as.numeric(DATA[,idx]), as.numeric(DATA[,(25:30)[eventIDX]]))$estimate
rownames(cors) <- colnames(DATA)[3:14]
colnames(cors) <- colnames(DATA)[25:30]

p = ggplot(DATA_CXR_melt, aes(x=Event, y=featureValue, fill=Event)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=.2, size=.1) + facet_wrap(~variable + CXR,ncol=10) + stat_compare_means(comparisons = list(c("FALSE","TRUE")), method = "wilcox.test") + theme(legend.position = "none")
png("Boxplot_CXR.png", res=150, height=2500, width=2500)
print(p)
dev.off()


##########################################################
# RandomForestSRC based MultiCentre Split Validation for All Evenets #
##########################################################

# Make 4-types of data -> Clinical only, Laboratory data only, Validator (CXR) data only, and fianlly multi-modal data
Clinical <- Labdata <- Validator <- AllData <- list()
for(Event_idx in 1:6){
    Clinical [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(3:14)])
    Labdata  [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(15:18)])
    Validator[[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(19:28)])
    AllData  [[Event_idx]] <- data.frame(Label = factor(DATAList[[ Event_idx ]]$Event), DATAList[[ Event_idx ]][, c(3:28)])
}

# Geometical split validation (External validation)
# Training RFSRC model
SPLIT = !DATA$Cohort %in% c("01", "04", "06", "10", "11")
RF_clinical <- RF_labdata <- RF_validator <- RF_allData <- list()
for(Event_idx in 1:6){
    RF_clinical [[Event_idx]] <- rfsrc(Label ~ ., data=Clinical [[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
    RF_labdata  [[Event_idx]] <- rfsrc(Label ~ ., data=Labdata[[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
    RF_validator[[Event_idx]] <- rfsrc(Label ~ ., data=Validator[[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
    RF_allData  [[Event_idx]] <- rfsrc(Label ~ ., data=AllData  [[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
}

# Performance estimation
RF_res <- list()
for(Event_idx in 1:6) RF_res[[Event_idx]] <- list(
    predict(RF_clinical [[Event_idx]], newdata=Clinical [[Event_idx]][which(!SPLIT),])$predicted,
    predict(RF_labdata  [[Event_idx]], newdata=Labdata  [[Event_idx]][which(!SPLIT),])$predicted,
    predict(RF_validator[[Event_idx]], newdata=Validator[[Event_idx]][which(!SPLIT),])$predicted,
    predict(RF_allData  [[Event_idx]], newdata=AllData  [[Event_idx]][which(!SPLIT),])$predicted
)

NA_IDX <- list(which(rowSums(is.na(Clinical [[Event_idx]][which(!SPLIT),])) == 0),
               which(rowSums(is.na(Labdata  [[Event_idx]][which(!SPLIT),])) == 0),
               which(rowSums(is.na(Validator[[Event_idx]][which(!SPLIT),])) == 0),
               which(rowSums(is.na(AllData  [[Event_idx]][which(!SPLIT),])) == 0))

rocTable <- list()
for(Event_idx in 1:6) rocTable[[Event_idx]] <- c(
    roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2])$auc,
    roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2])$auc,
    roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2])$auc,
    roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2])$auc)
rocTable

# get AUROC and Confidence interval by Delong methods
rocCITable <- list()
for(Event_idx in 1:6) rocCITable[[Event_idx]] <- c(
    paste0(round(ci(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2])), 3), collapse="-"),
    paste0(round(ci(roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2])), 3), collapse="-"),
    paste0(round(ci(roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2])), 3), collapse="-"),
    paste0(round(ci(roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2])), 3), collapse="-"))
rocCITable


# Get AUROC datas as table
rocTables <- list()
for(currentCohort in c("01", "04", "06", "10", "11")){
    for(Event_idx in 1:6){
        currocTable[[Event_idx]] <- rep(NA, 4)
        try(currocTable[[Event_idx]][1] <- roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]][DATA$Cohort[which(!SPLIT)[NA_IDX[[1]]]] == currentCohort] ~ RF_res[[Event_idx]][[1]][DATA$Cohort[which(!SPLIT)[NA_IDX[[1]]]] == currentCohort,2])$auc)
        try(currocTable[[Event_idx]][2] <- roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]][DATA$Cohort[which(!SPLIT)[NA_IDX[[2]]]] == currentCohort] ~ RF_res[[Event_idx]][[2]][DATA$Cohort[which(!SPLIT)[NA_IDX[[2]]]] == currentCohort,2])$auc)
        try(currocTable[[Event_idx]][3] <- roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]][DATA$Cohort[which(!SPLIT)[NA_IDX[[3]]]] == currentCohort] ~ RF_res[[Event_idx]][[3]][DATA$Cohort[which(!SPLIT)[NA_IDX[[3]]]] == currentCohort,2])$auc)
        try(currocTable[[Event_idx]][4] <- roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]][DATA$Cohort[which(!SPLIT)[NA_IDX[[4]]]] == currentCohort] ~ RF_res[[Event_idx]][[4]][DATA$Cohort[which(!SPLIT)[NA_IDX[[4]]]] == currentCohort,2])$auc)
    }
    rocTables[[ length(rocTables) + 1]] <- currocTable
}

# Final AUROC table
lapply(rocTables, function(x) do.call(rbind, x))

# Statistical test of each single modal approachses
Event_idx=6
roc.test(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2]),
         roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2]))
roc.test(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2]),
         roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2]))
roc.test(roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2]),
         roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2]))

roc.test(roc(AllData [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2]),
         roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2]))
roc.test(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2]),
         roc(AllData[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2]))
roc.test(roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2]),
         roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2]))


# Plot ROC Curves according to the data type
RF_allData
Event_idx = 6
plot(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[1]]] ~ RF_res[[Event_idx]][[1]][,2]))
plot(roc(Labdata  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[2]]] ~ RF_res[[Event_idx]][[2]][,2]))
plot(roc(Validator[[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[3]]] ~ RF_res[[Event_idx]][[3]][,2]))
plot(roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][NA_IDX[[4]]] ~ RF_res[[Event_idx]][[4]][,2]))


# Get relative importance of each variables, according to the Model
ImpDFs <- list()
for(Event_idx in 1:6){
    Imp = RF_allData[[Event_idx]]$importance
    ImpDF = data.frame( Features = rownames(Imp), Importance = Imp[,1], Modality = factor(c(rep("Clinical", 12), rep("Laboratory", 4), rep("CXR", 10))))
    ImpDFs [[ Event_idx ]] <- ImpDF[ order(ImpDF[,2], decreasing=T), ][1:6,]
}
ImpDF <- do.call(rbind, ImpDFs)
ImpDF$Model <- as.integer((0:35) / 6)
ggplot(ImpDF, aes(x=reorder(Features,Importance), y=Importance, fill=Modality)) + geom_bar(stat='identity') + facet_wrap(~Model, ncol=2, scale="free") + coord_flip() + theme_bw() + xlab("Features") + ylab("Relative feature importance")

#save.image(file = "2021_10_14.RData")
#setwd("C:\\Users\\sosal\\Documents\\Lunit\\COVID\\인수인계\\")
#load(file = "2021_10_14.RData")

# ROC test by delong's method.
ROC_comparison <- list()
for(Event_idx in 1:6){
    ROC_comparison[[Event_idx]] <- list(c(), c(), c())
    for(CV in 1:5){
        ROC_comparison[[Event_idx]][[1]] <- rbind(ROC_comparison[[Event_idx]][[1]], cbind(Clinical [[Event_idx]]$Label[CV_IDX == CV], RF_res[[CV]][[Event_idx]][[1]][,2]))
        ROC_comparison[[Event_idx]][[2]] <- rbind(ROC_comparison[[Event_idx]][[2]], cbind(Clinical [[Event_idx]]$Label[CV_IDX == CV], RF_res[[CV]][[Event_idx]][[2]][,2]))
        ROC_comparison[[Event_idx]][[3]] <- rbind(ROC_comparison[[Event_idx]][[3]], cbind(Clinical [[Event_idx]]$Label[CV_IDX == CV], RF_res[[CV]][[Event_idx]][[3]][,2]))
        ROC_comparison[[Event_idx]][[3]] <- rbind(ROC_comparison[[Event_idx]][[4]], cbind(Clinical [[Event_idx]]$Label[CV_IDX == CV], RF_res[[CV]][[Event_idx]][[4]][,2]))
    }
}

ROC_DF <- list()
for(Event_idx in 1:6) ROC_DF[[Event_idx]] <- data.frame(
    Label = as.numeric(Clinical[[Event_idx]][which(!SPLIT),"Label"])-1,
    CRF   = RF_res[[Event_idx]][[1]][,2],
    CXR   = RF_res[[Event_idx]][[2]][,2],
    ALL   = RF_res[[Event_idx]][[3]][,2])

ROC_DF_melt_List <- lapply(ROC_DF, function(x) melt(x, id.vars="Label"))
for(Event_idx in 1:6) colnames(ROC_DF_melt_List[[Event_idx]]) <- c("Label", "Type", "Pred")
p = lapply(ROC_DF_melt_List, function(x) ggplot(x, aes(d = Label, m = Pred, color = Type)) + geom_roc() + style_roc())


CohortrocTable <- list()
for(CohortLevel_IDX in 1:5){
    CohortLevel=which(DATA$Cohort[!SPLIT]==unique(DATA$Cohort[!SPLIT])[CohortLevel_IDX])
    rocTable <- list()
    for(Event_idx in 1:6){ 
        a1 <- a2 <- a3 <- NA
        try({a1 = roc(Clinical [[Event_idx]]$Label[which(!SPLIT)][CohortLevel] ~ RF_res[[Event_idx]][[1]][CohortLevel,2])$auc})
        try({a2 = roc(Validator[[Event_idx]]$Label[which(!SPLIT)][CohortLevel] ~ RF_res[[Event_idx]][[2]][CohortLevel,2])$auc})
        try({a3 = roc(AllData  [[Event_idx]]$Label[which(!SPLIT)][CohortLevel] ~ RF_res[[Event_idx]][[3]][CohortLevel,2])$auc})
        rocTable[[Event_idx]] <- c(a1, a2, a3)
    }
    CohortrocTable[[CohortLevel_IDX]] <- rocTable
}

rocTableMatrix <- do.call(rbind, rocTable)
colnames(rocTableMatrix) <- c("CRF", "CXR", "Combined")
rownames(rocTableMatrix) <- c(colnames(Event), "All")


rocList <- list(rocTableMatrix)
for(rocTable in lapply(CohortrocTable, function(x) do.call(rbind, x))){
    rocTableMatrix[] <- rocTable
    rocList[[ length(rocList) + 1]] <- rocTableMatrix
}


##################################
# Delopyment of Final model      #
##################################

SPLIT = DATA$Cohort %in% c("12", "13", "14", "04", "01")
RF_clinical <- RF_validator <- RF_allData <- list()
CRF_Variables  <- c("Label", "Age", "Fever", "Dyspnea")
CXR_Variables  <- c("Label", "Nodule", "Consolidation", "Pneumothorax", "Pleural_effusion", "Cardiomegaly", "Fibrosis", "Mediastinal_widening", "Atelectasis")
FinalVariables <- c("Label", "Age", "Fever", "Dyspnea", "Nodule", "Consolidation", "Pneumothorax", "Pleural_effusion", "Cardiomegaly", "Fibrosis", "Mediastinal_widening", "Atelectasis")
for(Event_idx in 1:6){
    RF_clinical [[Event_idx]] <- rfsrc(Label ~ ., data=Clinical [[Event_idx]][which(SPLIT),CRF_Variables], ntree=400, importance = TRUE)
    RF_validator[[Event_idx]] <- rfsrc(Label ~ ., data=Validator[[Event_idx]][which(SPLIT),CXR_Variables], ntree=400, importance = TRUE)
    RF_allData  [[Event_idx]] <- rfsrc(Label ~ ., data=AllData  [[Event_idx]][which(SPLIT),FinalVariables], ntree=400, importance = TRUE)
}
RF_res <- list()
for(Event_idx in 1:6) RF_res[[Event_idx]] <- list(
    predict(RF_clinical [[Event_idx]], newdata=Clinical [[Event_idx]][which(!SPLIT),])$predicted,
    predict(RF_validator[[Event_idx]], newdata=Validator[[Event_idx]][which(!SPLIT),])$predicted,
    predict(RF_allData  [[Event_idx]], newdata=AllData  [[Event_idx]][which(!SPLIT),])$predicted
)

# Performance estimation
rocTable <- list()
for(Event_idx in 1:6) rocTable[[Event_idx]] <- c(
    roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[1]][,2])$auc,
    roc(Validator[[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[2]][,2])$auc,
    roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2])$auc)

print(rocTable)


# Save Final model
write_rds(RF_allData, "MultiModalModel.rds", "xz", compression = 9L)
packageVersion("randomForestSRC")


# Statistical test of Models according to the data type.
roc.test(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[1]][,2]),
         roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2]))
roc.test(roc(Validator[[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[2]][,2]),
         roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2]))

RF_clinical [[Event_idx]] <- rfsrc(Label ~ ., data=Clinical [[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
RF_validator[[Event_idx]] <- rfsrc(Label ~ ., data=Validator[[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)
RF_allData  [[Event_idx]] <- rfsrc(Label ~ ., data=AllData  [[Event_idx]][which(SPLIT),], ntree=400, importance = TRUE)

ci(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[1]][,2]))
ci(roc(Validator[[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[2]][,2]))
ci(roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2]))

(roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[1]][,2]))
(roc(Validator[[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[2]][,2]))
(roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2]))


roc.test(
    roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[1]][,2]),
    roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2])
)

roc.test(
    roc(Clinical [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[2]][,2]),
    roc(AllData  [[Event_idx]]$Label[which(!SPLIT)] ~ RF_res[[Event_idx]][[3]][,2])
)


rocTable_MultiCenter <- list()
for(Centers in 1:4){
    rocTable_MultiCenter[[Centers]] <- list()
    AllDataIDX <- DATA$Cohort %in% c("02", "06", "09", "10")[Centers]
    ValidatIDX <- DATA$Cohort[which(!SPLIT)] %in% c("02", "06", "09", "10")[Centers]
    for(Event_idx in 1:6){
        rocTable_MultiCenter[[Centers]][[Event_idx]] <- c(NA,NA,NA)
        try(rocTable_MultiCenter[[Centers]][[Event_idx]][1] <- roc(Clinical [[Event_idx]]$Label[AllDataIDX] ~ RF_res[[Event_idx]][[1]][ValidatIDX,2])$auc)
        try(rocTable_MultiCenter[[Centers]][[Event_idx]][2] <- roc(Validator[[Event_idx]]$Label[AllDataIDX] ~ RF_res[[Event_idx]][[2]][ValidatIDX,2])$auc)
        try(rocTable_MultiCenter[[Centers]][[Event_idx]][3] <- roc(AllData  [[Event_idx]]$Label[AllDataIDX] ~ RF_res[[Event_idx]][[3]][ValidatIDX,2])$auc)
    }
}

DATA$Cohort[which(!SPLIT)]
DATA$Cohort %in% c("02", "06", "09", "10")[Centers]

ROC_comparison <- list()
for(CV in 1:5){
    ROC_comparison[[Event_idx]][[1]] <- rbind(ROC_comparison[[Event_idx]][[1]], cbind(Clinical [[Event_idx]]$Label[which(SPLIT)], RF_res[[CV]][[Event_idx]][[1]][,2]))
    ROC_comparison[[Event_idx]][[2]] <- rbind(ROC_comparison[[Event_idx]][[2]], cbind(Clinical [[Event_idx]]$Label[which(SPLIT)], RF_res[[CV]][[Event_idx]][[2]][,2]))
    ROC_comparison[[Event_idx]][[3]] <- rbind(ROC_comparison[[Event_idx]][[3]], cbind(Clinical [[Event_idx]]$Label[which(SPLIT)], RF_res[[CV]][[Event_idx]][[3]][,2]))
}

ROC_DF <- list()
for(Event_idx in 1:6) ROC_DF[[Event_idx]] <- data.frame(
    Label = as.numeric(Clinical[[Event_idx]][which(!SPLIT),"Label"])-1,
    CRF   = RF_res[[Event_idx]][[1]][,2],
    CXR   = RF_res[[Event_idx]][[2]][,2],
    ALL   = RF_res[[Event_idx]][[3]][,2])

ROC_DF_melt_List <- lapply(ROC_DF, function(x) melt(x, id.vars="Label"))
for(Event_idx in 1:6) colnames(ROC_DF_melt_List[[Event_idx]]) <- c("Label", "Type", "Pred")
p = lapply(ROC_DF_melt_List, function(x) ggplot(x, aes(d = Label, m = Pred, color = Type)) + geom_roc() + style_roc())

do.call(rbind, rocTable)
do.call("grid.arrange", c(p, ncol=3))

Events <- c("O2supply", "Mechnical_Ventilation", "ECMO", "ICU_admission", "Mortality")
Epi::ROC(DATA[which(!SPLIT), "Consolidation"], DATA[which(!SPLIT),Events[i]])


lapply(RF_allData, function(x){ 
    tmp = x$importance
    tmp[order(tmp[,3], decreasing=T),]
})


# Comparison our model vs Consolidation model
Consolidation_AUC <- c()
for(Event_idx in 1:6) Consolidation_AUC[Event_idx] <- roc(Validator[[Event_idx]]$Label[which(!SPLIT)] ~ Validator[[Event_idx]]$Consolidation[which(!SPLIT)])$auc



##################################################
# XGBoost based Cross Validation for All Evenets #
##################################################

# However, RFSRC is better than XGBoost

set.seed(8)
CV_IDX <- sample(rep(1:5, 200)[1:nrow(DATAList[[1]])])

CV_distribution_Table <- c()
for(CV in 1:5) CV_distribution_Table <- rbind(CV_distribution_Table, colSums(DATA[CV_IDX == CV,24:28]))
print(CV_distribution_Table)

params <- list(booster="gbtree",eta=0.001,max_depth=5,gamma=3,subsample=0.75,colsample_bytree=1,objective="binary:logistic",eval_metric="auc",num_class=2)
params <- list(booster="gbtree",objective="binary:logistic", eta=0.001, gamma=3, max_depth=10, min_child_weight=1, subsample=1, colsample_bytree=1, eval_metric="auc")

XG_clinical <- XG_validator <- XG_allData <- list()
for(CV in 1:5){
    XG_clinical[[CV]] <- XG_validator[[CV]] <- XG_allData[[CV]] <- list()
    for(Event_idx in 1:6){
        
        xgb.DMatrix(data=as.matrix(Clinical [[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(Clinical [[Event_idx]][CV_IDX != CV,1])-1)
        xgb.DMatrix(data=as.matrix(Validator[[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(Validator[[Event_idx]][CV_IDX != CV,1])-1)
        xgb.DMatrix(data=as.matrix(AllData  [[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(AllData  [[Event_idx]][CV_IDX != CV,1])-1)
        
        XG_clinical [[CV]][[Event_idx]] <- xgb.train(params=params,nrounds=10000,early_stopping_rounds=10,
            data=xgb.DMatrix(data=as.matrix(Clinical [[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(Clinical [[Event_idx]][CV_IDX != CV,1])-1),watchlist=list(
            val1=xgb.DMatrix(data=as.matrix(Clinical [[Event_idx]][CV_IDX == CV,-1]), label=as.numeric(Clinical [[Event_idx]][CV_IDX == CV,1])-1)
        ),verbose=1)
        XG_validator[[CV]][[Event_idx]] <- xgb.train(params=params,nrounds=10000,early_stopping_rounds=10,
            data=xgb.DMatrix(data=as.matrix(Validator[[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(Validator[[Event_idx]][CV_IDX != CV,1])-1),watchlist=list(
            val1=xgb.DMatrix(data=as.matrix(Validator[[Event_idx]][CV_IDX == CV,-1]), label=as.numeric(Validator[[Event_idx]][CV_IDX == CV,1])-1)
        ),verbose=1)
        XG_allData  [[CV]][[Event_idx]] <- xgb.train(params=params,nrounds=10000,early_stopping_rounds=10,
            data=xgb.DMatrix(data=as.matrix(AllData  [[Event_idx]][CV_IDX != CV,-1]), label=as.numeric(AllData  [[Event_idx]][CV_IDX != CV,1])-1),watchlist=list(
            val1=xgb.DMatrix(data=as.matrix(AllData  [[Event_idx]][CV_IDX == CV,-1]), label=as.numeric(AllData  [[Event_idx]][CV_IDX == CV,1])-1)
        ),verbose=1)

    }
}