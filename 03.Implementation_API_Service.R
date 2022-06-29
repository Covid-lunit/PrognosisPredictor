setwd("C:\\Users\\sosal\\Documents\\Lunit\\Radiology\\COVID\\VALIDATION")
library(randomForestSRC)
library(readr)

#packageVersion("randomForestSRC")   # 2.9.3
#packageVersion("readr")             # 1.3.1

# Load Best model
Model <- read_rds('MultiModalModel.rds')

#* Echo back the input
#* @param Clinical variables: Age, Fever, Dyspnea,
#* @param DLAD-10 variables: Nodule, Consolidation, Pneumothorax, Pleural effusion, Cardiomegaly, Fibrosis, Pneumoperitoneum, MEdiastinal widening, Calcification, Atelectasis
#* @get /Prognosis
function(Age=70,Fev=1,Dys=1,Nod=9.07,Con=99.18,PneTho=0.81,Ple=1.61,Car=3.67,Fib=2.23,PnePer=0.23,Med=0.64,Cal=58.4,Ate=1.06){

    Data = data.frame(
        Age = as.numeric(Age),
        Fever = as.numeric(Fev),
        Dyspnea = as.numeric(Dys),
        Nodule= as.numeric(Nod),
        Consolidation= as.numeric(Con),
        Pneumothorax=as.numeric(PneTho),
        Pleural_effusion=as.numeric(Ple),
        Cardiomegaly=as.numeric(Car),
        Fibrosis=as.numeric(Fib),
        Pneumoperitoneum=as.numeric(PnePer),
        Mediastinal_widening=as.numeric(Med),
        Calcification=as.numeric(Cal),
        Atelectasis=as.numeric(Ate)
    )
    
    Result <- unlist(lapply(Model, function(Model) predict(Model, Data)$predicted[,2]))
    names(Result) <- c("O2Supply", "MV", "ECMO", "ICU", "Mortality", "All")
    return(Result)
}


################### Run API
# plumber (run by Rscript) port 8000
# setwd("C:\\Users\\sosal\\Documents\\Lunit\\Radiology\\COVID\\VALIDATION")
# library(plumber)
# pr("4.Prediction.R") %>% pr_run(port=8000)
# curl "http://localhost:8000/Prognosis?Age=70&Fev=1&Dys=1&Nod=9.07&Con=99.18&Pne=0.81&Ple=1.61&Car=3.67&Fib=2.23&Pne=0.23&Med=0.64&Cal=58.4&Ate=1.06"

#Age=70,Fev=1,Dys=1,Nod=9.07,Con=99.18,PneTho=0.81,Ple=1.61,Car=3.67,Fib=2.23,PnePer=0.23,Med=0.64,Cal=58.4,Ate=1.06

