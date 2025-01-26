
registerDoParallel(cores=4)
import::from("numDeriv", "hessian")
source("Scripts/Paper1_func.R")
original_data<- read.csv("Data/Azith.csv")
ex_data<-original_data[-251,]
#Modify the study # The study is followed for 380 days(maxi), but 24 patients who have death recorded after the end of study
ex_data$Days_In_Study <- ifelse(is.na(ex_data$dayofdeath), ex_data$Days_In_Study, ex_data$dayofdeath)
ex_data$Days_In_Study <- ifelse(ex_data$Days_In_Study>380, 380, ex_data$Days_In_Study)
#other covariates
ex_data$trt <- as.numeric(ex_data$trtgroup == 1)
ex_data$latitude<-factor(ex_data$latitude)
ex_data$age<-ex_data$age/10
ex_data$fev1pp_00<-ex_data$fev1pp_00/10
ex_data$age_10_new<-ex_data$age-mean(ex_data$age)
ex_data$fev1_10_new<-ex_data$fev1pp_00-mean(ex_data$fev1pp_00)
ex_data$trt_smk <- ex_data$trt*ex_data$nowsmk
ex_data$trt_age <- ex_data$trt*ex_data$age_10_new
ex_data$trt_fev1 <- ex_data$trt*ex_data$fev1_10_new
ex_data$trt_gender <- ex_data$trt*ex_data$gender

ex_data <- fastDummies::dummy_cols(ex_data, select_columns = "latitude")

#####set initial model parameters for estimation
initial_para_est <- c(0,0,0.1,0,0.2,0,0,0,0,0,0,0,0,0.1,0,0,0,0,-0.5,0,0,0,2)

####################################
###MI
set.seed(1234)

ID<-ex_data$ID
T_s<-ex_data$Days_To_Onset1
X<-ex_data$Days_To_Onset1
X[is.na(X)]<-ex_data$Days_In_Study[which(is.na(X))]
Delta<-ifelse(!is.na(T_s),1,0)
Z_mu <- dplyr::select(ex_data, trt,age_10_new, gender,fev1_10_new,nowsmk,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44)
Z_pi <- dplyr::select(ex_data,trt,age_10_new,gender, fev1_10_new,nowsmk,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44) #nowsmk,oxygen,ster1yr
Data_cens<-data.frame(ID,T_s,X,Delta,Z_mu,Z_pi)
Z_mu<-data.matrix(Z_mu)
Z_pi<-data.matrix(Z_pi)
nsubj<-nrow(Z_mu)

row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)

#######calculate population mean
#rMST for averaged population
#average 
gender_mean<-mean(Z_mu[,3])
nowsmk_mean <- mean(Z_mu[,5])
center2_mean<-mean(Z_mu[,6])
center3_mean<-mean(Z_mu[,7])
center4_mean<-mean(Z_mu[,8])
center5_mean<-mean(Z_mu[,9])
center6_mean<-mean(Z_mu[,10])
