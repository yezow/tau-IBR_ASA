
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
initial_para_est <- c(0,0,0.1,0,0.2,0,0,0,0,0,0,0,0,0,0.1,0,0,0,0,0,-0.5,0,0,0,2)

####################################
###MI
set.seed(1234)
tau=365
ID<-ex_data$ID
T_s<-ex_data$Days_To_Onset1
X<-ex_data$Days_To_Onset1
X[is.na(X)]<-ex_data$Days_In_Study[which(is.na(X))]
Delta<-ifelse(!is.na(T_s),1,0)

##trt*age
Z_mu <- dplyr::select(ex_data, trt,age_10_new, gender,fev1_10_new,nowsmk,trt_age,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44)
Z_pi <- dplyr::select(ex_data,trt,age_10_new,gender, fev1_10_new,nowsmk,trt_age,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44) #nowsmk,oxygen,ster1yr
formula_fit<-as.formula("pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+trt_age+
                   latitude_38+latitude_39+latitude_40+latitude_42+latitude_44")

Data_cens<-data.frame(ID,T_s,X,Delta,Z_mu,Z_pi)
Z_mu<-data.matrix(Z_mu)
Z_pi<-data.matrix(Z_pi)
nsubj<-nrow(Z_mu)
  
row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
source("Scripts/MI_function.R") #run MI imputation and get MI_result
coef_MI<-MI_result[[1]]
var_MI<-MI_result[[2]]
  
EM_result<-MI_impute_converge()
w <- EM_result[[3]]
coef_EM <-EM_result[[1]]
var_EM <- EM_var(w,coef_EM)
  
pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
pseudo_fit <- lm(formula_fit,data=Data_cens) 
  
p_trt_age<-rbind(2*(1-pnorm(abs(coef_MI[7]),mean=0,sd=sqrt(var_MI[7,7]))),
                   2*(1-pnorm(abs(coef_MI[19]),mean=0,sd=sqrt(var_MI[19,19]))),
                   2*(1-pnorm(abs(coef_EM[7]),mean=0,sd=sqrt(var_EM[7,7]))),
                   2*(1-pnorm(abs(coef_EM[19]),mean=0,sd=sqrt(var_EM[19,19]))),
                   as.numeric(summary(pseudo_fit)$coefficients[,4][7]))

##trt*gender
Z_mu <- dplyr::select(ex_data, trt,age_10_new, gender,fev1_10_new,nowsmk,trt_gender,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44)
Z_pi <- dplyr::select(ex_data,trt,age_10_new,gender, fev1_10_new,nowsmk,trt_gender,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44) #nowsmk,oxygen,ster1yr
formula_fit<-as.formula("pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+trt_gender+
                   latitude_38+latitude_39+latitude_40+latitude_42+latitude_44")

Data_cens<-data.frame(ID,T_s,X,Delta,Z_mu,Z_pi)
Z_mu<-data.matrix(Z_mu)
Z_pi<-data.matrix(Z_pi)
nsubj<-nrow(Z_mu)

row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
source("Scripts/MI_function.R") #run MI imputation and get MI_result
coef_MI<-MI_result[[1]]
var_MI<-MI_result[[2]]

EM_result<-MI_impute_converge()
w <- EM_result[[3]]
coef_EM <-EM_result[[1]]
var_EM <- EM_var(w,coef_EM)

pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
pseudo_fit <- lm(formula_fit,data=Data_cens) 

p_trt_gender<-rbind(2*(1-pnorm(abs(coef_MI[7]),mean=0,sd=sqrt(var_MI[7,7]))),
                 2*(1-pnorm(abs(coef_MI[19]),mean=0,sd=sqrt(var_MI[19,19]))),
                 2*(1-pnorm(abs(coef_EM[7]),mean=0,sd=sqrt(var_EM[7,7]))),
                 2*(1-pnorm(abs(coef_EM[19]),mean=0,sd=sqrt(var_EM[19,19]))),
                 as.numeric(summary(pseudo_fit)$coefficients[,4][7]))

##trt*gev1
Z_mu <- dplyr::select(ex_data, trt,age_10_new, gender,fev1_10_new,nowsmk,trt_fev1,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44)
Z_pi <- dplyr::select(ex_data,trt,age_10_new,gender, fev1_10_new,nowsmk,trt_fev1,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44) #nowsmk,oxygen,ster1yr
formula_fit<-as.formula("pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+trt_fev1+
                   latitude_38+latitude_39+latitude_40+latitude_42+latitude_44")

Data_cens<-data.frame(ID,T_s,X,Delta,Z_mu,Z_pi)
Z_mu<-data.matrix(Z_mu)
Z_pi<-data.matrix(Z_pi)
nsubj<-nrow(Z_mu)

row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
source("Scripts/MI_function.R") #run MI imputation and get MI_result
coef_MI<-MI_result[[1]]
var_MI<-MI_result[[2]]

EM_result<-MI_impute_converge()
w <- EM_result[[3]]
coef_EM <-EM_result[[1]]
var_EM <- EM_var(w,coef_EM)

pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
pseudo_fit <- lm(formula_fit,data=Data_cens) 

p_trt_fev1<-rbind(2*(1-pnorm(abs(coef_MI[7]),mean=0,sd=sqrt(var_MI[7,7]))),
                 2*(1-pnorm(abs(coef_MI[19]),mean=0,sd=sqrt(var_MI[19,19]))),
                 2*(1-pnorm(abs(coef_EM[7]),mean=0,sd=sqrt(var_EM[7,7]))),
                 2*(1-pnorm(abs(coef_EM[19]),mean=0,sd=sqrt(var_EM[19,19]))),
                 as.numeric(summary(pseudo_fit)$coefficients[,4][7]))

##trt*smk
Z_mu <- dplyr::select(ex_data, trt,age_10_new, gender,fev1_10_new,nowsmk,trt_smk,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44)
Z_pi <- dplyr::select(ex_data,trt,age_10_new,gender, fev1_10_new,nowsmk,trt_smk,latitude_38,latitude_39,latitude_40,latitude_42,latitude_44) #nowsmk,oxygen,ster1yr
formula_fit<-as.formula("pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+trt_smk+
                   latitude_38+latitude_39+latitude_40+latitude_42+latitude_44")

Data_cens<-data.frame(ID,T_s,X,Delta,Z_mu,Z_pi)
Z_mu<-data.matrix(Z_mu)
Z_pi<-data.matrix(Z_pi)
nsubj<-nrow(Z_mu)

row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
source("Scripts/MI_function.R") #run MI imputation and get MI_result
coef_MI<-MI_result[[1]]
var_MI<-MI_result[[2]]

EM_result<-MI_impute_converge()
w <- EM_result[[3]]
coef_EM <-EM_result[[1]]
var_EM <- EM_var(w,coef_EM)

pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
pseudo_fit <- lm(formula_fit,data=Data_cens) 

p_trt_smk<-rbind(2*(1-pnorm(abs(coef_MI[7]),mean=0,sd=sqrt(var_MI[7,7]))),
                 2*(1-pnorm(abs(coef_MI[19]),mean=0,sd=sqrt(var_MI[19,19]))),
                 2*(1-pnorm(abs(coef_EM[7]),mean=0,sd=sqrt(var_EM[7,7]))),
                 2*(1-pnorm(abs(coef_EM[19]),mean=0,sd=sqrt(var_EM[19,19]))),
                 as.numeric(summary(pseudo_fit)$coefficients[,4][7]))
table5_result<-rbind(p_trt_age,p_trt_gender,p_trt_fev1,p_trt_smk)
write.csv(table5_result,file = "Results/Table5_result.csv")