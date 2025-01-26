
source("Scripts/Paper1_func.R")

# ########################################################
# Get simulated data for figure1 scatter plots in paper1)
# #######################################################
#pi model parameters
beta_0_pi<- -1
beta_1_pi<- c(1,2,-1.5) 
#mu model parameters
alpha_0<- -2
alpha_1<-c(1.2,2)

#one para in beta dist
nu <-3
#paraset
para_set<-c(alpha_0,alpha_1,beta_0_pi,beta_1_pi,nu)
#initial parameter sets used in estimation
initial_para_est <- c(-1,1,1,-1,1,1,0,2)

tau=30

Z1<-runif(nsubj,0,1)
Z2<-rbinom(n=nsubj, size=1, prob=0.7)
Z3<-runif(nsubj,0,1)
Z_mu<-cbind(Z1,Z2)
Z_pi<-cbind(Z1,Z2,Z3)
mu <- inv.logit(alpha_0+Z_mu%*%alpha_1)
logit_mu <- logit(mu)
alpha_b <- mu*nu
beta_b <-(1-mu)*nu
logit_pi <- beta_0_pi+Z_pi%*%beta_1_pi
pi <- inv.logit(logit_pi)
mu<-as.vector(mu)
pi<-as.vector(pi)
Mean_tau_T <- mu*tau*(1-pi)+tau*pi
Y<-Generate_T()

##########################
#for censored data
#EM
Data_cens <- generate_cens_dataset(Y)
row_MI <- which(Data_cens$X<tau & Data_cens$Delta==0)
EM_result<-MI_impute_converge()
w <- EM_result[[3]]
coef_EM <-EM_result[[1]]
var_EM <- EM_var(w,coef_EM)
se_EM <- sqrt(diag(var_EM))

#calcute estimated restricted mean
mu_est_EM <- inv.logit(EM_result[[1]][1]+Z_mu%*%EM_result[[1]][2:(ncol(Z_mu)+1)])
pi_est_EM <- inv.logit(EM_result[[1]][(ncol(Z_mu)+2)]+Z_pi%*%EM_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
Mean_tau_T_est_EM <- mu_est_EM*tau*(1-pi_est_EM)+tau*pi_est_EM

#calculate variance of restricted mean and CI (MI)
var_restricted_mean_EM_est<-Var_restricted_mean_func(var_EM,EM_result[[1]])[[1]]
restricted_mean_EM_lower <- Mean_tau_T_est_EM-1.964739*sqrt(var_restricted_mean_EM_est)
restricted_mean_EM_upper <- Mean_tau_T_est_EM+1.964739*sqrt(var_restricted_mean_EM_est)
restricted_mean_CI_EM <- 0
which_out_CI <- numeric()
for(i in 1:nrow(Z_mu)){
  if(Mean_tau_T[i]>restricted_mean_EM_lower[i] & Mean_tau_T[i]<restricted_mean_EM_upper[i]){
    restricted_mean_CI_EM <- restricted_mean_CI_EM+1
  } else{
    which_out_CI <- append(which_out_CI,i)
  }
}
restricted_mean_CI_EM_empirical_indep <- restricted_mean_CI_EM/nsubj

#######################################
#MI
source("Scripts/MI_function.R") #run MI imputation and get MI_result
coef_MI <-MI_result[[1]]
var_MI<-MI_result[[2]]
se_MI<- sqrt(diag(MI_result[[2]]))

#calcute estimated restricted mean
mu_est_MI <- inv.logit(MI_result[[1]][1]+Z_mu%*%MI_result[[1]][2:(ncol(Z_mu)+1)])
pi_est_MI <- inv.logit(MI_result[[1]][(ncol(Z_mu)+2)]+Z_pi%*%MI_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI

#calculate variance of restricted mean and CI (MI)
var_restricted_mean_MI_est<-Var_restricted_mean_func(var_MI,MI_result[[1]])[[1]]
restricted_mean_MI_lower <- Mean_tau_T_est_MI-1.964739*sqrt(var_restricted_mean_MI_est)
restricted_mean_MI_upper <- Mean_tau_T_est_MI+1.964739*sqrt(var_restricted_mean_MI_est)
restricted_mean_CI_MI <- 0
which_out_CI_MI <- numeric()
for(i in 1:nrow(Z_mu)){
  if(Mean_tau_T[i]>restricted_mean_MI_lower[i] & Mean_tau_T[i]<restricted_mean_MI_upper[i]){
    restricted_mean_CI_MI <- restricted_mean_CI_MI+1
  } else{
    which_out_CI_MI <- append(which_out_CI_MI,i)
  }
}
restricted_mean_CI_MI_empirical_indep <- restricted_mean_CI_MI/nsubj

######################################
#restricted mean model for censored data
pseudo <- pseudomean(Data_cens$X,Data_cens$Delta)
pseudo_fit <- lm(pseudo ~ Z_pi )
pseudo_RMST_coef<-pseudo_fit$coefficients
Pseudo_RMST_mean_tau_t_est <- pseudo_RMST_coef[1]+Z_pi%*%pseudo_RMST_coef[2:(ncol(Z_pi)+1)]
Pseudo_RMST_diff <- sum(Pseudo_RMST_mean_tau_t_est - Mean_tau_T)/length(Mean_tau_T)
Pseudo_RMST_MSE <- sum((Pseudo_RMST_mean_tau_t_est - Mean_tau_T)^2)/length(Mean_tau_T)

#estimates of the variance of restricted mean in RMST (when censoring is present)
RMST_prediction_pseudo <- predict.lm(pseudo_fit,se.fit = TRUE,interval = "confidence")
restricted_mean_CI_RMST_pseudo <- 0
RMST_which_out_CI <- numeric()
for(i in 1:nrow(Z_mu)){
  if(Mean_tau_T[i]>RMST_prediction_pseudo$fit[i,2] & Mean_tau_T[i]<RMST_prediction_pseudo$fit[i,3]){
    restricted_mean_CI_RMST_pseudo <- restricted_mean_CI_RMST_pseudo+1
  } else{
    RMST_which_out_CI <- append(RMST_which_out_CI,i)
  }
}
restricted_mean_CI_RMST_pseudo_empirical <- restricted_mean_CI_RMST_pseudo/nsubj
