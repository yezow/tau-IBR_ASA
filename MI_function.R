##################################################################
#Generate MI outputs
##################################################################

  coef_MI <- MI_impute_converge()[[1]] #obtain converged para estimates
  
  #step12: repeat bernoulli and uniform steps and get para estimators using 10 complete datasets
  Risk_set_1 <- MI_impute_risk_set(coef_MI) #fixed candidate set for imputation for I_i when coef are fixed
  coef_MI_para_10 <- list()
  var_MI_para_10 <- list()
  M_imputed_results <-foreach (i = 1:10, .combine=c,.packages = "boot") %dopar% {
   # source("call_libraries.R")
    set.seed(i)
    Impute_result_1 <- MI_impute_I_Y(Risk_set_1)
    I_i_impute_1 <- Impute_result_1[[1]]
    Y_i_impute_1 <- Impute_result_1[[2]]
    Y <- Data_cens$T_s
    Y[Data_cens$X>=tau & Data_cens$Delta==0]<-tau
    Y[row_MI][I_i_impute_1==1]<- tau
    Y[row_MI][I_i_impute_1==0]<-Y_i_impute_1[!is.na(Y_i_impute_1)]
    #re-estimate parameters
    z_mu <- Z_mu
    z_pi <- Z_pi
    #est_MI_10 <- nlminb(objective=nll_1,start = c(0,0,0.1,0,0.2,0,0,0,0,0,0,0,0,0.1,0,0,0,-0.5,0,0,0,0,2),control = list(iter.max=200),z_mu=z_mu,z_pi=z_pi,Y=Y)
    est_MI_10<-nlminb(objective=nll_1,start = initial_para_est,control = list(iter.max=200),z_mu=z_mu,z_pi=z_pi,Y=Y)
    list(est_MI_10$par[1:(2+ncol(Z_mu)+ncol(Z_pi))],Var_MI(est_MI_10$par,Y)[1:(2+ncol(Z_mu)+ncol(Z_pi)),1:(2+ncol(Z_mu)+ncol(Z_pi))],Y)
  }
  coef_MI_para_10<-M_imputed_results[c(1,4,7,10,13,16,19,22,25,28)]
  var_MI_para_10<-M_imputed_results[c(2,5,8,11,14,17,20,23,26,29)]
  MI_para_mean<-MI_combine(coef_MI_para_10,var_MI_para_10,10)[[1]] #combine 10 estimates of coefficients
  vcov_mean <- MI_combine(coef_MI_para_10,var_MI_para_10,10)[[2]]
  imputed_data<-M_imputed_results[c(3,6,9,12,15,18,21,24,27,30)]
  
  MI_result <-list(MI_para_mean,vcov_mean,imputed_data)
