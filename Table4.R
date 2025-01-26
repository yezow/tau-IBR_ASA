
############Fold change, CI and P values for beta regression using MI or EM
mu_model_results<-function(coef,var){
  fold_change_trt <- fold_change(coef[1:2])
  var_fold_change_trt <- fold_change_var(coef[1:2],var[1:2,1:2])
  fold_change_CI_lower<-fold_change_trt-1.96*sqrt(var_fold_change_trt)
  fold_change_CI_upper<-fold_change_trt+1.96*sqrt(var_fold_change_trt)
  P_mu<-2*(1-pnorm(abs(coef[2]),mean=0,sd=sqrt(var[2,2])))
  beta_model_results<-rbind(fold_change_trt,
                            fold_change_CI_lower,
                            fold_change_CI_upper,
                            P_mu)
  return(beta_model_results)
}


##########################Odds ratio, CI and P values for logistic regression using MI or EM
pi_model_results<-function(coef,var){
  odds_ratio_est<-exp(coef[13])
  odds_ratio_lower <- exp(coef[13]-1.96*sqrt(var[13,13]))
  odds_ratio_upper <- exp(coef[13]+1.96*sqrt(var[13,13]))
  P_pi<-2*(1-pnorm(abs(coef[13]),mean=0,sd=sqrt(var[13,13])))
  pi_model_result<-rbind(odds_ratio_est,
                         odds_ratio_lower,
                         odds_ratio_upper,
                         P_pi)
  return(pi_model_result)
}


###########Fit model
results_table4<-vector()
for (tau in c(90,180,270,365)){
    set.seed(1234)
    source("Scripts/EX_figure_source.R")
    source("Scripts/MI_function.R") #run MI imputation and get MI_result
    #MI_result <- MI_impute()
    MI_result_365<-MI_result
    coef_MI<-MI_result_365[[1]]
    var_MI<-MI_result_365[[2]]
    
    EM_result_365<-MI_impute_converge()
    w <- EM_result_365[[3]]
    coef_EM_365 <-EM_result_365[[1]]
    var_EM_365 <- EM_var(w,coef_EM_365)
    
    beta_model_results_MI<-mu_model_results(coef_MI,var_MI) ##beta regression results using MI
    beta_model_results_EM<-mu_model_results(coef_EM_365,var_EM_365) ##beta regression results using EM
    
    pi_model_results_MI<-pi_model_results(coef_MI,var_MI) ##logistic regression results using MI
    pi_model_results_EM<-pi_model_results(coef_EM_365,var_EM_365) ##logistic regression results using EM
    
    ###############################pseudo restricted mean model
    pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
    pseudo_fit <- lm(pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+
                       latitude_38+latitude_39+latitude_40+latitude_42+latitude_44,data=Data_cens) 
    CI_RMST_lower<-pseudo_fit$coefficients-1.96*coef(summary(pseudo_fit))[, 2]
    CI_RMST_upper<-pseudo_fit$coefficients+1.96*coef(summary(pseudo_fit))[, 2]
    RMST_model_results<-rbind((pseudo_fit$coefficients/tau)[2],
                              (CI_RMST_lower/tau)[2],
                              (CI_RMST_upper/tau)[2],
                              as.numeric(summary(pseudo_fit)$coefficients[,4][2]) )
    RMST_model_results ##RMST linear regression results using Pseudo obs
    
    results_table4<-cbind(results_table4,rbind(beta_model_results_MI,pi_model_results_MI,beta_model_results_EM,pi_model_results_EM,RMST_model_results))
}
colnames(results_table4)<-c("3 months","6 months","9 months","12 months")
write.csv(results_table4,file = "Results/Table4_result.csv")
