
###########Fit model
set.seed(1234)
tau=365
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R")
MI_result_365<-MI_result
coef_MI<-MI_result_365[[1]]
var_MI<-MI_result_365[[2]]

EM_result_365<-MI_impute_converge()
w <- EM_result_365[[3]]
coef_EM_365 <-EM_result_365[[1]]
var_EM_365 <- EM_var(w,coef_EM_365)


############Fold change, CI and P values for beta regression using MI or EM
mu_model_results<-function(coef,var){
    fold_change_trt <- fold_change(coef[1:2])
    var_fold_change_trt <- fold_change_var(coef[1:2],var[1:2,1:2])
    fold_change_age <- fold_change(coef[c(1,3)])
    var_fold_change_age <- fold_change_var(coef[c(1,3)],var[c(1,3),c(1,3)])
    fold_change_gender <- fold_change(coef[c(1,4)])
    var_fold_change_gender <- fold_change_var(coef[c(1,4)],var[c(1,4),c(1,4)])
    fold_change_fev1 <- fold_change(coef[c(1,5)])
    var_fold_change_fev1 <- fold_change_var(coef[c(1,5)],var[c(1,5),c(1,5)])
    fold_change_smk <- fold_change(coef[c(1,6)])
    var_fold_change_smk <- fold_change_var(coef[c(1,6)],var[c(1,6),c(1,6)])
    fold_change_CI_lower<-c(fold_change_trt-1.96*sqrt(var_fold_change_trt),fold_change_age-1.96*sqrt(var_fold_change_age),
                            fold_change_gender-1.96*sqrt(var_fold_change_gender),fold_change_fev1-1.96*sqrt(var_fold_change_fev1),
                            fold_change_smk-1.96*sqrt(var_fold_change_smk))
    fold_change_CI_upper<-c(fold_change_trt+1.96*sqrt(var_fold_change_trt),fold_change_age+1.96*sqrt(var_fold_change_age),
                            fold_change_gender+1.96*sqrt(var_fold_change_gender),fold_change_fev1+1.96*sqrt(var_fold_change_fev1),
                            fold_change_smk+1.96*sqrt(var_fold_change_smk))
    fold_change_est <- c(fold_change_trt,fold_change_age,fold_change_gender,fold_change_fev1,fold_change_smk)
    P_mu<-numeric()
    for (i in c(2:6)){
      P_mu<-append(P_mu,(2*(1-pnorm(abs(coef[i]),mean=0,sd=sqrt(var[i,i])))))
    }

    beta_model_results<-rbind(fold_change_est,
                             fold_change_CI_lower,
                             fold_change_CI_upper,
                             P_mu)
    return(beta_model_results)
}
beta_model_results_MI<-mu_model_results(coef_MI,var_MI) ##beta regression results using MI
beta_model_results_EM<-mu_model_results(coef_EM_365,var_EM_365) ##beta regression results using EM


#####################################Odds ratio, CI and P values for logistic regression using MI or EM
pi_model_results<-function(coef,var){
    odds_ratio_est<-exp(coef[13:17])
    odds_ratio_lower <- numeric()
    odds_ratio_upper <- numeric()
    for (j in 13:17){
      odds_ratio_lower <- append(odds_ratio_lower,exp(coef[j]-1.96*sqrt(var[j,j])))
      odds_ratio_upper <- append(odds_ratio_upper,exp(coef[j]+1.96*sqrt(var[j,j])))
    }
    P_pi<-numeric()
    for (i in c(13:17)){
      P_pi<-append(P_pi,(2*(1-pnorm(abs(coef[i]),mean=0,sd=sqrt(var[i,i])))))
    }
    pi_model_result<-rbind(odds_ratio_est,
                            odds_ratio_lower,
                            odds_ratio_upper,
                              P_pi)
    return(pi_model_result)
}
pi_model_results_MI<-pi_model_results(coef_MI,var_MI) ##logistic regression results using MI
pi_model_results_EM<-pi_model_results(coef_EM_365,var_EM_365) ##logistic regression results using EM

###############################pseudo restriced mean model
pseudo <- pseudomean(Data_cens$X,Data_cens$Delta,tmax = tau)
pseudo_fit <- lm(pseudo~ trt+age_10_new+gender+fev1_10_new+nowsmk+
                   latitude_38+latitude_39+latitude_40+latitude_42+latitude_44,data=Data_cens) 
summary(pseudo_fit)
pseudo_fit$coefficients/tau
CI_RMST_lower<-pseudo_fit$coefficients-1.96*coef(summary(pseudo_fit))[, 2]
CI_RMST_upper<-pseudo_fit$coefficients+1.96*coef(summary(pseudo_fit))[, 2]
RMST_model_results<-rbind((pseudo_fit$coefficients/tau)[2:6],
                          (CI_RMST_lower/tau)[2:6],
                          (CI_RMST_upper/tau)[2:6],
                          as.numeric(summary(pseudo_fit)$coefficients[,4][2:6]) )
RMST_model_results ##RMST linear regression results using Pseudo obs

table3_result<-data.frame(rbind(beta_model_results_EM,pi_model_results_EM,beta_model_results_MI,pi_model_results_MI,RMST_model_results))
write.csv(table3_result,file = "Results/Table3_result.csv")