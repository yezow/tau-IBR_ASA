
#################################
#difference b/w mu and pi
#use one KM curves to impute I and Y
#risk set in first interation is defined based on RMST model

############################################
#generate Ti for each pair of Z_mu and Z_pi
############################################
Generate_T <- function(){
  T <- numeric()
  for (i in 1:nrow(Z_mu)){
    b<-rbern(1,p=pi[i])
    if (b==1){
      T[i] <- tau
    }
    else if (b==0){
      T[i] <- rbeta(1,alpha_b[i],beta_b[i])
    }
  }
  T[T<1] <- tau*(T[T<1])
  return(T)
}

########################################################
# generate censoring variable
########################################################
Generate_censoring <- function(){
  B_c <- rbinom(nsubj,1,0.56) #0.56 for 30% censoring
  U_c <- runif(nsubj,0,tau)
  
  C <- ifelse(B_c==0,U_c,100)
  return(C)
}

########################################################
# generate dependent censoring variable
########################################################
Generate_de_censoring <- function(){
  C<-numeric(length = nsubj)
  cens_id <- sample(which(Z2==1),nsubj*percent)
  C[cens_id]<-runif(length(cens_id),0,tau)
  C[C==0]<-100
  return(as.numeric(C))
}


Generate_de_censoring_weibull <- function(){
  C<- sapply(1:nsubj, function(x) rweibull(1,scale=2.2/(lambda_0*(sum(Z_pi[x,][1:ncol(Z_pi)]))),shape=2))
  return(as.numeric(C))
}


########################################################
# generate a new dataset with censoring
########################################################
generate_de_cens_dataset<- function(Y){
  T_s <- Y #generate survival time
  Censor <- Generate_de_censoring() #generate censoring time 
  Delta <- ifelse(Censor<T_s,0,1)
  ID <- seq(1:nsubj)
  Data_cens_original <- data.frame(ID,T_s,Censor,Delta,Z_mu,Z_pi)
  Data_cens_original$X <- apply(Data_cens_original[,c(2,3)],1,FUN=min) #generate obs time
  return(Data_cens_original)
}

generate_de_cens_dataset_weibull<- function(Y){
  T_s <- Y #generate survival time
  Censor <- Generate_de_censoring_weibull() #generate censoring time 
  Delta <- ifelse(Censor<T_s,0,1)
  ID <- seq(1:nsubj)
  Data_cens_original <- data.frame(ID,T_s,Censor,Delta,Z_mu,Z_pi)
  Data_cens_original$X <- apply(Data_cens_original[,c(2,3)],1,FUN=min) #generate obs time
  return(Data_cens_original)
}

generate_cens_dataset<- function(Y){
  T_s <- Y #generate survival time
  Censor <- Generate_censoring() #generate censoring time 
  Delta <- ifelse(Censor<T_s,0,1)
  ID <- seq(1:nsubj)
  Data_cens_original <- data.frame(ID,T_s,Censor,Delta,Z_mu,Z_pi)
  Data_cens_original$X <- apply(Data_cens_original[,c(2,3)],1,FUN=min) #generate obs time
  return(Data_cens_original)
}

############################################
#estimation
############################################
#specify Y,z_mu,z_pi before estimation

nll_1 <- function(para,z_mu,z_pi,Y) {

  LINP <- para[1]+z_mu%*%para[2:(ncol(z_mu)+1)]
  MU <- inv.logit(LINP)
  Alpha_b <- MU*para[ncol(z_mu)+ncol(z_pi)+3]
  Beta_b <- (1-MU)*para[ncol(z_mu)+ncol(z_pi)+3]
  
  #model2 for T_i>tau
  LINP1 <- para[(ncol(z_mu)+2)] + z_pi%*%para[(ncol(z_mu)+3):(ncol(z_mu)+ncol(z_pi)+2)]
  PI <- inv.logit(LINP1)
  
  
  #specify log likelihood
  LOGL <- numeric()
  for (i in 1:length(Y)){
    if (Y[i]<tau){
      X2 <- Y[i]/tau
      LOGL[i] <- log(1-PI[i])+log(1/tau)+dbeta(X2,Alpha_b[i],Beta_b[i],log=TRUE)
    }
    else if(Y[i]>=tau) {
      LOGL[i] <- log(PI[i])
    }
  }
  -sum(LOGL)
}

#################################
#hessian matrix of MI
#################################
Var_MI<-function(coef_MI,Y){
  mu_est<- inv.logit(coef_MI[1]+Z_mu%*%coef_MI[2:(ncol(Z_mu)+1)])
  pi_est<- inv.logit(coef_MI[(ncol(Z_mu)+2)]+Z_pi%*%coef_MI[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  nu_est<-coef_MI[length(coef_MI)]
  alpha_b <- mu_est*coef_MI[length(coef_MI)]
  beta_b <- (1-mu_est)*coef_MI[length(coef_MI)]
  b<-ifelse(Y>=tau,1,0)
  
  J_b <- lapply(1:nsubj, function(x) pi_est[x]*(pi_est[x]-1)*matrix(c(1,Z_pi[x,]))%*%c(1,Z_pi[x,]))
  J_b <- Reduce(`+`, J_b)
  J_a<-lapply(1:nsubj, function(x) (1-b[x])*
                  (nu_est*(logit(Y[x]/tau)-(digamma(mu_est[x]*nu_est)-digamma((1-mu_est[x])*nu_est)))
                   *(1-2*mu_est[x])*(mu_est[x]*(1-mu_est[x]))-
                     nu_est^2*mu_est[x]^2*(1-mu_est[x])^2*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est)))*
                  matrix(c(1,Z_mu[x,]))%*%c(1,Z_mu[x,]))
  J_a<-Reduce(`+`, J_a[-which(b==1)])
  J_nu<-sapply(1:nsubj, function(x) (1-b[x])*(-mu_est[x]^2*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est))+
                                                 2*mu_est[x]*trigamma((1-mu_est[x])*nu_est)+
                                                 trigamma(nu_est)-trigamma((1-mu_est[x])*nu_est)))
  J_nu<-sum(J_nu) 
  J_alpha_nu<-lapply(1:nsubj, function(x) (1-b[x])*(logit(Y[x]/tau)-(digamma(mu_est[x]*nu_est)-digamma((1-mu_est[x])*nu_est))-
                                                       nu_est*mu_est[x]*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est))+
                                                       nu_est*trigamma((1-mu_est[x])*nu_est))*mu_est[x]*(1-mu_est[x])*c(1,Z_mu[x,]))
  J_alpha_nu <- Reduce(`+`, J_alpha_nu[-which(b==1)])
  
  #######construct final observed information matrix             
  J_H<-matrix(0L,nrow=ncol(Z_mu)+ncol(Z_pi)+3,ncol=ncol(Z_mu)+ncol(Z_pi)+3)
  J_H[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]<-J_a
  J_H[(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi)),(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi))]<-J_b
  J_H[1:(ncol(Z_mu)+1),ncol(Z_mu)+ncol(Z_pi)+3]<-J_alpha_nu
  J_H[ncol(Z_mu)+ncol(Z_pi)+3,1:(ncol(Z_mu)+1)]<-J_alpha_nu
  J_H[ncol(Z_mu)+ncol(Z_pi)+3,ncol(Z_mu)+ncol(Z_pi)+3]<-J_nu
  
  J_final<- -(J_H)
  Var_MI <- solve(J_final) #beta, alpha, v
  return(Var_MI)
  }

########################################################
# Multiple Imputation functions
########################################################
##################################################################
#MI: step 6 to 10 ( bernoulli b_i and uniform u_i)
##################################################################
MI_impute_I_Y<- function(Risk_set){#,pi_new
  ##################
  #imputation for and I_i and Y_i
  ##################
  #step5: use KM plot to updata new pi
  I_i_impute<-numeric()
  Y_j_impute <- numeric()
  for (i in 1:length(row_MI)){ # 1:length(row_MI_1)
    risk_set_pi <- Risk_set[[i]]
    data_impute <- Data_cens[risk_set_pi,]
    KM <- survival::survfit(survival::Surv(X, Delta) ~ 1,  type="kaplan-meier", conf.type="log", data=data_impute)
    #KM <- survfit(Surv(X, Delta) ~ 1,  data=data_impute)
    
    #inverse transform
    v<-0
    while (v<=Data_cens$X[row_MI[i]]){
      u <- runif(1,0,1)
      if (!any(KM$surv<=u)){ ############# for example
        if(KM$time[length(KM$time)]>=tau){
          I_i_impute[i] <- 1
        } else{
          I_i_impute[i] <- 0
          Y_j_impute<-append(Y_j_impute,runif(1,KM$time[length(KM$time)],tau))
        }
        v<-1000
        #print ("<u")
      } else { #################for example
        v <- KM$time[KM$surv<=u][1]
        if(v>=tau){
          I_i_impute[i]<-1
          Y_j_impute[i]<-NA
        }else{
          I_i_impute[i]<-0
          Y_j_impute[i]<-v
        }
      }
    }
  }
  return(list(I_i_impute,Y_j_impute))
}

#############################
#get candidate set for I_i imputation
#####################################
MI_impute_risk_set <- function(coef_MI){
  #step3: impute initial value of pi and mu for each ith using estimated paras above
  LINP_MI <- coef_MI[1]+Z_mu%*%coef_MI[2:(1+ncol(Z_mu))]
  LINP1_MI <- coef_MI[ncol(Z_mu)+2]  + Z_pi%*%coef_MI[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
  PI_MI <- inv.logit(LINP1_MI)
  MU_MI <- inv.logit(LINP_MI)
  
  #Step4: Identify risk set using distance between pi and mu
  loop_num <- c(1:nrow(Z_mu))
  Risk_set <- list()
  m <- 0
  
  for (i in row_MI){
    m = m+1
    X_i <- Data_cens$X[i]
    risk_row <- numeric()
    epsilon <- 0.01
    set1 <- i
    #first restriction: the sample size in the risk set >15 and epsilon <0.5
    while (length(risk_row) < 10 & epsilon <=0.2){
      for (j in loop_num[-set1] ){
        abs_diff <- max(abs(PI_MI[j]-PI_MI[i]),abs(MU_MI[j]-MU_MI[i]))
        if (abs_diff <= epsilon&Data_cens$X[j]>X_i){
          risk_row <- append(risk_row,j)  #the row num of subj in risk set of ith
        }
      }
      epsilon <- epsilon+0.005
      set1<-append(i,risk_row)
    }
    
    while (length(risk_row)==0){
      epsilon <- epsilon+0.005
      for (j in loop_num[-set1] ){
        abs_diff <- max(abs(PI_MI[j]-PI_MI[i]),abs(MU_MI[j]-MU_MI[i]))
        if (abs_diff <= epsilon&Data_cens$X[j]>X_i){
          risk_row <- append(risk_row,j)  #the row num of subj in risk set of ith
        }
      }
      epsilon <- epsilon+0.005
      set1<-append(i,risk_row)
    }
    
    ##Second restriction: the last event time should be failure or attain tau, otherwise increase sample size
    last_event_time_in_R <- max(Data_cens$X[risk_row]) #get the last observed time in the risk set
    last_censor_in_R <- Data_cens$Delta[risk_row][which(Data_cens$X[risk_row]==max(Data_cens$X[risk_row]))] #get the last censoring status in risk set
    while(any(last_censor_in_R==0) &last_event_time_in_R<tau){
      #print(i)
      epsilon <- epsilon+0.005
      for (j in loop_num[-set1] ){
        abs_diff <- max(abs(PI_MI[j]-PI_MI[i]),abs(MU_MI[j]-MU_MI[i]))
        if (abs_diff <= epsilon&Data_cens$X[j]>X_i){
          risk_row <- append(risk_row,j)  #the row num of subj in risk set of ith
        }
      }
        last_event_time_in_R <- max(Data_cens$X[risk_row]) #get the last observed time in the risk set
        last_censor_in_R <- Data_cens$Delta[risk_row][which(Data_cens$X[risk_row]==max(Data_cens$X[risk_row]))] #get the last censoring status in risk set
        set1 <- append(i,risk_row)
        print("add more subjs in risk set")
    }
    Risk_set[[m]] <- risk_row
  }
  return(Risk_set)
}



########################################################
# Get initial para estimates using KM curve
MI_impute_0 <- function(){

  KM <- survival::survfit(Surv(X, Delta) ~ 1,  type="kaplan-meier", conf.type="log", data=Data_cens)

  I_i_impute_initial <- numeric()
  Y_j_impute_initial <- numeric()
  #inverse transform
  for (i in 1:length(row_MI)){
    v=0
    while (v<=Data_cens$X[row_MI[i]]){
      u <- runif(1,0,1)
      if (!any(KM$surv<=u)){ ############# for example
        if(KM$time[length(KM$time)]>=tau){
          I_i_impute_initial[i] <- 1
          Y_j_impute_initial[i]<-NA
        } else{
          I_i_impute_initial[i] <- 0
          Y_j_impute_initial[i]<-runif(1,KM$time[length(KM$time)],tau)
        }
        v<-1000
        #print ("<u")
      } else { 
        v <- KM$time[KM$surv<=u][1]
        if(v>=tau){
          I_i_impute_initial[i]<-1
          Y_j_impute_initial[i]<-NA
        }else{
          I_i_impute_initial[i]<-0
          Y_j_impute_initial[i]<-v
        }
      }
    }
  }
  
  Y <- Data_cens$T_s
  Y[Data_cens$X>=tau &Data_cens$Delta==0]<-tau
  Y[row_MI][which(I_i_impute_initial==1)]<- tau
  Y[row_MI][which(I_i_impute_initial==0)]<-Y_j_impute_initial[!is.na(Y_j_impute_initial)]
  #step11-2: re-estimate parameters
  z_mu <- Z_mu
  z_pi <- Z_pi
  #est_MI <- nlminb(objective=nll_1,start = c(0,0,0.1,0,0.2,0,0,0,0,0,0,0,0,0.1,0,0,0,0,-0.5,0,0,0,2),control = list(iter.max=200),z_mu=z_mu,z_pi=z_pi,Y=Y)
  est_MI<-nlminb(objective=nll_1,start = initial_para_est,control = list(iter.max=200),z_mu=Z_mu,z_pi=Z_pi,Y=Y)
  coef_MI_initial<- est_MI$par
  coef_var_initial<-solve(hessian(nll_1, est_MI$par,z_mu=Z_mu,z_pi=Z_pi,Y=Y))[1:(2+ncol(Z_mu)+ncol(Z_pi)),1:(2+ncol(Z_mu)+ncol(Z_pi))]
  #print(coef_MI_initial)
  return(list(coef_MI_initial,coef_var_initial))
}



########################################################
# Get converged para estimates using "EM" algorithm based on log likelihood
nll_EM <- function(para,z_mu,z_pi,w,Delta,case) {
  # para<-c(0,0,0.1,0,0,0.2,0,0,0,0,0,0,0.1,0,0,0,0,-0.5,0,0,0,0,2)
  # z_mu<-Z_mu
  # z_pi<-Z_pi
  # w<-w
  # Delta <- Data_cens$Delta
  # case <-  Data_cens$case
  # para=initial_para
  
  LINP <- para[1]+z_mu%*%para[2:(ncol(z_mu)+1)]
  MU <- inv.logit(LINP)
  Alpha_b <- MU*para[ncol(z_mu)+ncol(z_pi)+3]
  Beta_b <- (1-MU)*para[ncol(z_mu)+ncol(z_pi)+3]
  LINP1 <- para[(ncol(z_mu)+2)] + z_pi%*%para[(ncol(z_mu)+3):(ncol(z_mu)+ncol(z_pi)+2)]
  PI <- inv.logit(LINP1)
  
  #specify log likelihood
  LOGL <- numeric()
  for (i in 1:nrow(z_mu)){
    if (case[i]==1){
      LOGL[i]<-log(1-PI[i])+dbeta(Data_cens$X[i]/tau,Alpha_b[i],Beta_b[i],log=TRUE)
      } else if (case[i]==2){
      LOGL[i]<-log(PI[i])
      } else if (case[i]==3){
      LOGL[i]<-(1-w[i])*log((1-PI[i])*(1-pbeta(Data_cens$X[i]/tau,Alpha_b[i],Beta_b[i])))+w[i]*log(PI[i])}
    }
  -sum(LOGL)
}


MI_impute_converge <- function(){
  initial_result <- MI_impute_0() #get para estimates and var estimates based on risk set using KM only
  initial_para <- initial_result[[1]]
  mu_est<- inv.logit(initial_para[1]+Z_mu%*%initial_para[2:(ncol(Z_mu)+1)])
  pi_est<- inv.logit(initial_para[(ncol(Z_mu)+2)]+Z_pi%*%initial_para[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  alpha_b <- mu_est*initial_para[ncol(Z_mu)+ncol(Z_pi)+3]
  beta_b <- (1-mu_est)*initial_para[ncol(Z_mu)+ncol(Z_pi)+3]
  
  #define the case 1 for subjs who are observed to experience the event
  # case 2: observed to attain tau
  #case 3:censored before tau
  Data_cens$case<-ifelse(Data_cens$X>=tau,2,3)
  Data_cens$case[Data_cens$Delta==1&Data_cens$X<tau]<-1
   
  diff_coef_converge<-3
  #converge_time<-0  
  while (diff_coef_converge>1e-4 ) {
    E_b_i<- as.numeric(lapply(row_MI,function(x) pi_est[x]/(pi_est[x]+(1-pi_est[x])*(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x])))))
    w<-as.numeric(lapply(1:nrow(Z_mu), function(x) ifelse(Data_cens$X[x]>=tau,1,0)))
    w[row_MI]<-E_b_i
    est_EM<-nlminb(objective=nll_EM,start = initial_para_est,control = list(iter.max=200),z_mu=Z_mu,z_pi=Z_pi,w=w,Delta=Data_cens$Delta,case=Data_cens$case)
    #est_EM<-nlminb(objective=nll_EM,start = c(0,0,0.1,0,0.2,0,0,0,0,0,0,0,0,0.1,0,0,0,0,-0.5,0,0,0,2),control = list(iter.max=200),z_mu=Z_mu,z_pi=Z_pi,w=w,Delta=Data_cens$Delta,case=Data_cens$case)
    coef_converge<- est_EM$par
    diff_coef_converge <- max(abs(coef_converge[1:(2+ncol(Z_mu)+ncol(Z_pi))]-initial_para[1:(2+ncol(Z_mu)+ncol(Z_pi))]))
    diff_coef_coverge_pi <- max(abs(coef_converge[(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi))]-initial_para[(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi))]))
    diff_coef_coverge_mu <- max(abs(coef_converge[1:(1+ncol(Z_mu))]-initial_para[1:(1+ncol(Z_mu))]))
    #print(diff_coef_converge)
    initial_para <- coef_converge
    mu_est<- inv.logit(initial_para[1]+Z_mu%*%initial_para[2:(ncol(Z_mu)+1)])
    pi_est<- inv.logit(initial_para[(ncol(Z_mu)+2)]+Z_pi%*%initial_para[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
    alpha_b <- mu_est*initial_para[length(initial_para)]
    beta_b <- (1-mu_est)*initial_para[length(initial_para)]
    print(diff_coef_converge)
  }
  #update w using the last EM para estimates
  E_b_i<- as.numeric(lapply(row_MI,function(x) pi_est[x]/(pi_est[x]+(1-pi_est[x])*(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x])))))
  w<-as.numeric(lapply(1:nrow(Z_mu), function(x) ifelse(Data_cens$X[x]>=tau,1,0)))
  w[row_MI]<-E_b_i
  var_converge<-solve(hessian(nll_EM, coef_converge,z_mu=Z_mu,z_pi=Z_pi,w=w,Delta=Data_cens$Delta,case=Data_cens$case))[1:(2+ncol(Z_mu)+ncol(Z_pi)),1:(2+ncol(Z_mu)+ncol(Z_pi))]
  return(list(coef_converge,var_converge,w))
}


##################################################################
#MI: combine multiple imputation results
##################################################################
MI_combine <- function(coef_list,vcov_list,m){
  para_data<-data.frame(coef_list)
  para_bar<-rowMeans(para_data)
  vcov_bar <- vcov_list[[1]]
  B <- (coef_list[[1]]-para_bar)%*%t(coef_list[[1]]-para_bar) #calculate the variance between the MI data sets
  for (i in 2:m){
    vcov_bar <- vcov_bar+vcov_list[[i]]
    B <- B+(coef_list[[i]]-para_bar)%*%t(coef_list[[i]]-para_bar)
  }
  vcov_bar <- vcov_bar/m
  B <-B/(m-1)
  
  #calculate the total variance
  comb.vcov <- vcov_bar+(1+1/m)*B
  return(list(para_bar,comb.vcov))
}


##########################################
#calculate the variance of restricted mean
##########################################
Var_restricted_mean_func <- function(coef_vcov,coef){
  beta_vcov <- coef_vcov[(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2),(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2)]
  alpha_vcov<- coef_vcov[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]
  var_restricted_mean <- numeric()
  var_pi <- numeric()
  var_mu <- numeric()
  for (i in 1:nrow(Z_pi)){#length(Z_pi)
    #calculate variance of pi
    logit_mu_est <- coef[1]+Z_mu[i,]%*%coef[2:(ncol(Z_mu)+1)]
    logit_pi_est <- coef[(ncol(Z_mu)+2)]+Z_pi[i,]%*%coef[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
    X_pi <- matrix(c(1,Z_pi[i,]), nrow = 1+ncol(Z_pi))
    var_logit_pi <- t(X_pi) %*% beta_vcov %*% X_pi
    var_pi[i] <- var_logit_pi*(exp(-logit_pi_est)^2)/(1+exp(-logit_pi_est))^4
    
    #calculate the first term in variance of restricted mean
    var_term_1<-tau^2*(1-inv.logit(logit_mu_est))^2*var_pi[i]
    
    #calculate variance of logit of pi
    X_mu <- matrix(c(1,Z_mu[i,]),nrow=1+ncol(Z_mu))
    var_logit_mu <- t(X_mu) %*% alpha_vcov %*% X_mu
    var_mu[i] <- var_logit_mu*(exp(-logit_mu_est)^2)/(1+exp(-logit_mu_est))^4
    
    #calculate the derivative of transformation function
    g_deriv <- tau*(1-1/(1+exp(-logit_pi_est)))
    
    #calculate the second term in variance of restricted mean
    var_term_2<-g_deriv^2*var_mu[i]
    
    var_restricted_mean[i] <- var_term_1+var_term_2
  }
  return(list(var_restricted_mean,var_pi,var_mu))
}


################################################
#fold change variance computation
fold_change <- function(coef){
  beta<-coef[2]
  exp_sum <- exp(sum(coef))
  fold<-(exp(beta)+exp_sum)/(1+exp_sum)
  return(fold)
}
fold_change_var<-function(coef,coef_vcov){
  exp_sum <- exp(sum(coef))
  beta<-coef[2]
  gradient_beta <- (exp(beta)+exp_sum)/(1+exp_sum)-(exp(beta)+exp_sum)*exp_sum/(1+exp_sum)^2
  gradient_intercept<-exp_sum/(1+exp_sum)-(exp(beta)+exp_sum)*exp_sum/(1+exp_sum)^2
  gradient <- matrix(c(gradient_intercept,gradient_beta),nrow=1)
  var <- gradient%*%coef_vcov%*%t(gradient)
  return(var)
}

############################################
# EM variance (Louis method)
############################################
EM_var<-function(w,coef_EM){
  ######calculate the var(U_beta) w.r.t unkown B_i
  Var_U_b <- lapply(row_MI, function(x) matrix(c(1,Z_pi[x,]))%*%c(1,Z_pi[x,])*w[x]*(1-w[x]))
  Var_U_b <- Reduce(`+`, Var_U_b)
  
  ######calculate the var(U_alpha) w.r.t unkown B_i
  mu_est<- inv.logit(coef_EM[1]+Z_mu%*%coef_EM[2:(ncol(Z_mu)+1)])
  pi_est<- inv.logit(coef_EM[(ncol(Z_mu)+2)]+Z_pi%*%coef_EM[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  nu_est<-coef_EM[length(coef_EM)]
  alpha_b <- mu_est*coef_EM[length(coef_EM)]
  beta_b <- (1-mu_est)*coef_EM[length(coef_EM)]
  
  #calculate the derivative of cdf of beta distribution w.r.t alpha
  deriv_alpha<-list()
  for (i in 1:nsubj){
    func_alpha<-function(alpha){pbeta(Data_cens$X[i]/tau, 
                                      inv.logit(append(1,Z_mu[i,])%*%alpha)*coef_EM[length(coef_EM)],
                                      (1-inv.logit(append(1,Z_mu[i,])%*%alpha))*coef_EM[length(coef_EM)])}
    alpha<-coef_EM[1:(ncol(Z_mu)+1)]
    deriv_alpha[[i]]<-grad(func_alpha,alpha,method="Richardson")
  } #all 0s are < 1e-22
  
  Var_U_a <- lapply(row_MI, function(x) (w[x]*(1-w[x]))/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2*deriv_alpha[[x]]%*%t(deriv_alpha[[x]]))
  Var_U_a <- Reduce(`+`, Var_U_a)
  
  ######calculate the var(U_nu) w.r.t unkown B_i
  deriv_nu<-numeric()
  for (i in 1:nsubj){
    func_nu<-function(nu){pbeta(Data_cens$X[i]/tau,mu_est[i]*nu,(1-mu_est[i])*nu)}
    nu<-coef_EM[length(coef_EM)]
    deriv_nu<-append(deriv_nu,grad(func_nu,nu,method="Richardson"))
  }
  
  Var_U_nu <- sum(sapply(row_MI, function(x) (w[x]*(1-w[x]))/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2*deriv_nu[x]^2))
  
  ######calculate the cov(U_alpha,U_nu) w.r.t unkown B_i
  Var_U_alpha_nu <- lapply(row_MI, function(x) (w[x]*(1-w[x]))/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2*deriv_alpha[[x]]*deriv_nu[x])
  Var_U_alpha_nu <-Reduce(`+`, Var_U_alpha_nu)
  
  ######calculate the J_beta
  J_b<-lapply(1:nsubj, function(x) pi_est[x]*(pi_est[x]-1)*matrix(c(1,Z_pi[x,]))%*%c(1,Z_pi[x,]))
  J_b <- Reduce(`+`, J_b)
  
  ######calculate the J_alpha
  row_n1<-1:nsubj
  row_n1 <- row_n1[-row_MI]
  #calculate hessian matrix of cdf of alpha
  Hessian_alpha<-list()
  for (i in 1:nsubj){
    func_alpha<-function(alpha){pbeta(Data_cens$X[i]/tau, 
                                      inv.logit(append(1,Z_mu[i,])%*%alpha)*coef_EM[length(coef_EM)],
                                      (1-inv.logit(append(1,Z_mu[i,])%*%alpha))*coef_EM[length(coef_EM)])}
    alpha<-coef_EM[1:(ncol(Z_mu)+1)]
    Hessian_alpha[[i]]<-hessian(func_alpha,alpha,method="Richardson")
  }
  
  J_a_1<-lapply(row_n1, function(x) (1-w[x])*
                  (nu_est*(logit(Data_cens$X[x]/tau)-(digamma(mu_est[x]*nu_est)-digamma((1-mu_est[x])*nu_est)))
                   *(1-2*mu_est[x])*(mu_est[x]*(1-mu_est[x]))-
                     nu_est^2*mu_est[x]^2*(1-mu_est[x])^2*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est)))*
                  matrix(c(1,Z_mu[x,]))%*%c(1,Z_mu[x,]))
  J_a_1<-Reduce(`+`, J_a_1[-which(w[row_n1]==1)])
  J_a_2<-lapply(row_MI, function(x) (w[x]-1)/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))*
                  Hessian_alpha[[x]]+(w[x]-1)/((1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2)*
                  deriv_alpha[[x]]%*%t(deriv_alpha[[x]]))
  J_a_2<-Reduce(`+`, J_a_2) 
  J_a <-J_a_1+J_a_2 ########add 1 into Z_mu and Z_pi
  
  ######calculate the J_nu
  J_nu_1<-sapply(row_n1, function(x) (1-w[x])*(-mu_est[x]^2*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est))+
                                                 2*mu_est[x]*trigamma((1-mu_est[x])*nu_est)+
                                                 trigamma(nu_est)-trigamma((1-mu_est[x])*nu_est)))
  J_nu_1<-sum(J_nu_1)     
  
  #calculate second derivative of cdf of nu
  Hessian_nu<-numeric()
  for (i in 1:nsubj){
    func_nu<-function(nu){pbeta(Data_cens$X[i]/tau,mu_est[i]*nu,(1-mu_est[i])*nu)}
    nu<-coef_EM[length(coef_EM)]
    Hessian_nu<-append(Hessian_nu,hessian(func_nu,nu,method="Richardson"))
  }
  J_nu_2<-sapply(row_MI, function(x) (w[x]-1)/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))*
                   Hessian_nu[x]+(w[x]-1)/((1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2)*
                   deriv_nu[x]^2)              
  J_nu<-J_nu_1+sum(J_nu_2)            
  
  #########stop here
  ######calculate the J_alpha_nu
  J_alpha_nu_1<-lapply(row_n1, function(x) (1-w[x])*(logit(Data_cens$X[x]/tau)-(digamma(mu_est[x]*nu_est)-digamma((1-mu_est[x])*nu_est))-
                                                       nu_est*mu_est[x]*(trigamma(mu_est[x]*nu_est)+trigamma((1-mu_est[x])*nu_est))+
                                                       nu_est*trigamma((1-mu_est[x])*nu_est))*mu_est[x]*(1-mu_est[x])*c(1,Z_mu[x,]) )
  #calculate the hessian matrix of alpha&nu
  Hessian_alpha_nu<-list()
  for (i in 1:nsubj){
    func_alpha_nu<-function(vec_alpha_nu){pbeta(Data_cens$X[i]/tau, 
                                                inv.logit(append(1,Z_mu[i,])%*%vec_alpha_nu[1:(1+ncol(Z_mu))])*vec_alpha_nu[length(vec_alpha_nu)],
                                                (1-inv.logit(append(1,Z_mu[i,])%*%vec_alpha_nu[1:(1+ncol(Z_mu))]))*vec_alpha_nu[length(vec_alpha_nu)])}
    vec_alpha_nu<-coef_EM[c(1:(1+ncol(Z_mu)),length(coef_EM))]
    Hessian_alpha_nu[[i]]<-hessian(func_alpha_nu,vec_alpha_nu,method="Richardson")[,length(vec_alpha_nu)][-length(vec_alpha_nu)]
  }
  
  J_alpha_nu_2<-lapply(row_MI,function(x) (w[x]-1)/(1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))*
                         Hessian_alpha_nu[[x]]+(w[x]-1)/((1-pbeta(Data_cens$X[x]/tau, alpha_b[x],beta_b[x]))^2)*
                         deriv_alpha[[x]]*deriv_nu[[x]])   
  J_alpha_nu <- Reduce(`+`, J_alpha_nu_1[-which(w[row_n1]==1)])+Reduce(`+`,J_alpha_nu_2)
  
  #######construct final observed information matrix             
  COV_H<-matrix(0L,nrow=ncol(Z_mu)+ncol(Z_pi)+3,ncol=ncol(Z_mu)+ncol(Z_pi)+3)
  COV_H[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]<-Var_U_a
  COV_H[(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi)),(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi))]<-Var_U_b
  COV_H[1:(ncol(Z_mu)+1),ncol(Z_mu)+ncol(Z_pi)+3]<-Var_U_alpha_nu
  COV_H[ncol(Z_mu)+ncol(Z_pi)+3,1:(ncol(Z_mu)+1)]<-Var_U_alpha_nu
  COV_H[ncol(Z_mu)+ncol(Z_pi)+3,ncol(Z_mu)+ncol(Z_pi)+3]<-Var_U_nu
  
  J_H<-matrix(0L,nrow=ncol(Z_mu)+ncol(Z_pi)+3,ncol=ncol(Z_mu)+ncol(Z_pi)+3)
  J_H[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]<-J_a
  J_H[(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi)),(ncol(Z_mu)+2):(ncol(Z_mu)+2+ncol(Z_pi))]<-J_b
  J_H[1:(ncol(Z_mu)+1),ncol(Z_mu)+ncol(Z_pi)+3]<-J_alpha_nu
  J_H[ncol(Z_mu)+ncol(Z_pi)+3,1:(ncol(Z_mu)+1)]<-J_alpha_nu
  J_H[ncol(Z_mu)+ncol(Z_pi)+3,ncol(Z_mu)+ncol(Z_pi)+3]<-J_nu
  
  J_final<- -(COV_H+J_H)
  Var_EM <- solve(J_final) #beta, alpha, v
  return(Var_EM)
}








