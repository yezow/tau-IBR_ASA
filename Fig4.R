######Estimate RMST for trt/placebo group for patients with population mean
RMST_pi_func<-function(MI_result,new_z_mu,tau){
      mu_est_MI <- inv.logit(MI_result[[1]][1]+new_z_mu%*%c(MI_result[[1]][2:(ncol(new_z_mu)+1)]))
      pi_est_MI <- inv.logit(MI_result[[1]][ncol(new_z_mu)+2]+new_z_mu%*%c(MI_result[[1]][(ncol(new_z_mu)+3):(ncol(new_z_mu)+ncol(new_z_mu)+2)]))
      Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
      
      beta_vcov<-MI_result[[2]][(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2),(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2)]
      X_pi <- matrix(c(1,Z_pi), nrow = 1+ncol(Z_pi))
      var_logit_pi <- t(X_pi) %*% beta_vcov %*% X_pi
      sd_pi<- sqrt(var_logit_pi*(exp(-logit(pi_est_MI))^2)/(1+exp(-logit(pi_est_MI)))^4)
      pi_lower<-as.numeric(pi_est_MI)-1.96*as.numeric(sd_pi)
      pi_higher<-as.numeric(pi_est_MI)+1.96*as.numeric(sd_pi)
      
      #calculate variance of logit of pi
      alpha_vcov<- MI_result[[2]][1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]
      X_mu <- matrix(c(1,Z_mu), nrow = 1+ncol(Z_mu))
      var_logit_mu <- t(X_mu) %*% alpha_vcov %*% X_mu
      sd_mu <- sqrt(var_logit_mu*(exp(-logit(mu_est_MI))^2)/(1+exp(-logit(mu_est_MI)))^4)
      mu_lower<-as.numeric(mu_est_MI)-1.96*as.numeric(sd_mu)
      mu_higher<-as.numeric(mu_est_MI)+1.96*as.numeric(sd_mu)
      
      sd_restricted_mean_MI_est<-sqrt(Var_restricted_mean_func(MI_result[[2]],MI_result[[1]])[[1]])
      Mean_tau_T_est_MI_lower<-as.numeric(Mean_tau_T_est_MI)-1.96*as.numeric(sd_restricted_mean_MI_est)
      Mean_tau_T_est_MI_higher<-  as.numeric(Mean_tau_T_est_MI)+1.96*as.numeric(sd_restricted_mean_MI_est)
     return(list(as.numeric(Mean_tau_T_est_MI),
                 as.numeric(Mean_tau_T_est_MI_lower),
                 as.numeric(Mean_tau_T_est_MI_higher),
                 as.numeric(pi_est_MI),
                 as.numeric(pi_lower),
                as.numeric(pi_higher),
                as.numeric(sd_restricted_mean_MI_est),
                as.numeric(sd_pi),
                as.numeric(mu_est_MI),
                as.numeric(mu_lower),
                as.numeric(mu_higher),
                as.numeric(sd_mu)))     
}

######Fit model
set.seed(1234)
tau=365
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R")
MI_result_365<-MI_result
new_z_mu<-matrix(c(1,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_trt_365<-unlist(RMST_pi_func(MI_result_365,new_z_mu,tau))
new_z_mu<-matrix(c(0,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_placebo_365<-unlist(RMST_pi_func(MI_result_365,new_z_mu,tau))


tau=90
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R")
MI_result_90<-MI_result
new_z_mu<-matrix(c(1,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_trt_90<-unlist(RMST_pi_func(MI_result_90,new_z_mu,tau))
new_z_mu<-matrix(c(0,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_placebo_90<-unlist(RMST_pi_func(MI_result_90,new_z_mu,tau))

tau=180
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R")
MI_result_180<-MI_result
new_z_mu<-matrix(c(1,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_trt_180<-unlist(RMST_pi_func(MI_result_180,new_z_mu,tau))
new_z_mu<-matrix(c(0,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_placebo_180<-unlist(RMST_pi_func(MI_result_180,new_z_mu,tau))

tau=270
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R")
MI_result_270<-MI_result
new_z_mu<-matrix(c(1,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_trt_270<-unlist(RMST_pi_func(MI_result_270,new_z_mu,tau))
new_z_mu<-matrix(c(0,0,gender_mean,0,nowsmk_mean,center2_mean,center3_mean,center4_mean,center5_mean,center6_mean),nrow=1)
Z_mu<-new_z_mu
Z_pi<-new_z_mu
results_placebo_270<-unlist(RMST_pi_func(MI_result_270,new_z_mu,tau))

####################
#calculate p values for each comparison
2*pnorm(-abs((results_trt_90[1]-results_placebo_90[1])/sqrt(results_trt_90[7]^2+results_placebo_90[7]^2)))
2*pnorm(-abs((results_trt_180[1]-results_placebo_180[1])/sqrt(results_trt_180[7]^2+results_placebo_180[7]^2)))
2*pnorm(-abs((results_trt_270[1]-results_placebo_270[1])/sqrt(results_trt_270[7]^2+results_placebo_270[7]^2)))
2*pnorm(-abs((results_trt_365[1]-results_placebo_365[1])/sqrt(results_trt_365[7]^2+results_placebo_365[7]^2)))

2*pnorm(-abs((results_trt_90[4]-results_placebo_90[4])/sqrt(results_trt_90[8]^2+results_placebo_90[8]^2)))
2*pnorm(-abs((results_trt_180[4]-results_placebo_180[4])/sqrt(results_trt_180[8]^2+results_placebo_180[8]^2)))
2*pnorm(-abs((results_trt_270[4]-results_placebo_270[4])/sqrt(results_trt_270[8]^2+results_placebo_270[8]^2)))
2*pnorm(-abs((results_trt_365[4]-results_placebo_365[4])/sqrt(results_trt_365[8]^2+results_placebo_365[8]^2)))

2*pnorm(-abs((results_trt_90[9]-results_placebo_90[9])/sqrt(results_trt_90[12]^2+results_placebo_90[12]^2)))
2*pnorm(-abs((results_trt_180[9]-results_placebo_180[9])/sqrt(results_trt_180[12]^2+results_placebo_180[12]^2)))
2*pnorm(-abs((results_trt_270[9]-results_placebo_270[9])/sqrt(results_trt_270[12]^2+results_placebo_270[12]^2)))
2*pnorm(-abs((results_trt_365[9]-results_placebo_365[9])/sqrt(results_trt_365[12]^2+results_placebo_365[12]^2)))
####################
#RMST curve
# library(ggplot2)
# library(tikzDevice)
# library('latex2exp')
# library(ggpubr)
# library(gridExtra)
# library(grid)
x_axis<- rep(c(3,6,9,12),each=2)
group<-rep(c("t-IBR (MI) Azithromycin","t-IBR (MI) Placebo"),4)
mean <- c(results_trt_90[1],results_placebo_90[1],results_trt_180[1],results_placebo_180[1],results_trt_270[1],results_placebo_270[1],results_trt_365[1],results_placebo_365[1])
c_lower<-c(results_trt_90[2],results_placebo_90[2],results_trt_180[2],results_placebo_180[2],results_trt_270[2],results_placebo_270[2],results_trt_365[2],results_placebo_365[2])
c_upper<-c(results_trt_90[3],results_placebo_90[3],results_trt_180[3],results_placebo_180[3],results_trt_270[3],results_placebo_270[3],results_trt_365[3],results_placebo_365[3])
mean<-mean/30
c_lower <- c_lower/30 #change it to months
c_upper <- c_upper/30
data_plot<- data.frame(x_axis,group,mean,c_lower,c_upper)
data_plot$x_axis<-as.factor(data_plot$x_axis)
data_plot$group <- factor(data_plot$group)
#plot1
options(repr.plot.width = 14, repr.plot.height = 8)
pd <- position_dodge(0.1) # move them .05 to the left and right
plot1<-ggplot(data_plot, aes(x=x_axis, y=mean, group=group,linetype = group,color=group)) + 
  geom_errorbar(aes(ymin=c_lower, ymax=c_upper), width=.3,position=pd,size=0.9) +
  geom_line(position=pd,size=0.9) +
  ylim(2,9)+
  geom_point(position=pd,size=0.7)+ 
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c("black","blue"))+
  labs(x =TeX("Varying Choices of $\\tau$ (months)"), y = TeX("Restricted Mean Survival Time (months)"))+
  theme(legend.position="top",legend.title = element_blank())
plot1

##########pi curve
mean <- c(results_trt_90[4],results_placebo_90[4],results_trt_180[4],results_placebo_180[4],results_trt_270[4],results_placebo_270[4],results_trt_365[4],results_placebo_365[4])
c_lower<-c(results_trt_90[5],results_placebo_90[5],results_trt_180[5],results_placebo_180[5],results_trt_270[5],results_placebo_270[5],results_trt_365[5],results_placebo_365[5])
c_upper<-c(results_trt_90[6],results_placebo_90[6],results_trt_180[6],results_placebo_180[6],results_trt_270[6],results_placebo_270[6],results_trt_365[6],results_placebo_365[6])
data_plot<- data.frame(x_axis,group,mean,c_lower,c_upper)
data_plot$x_axis<-as.factor(data_plot$x_axis)
data_plot$group <- factor(data_plot$group)
#plot1
options(repr.plot.width = 14, repr.plot.height = 8)
pd <- position_dodge(0.1) # move them .05 to the left and right
plot2<-ggplot(data_plot, aes(x=x_axis, y=mean, group=group,linetype = group,color=group)) + 
  geom_errorbar(aes(ymin=c_lower, ymax=c_upper), width=.3,position=pd,size=0.9) +
  geom_line(position=pd,size=0.9) +
  ylim(0,1)+
  geom_point(position=pd,size=0.7)+ 
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c("black","blue"))+
  labs(x =TeX("Varying Choices of $\\tau$ (months)"), y = TeX("Probability of Being Event-free at $\\tau$"))+
  theme(legend.position="top",legend.title = element_blank())
plot2

##########mu curve
mean <- c(results_trt_90[9],results_placebo_90[9],results_trt_180[9],results_placebo_180[9],results_trt_270[9],results_placebo_270[9],results_trt_365[9],results_placebo_365[9])
c_lower<-c(results_trt_90[10],results_placebo_90[10],results_trt_180[10],results_placebo_180[10],results_trt_270[10],results_placebo_270[10],results_trt_365[10],results_placebo_365[10])
c_upper<-c(results_trt_90[11],results_placebo_90[11],results_trt_180[11],results_placebo_180[11],results_trt_270[11],results_placebo_270[11],results_trt_365[11],results_placebo_365[11])
data_plot<- data.frame(x_axis,group,mean,c_lower,c_upper)
data_plot$x_axis<-as.factor(data_plot$x_axis)
data_plot$group <- factor(data_plot$group)
#plot1
options(repr.plot.width = 14, repr.plot.height = 8)
pd <- position_dodge(0.1) # move them .05 to the left and right
plot3<-ggplot(data_plot, aes(x=x_axis, y=mean, group=group,linetype = group,color=group)) + 
  geom_errorbar(aes(ymin=c_lower, ymax=c_upper), width=.3,position=pd,size=0.9) +
  geom_line(position=pd,size=0.9) +
  ylim(0,1)+
  geom_point(position=pd,size=0.7)+ 
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_color_manual(values=c("black","blue"))+
  labs(x =TeX("Varying Choices of $\\tau$ (months)"), y = TeX("Mean Survival Time Divided by $\\tau$ When Below $\\tau$"))+
  theme(legend.position="top",legend.title = element_blank())
plot3
pdf("Results/fig4.pdf", width = 12, height = 5) # Open a new pdf file       # smaller font size
grid.arrange(plot1, plot2,plot3, ncol=3)
#14: 4.5