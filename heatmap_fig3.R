
###########Fit model
set.seed(1234)
tau=365
source("Scripts/EX_figure_source.R")
source("Scripts/MI_function.R") #run MI imputation and get MI_result
#########heatmap
mu_est_MI <- inv.logit(MI_result[[1]][1]+Z_mu%*%c(MI_result[[1]][2:(ncol(Z_mu)+1)]))
pi_est_MI <- inv.logit(MI_result[[1]][ncol(Z_mu)+2]+Z_pi%*%c(MI_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]))
Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI

Imputed_data<-matrix(unlist(MI_result[[3]]),nrow=1111)
mean_RMST<-rowMeans(Imputed_data)
heatmap_data<-dplyr::select(Data_cens,trt,age_10_new,gender,fev1_10_new,nowsmk)
heatmap_data$mean_RMST_365<-mean_RMST
heatmap_data$age<-(heatmap_data$age_10_new+mean(ex_data$age))*10
heatmap_data$fev1<-(heatmap_data$fev1_10_new+mean(ex_data$fev1pp_00))*10
heatmap_data$age_category[heatmap_data$age<60]<-"age<60"
heatmap_data$age_category[heatmap_data$age>=60&heatmap_data$age<70]<-"60<age<70"
heatmap_data$age_category[heatmap_data$age>=70]<-"age>70"
heatmap_data$FEV1_category[heatmap_data$fev1<30]<-"FEV1<30"
heatmap_data$FEV1_category[heatmap_data$fev1>=30 & heatmap_data$fev1<50]<-"30<FEV1<50"
heatmap_data$FEV1_category[heatmap_data$fev1>=50 & heatmap_data$fev1<80]<-"FEV1>50"
heatmap_data$gender_1<-ifelse(heatmap_data$gender==1,"Male","Female")
heatmap_data$predicted_RMST<-as.vector(Mean_tau_T_est_MI)

#heatmap_list<-cbind(mean_RMST,as.vector(Mean_tau_T_est_MI))
heatmap_list<-as.matrix(Mean_tau_T_est_MI)
rownames(heatmap_list) <- seq(1,1111,1)
annotation_data<-data.frame(Treatment=ifelse(Data_cens$trt==1,"Azithromycin","Placebo"),
                            Age=heatmap_data$age_category,
                            FEV1=heatmap_data$FEV1_category,
                            Gender=heatmap_data$gender_1,
                            "Still_smoking"=ifelse(Data_cens$nowsmk,"Yes","No")
)
rownames(annotation_data)<-rownames(heatmap_list)
#library(RColorBrewer)
mycolors <- c("deepskyblue","deeppink")
#mycolors <- c("deepskyblue","darkkhaki","cornsilk","coral","palegreen1","deeppink","green","yellow")
names(mycolors) <- unique(annotation_data$Treatment)
mycolors1 <- c("cyan1","blue4","coral")
names(mycolors1) <- unique(annotation_data$Age)
mycolors2<-c("plum3","purple","orange")
names(mycolors2) <- unique(annotation_data$FEV1)
#mycolors3 <- list(Trt = mycolors,Age=mycolors1,FEV1=mycolors2)
mycolors3<-c("yellow","blue")
names(mycolors3) <- unique(annotation_data$Gender)
mycolors4<-c("orchid1","royalblue1")
names(mycolors4)<-unique(annotation_data$Still_smoking)
mycolors5 <- list(Trt = mycolors,Age=mycolors1,FEV1=mycolors2,Gender=mycolors3,Still_smoking=mycolors4)

mat_breaks <- seq(min(heatmap_list,na.rm = TRUE), max(heatmap_list,na.rm = TRUE), length.out = 100)
#library(viridis)
png("Results/heatmap_fig3.png",    # create PNG for the heat map
    width = 2.2*500,        # 5 x 300 pixels
    height = 3*500,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
par(oma=c(0,10,0,10))
#library(pheatmap)
pheatmap(heatmap_list,breaks= mat_breaks,
         col=plasma(100),
         annotation_row=annotation_data,
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames=FALSE,
         fontsize=6,
         angle_col = "45",
         annotation_colors = mycolors5)
