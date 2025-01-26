#Manuscript Title: “tau-Inflated Beta Regression Model for Estimating $\tau$-Restricted Means and Event-free Probabilities for Censored Time-to-event Data"
#Authors: Yizhuo Wang, Susan Murray
#Code was written by Yizhuo Wang.
#In case of questions or comments please contact yizhuow@umich.edu
library(foreach)
library(doParallel)
library(parallel)
library(boot)
library(purrr)
library(survival)
library(survminer)
library(data.table)
library(survRM2)
library(pseudo)
library(mvtnorm)
library(coda)
library(numDeriv)
library(matrixStats)
library(betareg)
library(stats)
library(simcausal)
library(latex2exp)
library(dplyr)
library(base)
library(eoffice)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(tikzDevice)
library(ggpubr)
library(grid)

options(warn=0) #keep running the code with warnings
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

source("Scripts/Table1andTable2.R") #takes few days to get full simulation results; could do dry-run for 100 subjects first
source("Scripts/Table3.R")
source("Scripts/Table4.R")
source("Scripts/Table5.R")
source("Scripts/fig2_biasplot.R")
dev.off() #run this line to open the saved plot in the folder
source("Scripts/heatmap_fig3.R")
dev.off() #run this line to open the saved plot in the folder
source("Scripts/Fig4.R")
dev.off() #run this line to open the saved plot in the folder
