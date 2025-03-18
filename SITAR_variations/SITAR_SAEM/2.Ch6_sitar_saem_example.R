#########################
# library
#########################

library(sitar)
library(saemix)
library(dplyr)
library(splines)
library(ggplot2)
library(MASS)

source("SITAR_variations/SITAR_SAEM/1.sitar_saem_help_functions.R")

##########################
# Data
##########################
data <- readRDS("SITAR_variations/SITAR_SAEM/Data/sim_data.rds")

data=data$data
data <- data %>%rename(age=x,height=y,id=id) %>%  arrange(id,age)

x <- data$age
y <- data$height
id <- as.factor(data$id)

# ##############################################
# # Fit sitar model LB
# ##############################################

### find the df based on BIC
sitar_model=sitar(df=6,x = age, y = height, id = id, data = data)
df_val=dfpower(sitar_model,df = 3:8, fixed = c( 'a+b+c'),
               xpowers = 1, ypowers = 1)
df_sel=as.numeric(names(df_val[df_val>0][which.min(df_val[df_val>0])]))


sitar_model=sitar(df=df_sel,x = age, y = height, id = id, data = data)

# ##############################################
# # Fit sitar model with SAEM
# ##############################################

df =  df_sel

fulldata <- data.frame(x=data$age-mean(data$age), y=data$height, id=data$id)

## saemix data type
pb.data <- saemixData(name.data = fulldata, header = TRUE,
                      sep = " ", na = NA, name.group = c("id"),
                      name.predictors = c("x"), name.response = c("y"),
                      units = list(x = "years", y = "cm"),
                      name.X = "x", verbose = F)
#--------------------------------------------------------------------------------------
# SAEM Model LM
saemix.model_lm=sitar_saemixModel(df=df,data=fulldata,transform.par=rep(0,df+3),
                                  fixed.estim=rep(1,df+3),
                                  covariance.model=matrix(c(rep(c(rep(1,3),rep(0,df)),3), 
                                                            rep(rep(0,df+3),df)),ncol=df+3,byrow=TRUE), verbose = F)


#####################
# stop_criteria
#####################

k1_lm=number_iter_k1(x, y, id, df, k1_initial= 200,seed = 42,  error_abs=0.01, error_rel=0.00085,start=attr(saemix.model_lm, "start"))
k1_sel_lm=k1_lm$k1_stop$k1_rel

#####################
# control values 
#####################
opt_lm <- list(seed = 123 ,save = FALSE, save.graphs = FALSE, print=FALSE,nbiter.saemix = c(k1_sel_lm, 200))
opt_default <- list(seed = 123 ,save = FALSE, save.graphs = FALSE, print=FALSE)

##############################
###### run sitar-saem
##############################

# Run SAEM for LM
start_time<-proc.time()[3]
saemix.fit_lm<-try(saemix(saemix.model_lm,pb.data,opt_default),silent = T)
end<-proc.time()[3]
time_saem.lm_defaul<- end-start_time

##################### SAEM TUN

# Run SAEM for LM
start_time<-proc.time()[3]
saemix.fit_lm_tun<-try(saemix(saemix.model_lm,pb.data,opt_lm),silent = T)
end<-proc.time()[3]
time_saem.lm<- end-start_time


results_saem= list(saemix.fit_lm=saemix.fit_lm,
                    saemix.fit_lm_tun=saemix.fit_lm_tun)

resuts_sitar_all = list(k1_lm=k1_lm,  results_saem=results_saem )
  
AIC_BIC= lapply(results_saem, function(x) try(x@results@bic.is,silent = T))


##########################BIC tab######################3
BIC_tab=bind_cols(lapply(AIC_BIC, bind_cols))
names(BIC_tab)=names(AIC_BIC)
BIC_tab$LB=BIC(sitar_model)
BIC_tab$df= df_sel

##########################model selected sitar y saem based in BIC ######################3

saem_model=results_saem$saemix.fit_lm
sitar_model = sitar_model

##########################Population parameter estimation######################

Fixed_effect_tab= data.frame( SAEM=saem_model@results@fixed.effects,
                              LB=c(sitar_model$coefficients$fixed[c((df+1):(df+3),1:(df))]), 
                              coef=c("a", "b", "c", paste0("s", 1:df)))

Variance_tab = data.frame(SAEM = c(diag(saem_model@results@omega[1:3,1:3]),(saem_model@results@respar[1])),
                          LB = as.numeric(VarCorr(sitar_model)[,1]),
                          coef= c("sigma_a^2", "sigma_b^2", "sigma_c^2", "sigma_eps^2"))  

Variance_tab[4,1] = Variance_tab[4,1] **2

tab_pop_param= rbind(Fixed_effect_tab,Variance_tab)

########### fit model ################################

fit_res=predict(sitar_model, level=0:1)
fit_res$age=data$age
fit_res$height=data$height
fit_res=fit_res %>% arrange(id, age)
fit_res$ppred=saem_model@results@ppred
fit_res$ipred=saem_model@results@ipred
fit_res=as.data.frame(fit_res)

###########################################
### plots
###########################################

saem_fit=fit_res %>% ggplot(aes(x=age, y=ipred, group=id))+geom_line()+
  labs(y = "height")+
  theme_bw()

LB_fit=fit_res %>% ggplot(aes(x=age, y=predict.id, group=id))+geom_line()+
  labs(y = "height")+
  theme_bw()

cowplot::plot_grid(LB_fit,saem_fit)

# Convergence
plot(saem_model,plot.type="convergence")

# Individual plot for subject 1, smoothed
plot(saem_model,plot.type="individual.fit",ilist=1:9,smooth=TRUE)

# Diagnostic plot: observations versus population predictions
par(mfrow=c(1,2))
plot(saem_model,plot.type="observations.vs.predictions",level=0:1,new=FALSE)

par(mfrow=c(1,1))
# Scatter plot of residuals

plot(saem_model,plot.type="residuals.scatter")
npde.obj<-npdeSaemix(saem_model)
plot(npde.obj)

######################################################################
############### Estimate the Derivative, PHV and aPHV ########################
#######################################################################

deriv_saem= deriv_sitar_saemixModel(object=saem_model, xoffset=mean(data$age) ,df=df_sel,
                                   n.smooth=100,smooth.deriv=T)

phv_res_saem=phv_estimates(x_grid=deriv_saem$deriv_smooth$x,
                           id=deriv_saem$deriv_smooth$id, 
                           y_grid=deriv_saem$deriv_smooth$y)

deriv_saem$deriv_smooth$LB=predict(sitar_model,  deriv=1, level=1,
                                   newdata=data.frame(age=deriv_saem$deriv_smooth$x,
                                                      id=deriv_saem$deriv_smooth$id))

phv_sitar=plot(sitar_model, opt="V", apv=T )
phv_res=merge(phv_sitar$apv,phv_res_saem, by="id", suffixes = c(".sitar",".saem")) 

#############################################################
############################# plot deriv
#############################################################

par(mfrow=c(1,1))

deriv_saem_plot=deriv_saem$deriv_smooth %>% ggplot(aes(x=x, y=y))+
  labs(x="Age(years)", y="Height Velocity (cm/yrs)", colour="Method")+
  geom_line(aes(group= id, col="SAEM"))+ 
  geom_vline(data=phv_res, aes(xintercept=apv.saem), lty=2, alpha=0.4, col="grey")+
  theme_bw()+ scale_color_manual(values ="#00BFC4")

deriv_LB_plot=deriv_saem$deriv_smooth %>% ggplot(aes(x=x, y=LB))+
  labs(x="Age(years)", y="Height Velocity (cm/yrs)", colour="Method")+
  geom_line(aes(group= id, col="LB"))+
  geom_vline(data=phv_res, aes(xintercept=apv.sitar), lty=2, alpha=0.4, col="grey")+
  theme_bw()+ scale_color_manual(values = "#F8766D")

cowplot::plot_grid(deriv_LB_plot,deriv_saem_plot)
