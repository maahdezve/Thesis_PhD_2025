#########################.
# Libraries ----
#########################.

library(dplyr)
library(nlme)
library(splines)
library(cowplot)
library(dplyr)
library(ggplot2)

#########################.
# load functions  ----
#########################.
source("SPNLME_models/1.SPNLME_help_functions.R")

#########################.
# Dataset-----
#########################.

# --------------------- Preece-Baines growth dataset -------------

Data <- readRDS(file = "SPNLME_models/Data/sim_data.rds")
names(Data$data)

sim_data =Data$data

# Growth curves
ggplot(data = Data$data.grid, aes(x = x, y = f, col = as.factor(id), group = id)) +
  geom_line(show.legend = F)+ theme_bw()

# Derivatives
ggplot(data =  Data$data.grid, aes(x = x, y = df, col = as.factor(id), group = id)) + 
  geom_line(show.legend = F)+ theme_bw()

#########################.
# Apply Method-----
#########################.
#---------------------------------------------------------------------------------.
#  (I) P-spline features 
#---------------------------------------------------------------------------------.
nj=round(sim_data %>% count(id) %>% summarise(nj=mean(n)) %>% pull(nj))
nseg <-  ceiling(nj/2)     # Number of segments
bdeg <- 3     # polynomial degree q
pord <- 2    # Penalty order p

#---------------------------------------------------------------------------------.
# (II) Fitting Subject specific curves Models with semiparametric methods--------------
#---------------------------------------------------------------------------------.

x <- sim_data$x
y <- sim_data$y
id <- sim_data$id

# PB
fm_pb<- semiparametric_model(x=x, y=y, id=id,pord=2,bdeg=3, 
                    start =c(a_i = 170, b_i = 12.5, c_i = 162, s0 = 0.2, s1 = 1),
                    PB_model = 1)

start.pb=fm_pb$model.out$coefficients$fixed

# SP1
fm_SP1 <-semiparametric_model(x=x, y=y, id=id, pord=2,bdeg=3, 
                    start=start.pb ,ndx=nseg, PB_model = 2)

#SP2

fm_SP2<-semiparametric_model(x=x, y=y, id=id,pord=2,bdeg=3, 
                            start=start.pb ,ndx=nseg,  PB_model = 3)
#SP3

fm_SP3<- semiparametric_model(x=x, y=y, id=id, pord=2,bdeg=3, 
                     start = start.pb,ndx=nseg,  PB_model = 4)

# #---------------------------------------------------------------------------------.
# # (III) Estimate the  derivative, PHV and APHV --------------
# #---------------------------------------------------------------------------------.
# 

deriv_PB= Deriv_sp(model = fm_pb,x.grid=seq(min(sim_data$x), max(sim_data$x), length.out = 100))
deriv_SP1= Deriv_sp(model = fm_SP1,x.grid=seq(min(sim_data$x), max(sim_data$x), length.out = 100))
deriv_SP2= Deriv_sp(model = fm_SP2, x.grid=seq(min(sim_data$x), max(sim_data$x), length.out = 100))
deriv_SP3= Deriv_sp(model = fm_SP3, x.grid=seq(min(sim_data$x), max(sim_data$x), length.out = 100))


deriv=bind_rows(list(deriv_PB$res,deriv_SP1$res,deriv_SP2$res,deriv_SP3$res))
deriv_apv=bind_rows(list(deriv_PB$aphv_est,deriv_SP1$aphv_est,deriv_SP2$aphv_est,deriv_SP3$aphv_est))

# #---------------------------------------------------------------------------------.
# # (IV) plots --------------
# #---------------------------------------------------------------------------------.
# 
sim_data$PB=(predict(fm_pb$model.out, level= 0:1))$predict.id 
sim_data$SP1=predict(fm_SP1$model.out, level= 0:2)$predict.id 
sim_data$SP2=predict(fm_SP2$model.out, level= 0:1)$predict.id 
sim_data$SP3=predict(fm_SP3$model.out, level= 0:2)$predict.id 

id_sel = c(1:3,84:86)

#-------------- (IV-a): Fitted function ------------
sim_data %>% filter(id%in%id_sel) %>% ggplot(aes(x=x, y=y, group=id)) + 
  geom_point() + 
  geom_line(aes(y=PB, col="PB", lty="PB")) +
  geom_line(aes(y=SP1, col="SP1", lty="SP1")) +
  geom_line(aes(y=SP2, col="SP2", lty="SP2")) +
  geom_line(aes(y=SP3, col="SP3", lty="SP3")) +
  labs(x = "Age(years)", y = "Height (cm)", colour = "Method", linetype="Method") +
  facet_wrap(~id, labeller = label_both)+
  theme_bw()

#-------------- (IV-b): first derivative estimation ------------

deriv%>% rename(id=id) %>%  filter(id%in%id_sel)  %>% 
  ggplot(aes(x=x,y = df_all)) + 
  geom_line(aes(group=interaction(PB_model,id),col=as.factor(PB_model)))+
  labs(x = "Age(years)", y = "Height Velocity (cm/yrs)", colour = "Method") +
  geom_vline(data=deriv_apv%>% filter(id%in%id_sel) ,aes(xintercept=apv, col=as.factor(str)))+
  facet_wrap(~id, labeller = label_both)+
  theme_bw()

#-------------- (IV-c): Parts semiparametric model  ------------


##------------------- parts 

parts_SP3=semiprametric_parts(model=fm_SP3, ndx=nseg) %>%
  dplyr::select(age,id,PB_all_f:fi_ale)%>% mutate(Model="SP3")


########## SP3

ffi_plot_fix=parts_SP3%>%  ggplot(aes(x=age, y=f_fixed, group=interaction(Model), col=Model))+geom_line(aes(lty=Model))+
  labs(x="Age", y=expression(f:X*beta+Zb))+theme_bw()+
  theme(legend.position = "none")+facet_wrap(~Model, labeller = label_both)+
  scale_color_manual(values = c("#C77CFF"))


PB_ffi_plot=parts_SP3 %>% ggplot(aes(x=age, y=PB_all_f, col=Model))+geom_line(aes(lty=Model,group=interaction(id,Model)))+
  geom_line(aes(y=PB_pop_f), col=1, lty=2)+
  labs( x="Age", y=expression(PB(x)))+theme_bw()+
  theme(legend.position = "none")+facet_wrap(~Model, labeller = label_both)+
  scale_color_manual(values = c("#C77CFF"))

ffi_plot_ale=parts_SP3%>% ggplot(aes(x=age, y=fi_ale, group=interaction(id, Model), col=Model))+geom_line(aes(lty=Model))+
  labs( x="Age", y=expression(fi:X*beta[i]+Zb[i]))+
  geom_hline(yintercept=0, col=1)+theme_bw()+
  theme(legend.position = "bottom")+facet_wrap(~Model, labeller = label_both)+
  scale_color_manual(values = c("#C77CFF"))

# when position = "bottom," legend is the third grob in the list
shared_legend_5 <- ggpubr::get_legend(ffi_plot_ale)
merge_5=plot_grid(PB_ffi_plot,ffi_plot_fix,ffi_plot_ale+theme(legend.position = "none"), ncol=3)
merge_5
#
