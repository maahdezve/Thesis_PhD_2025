#########################.
# Libraries ----
#########################.

library(dplyr)
library(nlme)
library(splines)
library(cowplot)
library(dplyr)
library(ggplot2)
library(saemix)
library(VCA)

#########################.
# functions  ----
#########################.

source("SPNLME_models/1.SPNLME_help_functions.R")

#########################.
# Dataset-----
#########################.


Data <- readRDS(file = "SPNLME_models/Data/data_simul_mix.rds")
names(Data$data)

sim_data=Data$data %>%  filter(n_j==20)

# Growth curves
ggplot(data = sim_data, aes(x = x, y = y, col = group, group = id)) +
  geom_line(aes(y=f),show.legend = F)+
  geom_point(show.legend = F)+ 
  theme_bw()

# Derivatives
ggplot(data = sim_data, aes(x = x, y = df, col = group, group = id)) + 
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
# (II) Fitting Models with semiparametric methods by groups--------------
#---------------------------------------------------------------------------------.

x <- sim_data$x
y <- sim_data$y
id <- sim_data$id

x_min_v=  min(x)
x_max_v=  max(x)

################## Select the SPNLME model

model_sel <- 3  #  1=PB, 2=SP1, 3=SP2,  4=SP3

#####################################################3
# g0
######################################################

data_sel_g0= sim_data%>% filter(group == "0", type=="Training" )
data_sel_g1= sim_data%>% filter(group == "1", type=="Training"  )

fm_g0 <-semiparametric_model(
  x = data_sel_g0$x, y = data_sel_g0$y, id = data_sel_g0$id,
  pord = 2, bdeg = 3, start = c(a_i = 175, b_i = 13.5, c_i = 162, s0 = 0.1, s1 = 1), 
  ndx = nseg, PB_model = model_sel)

fm_g0$model.out

#####################################################3
# g1
######################################################

fm_g1 <-semiparametric_model(
  x = data_sel_g1$x, y = data_sel_g1$y, id = data_sel_g1$id, 
  pord = 2, bdeg = 3, start = c(a_i =176 , b_i = 14.2, c_i = 164, s0 = 0.1, s1 = 0.7), ndx = nseg, 
  PB_model = model_sel)

fm_g1$model.out

# #---------------------------------------------------------------------------------.
# # (III) Compute lassification rule --------------
# #---------------------------------------------------------------------------------.
# 

pi0 <- n_distinct(sim_data %>% filter(group == "0") %>% pull(id)) / (n_distinct(sim_data %>% pull(id)))

########### Data from g0 group with parameters g0 group

dens_data <- pi.is_cond(data = fm_g0$data_sp, model = fm_g0$model.out, 
                               spline_struct = model_sel, ndx = nseg, pord = pord, 
                               bdeg = bdeg, v=15, x_min=x_min_v,x_max= x_max_v)

########### Data from g0 group with parameters g1 group

dens_alt <- pi.is_cond(data = fm_g0$data_sp, model = fm_g1$model.out, 
                              spline_struct = model_sel, ndx = nseg, pord = pord, 
                              bdeg = bdeg, v=15, x_min=x_min_v,x_max= x_max_v)

#------------------------------------------------
########### Bayes rule for data from g0 group
#------------------------------------------------

class_g0=data.frame(p_clas = pi0 * dens_data / (pi0 * dens_data + (1 - pi0) * dens_alt), 
                    id = unique( fm_g0$data_sp$id), 
                    AIC = AIC(fm_g0$model.out), BIC = BIC(fm_g0$model.out))%>% 
  mutate(pred_class = ifelse(p_clas >= 0.5, 0, 1),nj=nj, ndx=nseg, PB_model=model_sel, orig_group="0")


class_g0_id=class_g0 %>% filter(pred_class==orig_group) %>% pull(id)


########### Data from g1 group with parameters g1 group
dens_data <- pi.is_cond(data = fm_g1$data_sp, model = fm_g1$model.out, 
                               spline_struct = model_sel,  ndx = nseg, pord = pord, 
                               bdeg = bdeg, v=15, x_min=x_min_v,x_max= x_max_v)
########### Data from g1 group with parameters g0 group
dens_alt <- pi.is_cond(data = fm_g1$data_sp, model = fm_g0$model.out, 
                              spline_struct = model_sel, ndx = nseg, pord = pord, 
                              bdeg = bdeg, v=15, x_min=x_min_v,x_max= x_max_v)

#------------------------------------------------
########### Bayes rule for data from g1 group
#------------------------------------------------

class_g1=data.frame(p_clas = (1-pi0)* dens_data / ((1-pi0) * dens_data + (pi0) * dens_alt),
                    id = unique( fm_g1$data_sp$id), 
                    AIC = AIC(fm_g1$model.out), BIC = BIC(fm_g1$model.out))%>% 
  mutate(pred_class = ifelse(p_clas >= 0.5, 1, 0), nj=nj, ndx=nseg, PB_model=model_sel, orig_group="1")

class_g1_id=class_g1 %>% filter(pred_class==orig_group) %>% pull(id)

# #---------------------------------------------------------------------------------.
# # (IV) confussion matrix  --------------
# #---------------------------------------------------------------------------------.
# 

conf_class = bind_rows(class_g0, class_g1) %>% arrange(id, nj) %>% 
  mutate(PB_model=as.factor(PB_model),
         PB_model_lab= as.factor(ifelse(PB_model==1,"PB",
                                        ifelse(PB_model==2,"SP1", 
                                               ifelse(PB_model==3,"SP2", "SP3")))))

(conf_mat=table(conf_class$orig_group,conf_class$pred_class))


print(paste0("Clasification in g0 group: ", diag(conf_mat)[1] , " individuals"))
print(paste0("Clasification in g1 group: ", diag(conf_mat)[2], " individuals"))



# #---------------------------------------------------------------------------------.
# # (V) plot  --------------
# #---------------------------------------------------------------------------------.
# 
ggplot(data= data_sel_g0,aes(x=x, y=y, group=id))+ geom_line(col="grey")+
  geom_line(data = data_sel_g0 %>% filter(!id %in% class_g0_id),aes(group=id), col="red")+
  labs(subtitle = "results_g0")+theme_bw()


ggplot(data= data_sel_g1,aes(x=x, y=y, group=id))+ geom_line(col="grey")+
  geom_line(data = data_sel_g1 %>% filter(!id %in% class_g1_id),aes(group=id), col="red")+
  labs(subtitle = "results_g1")+theme_bw()



