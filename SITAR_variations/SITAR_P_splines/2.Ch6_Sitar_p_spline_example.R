#####################
# Library
#####################

library(nlme)
library(sitar)
library(splines)
library(tidyverse)
library(glue)
library(cowplot)

source("SITAR_variations/SITAR_P_splines/1.Sitar_p_spline_help_function.R")

##########################
# Data
##########################
data=readRDS("SITAR_variations/SITAR_P_splines/Data/sim_data.rds")

data=data$data

data <- data %>%rename(age=x,height=y,id=id) %>%  arrange(id,age) %>% mutate(id = as.factor(id))
nind = n_distinct(data$id)

#####################
# Fit sitar models
#####################

m1 <- sitar(x = age, y = height, id = id, data = data, df = 5)
m1_ops=dfpower(m1, df = 4:15, fixed = c('a+b+c'),
        xpowers = 1, ypowers = 1, maxIter = 8,FUN=BIC)

#----------------- classic sitar ncs
sitar_ncs = sitar(x = age, y = height, id = id, data = data, df = 4)

#-----------------  sitar B-spline
sitar_bs = sitar_B_spline(data=data, df =4,deriv =T)
  
#----------------- sitar P-spline
sitar_ps = sitar_P_spline(data=data, df =4,deriv=T)


#----------------- BIC
BIC(sitar_ncs)
BIC(sitar_bs$model)
BIC(sitar_ps$model)

#################################################
# fitted curves Plots
#################################################

p_sp=ggplot(data = data, aes(x = age+mean(age),y = height), show.legend = FALSE) + 
  geom_point(size = 1, shape = 21, fill = "white", show.legend = FALSE) +
  geom_line(aes(y = fitted(sitar_ps$model), col = id, alpha = 0.5, group = id), show.legend = FALSE)+
  geom_line(aes(y = predict(sitar_ps$model, level = 1)), colour = "black", linetype = "dashed", linewidth = 1, show.legend = FALSE)+
  theme_bw() +  labs(x="age")+
  ggtitle("P-splines")


sitar=ggplot(data = data, aes(x = age+mean(age),y = height), show.legend = FALSE) + 
  geom_point(size = 1, shape = 21, fill = "white", show.legend = FALSE) +
  geom_line(aes(y = fitted(sitar_ncs), col = id, alpha = 0.5, group = id), show.legend = FALSE)+
  geom_line(aes(y = predict(sitar_ncs, level = 0)), colour = "black", linetype = "dashed", linewidth = 1,show.legend = FALSE)+
  theme_bw() +  labs(x="age")+
  ggtitle("Classic sitar")


B_sp=ggplot(data = data, aes(x = age+mean(age),y = height), show.legend = FALSE) + 
  geom_point(size = 1, shape = 21, fill = "white", show.legend = FALSE) +
  geom_line(aes(y = fitted(sitar_bs$model), col = id, alpha = 0.5, group = id), show.legend = FALSE)+
  geom_line(aes(y = predict(sitar_bs$model, level = 0)), colour = "black", linetype = "dashed", linewidth = 1, show.legend = FALSE)+
  theme_bw() +  labs(x="age")+
  ggtitle("B-splines")


plot_grid(sitar,B_sp,p_sp, ncol=3)


############################################3
### Derivative plot
###############################################
plot(sitar_ncs, opt = 'V', apv=T,col = 4, main="sitar ncs")

sitar_bs$deriv %>%  ggplot(aes(x=age, y=der_bs, group=id, col=as.factor(id)))+
  geom_line()+ labs(x="Age",y="Height derivative (cm/years)")+geom_hline(yintercept = 0)+
  geom_vline(data=sitar_bs$PHV,aes(xintercept=apv, col=as.factor(id)), lty=2, alpha=0.4)+
  labs(subtitle = "sitar B-spline")+
  theme_bw()+
  theme(legend.position = "none")

sitar_ps$deriv  %>%  ggplot(aes(x=age, y=der_ps, group=id, col=as.factor(id)))+geom_line()+
  geom_line()+ labs(x="Age",y="Height derivative (cm/years)")+geom_hline(yintercept = 0)+
  geom_vline(data=sitar_ps$PHV,aes(xintercept=apv, col=as.factor(id)), lty=2, alpha=0.5)+
  labs(subtitle = "sitar P-spline")+ 
  theme_bw()+
   theme(legend.position = "none")

 