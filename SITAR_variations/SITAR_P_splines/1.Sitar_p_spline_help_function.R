########################
# libraries
########################

library(splines)
library(glue)


############################
# Bspline
############################

predict.bbase <- function(knots, bdeg=3, newx, deriv = 0, intercept = FALSE) {
  B <- spline.des(knots = knots, x = newx, derivs = deriv,
                  ord = bdeg + 1, outer.ok = TRUE)$design
  if(!intercept) {
    B <- B[,-1] # Remove intercept
  }
  B
} 

#########################################
# peak function
#########################################

peak <-function (x) {
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 1
  if (length(pks)!=0) {
    pks[which.max(x[pks])]
  }
  else{
    which.max(x)
  } 
}  

#########################################
# apv/phv function
#########################################

phv_estimates<-function(x_grid, id, y_grid, data.grid=NULL) {
  phv_info<-as.data.frame(matrix(NA, ncol=2, nrow=length(unique(id))))
  
  data.pred=data.frame(x_grid=x_grid, id=id, y_grid=y_grid)
  
  id_sel=unique(id)
  for(i in 1:n_distinct(data.pred$id)){
    age_int<-data.pred$x_grid[data.pred$id==id_sel[i]] 
    phv_info[i,1] = spline(age_int,data.pred$y_grid[data.pred$id==id_sel[i]], method="natural")$x[peak(spline(age_int,data.pred$y_grid[data.pred$id==id_sel[i]], method="natural")$y)] # aphv
    phv_info [i,2]=spline(age_int,data.pred$y_grid[data.pred$id==id_sel[i]], method="natural")$y[peak(spline(age_int,data.pred$y_grid[data.pred$id==id_sel[i]], method="natural")$y)] #phv
    
  }
  
  colnames(phv_info)<-c("apv", "pv")
  phv_info$id=unique(data.pred$id)
  return(phv_info)
}



# bs_f <- function(x, knots, bdeg, intercept = FALSE) {
#   res <- spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
#   if(!intercept) {
#     res <- res[,-1] # Remove intercept
#   }
#   res
# }


##############################################
# B-splines sitar
##############################################
sitar_B_spline=function(data, df,vi_b=0.1,vi_c=0.1,bdeg =3, deriv=F,
                        control = nlmeControl(msMaxIter = 100, returnObject = TRUE) ){

bs_f <- function(x, knots, bdeg, intercept = FALSE) {
  res <- spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
  if(!intercept) {
    res <- res[,-1] # Remove intercept
  }
  
  res
}

x  <- data$age
id <- data$id
y  <- data$height

xl <- min(x)
xr <- max(x)

nseg <- df

dx <- (xr - xl) / nseg

knots_aux <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)

#----------------------------
# data transformation
#----------------------------

xoffset <- mean(x)
bstart <- mean(x) - xoffset
x <- x - xoffset
knots <- knots_aux - xoffset

fulldata <- data.frame(x, y, id)


##----------------------------
# funtion lm: initial values
##----------------------------

spline.lm <- lm(y ~ bs_f(x, knots = knots, bdeg = bdeg))
start_p1 <- coef(spline.lm)[c(2:(nseg+bdeg), 1)]              ### start values for ns parameters and a
start_p1 <- c(start_p1, bstart+vi_b,vi_c)

##----------------------------
## sitar_nlme_BS
##----------------------------

## define parameters
ss<-paste0("s",1:(nseg+bdeg-1))
fixed<- c(ss,"a","b","c")
pars<-paste(c("x",fixed), collapse = ' , ')

# 1. define spline
Bspline_prod<-glue('rowSums((cbind(', paste(ss, collapse = ' , '),')*bs_f((x-(b))*exp(c), knots = knots, bdeg = bdeg)))')

#2. define fixed form
fixed<-paste(fixed, collapse = ' + ')
fixed <- glue('{fixed} ~ 1')

#3. define fitnlme_p function and apply nlme
fitnlme_p1<-glue("fitenv <- new.env()\n",
                 "fitenv$fitnlme_p1 <- function(<<pars>>) {\n",
                 "a+<<Bspline_prod>>}\n", "on.exit(detach(fitenv))\n",
                 "attach(fitenv)\n", "nlme(y ~ fitnlme_p1(<<pars>>),\n",
                 "fixed = <<fixed>>,\n",
                 "random = a + b + c ~ 1 | id, \n",
                 "control = control, \n",
                 "data = fulldata,  start = start_p1)",
                 .open = "<<",  .close = ">>"
)


time_saem_1=proc.time()[3]
m1_p1 <- try(eval(parse(text = fitnlme_p1)), silent = T)
time_models = proc.time()[3]-time_saem_1

if(deriv){

  res_der_bs=list()
  
  for( i in 1:nind){
    data_sel =data %>% filter(id == i )
    new_x=seq(min(data_sel$age), max(data_sel$age), l=100)
  
    b_i= m1_p1$coefficients$fixed[c("b")]+m1_p1$coefficients$random$id[i,c("b")]
    c_i=m1_p1$coefficients$fixed[c("c")]+m1_p1$coefficients$random$id[i,c("c")]
    
    new_x_trans = (new_x-xoffset-b_i)/ exp(-c_i)
    
    B_der=predict.bbase(knots = knots ,newx=new_x_trans, deriv = 1)
    
    spline_coef=m1_p1$coefficients$fixed[c(1: (length(m1_p1$coefficients$fixed) -3))]
    
    Der = data.frame(der_bs= B_der%*%spline_coef, id=i,age=new_x)
    
    res_der_bs[[i]]= Der %>%  mutate(der_bs= der_bs*(1)/exp(-c_i))
    
  }
  
  res_der_bs = do.call("rbind", res_der_bs)
  
  peak_val_bs = phv_estimates(res_der_bs$age, res_der_bs$id, res_der_bs$der_bs)
  
  
  return(list(model=m1_p1,deriv=res_der_bs, PHV=peak_val_bs))
  
} else(return(m1_p1))

}

##############################################
# P-splines sitar (transformation)
##############################################

sitar_P_spline=function(data, df,vi_b=0.1,vi_c=0.1,  bdeg =3,  pord = 2, deriv=F,
                        control = nlmeControl(msMaxIter = 100, returnObject = TRUE) ){
  
bs_f <- function(x, knots, bdeg, intercept = FALSE) {
  res <- spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
  if(!intercept) {
    res <- res[,-1] # Remove intercept
  }
  
  res
}

x  <- data$age
id <- data$id
y  <- data$height

 nseg <- df

xl <- min(x)
xr <- max(x)

dx <- (xr - xl) / nseg

knots_aux <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)

#----------------------------
# data transformation
##----------------------------


xoffset <- mean(x)
bstart <- mean(x) - xoffset
x <- x - xoffset
knots <- knots_aux - xoffset

ID <- factor(rep(1, length(x)))
fulldata <- data.frame(x, y, id, ID)

##----------------------------
# B-spline transformation
##----------------------------

pord <- 2
m <- nseg + bdeg
D <- diff(diag(m), differences = pord)
P.svd <- svd(crossprod(D))
Ts <- (P.svd$u)[,1:(m-pord)]%*%diag(sqrt(1/(P.svd$d)[1:(m-pord)]))

# Wood
B <- bs_f(x, knots = knots, bdeg = bdeg, intercept = TRUE)
Z <- B%*%Ts
X <- B%*%((P.svd$u)[,-(1:(m-pord))])
D.temp <- sweep(X, 2, colMeans(X))
Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
X <- X%*%Xf
Tn <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf


#----------------------------
# funtion lm: initial values
##----------------------------

B <- bs_f(x, knots = knots, bdeg = bdeg, intercept = TRUE)
aux1 <- (B%*%Tn)[,-1] # Exclude intercept
aux2 <- B%*%Ts
spline.lm <- lm(y ~ aux1)
start_p3 <- coef(spline.lm)[c(2, 1)]              ### start values for ns parameters and a
start_p3 <- c(start_p3, bstart+vi_b,vi_c)#+0.01, 0.01)
names(start_p3)=c(paste0("s",1:1),"a","b","c")

##----------------------------
## sitar_nlme
##----------------------------

## define parameters
ss<-paste0("s",1:(nseg+bdeg-1))
fixed<- c(ss[1],"a","b","c")
pars<-paste(c("x",ss,fixed[-1]), collapse = ' , ')
random = glue('list(ID = pdIdent(', paste(ss[2:(nseg+bdeg-1)], collapse = ' + '),' ~ 1),\n',
              'id = pdSymm(a + b + c ~ 1))')

# 1. define spline

Transf = glue(
  "B <- bs_f((x-(b))*exp(c), knots = knots, bdeg = bdeg, intercept = TRUE) \n",
  "aux1 <- (B%*%Tn)[,-1, drop = FALSE] \n",
  "aux2 <- B%*%Ts \n")

spline_trans<-glue('rowSums((cbind(', paste(ss[1], collapse = ' , '),')*aux1)) +  rowSums((cbind(', paste(ss[2:(nseg+bdeg-1)], collapse = ' , '),')*aux2))')

#2. define fixed form
fixed<-paste(fixed, collapse = ' + ')
fixed <- glue('{fixed} ~ 1')

#3. define fitnlme_p function and apply nlme
fitnlme_p3<-glue("fitenv <- new.env()\n",
                 "fitenv$fitnlme_p3 <- function(<<pars>>) {\n",
                 "<<Transf>>\n",
                 "a+<<spline_trans>>}\n", "on.exit(detach(fitenv))\n",
                 "attach(fitenv)\n", "nlme(y ~ fitnlme_p3(<<pars>>),\n",
                 "fixed = <<fixed>>,\n",
                 "random = <<random>>, \n",
                 "control = control, \n",
                 "data = fulldata,  start = start_p3)",
                 .open = "<<",  .close = ">>"
)


time_saem_3=proc.time()[3]
m1_p3<- try(eval(parse(text = fitnlme_p3)), silent = T)
time_models = proc.time()[3]-time_saem_3
if(deriv){

  m <- df + bdeg
  D <- diff(diag(m), differences = pord)
  P.svd <- svd(crossprod(D))
  Ts <- (P.svd$u)[,1:(m-pord)]%*%diag(sqrt(1/(P.svd$d)[1:(m-pord)]))
  
  # Wood
  B <- bs_f(x, knots = knots, bdeg = bdeg, intercept = TRUE)
  Z <- B%*%Ts
  X <- B%*%((P.svd$u)[,-(1:(m-pord))])
  D.temp <- sweep(X, 2, colMeans(X))
  Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
  X <- X%*%Xf
  Tn <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
  
  
  res_der_ps=list()
  
  for( i in 1:nind){
    data_sel = data %>% filter(id == i )
    new_x=seq(min(data_sel$age), max(data_sel$age), l=100)
    
    b_i= m1_p3$coefficients$fixed[c("b")]+m1_p3$coefficients$random$id[i,c("b")]
    c_i=m1_p3$coefficients$fixed[c("c")]+m1_p3$coefficients$random$id[i,c("c")]
    
    new_x_trans = (new_x-xoffset-b_i)/ exp(-c_i)
    
    B_der=predict.bbase(knots = knots ,newx=new_x_trans, deriv = 1, intercept = TRUE)
    
    spline_coef_x=m1_p3$coefficients$fixed
    spline_coef_z=m1_p3$coefficients$random$ID
    
    X_der=(B_der%*%Tn)[,-1, drop = FALSE]
    Z_der=B_der%*%Ts
    
    Der = data.frame(der_ps= X_der%*%spline_coef_x["s1"] + Z_der%*%c(spline_coef_z) , id=i,age=new_x)
    
    res_der_ps[[i]]= Der %>%  mutate(der_ps= der_ps*(1)/exp(-c_i))
    
  }
  
  res_der_ps = do.call("rbind", res_der_ps)
  
  
  peak_val_ps = phv_estimates(res_der_ps$age, res_der_ps$id, res_der_ps$der_ps)

  
  return(list(model=m1_p3,deriv=res_der_ps, PHV=peak_val_ps))
  
} else(return(m1_p3))

}
