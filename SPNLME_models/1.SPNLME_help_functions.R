############################
# LIBRARIES
############################
library(dplyr)
library(nlme)
library(splines)
library(saemix)
library(glue)
library(VCA)

############################
# Bspline
############################
bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3, eps = 1e-15) {
  dx <- (xr - xl) / nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  # package spline
  B <- spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
  res <- list(B = B, knots = knots,bdeg=bdeg)
  res 
}
#-----------------------------------------------------------.
# B-spline basis (or its derivative) at new values ---------
#-----------------------------------------------------------.
predict.bbase <- function(object, newx, deriv = 0) {
  B <- spline.des(knots = object$knots, x = newx, derivs = deriv,
                  ord = object$bdeg + 1, outer.ok = TRUE)$design
  B
} 
################################
# transformation to mixed-model
################################
MM.basis <- function (x, xl, xr, nseg, bdeg, pord,param=1) {
  Bb <- bbase(x, xl, xr, nseg, bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  n <- nrow(B)
  D <- diff(diag(m), differences = pord)
  P.svd <- svd(crossprod(D))
  
  # Al multiplicar por diag(sqrt(1/(P.svd$d)[1:(m-pord)])) nos aseguramos que la matrix de var-cov es la diagonal
  Ts <- (P.svd$u)[,1:(m-pord)]%*%diag(sqrt(1/(P.svd$d)[1:(m-pord)]))
  Z <- B%*%Ts
  
  
  
  if (param == 1){
    X <- B%*%((P.svd$u)[,-(1:(m-pord))])
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    Tn <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
  } else    if (param == 2){ # Knots
    X <- NULL
    for(i in 0:(pord-1)){
      X <- cbind(X,x^i)
    }
    Tn <- NULL
    for(i in 0:(pord-1)){
      Tn <- cbind(Tn, knots[-c((1:(bdeg - 1)),
                               (length(knots)- (bdeg - 1) + 1):length(knots))]^i)
    }
    # Note that B%*%U.X == X
  } 
  
  list(X = X, Z = Z, Bb = Bb, m = m, D = D, knots = knots, Tn = Tn, Ts = Ts)
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

#########################################
# SPNLME main function
#########################################


semiparametric_model<-function(x,y,id,start,PB_model=1,pord=2,bdeg=3, ndx=4){
  
  ########################
  # x= vector covariate
  # y = vector obs
  # id=vector id
  # start = vector with start values with names
  # PB_model=1 (PB), 2 (PB_f), 3 (PB_fi), 4 (PB_ffi)
  ##########################
  
  # Transformation

  Matrix_trans <- MM.basis(x,min(x), max(x), nseg = ndx, bdeg = 3, pord = 2)
  
  X<- Matrix_trans$X
  Z<- Matrix_trans$Z
  
  DX <- dim(X)
  DZ <- dim(Z)
  
  n   <- length(y)
  Id2 <- factor(rep(1, length(y)))
  
  data_sp <-data.frame(data.frame(x,y,id),Id2,X,Z)
  colnames(data_sp)[5:(dim(data_sp)[2])] <- c(paste0("X",1:(DX[2])),paste0("Z",1:(DZ[2])))
  start_org=start
  
  
  #####################################################3
  # Fit Model
  ######################################################
  #
  ## define parameters
  X_i = paste0("X",2:DX[2])
  Z_i = paste0("Z",1:DZ[2])
  fixed<- c("a_i", "b_i", "c_i", "s0", "s1")
  
  
  pars<-paste(c("x",fixed), collapse = ',')
  
  fixed<-paste(fixed, collapse = '+')
  random<-paste(c("a_i", "b_i", "c_i"), collapse = '+')
  
  fixed <- glue('{fixed} ~ 1')
  random <- glue('{random} ~ 1| id')
  
  fit_model<-glue(
    "fitenv <- new.env()\n", "fitenv$SSpb <- function(<<pars>>) {\n",
    "a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i)) + exp(s1*(x-b_i)))}\n", "on.exit(detach(fitenv))\n", 
    "attach(fitenv)\n", "nlme(y ~ SSpb(<<pars>>),\n",
    "fixed = <<fixed>>,\n",
    "random =  <<random>>, \n",
    "data = data_sp,  start = start)\n",
    .open = "<<",  .close = ">>"
  )
  
  time <- proc.time()[3]
  model.out2 <- try(eval(parse(text = fit_model)), silent = T)
  time_t <- proc.time()[3]-time
  
  
  if(class(model.out2)[1]!="try-error") start=model.out2$coefficients$fixed else start=start
  
  if(PB_model==2){
    
    fixed<- c("a_i", "b_i", "c_i", "s0", "s1")
    pars<-paste(c("x","X2",fixed,"X_2","f_1"), collapse = ',')
    
    fixed<-paste(c(fixed, "X_2"), collapse = '+')
    
    random1<-glue('pdIdent(f_1 ~', paste(Z_i, collapse = '+'),' - 1)')
    random2<-glue('pdSymm(a_i + b_i + c_i ~ 1)')
    
    
    fixed <- glue('{fixed} ~ 1')
    
    fit_model<-glue(
      "fitenv <- new.env()\n", "fitenv$SSpb_f <- function(<<pars>>) {\n",
      "a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i)) + exp(s1*(x-b_i))) +  (X_2*X2+f_1)}\n", "on.exit(detach(fitenv))\n", 
      "attach(fitenv)\n", "nlme(y ~ SSpb_f(<<pars>>),\n",
      "fixed = <<fixed>>,\n",
      "random =  list(Id2=<<random1>>, \n",
      "id=<<random2>>), \n",
      "data = data_sp,  start = c(start,X_2=0),\n",
      "control = nlmeControl(msMaxIter = 1e8, msMaxEval = 1e8,MaxIter = 100,pnlsTol=0.008))",
      .open = "<<",  .close = ">>"
    )
    time <- proc.time()[3]
    model.out <- try(eval(parse(text = fit_model)), silent = T)
    time_t <- proc.time()[3]-time
    
    if(class(model.out)[1]=="try-error"){
      
      start=start_org
      time <- proc.time()[3]
      model.out <- try(eval(parse(text = fit_model)), silent = T)
      time_t <- proc.time()[3]-time
      
    } 
    
    
  }else if(PB_model==3){
    
    fixed<- c("a_i", "b_i", "c_i", "s0", "s1")
    pars<-paste(c("x",fixed,"f_i"), collapse = ',')
    
    fixed<-paste(c(fixed), collapse = '+')
    random1<-glue('pdSymm(a_i + b_i + c_i ~ 1)')
    random2<-glue('pdSymm(f_i ~X2-1)')
    random3<-glue('pdIdent(f_i ~', paste(Z_i, collapse = '+'),' - 1)')
    
    
    fixed <- glue('{fixed} ~ 1')
    
    fit_model<-glue(
      "fitenv <- new.env()\n", "fitenv$SSpb_fi <- function(<<pars>>) {\n",
      "a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i)) + exp(s1*(x-b_i))) +  f_i}\n","on.exit(detach(fitenv))\n", 
      "attach(fitenv)\n", "nlme(y ~ SSpb_fi(<<pars>>),\n",
      "fixed = <<fixed>>,\n",
      "random =  list( id = pdBlocked(list(<<random1>>,<<random2>>,<<random3>>))), \n",
      "data = data_sp,  start = start,\n",
      "control = nlmeControl(msMaxIter = 1e8, msMaxEval = 1e8,MaxIter = 100,pnlsTol=0.008))",
      .open = "<<",  .close = ">>"
    )
    
    time <- proc.time()[3]
    model.out <- try(eval(parse(text = fit_model)), silent = T)
    time_t <- proc.time()[3]-time
    
    if(class(model.out)[1]=="try-error"){
      
      start=start_org
      time <- proc.time()[3]
      model.out <- try(eval(parse(text = fit_model)), silent = T)
      time_t <- proc.time()[3]-time
      
    } 
    
  }else if(PB_model==4){
    fixed<- c("a_i", "b_i", "c_i", "s0", "s1")
    pars<-paste(c("x","X2",fixed,"X_2","f_1","f_i"), collapse = ',')
    
    fixed<-paste(c(fixed,"X_2"), collapse = '+')
    
    random0<-glue('pdIdent(f_1 ~', paste(Z_i, collapse = '+'),' - 1)')
    random1<-glue('pdSymm(a_i + b_i + c_i ~ 1)')
    random2<-glue('pdSymm(f_i ~X2-1)')
    random3<-glue('pdIdent(f_i ~', paste(Z_i, collapse = '+'),' - 1)')
    
    
    
    fixed <- glue('{fixed} ~ 1')
    random <- glue('{random} ~ 1| id')
    
    fit_model<-glue(
      "fitenv <- new.env()\n", "fitenv$SSpb_ffi <- function(<<pars>>) {\n",
      "a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i)) + exp(s1*(x-b_i))) + (X_2*X2+f_1) + f_i}\n","on.exit(detach(fitenv))\n", 
      "attach(fitenv)\n", "nlme(y ~ SSpb_ffi(<<pars>>),\n",
      "fixed = <<fixed>>,\n",
      "random =  list( Id2 = <<random0>>,\n",
      "id = pdBlocked(list(<<random1>>,<<random2>>,<<random3>>))), \n",
      "data = data_sp,  start =  c(start,X_2=0),\n",
      "control = nlmeControl(msMaxIter = 1e8, msMaxEval = 1e8,MaxIter = 100,pnlsTol=0.008))",
      .open = "<<",  .close = ">>"
    )
    
    time <- proc.time()[3]
    model.out <- try(eval(parse(text = fit_model)), silent = T)
    time_t <- proc.time()[3]-time
    
    
    if(class(model.out)[1]=="try-error"){
      
      start=start_org
      time <- proc.time()[3]
      model.out <- try(eval(parse(text = fit_model)), silent = T)
      time_t <- proc.time()[3]-time
      
    } 

  }
  
  
  if(PB_model==1)model.out=model.out2
  return(list(model.out=model.out, data_sp=data_sp, PB_model=PB_model,Matrix_trans=Matrix_trans, time=time_t))
  
}

#########################################
# SPNLME deriv pred function
#########################################


Deriv_sp=function(model, x.grid=NULL ){
  
  if(class(model$model.out)[1]!="try-error"){
     data_sp=model$data_sp 
   
    deriv_pb=list()
    
    PB1 = expression(a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i))+exp(s1*(x-b_i))))
    dPB1<-D(PB1,"x")
    
    PB1_pop = expression(a_pop-(2*(a_pop-c_pop))/(exp(s0*(x-b_pop))+exp(s1*(x-b_pop))))
    dPB1_pop<-D(PB1_pop,"x")
    
    case_sel=unique(data_sp$id)
    
   if( is.null(x.grid)) {
      grid_age = rep(seq(min(data_sp$x), max(data_sp$x),length.out=100), times=n_distinct(data_sp$id))
    } else{
      grid_age = rep(x.grid, times=n_distinct(data_sp$id))
      
    }
    
    age_matrix=data.frame(id=rep(unique(data_sp$id), each=n_distinct(grid_age)),x=grid_age)
    
    for( i in 1:n_distinct(data_sp$id)){
      
      x=age_matrix%>% filter(id==case_sel[i]) %>% pull(x)
      a_pop = model$model.out$coefficients$fixed[1]
      b_pop = model$model.out$coefficients$fixed[2]
      c_pop = model$model.out$coefficients$fixed[3]
      
      a_i =  a_pop + model$model.out$coefficients$random$id[i,1]
      b_i = b_pop + model$model.out$coefficients$random$id[i,2]
      c_i = c_pop + model$model.out$coefficients$random$id[i,3]
      s0 = model$model.out$coefficients$fixed[4]
      s1 = model$model.out$coefficients$fixed[5]
      
      
      deriv_pb[[i]] = data.frame(id=case_sel[i],x=x,fit_df_pb=eval(dPB1),fit_df_pop_pb=eval(dPB1_pop))
    }
    
    deriv_pb=do.call("rbind",deriv_pb)
    
    if(model$PB_model!=1){
      
      
      Bp_deriv <- predict.bbase(object=model$Matrix_trans$Bb, newx = grid_age, deriv = 1)
      
      DX=dim(model$Matrix_trans$X)
      DZ=dim(model$Matrix_trans$Z)
      
      X_der_prim = Bp_deriv%*%model$Matrix_trans$Tn
      Z_der_prim = Bp_deriv%*%model$Matrix_trans$Ts
      
      if(model$PB_model %in% c(2,4))  df_fixed = X_der_prim[,2]*model$model.out$coefficients$fixed[6] + as.matrix(Z_der_prim)%*%c(model$model.out$coefficients$random$Id2)
      if(model$PB_model%in% c(3,4)){
        
        case_sel=unique(data_sp$id)
        dfi_ale=list()
        
        data_sp_2 <- cbind(age_matrix,X_der_prim,Z_der_prim)
        names(data_sp_2)=c( "id","x",paste0("dX",1:DX[2]),paste0("dZ",1:DZ[2]))
        
        for(i in 1:n_distinct(data_sp$id)){
          
          dat = data_sp_2 %>%  filter(id==case_sel[i])
          
          dfi_ale[[i]] = as.matrix(dat[,4])%*%c(model$model.out$coefficients$random$id[i,4])+
            as.matrix(dat[,c(5:(4+DZ[2]))])%*%c(model$model.out$coefficients$random$id[i,c(5:(4+DZ[2]))])
        }
        
        dfi_ale=do.call("rbind",dfi_ale)
        
      }
      
    }
    if(model$PB_model==1){
      res=cbind(deriv_pb, PB_model=model$PB_model, df_all=deriv_pb$fit_df_pb, df_all_pop=deriv_pb$fit_df_pop_pb)
      aphv_est=phv_estimates(x_grid=res$x, y_grid=res$fit_df_pb, id=res$id)
    }
    if(model$PB_model==2){
      res=cbind(deriv_pb,df_fixed=df_fixed, PB_model=model$PB_model, df_all=deriv_pb$fit_df_pb+df_fixed, df_all_pop=deriv_pb$fit_df_pop_pb+df_fixed)
      aphv_est=phv_estimates(x_grid=res$x, y_grid=res$fit_df_pb+res$df_fixed, id=res$id)}
    if(model$PB_model==3){
      res=cbind(deriv_pb,dfi_ale=dfi_ale, PB_model=model$PB_model, df_all=deriv_pb$fit_df_pb+dfi_ale, df_all_pop=deriv_pb$fit_df_pop_pb)
      aphv_est=phv_estimates(x_grid=res$x, y_grid=res$fit_df_pb+res$dfi_ale, id=res$id)}
    if(model$PB_model==4){
      res=cbind(deriv_pb,df_fixed=df_fixed,dfi_ale=dfi_ale, PB_model=model$PB_model, df_all=deriv_pb$fit_df_pb+dfi_ale+df_fixed, df_all_pop=deriv_pb$fit_df_pop_pb+df_fixed)
      aphv_est=phv_estimates(x_grid=res$x, y_grid=res$fit_df_pb+res$df_fixed+res$dfi_ale, id=res$id)}
    
    return(list(res=res,aphv_est=cbind(aphv_est,str=model$PB_model)))
  } else   return(NA)
  
}

##########################################
## Plot parts SPNLME model
#########################################

semiprametric_parts=function(model, ndx=4){
  # ajuste PB(x)

  DX=dim(model[["Matrix_trans"]][["X"]])
  DZ=dim(model[["Matrix_trans"]][["Z"]])
  
  data_sp =model$data_sp
  PB_pop = NULL
  PB_all = NULL
  
  id_uniqu = unique(data_sp$id)

  if(class(model$model.out)[1]!="try-error"){
    
    for( i in 1:n_distinct(data_sp$id) ){
      
      x = data_sp %>%  filter(id==id_uniqu[i]) %>% pull(x)
      a_pop = model$model.out$coefficients$fixed[1]
      b_pop = model$model.out$coefficients$fixed[2]
      c_pop = model$model.out$coefficients$fixed[3]
      
      a_tot =  a_pop + model$model.out$coefficients$random$id[i,1]
      b_tot = b_pop + model$model.out$coefficients$random$id[i,2]
      c_tot = c_pop + model$model.out$coefficients$random$id[i,3]
      s0 = model$model.out$coefficients$fixed[4]
      s1 = model$model.out$coefficients$fixed[5]
      
      PB_pop = c( PB_pop , a_pop - (2*(a_pop - c_pop))/(exp(s0*(x - b_pop)) + exp(s1*(x - b_pop))))
      PB_all = c( PB_all , a_tot - (2*(a_tot - c_tot))/(exp(s0*(x - b_tot)) + exp(s1*(x - b_tot))))
      
    }
    
    data_sp$PB_all_f=PB_all ### 
    data_sp$PB_pop_f=PB_pop
    
    ################### ajuste f(x)/fixed comp
    
    if(model$PB_model%in%c(2,4)){
      f_fixed = data_sp$X2*model$model.out$coefficients$fixed[6] +
        as.matrix(data_sp %>% dplyr::select(starts_with("Z")))%*%c(model$model.out$coefficients$random$Id2)
      
      data_sp$f_fixed=f_fixed
    }
    
    ################### ajuste fi(x)/random comp
    if(model$PB_model%in%c(3,4)){
      fi_ale=list()
      
      for(i in 1:n_distinct(data_sp$id)){
        
        dat = data_sp %>%  filter(id==id_uniqu[i]) %>% dplyr::select(paste0("X",2:DX[2]),paste0("Z",1:DZ[2]))
        
        fi_ale[[i]] = as.matrix(dat$X2)%*%c(model$model.out$coefficients$random$id[i,4])+
          as.matrix(dat %>% dplyr::select(starts_with("Z")))%*%c(model$model.out$coefficients$random$id[i,c(5:(4+DZ[2]))])
        
      }
      
      fi_ale=do.call("rbind",fi_ale)
      data_sp$fi_ale=fi_ale
    }
    data_sp= data_sp%>% rename(age=x)
  } else {
    data_sp=NA
  }
  
  
  return(data_sp)
}

###############################
# classification with SPNLME
###############################
pi.is_cond<-function (data, model ,spline_struct=1 , ndx, pord, bdeg, v=10, eps=1e-10, x_min, x_max)
{
  MM=100
  KM=round(MM/4,0)
  
  data=data %>% arrange(id)
  #### data info indiv
  yobs = data$y
  xind = data$x
  id = data$id
  n_id=n_distinct(id)
  
  #### data info spline
  Matrix_trans <- MM.basis(xind,x_min, x_max, nseg = ndx, bdeg = 3, pord = 2)
  
  X<- as.data.frame(Matrix_trans$X)
  colnames(X)=paste0("X",1:dim(X)[2])
  Z<- as.data.frame(Matrix_trans$Z)
  colnames(Z)=paste0("Z",1:dim(Z)[2])
  DX = dim(X)[2]-1
  DZ = dim(Z)[2]
  
  #### model info define random eff
  names_re_par = as.data.frame(model$coefficients$random$id) %>%
    dplyr::select(!(tidyselect::starts_with("f.") | tidyselect::starts_with("f_i."))) %>%
    colnames()
  
  n_random_effects_param = length(names_re_par)#dim(model$coefficients$random$id)[2] # number random effect parametric part
  
  names_re_nopar=as.data.frame(model$coefficients$random$id) %>%
    dplyr::select((tidyselect::starts_with("f.") | tidyselect::starts_with("f_i."))) %>%
    colnames()
  
  n_random_effects_no_param = length(names_re_nopar)
  
  n_ran_total=n_random_effects_param+n_random_effects_no_param
  
  n_fixed_effects_param = length(model$coefficients$fixed)
  
  sel_r = 1:n_random_effects_param
  
  model_str_param= expression("a_i-(2*(a_i-c_i))/(exp(s0*(x-b_i)) + exp(s1*(x-b_i)))")
  
  #### model Omega matrix RE
  if( spline_struct%in%c(1,2)){
    
    # Check the value of spline_struct
    if (spline_struct == 2) {
      # Extract variances for the specified range when spline_struct is 2
      sm1 <- matrix(as.numeric(VarCorr(model)[-(1:(ndx + bdeg - pord + 2)), 2]), ncol = 1, byrow = FALSE)
    } else {
      # Extract the first n_random_effects_param variances when spline_struct is not 2
      sm1 <- matrix(as.numeric(VarCorr(model)[1:n_random_effects_param, 2]), ncol = 1, byrow = FALSE)
    }
    
    # Check if there's only one random effect
    if (n_random_effects_param == 1) {
      Omega <- sm1[1, 1]
    } else {
      # For multiple random effects, initialize Omega as a diagonal matrix
      Omega <- diag(sm1[sel_r, 1])
      
      # Correlation matrix for the random effects
      corr.m <- as.matrix(corMatrix(model$modelStruct$reStruct$id)[sel_r, sel_r])
      
      # Update off-diagonal elements of Omega using the correlation matrix
      for (i in 2:n_random_effects_param) {
        for (j in 1:(n_random_effects_param - 1)) {
          if (i != j) {
            Omega[i, j] <- Omega[j, i] <- corr.m[i, j] * Omega[i, i] * Omega[j, j]
          }
        }
      }
      
      diag(Omega)=diag(diag(sm1[sel_r,1]))**2
    }
    
  }else  if( spline_struct==3){
    
    # Extract variance-covariance values into a numeric matrix
    sm1 <- matrix(as.numeric(VarCorr(model)[seq_len(n_ran_total), 2]),
                  ncol = 1,
                  byrow = FALSE
    )
    
    # Initialize Omega based on the number of random effects parameters
    if (n_random_effects_param == 1) {
      Omega <- sm1[1, 1]
    } else {
      # Construct diagonal matrix for Omega with random effects variances
      Omega <- diag(sm1[sel_r, 1])
      
      # Extract correlation matrix for random effects
      corr.m <- as.matrix(corMatrix(model$modelStruct$reStruct$id)[sel_r, sel_r])
      
      # Update off-diagonal elements of Omega based on correlations
      for (i in seq_len(n_random_effects_param)) {
        for (j in seq_len(i - 1)) {
          Omega[i, j] <- Omega[j, i] <- corr.m[i, j] * Omega[i, i] * Omega[j, j]
        }
      }
      
      # Square the diagonal elements for proper variance scaling
      diag(Omega) = diag(diag(sm1[sel_r,1]))**2
    }
    Omega_X_fi=as.matrix((sm1[n_random_effects_param+1,1]))**2
    Omega_Z_fi=diag(sm1[(n_random_effects_param+2):length(sm1),1])**2
    
    Omega=Matrix::bdiag(Omega, Omega_X_fi, Omega_Z_fi)
  }else  if( spline_struct==4){
    # Extract the variance-covariance components from VarCorr
    sm1 = matrix(as.numeric(VarCorr(model)[-c(1:(ndx + bdeg - pord + 2)), 2]), ncol = 1, byrow = FALSE)
    
    # Initialize Omega depending on the number of random effects parameters
    if (n_random_effects_param == 1) {
      Omega = sm1[1, 1]
    } else {
      Omega = diag(sm1[sel_r, 1])
      corr.m = as.matrix(corMatrix(model$modelStruct$reStruct$id)[sel_r, sel_r])
      
      # Compute the correlation matrix for the random effects
      for (i in 2:n_random_effects_param) {
        for (j in 1:(n_random_effects_param - 1)) {
          if (i != j) {
            Omega[i, j] = Omega[j, i] = corr.m[i, j] * Omega[i, i] * Omega[j, j]
          }
        }
      }
      # Correct diagonal values of Omega using sm1 values
      diag(Omega)=diag(diag(sm1[sel_r,1]))**2
    }
    
    # Extract Omega_X_fi and Omega_Z_fi
    Omega_X_fi = sm1[n_random_effects_param + 1, 1]**2
    Omega_Z_fi = diag(sm1[(n_random_effects_param + 2):(length(sm1) - 1), 1])**2
    
    # Combine all components into a block diagonal matrix
    Omega = Matrix::bdiag(Omega, Omega_X_fi, Omega_Z_fi)
  }
  
  sigma = model$sigma
  coef = model$coefficients
  
  ################3Random values
  
  Omega1 = as.matrix(Omega)
  IOmega.phi1 = MPinv(Omega1)
  det_Omega1 = prod(eigen(Omega1)$values[eigen(Omega1)$values>0+eps])
  
  ################ Prepare matrices
  yM <- rep(yobs, MM) # y observations
  XM <- do.call(c, rep(list(xind), MM)) # x value
  
  ################ Fixed effect asociated RE param_part
  index = rep(1:n_id, times = count(data, id)$n)
  
  fixed_eff = matrix(rep(c(coef$fixed[1:n_random_effects_param]), n_id),ncol=n_random_effects_param,nrow= n_id,byrow=TRUE) # repeat fixed effects as inds.
  #   
  
  ################### Fill io matrix
  io <- matrix(0, nrow = n_id, ncol = max(count(data, id)$n)) # matrix to save the etimation , how it is no balanced they took the ncol as max of nind.obs
  for (isuj in 1:n_id) io[isuj, 1:count(data, id)$n[isuj]] <- 1 # Put 1 if it has a measurement in the time j
  
  ioM <- matrix(rep(t(io), MM), ncol = dim(io)[2], byrow = TRUE)
  ind.ioM <- which(t(ioM) != 0) # to know were the id has measurements in the matrix io
  
  DYF <- matrix(0, nrow = dim(ioM)[2], ncol = dim(ioM)[1])
  
  #Repeat the fixed effects MM paramtric model
  
  mean.ind.param <- matrix(rep(t(fixed_eff), MM), byrow = TRUE,ncol = n_random_effects_param) #REpeat the fixed effects MM N*NN x random_eff
  
  
  #detrminate the random effects for the parametric model
  ind.param_fix <- matrix(rep(t(coef$random$id[,1:n_random_effects_param]),
                              MM), byrow = TRUE, ncol = n_random_effects_param) # REpeat the individual parameters MM
  
  
  # Repeat the standard deviation for the individual parameters ESIIMATED  MM times parametric model
  var.ind.param <- matrix(rep(t((diag(cov(as.matrix(coef$random$id))))), n_id), byrow = TRUE, ncol = n_ran_total)
  mean.var.ind.param <- matrix(rep(colMeans(var.ind.param),n_id),ncol=n_ran_total,nrow=n_id,byrow=TRUE) # similar to the    diag(Omega), they take de mean of the conditional variance
  sd.ind.param <- matrix(rep(t(sqrt(var.ind.param)), MM), byrow = TRUE,  # sqrt of the coditional varias of random eff N*MM x rand_eff
                         ncol = n_ran_total)
  
  
  if( spline_struct%in%c(3,4)){ # add the random for the f_i in the case of splines
    
    # REpeat the individual parameters for f_i ESIIMATED  MM times
    ind.param.X <- colMeans(matrix(rep(t(coef$random$id[,"f_i.X2"]), MM), byrow = TRUE, ncol = 1))
    ind.param.Z <- colMeans(matrix(rep(t(coef$random$id[,paste0("f_i.Z",1:(ndx+bdeg-pord))]), MM), byrow = TRUE, ncol = ndx+bdeg-pord)) # REpeat the individual parameters ESIIMATED  MM times
    
    mean.ind.param=cbind(mean.ind.param,matrix(0, ncol=n_random_effects_no_param, nrow=(n_id*MM)))
    ind.param=cbind(ind.param_fix,rep(ind.param.X,dim(coef$random$id)[1]*MM),
                    matrix(rep(ind.param.Z, each=dim(coef$random$id)[1]*MM),ncol=length(ind.param.Z), byrow=F))
    
  }else{ind.param=ind.param_fix}
  
  meana <- rep(0, n_id) # vector o with the # of id
  LL <- matrix(0, nrow = KM, ncol = 1)
  LL.i <- matrix(0, nrow = KM, ncol =  n_id)
  
  
  if(is.nan(log(det(Omega1)))){
    LL.i_nan = rep(NaN, n_id)
    return(LL.i_nan)
  }else{
    
    c2 <- log(det_Omega1) + dim(Omega1)[1] * log(2 * pi)
    c1 <- log(2 * pi)
    
    for (km in seq_len(KM)) {
      
      ################# likelihood random eff
      
      r = trnd.mlx(v,  n_id* MM, n_ran_total)
      
      new.ind.param <- matrix(rep(colMeans(ind.param),  n_id * MM),ncol= n_ran_total,nrow=  n_id  *  MM,byrow=TRUE) +
        matrix(rep(colMeans(sd.ind.param),  n_id  * MM),ncol= n_ran_total,nrow=  n_id * MM,byrow=TRUE) * r
      
      
      
      # new.ind.param <- matrix(rep(colMeans(ind.param),  n_id * MM),ncol= n_ran_total,nrow=  n_id  *  MM,byrow=TRUE) +
      #   matrix(rep(sqrt(diag(Omega1)),  n_id  * MM),ncol= n_ran_total,nrow=  n_id * MM,byrow=TRUE) * r
      # 
      
      new.ind.param =new.ind.param+mean.ind.param
      new.random.eff <- new.ind.param - mean.ind.param # Diference between (mean individual parameters+ sigma*t) and fixed effects
      
      
      d2 <- (-0.5) * (rowSums(new.random.eff * (new.random.eff %*% IOmega.phi1)) +   c2) #  density function random effects
      e2 <- matrix(d2, nrow = n_id, ncol = MM)
      
      ################# likelihood  t-student importance sampling
      
      d3 <- rowSums(log(tpdf.mlx(r, v)))
      e3 <- matrix(d3, nrow = n_id, ncol = MM) -
        matrix(rep(0.5 * rowSums(log(mean.var.ind.param)), MM), ncol = MM) ## REVIEW....
      
      
      ################  likelihood observ
      data.param <- vector("list", MM)
      # Create base data frame
      dat1 <- data.frame(x = xind, y = yobs, id = id)
      
      # Loop through iterations
      for (j in seq_len(MM)) {
        # Calculate indices for the current slice
        start_idx <- (j - 1) * n_id + 1
        end_idx <- j * n_id
        
        
        # Extract individual parameters and merge with data
        ind_par <- data.frame(new.ind.param[start_idx:end_idx,], id = unique(id))
        colnames(ind_par) <- c(colnames(coef$random$id), "id")
        data.param[[j]] <- merge(dat1, ind_par, by = "id")
      }
      
      # Combine all iterations into a single data frame
      data.param <- do.call(rbind, data.param)
      
      # Add fixed effect parameters
      data.param <- data.param %>%
        mutate(s0 = coef$fixed["s0"],
               s1 = coef$fixed["s1"])
      
      # Evaluate the model string to calculate f_par
      estim_model <- data.param %>%
        mutate(f_par = eval(parse(text = model_str_param)))
      
      estim_model_pop= estim_model %>% group_by(id,x) %>% mutate(f_par_mean=mean(f_par))
      
      # id_sel= unique(id)[2]
      # 
      # if(km==KM){
      #   new_data = data[c("id", "x")]
      #   curv_prom = data.frame(x= data_model$x,"a_i"=model$coefficients$fixed["a_i"],
      #                          "b_i"=model$coefficients$fixed["b_i"],
      #                          "c_i"=model$coefficients$fixed["c_i"],
      #                          "s0"=model$coefficients$fixed["s0"],
      #                          "s1"=model$coefficients$fixed["s1"])  %>%
      #     mutate(pred_pb = eval(parse(text = model_str_param)))
      # 
      #   p1= estim_model %>% filter(id %in% id_sel) %>% ggplot(aes(x=x,y=y))+
      #     geom_point(aes(y=f_par, group=id, col="fpar"))+labs(y="y_PB")+
      #     geom_line(data=curv_prom, aes(x=x ,y=pred_pb, col="pop_pb"), lty=2)+
      #     geom_line(data=estim_model_pop%>% filter(id %in% id_sel), 
      #               aes(x=x ,y=f_par_mean, group=id, col="mean_sim"))+
      #     geom_point(aes(col="data_id"))+labs(colour="model")+
      #     ggtitle(paste0("model", spline_struct))+
      #     labs(colour="model")+theme_bw()+
      #     theme(legend.position = "bottom")
      # }
      
      
      ################### ajuste f(x)/fixed comp
      
      # Calculate fixed component for specific spline structures
      if (spline_struct %in% c(2, 4)) {
        f_fixed <- as.matrix(X["X2"]) * coef$fixed["X_2"] +
          as.matrix(Z) %*% t(coef$random$Id2)
        
        estim_model$f_fixed <- rep(c(f_fixed), MM)
      }
      
      # Random component adjustment for specific spline structures
      if (spline_struct %in% c(3, 4)) {
        # Replicate and bind data for all MM iterations
        fi_ale <- do.call(rbind, replicate(MM, cbind(X[paste0("X",2:(DX+1))],Z),
                                           simplify = FALSE))
        
        # Extract and adjust relevant columns
        estim_model_ale <- estim_model %>% dplyr::select(tidyselect::starts_with("f_i."))
        estim_model_ale <- estim_model_ale* fi_ale 
        
        # Sum rows for the adjusted random effects
        estim_model$fi_ale <- rowSums(estim_model_ale)
      }
      
      # Combine components based on spline structure
      f <- c(estim_model$f_par)  # Initialize with base component
      if (spline_struct %in% c(2, 4)) f <- f + estim_model$f_fixed
      if (spline_struct %in% c(3, 4)) f <- f + estim_model$fi_ale
      
      estim_model$f =f
      
      
      # if(km==KM){
      # 
      # 
      #   # if(km==KM){
      #   if(spline_struct %in% c(1)){
      #     
      #     fit_model=data.frame(model$fitted) %>% rename(id_fit=id) %>% mutate(x=data_model$x)
      #     
      #     plot2_dat =estim_model %>% group_by(id,x,y) %>% summarise(f_fit=f,f=f_par,f_p=f_par, mean_f=mean(f_par), 
      #                                                               mean_all=mean(f_p)) %>% arrange(id,x)
      #     p2= plot2_dat %>% filter(id %in% c(id_sel)) %>%
      #       ggplot(aes(x=x,y=f))+labs(y="y_nopar", colour="model")+
      #       theme_bw()+theme(legend.position = "bottom")
      #     
      #   } else if(spline_struct %in% c(2)){
      #     
      #     fit_model=data.frame(model$fitted) %>% rename(id_fit=id, fixed2=fixed, fixed=Id2) %>%
      #       mutate(x=data_model$x)
      #     
      #     plot2_dat =estim_model %>% group_by(id,x,y) %>% summarise(f_fit=f,f=f_fixed,f_p=f_par+f_fixed, 
      #                                                               mean_f=mean(f_fixed), mean_all=mean(f_p)) %>% arrange(id,x)
      #  
      #     p2= plot2_dat %>% filter(id %in% c(id_sel)) %>%
      #       ggplot(aes(x=x,y=f))+
      #       geom_point(aes( col="f_nopar", group=id), alpha=0.5)+
      #       geom_line(aes(x=x ,y=mean_f,col="mean_sim"))+labs(y="y_nopar", colour="model")+
      #       # geom_line(data=estim_model%>% filter(id %in% c(id_sel)),aes( x=x, y=f_fixed,col="nonpar_model"), lty=2)+
      #       theme_bw()+theme(legend.position = "bottom") # el mean_sim y el pop_model son iguales pues este no varia
      #     
      #      } else if(spline_struct %in% c(3)){
      #     
      #     # fi_model=data.frame(model$coefficients$random$id) %>% mutate(id=rownames(model$coefficients$random$id))
      #     # 
      #     # fit_model=merge(data.frame(model$fitted)%>% rename(id_fit=id) %>%
      #     #                   mutate(x=data_model$x, id=id),fi_model[,-c(1:3)], by="id")          
      #     fit_model=data.frame(model$fitted)%>% rename(id_fit=id) %>%
      #                       mutate(x=data_model$x)
      #     # 
      #     # # Extract and adjust relevant columns
      #     # 
      #     # fit_model_f_i <- fit_model %>% dplyr::select(tidyselect::starts_with("f_i."))
      #     # fit_model_f_i <- fit_model_f_i* cbind( X[,2], Z)
      #     # fit_model$nonpar_fi <- rowSums(fit_model_f_i)
      #     
      #     
      #     plot2_dat =estim_model %>% group_by(id,x,y) %>% summarise(f_fit=f,f=fi_ale,f_p=f_par, mean_f=mean(fi_ale), 
      #                                                               mean_all=mean(f_p)) %>% arrange(id,x)
      #     
      #     p2= plot2_dat %>% filter(id %in% c(id_sel)) %>%
      #       ggplot(aes(x=x,y=f))+
      #       geom_point(aes( col="f_nopar", group=id), alpha=0.5)+
      #       geom_line(aes(x=x ,y=mean_f,col="mean_sim"))+labs(y="y_nopar", colour="model")+
      #       # geom_line(data=fit_model %>%filter(id %in% c(id_sel)) , aes( x=x, y=nonpar_fi,col="nonpar_model"), lty=2)+
      #       theme_bw()+theme(legend.position = "bottom")
      #     
      #   } else if(spline_struct %in% c(4)){
      #     
      #     f_fixed2 <- data_model["X2"]* coef$fixed["X_2"] +
      #       as.matrix(data_model[paste0("Z", 1:(ndx+bdeg-pord))]) %*% t(coef$random$Id2)
      #     
      #     fit_model=data.frame(model$fitted)%>% rename(id_fit=id, fixed2=fixed, fixed=Id2) %>%
      #                       mutate(x=data_model$x)
      #     
      #     fit_model$nonpar_fix=f_fixed2[,1]
      #     
      #     # Extract and adjust relevant columns
      #     
      #     plot2_dat =estim_model %>% group_by(id,x,y) %>% summarise(f_fit=f,f=fi_ale+f_fixed,f_p=f_par+f_fixed, mean_f=mean(fi_ale+f_fixed),
      #                                                               mean_all=mean(f_p)) %>% arrange(id,x)
      #  
      #     
      #     p2= plot2_dat %>% filter(id %in% c(id_sel)) %>%
      #       ggplot(aes(x=x,y=f))+
      #       geom_point(aes( col="f_nopar", group=id), alpha=0.5)+
      #       geom_line(aes(x=x ,y=mean_f,col="mean_sim"))+labs(y="y_nopar", colour="model")+
      #       # geom_line(data=fit_model , aes( x=x, y=nonpar_fix,col="nonpar_model"), lty=2)+
      #       theme_bw()+theme(legend.position = "bottom")
      #     
      #      }
      #  
      #   # }
      # 
      #   # if(km==KM){
      # 
      # 
      #   p3= plot2_dat %>% filter(id %in% c(id_sel)) %>%
      #     ggplot(aes(x=x,y=y))+
      #     geom_point(aes(x=x,y=f_fit, col="f_estim", group=id), alpha=0.5)+
      #     geom_line(aes(x=x ,y=mean_all, col="mean_sim" ))+labs(y="y_fit", colour="model")+
      #     geom_line(data=fit_model, aes(x=x ,y=fixed, col="pob_model"), lty=2)+
      #    geom_point(aes(col="data"))+ theme_bw()+
      #    theme(legend.position = "bottom")
      # 
      #   # }
      # 
      #   plot_all=cowplot::plot_grid(p1+theme(legend.position = "bottom"),
      #                               p2+theme(legend.position = "bottom"),
      #                               p3+theme(legend.position = "bottom"),
      #                               ncol=3)
      #   print(plot_all)
      # 
      #   # pdf(file = paste0("plots_compare.pdf"),   # The directory you want to save the file in
      #   #     width = 22, # The width of the plot in inches
      #   #     height = 10) # The height of the plot in inches
      #   # list(plot_all)
      #   # dev.off()
      # 
      # 
      # }
      
      
      g <- rep(sigma,l=length(f)) # repeat the sigma error measurement estimated
      
      DYF[ind.ioM] <- -0.5 * ((yM - f)/g)^2 - log(g) - 0.5 *  c1 # density of (y- hat(y))
      e1 <- matrix(colSums(DYF), nrow = n_id, ncol = MM)
      dim(e1)
      ###### compute all likelihood
      sume <- e1 + e2 - e3 # sum the densities because are normal distribution and the log are applied; dim= # ind, MM
      newa <- rowMeans(exp(sume), na.rm = TRUE)  # sacar la media por individuo sobre los 200 valores, vector of length number of ind
      meana <- meana + 1/km * (newa - meana) # review account the before estimation
      LL.i[km,]<-(cutoff(meana)) # value of density in each x
      
    }
    
    return(LL.i[KM,])
    
  }
  
  
  
}
