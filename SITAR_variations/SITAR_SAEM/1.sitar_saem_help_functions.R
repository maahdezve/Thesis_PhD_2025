##########################
# libraries
##########################

library(MASS)
library(dplyr)
library(nlme)
library(sitar)
library(saemix)
library(splines)
library(glue)
require(JM)


########################################################################
# Functions 
#########################################################################


#########################################
# Initial_values sitar_saem
#########################################

initial_values = function(x,y, knots, bounds, df, type="lm", start=NULL, It_rw=20000, var_rw = 0.2, burn_rw=0.1, seed=111){
  
  yoffset <- mean(y)
  xoffset <- mean(x)
  bstart <- mean(x) - xoffset
  x <- x - xoffset
  knots <- knots - xoffset
  bounds <- bounds - xoffset
  
  if(type=="lm"){
    
    ##########################################
    #############  lm
    ##########################################
    spline.lm <- lm(y ~ ns(x, knots = knots, Bound = bounds))
    start <- coef(spline.lm)[c(2:(df + 1), 1)]              ### start values for ns parameters and a
    start_var=diag(vcov(spline.lm))
    start_var = c(start_var[1],0.5,0.5,start_var[-1])
    start<-c(start, bstart,0)
    start=start[c((df+1):(df+3),1:(df))] 
    if(sum(is.na(start))!=0)start=c(mean(y),bstart,0,rep(runif(1,20,50),df))
    
    names(start) = names(start_var) = c("a","b","c",paste0("s",1:df))
    
    return(data.frame(start=start, start_var=start_var))
    
  }else if(type=="rw"){
    
    ##########################################
    ############# Random walk
    ##########################################
    
    numb<-seq(4,4+df-1,1)
    ss<-paste0("s",1:df)
    ss_spline_i_1<-glue:::glue_collapse(glue('{ss}=theta[w+1,{numb}]'), ",")
    ss_spline_i<-glue:::glue_collapse(glue('{ss}=theta[w,{numb}]'), ",")
    spline_prod2<-glue( 'rowSums((ns((x-(b))/exp(-c),k=knots,B=bounds))%*%rbind(', paste(ss, collapse = ','),'))')
    
    ## e= y-ax^2-bx, it has normal distribution, the the density is given by:
    par_f= c("a","b","c", ss)
    par_f=paste(par_f, collapse = ",")
    
    f_glue<-glue(
      "f_e<-function(<<par_f>>) {\n",
      "sit_f= a-yoffset+<<spline_prod2>>\n",
      "sum((-0.5*log(2*pi*1) -0.5*(y-yoffset-sit_f)^2))\n",
      "}",
      .open = "<<",  .close = ">>"
    )
    
    f_e<- eval(parse(text = f_glue))
    
    #############333 random walk
    
    set.seed(seed)
    It=It_rw
    var=var_rw# similar al aplha gradiente descendiente
    theta=matrix(NA, ncol=df+3, nrow=It+1)
    prob=matrix(NA, ncol=1, nrow=It+1)
    theta[1,]=c(rnorm(1,130-yoffset,10), rnorm(1,0,0.5), rnorm(1,0,0.5),rep(runif(1,20,50),df))#start.b
    count_v=0 ## muestras aceptadas
    
    for(w in 1:(It+1)){
      
      if(w==It+1)break
      
      theta[w+1,]=theta[w,]+var*c(rnorm(1,0,1),rnorm(1,0,0.5),rnorm(1,0,0.5),runif(df,-1,1))
      
      p1<-glue(
        "f_e(a=theta[w+1,1],\n",
        "b=theta[w+1,2],\n",
        "c=theta[w+1,3],\n",
        "<<ss_spline_i_1>>)",
        .open = "<<",  .close = ">>"
      )
      
      prob[w+1,]= eval(parse(text = p1))
      
      p2<-glue(
        "f_e(a=theta[w,1],\n",
        "b=theta[w,2],\n",
        "c=theta[w,3],\n",
        "<<ss_spline_i>>)\n",
        .open = "<<",  .close = ">>"
      )
      
      prob[w,]=eval(parse(text = p2))
      
      alpha= try(min(1,exp(prob[w+1,1]-prob[w,1])), silent = T)
      # if(is.nan(alpha)) next
      # alpha= min(0,prob[i+1,]-prob[i,] )
      # alpha= min(1,prob[i+1,]/prob[i,] )
      u =runif(1)
      
      if(u<alpha) {
        count_v=count_v+1
        theta[w+1,]=theta[w+1,]
        prob[w+1,]=prob[w+1,]
      }else{
        theta[w+1,]=theta[w,]
        prob[w+1,]=prob[w,]
      }
    }
    count_v
    ##########################################
    ############# burning  quitar el 10% de las muestras 
    ##########################################
    burn = round(It*burn_rw)
    theta=theta[-c(1:(burn+1)),]
    prob=prob[-c(1:(burn+1))]
    
    ###########################################
    # muestreo sistematico cada n pasos
    ############################################
    
    id_sis=seq(1,length(prob), by=10)
    
    theta_sel=theta[id_sis,]
    prob_sel=prob[id_sis]
    
    ###########################################
    # mvalores iniciales seleccionados
    ############################################
    
    # maximo apost initial value
    param_sel=theta_sel[which.max(prob_sel),] # maximo apost
    param_sel[2:3]=tanh(param_sel[2:3])
    
    # media condicional initial value
    param_sel_mean=colMeans(theta_sel)
    param_sel_mean[2:3]=tanh(param_sel_mean[2:3])
    
    # varainza cond inicial values
    param_sel_sd=apply(theta_sel,2, sd)
    
    star_map = param_sel
    start_cond= param_sel_mean
    start_sd_cond =param_sel_sd
    
    names(star_map) = names(start_cond) = names(start_sd_cond) = c("a","b","c",paste0("s",1:df))
    
    start = data.frame(star_map, start_cond, start_sd_cond)
    row.names(start) = c("a","b","c",paste0("s",1:df))
    attr(start, "count_v")=count_v
    
    return(start)
    
    
  }
  else{
    print("error in type, choose among the options lm, rw  or  free")
  }
  
  
}


##################################################
# star_criteria and convergence index sitar_saem
##################################################


number_iter_k1 =  function(x, y, id, df, k1_initial= 1000,seed = 789, error_abs=0.01, error_rel=0.001,start, m=2){
  
  ### modify saemix fnction to get historical values for random effects to calculate the convergence index
  
  saemix.mod<-function(model,data,control=list()) {
    
    # Convergence plots during fit (special function, not user-level)
    convplot.infit<-function(allpar,K1,niter=0) {
      # Convergence plots for all the fixed effects, random effects and residual variability
      oldpar <- par(no.readonly = TRUE)    # code line i
      on.exit(par(oldpar))            # code line i + 1 
      np<-dim(allpar)[2]
      K<-dim(allpar)[1]
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      if(n1>5 | n2>5) {n1<-3;n2<-4}
      if(niter==0) niter<-K
      par(mfrow=c(n1,n2))
      for(j in 1:np) {
        plot(1:niter,allpar[1:niter,j],type="l", xlab="Iteration", ylab=colnames(allpar)[j])
        abline(v=K1)
      }
    }
    if(!is(model,"SaemixModel")) {
      message("Please provide a valid model object (see the help page for SaemixModel)\n")
      return("Missing model")
    }
    if(!is(data,"SaemixData")) {
      message("Please provide a valid data object (see the help page for SaemixData)\n")
      return("Missing data")
    }
    
    saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
    #  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
    opt.warn<-getOption("warn")
    if(!saemixObject["options"]$warnings) options(warn=-1)
    
    saemix.options<-saemixObject["options"]
    saemix.model<-saemixObject["model"]
    saemix.data<-saemixObject["data"]
    saemix.data@ocov<-saemix.data@ocov[saemix.data@data[,"mdv"]==0,,drop=FALSE]
    saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
    saemix.data@ntot.obs<-dim(saemix.data@data)[1]
    #  showall(saemixObject)
    
    # Initialising random generator
    set.seed(saemix.options$seed)
    
    ############################################
    #  Main Algorithm
    ############################################
    
    # Initialisation - creating several lists with necessary information extracted (Uargs, Dargs, opt,varList, suffStat)
    xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
    saemix.model<-xinit$saemix.model
    Dargs<-xinit$Dargs
    Uargs<-xinit$Uargs
    varList<-xinit$varList
    phiM<-xinit$phiM
    mean.phi<-xinit$mean.phi
    DYF<-xinit$DYF
    opt<-xinit$opt
    betas<-betas.ini<-xinit$betas
    fixed.psi<-xinit$fixedpsi.ini
    var.eta<-varList$diag.omega
    initials_val=list(allpar=xinit$allpar0,mean.phi.l=mean.phi,phi.l=phiM,Omega.l=varList$omega)
    
    if (Dargs$modeltype=="structural"){
      theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
      parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=(Uargs$nb.parameters+length(Uargs$i1.omega2)+length(saemix.model["indx.res"])))
      colnames(parpop)<-c(saemix.model["name.modpar"], saemix.model["name.random"], saemix.model["name.sigma"][saemix.model["indx.res"]])
      allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(Uargs$nb.betas+length(Uargs$i1.omega2)+length(saemix.model["indx.res"])))
      colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"], saemix.model["name.sigma"][saemix.model["indx.res"]])
      phi.l =list() ##** added
      mean.phi.l =list() ##** added
      Omega.l=list() ##** added
      etaM.l=list() ##** added
      respar_vec=matrix(NA, ncol=2, nrow=saemix.options$nbiter.tot)
    } else{
      theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2])
      parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=(Uargs$nb.parameters+length(Uargs$i1.omega2)))
      colnames(parpop)<-c(saemix.model["name.modpar"], saemix.model["name.random"])
      allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(Uargs$nb.betas+length(Uargs$i1.omega2)))
      colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"])
    }
    
    parpop[1,]<-theta0
    allpar[1,]<-xinit$allpar0
    
    # using several Markov chains - only useful if passed back to main routine...
    # 	chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=saemix.options$nb.chains)
    # 	NM<-chdat["NM"]
    # 	IdM<-chdat["dataM"]$IdM
    # 	yM<-chdat["dataM"]$yM
    # 	XM<-chdat["dataM"][,saemix.data["name.predictors"],drop=FALSE]
    
    # List of sufficient statistics - change during call to stochasticApprox
    suffStat<-list(statphi1=0,statphi2=0,statphi3=0,statrese=0)
    phi<-array(data=0,dim=c(Dargs$N, Uargs$nb.parameters, saemix.options$nb.chains))
    
    # structural model, check nb of parameters
    structural.model<-saemix.model["model"]
    #  nb.parameters<-saemix.model["nb.parameters"]
    
    # Running the algorithm
    # hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...');
    if(saemix.options$displayProgress) par(ask=FALSE)
    if(saemix.options$warnings) cat("Running main SAEM algorithm\n")
    if(saemix.options$warnings) print(date())
    
    #-----------------------------------
    #
    #----------------------------------------
    
    
    for (kiter in 1:saemix.options$nbiter.tot) { # Iterative portion of algorithm
      
      # SAEM convergence plots
      if(kiter%%saemix.options$nbdisplay==0) {
        if(saemix.options$warnings) cat(".")
        if(saemix.options$displayProgress)
          try(convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2)))
      }
      # Burn-in - resetting sufficient statistics
      if(opt$flag.fmin && kiter==saemix.options$nbiter.sa) {
        Uargs$COV1<-Uargs$COV[,Uargs$ind.fix11]
        ind.prov<-!(varList$ind.eta %in% Uargs$i0.omega2)
        varList$domega2<-varList$domega2[ind.prov,ind.prov,drop=FALSE] # keep in domega2 only indices of parameters with IIV
        varList$ind0.eta<-Uargs$i0.omega2
        varList$ind.eta<-1:(Uargs$nb.parameters)
        if(length(varList$ind0.eta)>0) varList$ind.eta<-varList$ind.eta[!(varList$ind.eta %in% varList$ind0.eta)] # update ind.eta, now only parameters with IIV
        Uargs$nb.etas<-length(varList$ind.eta)
        suffStat$statphi1<-0
        suffStat$statphi2<-0
        suffStat$statphi3<-0
      }
      
      # E-step
      xmcmc<-estep(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM)
      varList<-xmcmc$varList
      DYF<-xmcmc$DYF
      phiM<-xmcmc$phiM
      psiM<-transphi(phiM,saemix.model["transform.par"])
      etaM.l[[kiter]] = xmcmc$etaM ##** added
      
      
      # M-step
      if(opt$stepsize[kiter]>0) {
        ############# Stochastic Approximation
        xstoch<-mstep(kiter, Uargs, Dargs, opt, structural.model, DYF, phiM, varList, phi, betas, suffStat)
        varList<-xstoch$varList
        mean.phi<-xstoch$mean.phi
        mean.phi.l[[kiter]]<-mean.phi ##** added
        phi<-xstoch$phi
        phi.l[[kiter]]<-phi ##** added
        Omega.l[[kiter]] =varList$omega ##** added
        respar_vec[kiter,]=varList$pres
        
        betas<-xstoch$betas
        suffStat<-xstoch$suffStat
        
        beta.I<-betas[Uargs$indx.betaI]
        fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
        betaC<-betas[Uargs$indx.betaC]
        var.eta<-mydiag(varList$omega)
        l1<-betas.ini
        l1[Uargs$indx.betaI]<-fixed.psi
        l1[Uargs$indx.betaC]<-betaC
        
        if(Dargs$modeltype=="structural") {
          allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
        } else{
          allpar[(kiter+1),]<-c(l1,var.eta[Uargs$i1.omega2])
        }
        
      } else { #end of loop on if(opt$stepsize[kiter]>0)
        allpar[(kiter+1),]<-allpar[kiter,]
        phi.l[[kiter]]=phiM[1:saemix.data["N"],] ##*added ****revisar
        Omega.l[[kiter]] =varList$omega ##** added
        respar_vec[kiter,]=varList$pres ##** added
        mean.phi.l[[kiter]]=matrix(rep(allpar[1,1:9], saemix.data["N"]), ncol=9, byrow = T) ##** added
      }
      if(Dargs$modeltype=="structural") {
        theta<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
      } else{
        theta<-c(fixed.psi,var.eta[Uargs$i1.omega2])
      }
      
      # End of loop on kiter
    }
    
    
    etaM<-xmcmc$etaM # only need etaM here (re-created in estep otherwise)
    if(saemix.options$warnings) cat("\n    Minimisation finished\n")
    if(saemix.options$warnings) print(date())
    
    ############# After end of iterations
    fixed.effects<-0*betas
    fixed.effects[Uargs$indx.betaI]<-fixed.psi
    fixed.effects[Uargs$indx.betaC]<-betaC
    varList$omega[Uargs$i0.omega2,]<-0
    varList$omega[,Uargs$i0.omega2]<-0
    
    ##### Compute the individual parameters (MAP)
    phi[,Uargs$i0.omega2,1:saemix.options$nb.chains]<-mean.phi[,Uargs$i0.omega2]
    phi.samp<-phi
    phi<-apply(phi,c(1,2),mean)
    
    ##### Conditional means and variances used for the estimation of the log-likelihood via Importance Sampling
    cond.mean.phi<-phi
    sphi1<-phi
    sphi1[,varList$ind.eta]<-suffStat$statphi1
    cond.mean.phi[,Uargs$i1.omega2]<-sphi1[,Uargs$i1.omega2]
    cond.var.phi<-array(data=0,dim=dim(phi))
    cond.var.phi[,Uargs$i1.omega2]<-suffStat$statphi3[,Uargs$i1.omega2]-cond.mean.phi[,Uargs$i1.omega2]**2
    cond.mean.psi<-transphi(cond.mean.phi,saemixObject["model"]["transform.par"])
    
    cond.mean.eta<-matrix(0,nrow=dim(etaM)[1],ncol=Uargs$nb.parameters)
    cond.mean.eta[,varList$ind.eta]<-etaM
    cond.mean.eta<-array(t(cond.mean.eta),dim=c(Uargs$nb.parameters, Dargs$N, saemix.options$nb.chains))
    cond.mean.eta<-t(apply(cond.mean.eta,c(1,2),mean))
    shrinkage<-100*(1-apply(cond.mean.eta,2,var)/mydiag(varList$omega))
    names(shrinkage)<-paste("Sh.",saemixObject["model"]["name.modpar"],".%",sep="")
    
    # Updating objects
    saemix.model["Mcovariates"]<-Uargs$Mcovariates
    saemix.model["indx.res"]<-Uargs$ind.res
    saemix.model["indx.fix"]<-Uargs$indx.betaI
    saemix.model["indx.cov"]<-Uargs$indx.betaC
    saemix.model["indx.omega"]<-Uargs$i1.omega2
    
    # Filling in result object
    if(Dargs$modeltype=="structural") {
      saemix.res<-new(Class="SaemixRes",status="fitted",modeltype=Dargs$modeltype,name.fixed=saemix.model["name.fixed"],
                      name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"], 
                      fixed.effects=c(fixed.effects),fixed.psi=c(fixed.psi),betas=betas,betaC=betaC,
                      omega=varList$omega,respar=varList$pres,cond.mean.phi=cond.mean.phi,cond.var.phi=cond.var.phi,
                      mean.phi=mean.phi, phi=phi,phi.samp=phi.samp,parpop=parpop,allpar=allpar,MCOV=varList$MCOV)
      ##** added------------
      setClass("eval_conv", representation(comp = "list"))
      eval_conv_object <- new("eval_conv")
      eval_conv_object@comp <- list(allpar=allpar,mean.phi.l=mean.phi.l,phi.l=phi.l, Omega.l=Omega.l,etaM.l=etaM.l,respar=respar_vec,initials_val=initials_val)
      ##**-----------------
    } else{
      saemix.res<-new(Class="SaemixRes",status="fitted",modeltype=Dargs$modeltype,name.fixed=saemix.model["name.fixed"], 
                      name.random=saemix.model["name.random"],name.sigma=saemix.model["name.sigma"], fixed.effects=c(fixed.effects),
                      fixed.psi=c(fixed.psi),betas=betas,betaC=betaC, omega=varList$omega,cond.mean.phi=cond.mean.phi,
                      cond.var.phi=cond.var.phi, mean.phi=mean.phi, phi=phi,phi.samp=phi.samp,parpop=parpop,allpar=allpar,MCOV=varList$MCOV)
    }
    saemix.res["indx.res"]<-Uargs$ind.res
    saemix.res["indx.fix"]<-Uargs$indx.betaI
    saemix.res["indx.cov"]<-Uargs$indx.betaC
    saemix.res["indx.omega"]<-Uargs$i1.omega2
    saemix.res["npar.est"]<-Uargs$nb.parest
    saemix.res["nbeta.random"]<- sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(diag(saemix.model["covariance.model"])))
    saemix.res["nbeta.fixed"]<-  sum(saemix.model["betaest.model"]%*%diag(saemix.model["fixed.estim"])%*%as.matrix(-1*diag(saemix.model["covariance.model"])+1))
    saemix.res["cond.mean.psi"]<-cond.mean.psi
    saemix.res["cond.mean.eta"]<-cond.mean.eta
    saemix.res["cond.shrinkage"]<- shrinkage
    
    # Updating elements of saemixObject
    saemixObject["model"]<-saemix.model
    saemixObject["results"]<-saemix.res
    saemixObject["options"]<-saemix.options
    
    #  saemixObject["rep.data"]<-chdat # Utile ? maybe remove rep.data
    
    # ECO TODO check
    # a la fin: mais verifier, pe pb de distribution ??? ie allpar sur l'echelle des betas et pas parpop ? a verifier
    # saemix.res["allpar"]<-allpar
    # saemix.res["parpop"]<-allpar[,-c(indx.betaC)]
    #### Final computations
    # Compute the MAP estimates of the PSI_i's 
    if(saemix.options$map) {
      x<-try(saemixObject<-map.saemix(saemixObject))
      if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the MAP parameters\n") else {
        if(!saemix.options$save.graphs) saemixObject<-saemix.predict(saemixObject, type=c("ipred","ppred")) # if no graphs, compute predictions all the same
      }
    }
    
    # Compute the Fisher Information Matrix & update saemix.res
    if(saemix.options$fim) {
      x<-try(saemixObject<-fim.saemix(saemixObject))
      if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the FIM\n")
    }
    
    # Estimate the log-likelihood via importance Sampling/Gaussian quadrature
    if(saemix.options$ll.is) {
      x<-try(saemixObject<-llis.saemix(saemixObject))
      if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the likelihood by IS\n")
    }
    if(saemix.options$ll.gq) {
      x<-try(saemixObject<-llgq.saemix(saemixObject))
      if(inherits(x,"try-error") & saemixObject@options$warnings) message("Problem estimating the likelihood by GQ\n")
    }
    
    #### Pretty printing the results (TODO finish in particular cov2cor)
    if(saemix.options$print) print(saemixObject,digits=2)
    
    #### Save the results to a file
    if(saemix.options$save | saemix.options$save.graphs) {
      # create directory to save the results
      if(saemix.options$directory!="") xsave<-dir.create(saemix.options$directory) else xsave<-TRUE
      if(!xsave) {
        # Check that we're not trying to create a directory with the same name as a file
        if(!file_test("-d",saemix.options$directory)) {
          if(saemix.options$warnings) message("Unable to create directory",saemix.options$directory)
          saemix.options$directory<-"newdir"
          dir.create(saemix.options$directory)         
          xsave<-file_test("-d",saemix.options$directory)
          if(!xsave) {
            saemix.options$directory<-""
            xsave<-TRUE
            if(saemix.options$warnings) message(", saving in current directory.\n")
          } else {if(saemix.options$warnings) message(", saving results in newdir instead.\n")}
        } else {
          xsave<-TRUE
          if(saemix.options$warnings) message("Overwriting files in directory",saemix.options$directory,"\n")
        }
      }
    }
    if(saemix.options$save) {
      namres<-ifelse(saemix.options$directory=="","pop_parameters.txt", file.path(saemix.options$directory,"pop_parameters.txt"))
      xtry<-try(sink(namres))
      if(!inherits(xtry,"try-error")) {
        print(saemixObject)
        sink()
        namres<-ifelse(saemix.options$directory=="","indiv_parameters.txt", file.path(saemix.options$directory,"indiv_parameters.txt"))
        if(length(saemixObject["results"]["map.psi"])>0) {
          tab<-cbind(Id=unique(saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]),saemixObject["results"]["map.psi"])
          colnames(tab)[1]<-saemixObject["data"]["name.group"]
          write.table(tab,namres,quote=FALSE, row.names=FALSE)
        }
      } else {
        message("Unable to save results, check writing permissions and/or path to directory.\n")
      }
    }
    # ECO TODO finish, adding all
    if(saemix.options$save.graphs) {
      saemixObject<-saemix.predict(saemixObject)
      if(saemix.options$directory=="") namgr<-"diagnostic_graphs.ps" else
        namgr<-file.path(saemix.options$directory,"diagnostic_graphs.ps")
      xtry<-try(postscript(namgr,horizontal=TRUE))
      if(!inherits(xtry,"try-error")) {
        oldpar <- par(no.readonly = TRUE)    # code line i
        on.exit(par(oldpar))            # code line i + 1 
        par(mfrow=c(1,1))
        try(plot(saemixObject,plot.type="data"))
        
        try(plot(saemixObject,plot.type="convergence"))
        
        if(length(saemixObject["results"]["ll.is"])>0) {
          par(mfrow=c(1,1))
          try(plot(saemixObject, plot.type="likelihood"))
        }
        
        try(plot(saemixObject,plot.type="observations.vs.predictions"))
        
        try(plot(saemixObject,plot.type="random.effects"))
        
        try(plot(saemixObject,plot.type="correlations"))
        
        # Note: can replace all this by:
        #    default.saemix.plots(saemixObject)
        
        dev.off()
        
        if(saemix.options$directory=="") namgr<-"individual_fits.ps" else
          namgr<-file.path(saemix.options$directory,"individual_fits.ps")
        postscript(namgr,horizontal=FALSE)
        try(plot(saemixObject,plot.type="individual.fit"))
        dev.off()
      } else {
        message("Unable to save results, check writing permissions and/or path to directory.\n")
      }
    }
    
    options(warn=opt.warn)
    return(list(saemixObject=saemixObject,eval_conv_object=eval_conv_object)) ##** modified
  }
  
  
  #--------------------------------------------------------
  ######  fit saem 
  #--------------------------------------------------------
  knots <-  quantile(x, (1:(df - 1)) / df)
  bounds <- range(x) + abs(0.04) * c(-1, 1) * diff(range(x))
  fulldata <- data.frame(x, y, id)
  
  ##### SAEM COMPONENTS
  
  start.b<-start
  
  transform.par=rep(0,df+3)
  fixed.estim=rep(1,df+3)
  covariance.model=matrix(c(rep(c(rep(1,3),rep(0,df)),3), rep(rep(0,df+3),df)),ncol=df+3,byrow=TRUE)
  omega.init=diag(1,(df+3))
  saemix.model.2<-sitar_saemixModel(df=df, data=fulldata,start=start.b,transform.par,fixed.estim,covariance.model,omega.init)
  

  
  #### Data_saemix
  pb.data <- saemixData(name.data = fulldata, header = TRUE,
                        sep = " ", na = NA, name.group = c("id"),
                        name.predictors = c("x"), name.response = c("y"),
                        units = list(x = "years", y = "cm"),
                        name.X = "x", verbose = F)
  
  #### Options saemix function
  opt <- list(seed = seed,print=FALSE,save = FALSE, save.graphs = FALSE, nbiter.saemix = c(k1_initial,10))
  saem_models<-try(saemix.mod(model=saemix.model.2,data=pb.data,control =opt), silent=T)
  
  
  #--------------------------------------------------------
  ########### fin fit saem 
  #--------------------------------------------------------
  
  saemixObject= saem_models$saemixObject
  random_eff_val_it= saem_models$eval_conv_object@comp
  
  saemix.model <- saemixObject["model"]
  saemix.data <- saemixObject["data"]
  saemix.res <- saemixObject["results"]
  yobs <- saemix.data["data"][, saemix.data["name.response"]]
  saemix.options <- saemixObject["options"]
  npred <- length(saemix.data["name.predictors"])
  xind <- saemix.data["data"][, c(2,3,5,6,7,8), drop = FALSE]
  i1.omega2 <- saemix.model["indx.omega"]
  
  # Omega <- saemix.res["omega"] ## Matrix varaince and covariance random effects
  # pres <- saemix.res["respar"]  ## sigma estimation
  ## phi: Individual parameter, it is defined as  Phi= C_imu+ eta_i ;  
  # C_i  vector of covariance, mu is an unknown vector of fixed effects of size;
  # eta_i is an unknown vector of normally distributed random effects
  
  mean.phi <- random_eff_val_it$mean.phi.l## Estimate of the corresponding mean individual φ 
  cond.mean.phi <-random_eff_val_it$phi.l## Estimate of the corresponding individual φ 
  etaM <-random_eff_val_it$etaM.l## Estimate of the corresponding random effects
  omegaM <-random_eff_val_it$Omega.l## Estimate of the corresponding Omega matrix. Variance and covariance matrix RE
  respar <-random_eff_val_it$respar ## sigma estimation
  #--------------------------------------------------------
  ########### Get log-verosimilitud completa
  #--------------------------------------------------------
  
  LL_c=c()
  fix_param=list()
  
  # Initialize iteration counter
  it <- 1
  
  # Extract parameters for readability and efficiency
  nbiter_burn <- saemix.options$nbiter.burn
  nbiter_tot <- saemix.options$nbiter.tot
  N <- saemix.data["N"]
  nind_obs <- saemix.data["nind.obs"]
  yM <- yobs
  XM <- xind
  nphi1 <- length(i1.omega2)
  # Precompute indices for observations
  io <- matrix(0, nrow = N, ncol = max(nind_obs))
  for (i_suj in 1:N) {
    io[i_suj, 1:nind_obs[i_suj]] <- 1
  }
  ind_ioM <- which(t(io) != 0)
  
  # Main loop
  for (k in (nbiter_burn + 1):nbiter_tot) {
    # Store fixed parameters
    fix_param[[it]] <- data.frame(
      fixed.effects = c("a", "b", "c", paste0("s", 1:df)),
      value = saemix.res["allpar"][k + 1, ][1:(df + 3)],
      k = k
    )
    
    # Extract current parameters
    Omega <- omegaM[[k]]
    pres <- respar[k, ]
    IOmega_phi1 <- solve(Omega[i1.omega2, i1.omega2])
    mean_phi1 <- mean.phi[[k]][, i1.omega2]
    mtild_phiM1 <- cond.mean.phi[[k]][, i1.omega2, 1]
    dphiM <- mtild_phiM1 - mean_phi1
    
    # Compute e2 vectorized
    log_det_Omega <- log(det(Omega[i1.omega2, i1.omega2, drop = FALSE]))
    e2 <- -0.5 * rowSums((dphiM %*% IOmega_phi1) * dphiM) - 0.5 * (log_det_Omega + nphi1 * log(2 * pi))
    
    # Compute predicted values and residuals
    phiM <- cond.mean.phi[[k]][, i1.omega2, 1]
    fit_y <- lapply(1:N, function(i) {
      tim <- XM[XM$id == i, c("x", "ytype")]
      ypred <- phiM[i, 1] + drop(ns((tim$x - phiM[i, 2]) * exp(phiM[i, 3]), k = knots, B = bounds) %*% mean.phi[[k]][i, 4:(df + 3)])
      data.frame(ypred = ypred, id = i, x = tim$x, ytype = tim$ytype)
    })
    fit_y <- do.call(rbind, fit_y)
    
    # Compute residuals and likelihood
    g <- error(fit_y$ypred, pres, XM$ytype)
    DYF <- matrix(0, nrow = max(nind_obs), ncol = N)
    DYF[ind_ioM] <- -0.5 * ((yM - fit_y$ypred) / g)^2 - log(g) - 0.5 * log(2 * pi)
    e1 <- colSums(DYF)
    
    # Compute complete log-likelihood
    LL_c[it] <- sum(e1 + e2)
    fix_param[[it]]$LL_c <- LL_c[it]
    
    it <- it + 1
  }
  
  ############################################################################################
  #  resultados log-verosimilitud completa y criterio de parada para seleccion de k1
  ############################################################################################
  
  fix_param_it <- do.call(rbind, fix_param)   
  LL_values= fix_param_it %>% group_by(it=k) %>% summarise(Comp_loglik=-2*unique(LL_c))
  
  LL_values_burn= LL_values %>%filter(it>=round(saemix.options$nbiter.tot*0.1) ) %>%
    mutate(Comp_loglik = as.numeric(Comp_loglik))
  
  
  LL_values_burn <- LL_values_burn %>%
    mutate(
      dif_abs = abs(diff(c(NA, Comp_loglik))),
      dif_rel = dif_abs / abs(Comp_loglik),
      it= round(saemix.options$nbiter.tot*0.13)+(1:dim(LL_values_burn)[1])
    ) 
  
  # Initialize vectors to store windowed means
  window_abs <- rep(NA_real_, nrow(LL_values_burn))
  window_rel <- rep(NA_real_, nrow(LL_values_burn))
  
  # Calculate windowed means
  for (i in seq_len(nrow(LL_values_burn) - m + 1)) {
    window_abs[i + m - 1] <- mean(LL_values_burn$dif_abs[i:(i + m - 1)], na.rm = TRUE)
    window_rel[i + m - 1] <- mean(LL_values_burn$dif_rel[i:(i + m - 1)], na.rm = TRUE)
  }
  
  # Add windowed means to the data frame
  LL_values_burn <- LL_values_burn %>%
    mutate(
      window_abs = window_abs,
      window_rel = window_rel
    )
  
  # Select the iteration where both absolute and relative differences are below thresholds
  selected_iteration_rel <- LL_values_burn %>%
    filter(window_rel < error_rel) %>%
    arrange(it) %>% 
    slice(1) 
  
  
  if(dim(selected_iteration_rel)[1]==0){
    selected_iteration_rel <- LL_values_burn %>%na.omit()%>%
      filter(dif_rel < error_abs) %>%
      arrange(it) %>% 
      slice(1) 
    
    if(dim(selected_iteration_rel)[1]==0){
      
      selected_iteration_rel <- LL_values_burn%>%na.omit()%>%
        slice_tail(n = 1)
      
    }
    
  }
  
  
  selected_iteration_abs <- LL_values_burn %>%
    filter(window_abs < error_abs) %>%
    arrange(it) %>% 
    slice(1) 
  
  if(dim(selected_iteration_abs)[1]==0){
    selected_iteration_abs <- LL_values_burn %>%na.omit()%>%
      filter(dif_abs < error_abs) %>%
      arrange(it) %>% 
      slice(1) 
    
    if(dim(selected_iteration_abs)[1]==0){
      
      selected_iteration_abs <- LL_values_burn%>%na.omit()%>%
        slice_tail(n = 1)
      
    }
    
  }
  
  mean_Comp_loglik= mean(LL_values_burn$Comp_loglik)
  sd_Comp_loglik= sd(LL_values_burn$Comp_loglik)
  
  # compl_log_lik= LL_values %>% ggplot(aes(x=it, y=Comp_loglik))+geom_line()+
  #   geom_hline(yintercept =mean_Comp_loglik, col="red" )+theme_bw()
  
  
  ############## plots
  
  compl_log_lik1= LL_values %>% ggplot(aes(x=it, y=Comp_loglik))+geom_line()+labs(y="-2LL_c", x="kiter")+
    geom_hline(yintercept =mean_Comp_loglik, col="red" )+theme_bw()+
    geom_vline(xintercept=selected_iteration_abs$it, col="blue")+
    geom_vline(xintercept=selected_iteration_rel$it, col="green")
  
  # compl_log_lik2= LL_values %>% filter(it%in%(iter_min_SE.abs-200):(iter_min_SE.abs+200))%>% ggplot(aes(x=it, y=Comp_loglik))+
  #   geom_line()+labs(y="-2LL_c", x="kiter")+
  #   geom_hline(yintercept =mean_Comp_loglik, col="red" )+theme_bw()+
  #   geom_vline(xintercept=iter_min_SE.abs, col="blue")
  # 
  # plot_all=cowplot::plot_grid(compl_log_lik1,compl_log_lik2, 
  #                    LL_values%>%filter(it>100) %>%  ggplot( aes(y = Comp_loglik, x="")) +
  #                      geom_boxplot()+theme_bw()+labs(y="-2LL_c"),labels = "AUTO", ncol=3)
  # 
  
  results_stop = list()
  results_stop$LL_values=LL_values
  results_stop$sel_inf=rbind(selected_iteration_abs,selected_iteration_rel)
  results_stop$k1_stop =  data.frame(k1_abs=selected_iteration_abs$it, 
                                     k1_rel=selected_iteration_rel$it)
  results_stop$plots= compl_log_lik1
  
  return(results_stop)
}


#########################################
# sitar-saem function
#########################################

sitar_saemixModel<-function(df,data,
                            transform.par,
                            fixed.estim,
                            covariance.model=NULL,
                            omega.init=NULL,
                             verbose=F, 
                            start=NULL){
  
  x <- data$x
  y <- data$y
  id <- as.factor(data$id)
  # Define knots and bounds
  knots <-  quantile(x, (1:(df - 1)) / df)
  bounds <- range(x) + abs(0.04) * c(-1, 1) * diff(range(x))
  
  # Adjust x and knots for centering
  # xoffset <- mean(x)
  # bstart <- mean(x) - xoffset
  # x <- x - xoffset
  # knots <- knots - xoffset
  # bounds <- bounds - xoffset
  
  fulldata <- data.frame(x, y, id)
  
  #######################
  # Valores iniciales lm
  #######################
  if(is.null(start)){

  spline.lm <- lm(y ~ ns(x, knots = knots, Bound = bounds))
  start.all <- coef(spline.lm)[c(2:(df + 1), 1)]              ### start values for ns parameters and a
  start_var=diag(vcov(spline.lm))
  start_var = c(start_var[1],0.5,0.5,start_var[-1])
  start.all<-c(start.all, 0,0)
  start.all=start.all[c((df+1):(df+3),1:(df))] 
  if(sum(is.na(start.all))!=0)start.all=c(mean(y),0,0,rep(runif(1,20,50),df))
  
  names(start.all) = names(start_var) = c("a","b","c",paste0("s",1:df))
  
  start_df=data.frame(start=start.all, start_var=start_var)
  
  start= start.all

  }  
  
  # 
  # start<-start[c((df+1):(df+3),1:df)]
  # transform.par=rep(0,df+3)
  # fixed.estim=rep(1,df+3)
  # covariance.model=matrix(c(rep(c(rep(1,3),rep(0,df)),3), rep(rep(0,df+3),df)),ncol=df+3,byrow=TRUE)
  # omega.init=matrix(c(1,0,0,rep(0,df),0,1,0,rep(0,df),0,0,1,rep(0,df), rep(rep(0,df+3),df)),ncol=df+3,byrow=TRUE)
  # 
  #---------------------------------
  ### function replace
  #---------------------------------
  
  ss1<-paste0("'s",1:df, "'")
  pars.saem<-paste(c("'a'","'b'","'c'",ss1), collapse = ',')
  numb<-seq(4,4+df-1,1)
  ss<-paste0("s",1:df)
  ss_v<-glue:::glue_collapse(glue('{ss} <-psi[id,{numb}]'), " \n")
  
  if(is.null(omega.init)){
  
  om1=diag(start_df$start_var)
  om1[1,2] = om1[2,1] = om1[2,3]= om1[3,2] = om1[1,3] = om1[3,1]=0.1
  omega.init=om1

  }
  #---------------------------------
  ## saemix_glue
  #---------------------------------
  print(paste0("glue saem",df))

  spline_prod2<-glue( 'rowSums((cbind(', paste(ss, collapse = ','),')*ns((tim-b)*exp(c),k=knots,B=bounds)))')
  
  saem.fitcode_p<-glue(
    "model.sm<-function(psi,id,xidep) {\n",
    "\n",
    "tim<-xidep[,1] \n",
    "\n",
    "a<-psi[id,1]\n",
    "b<-psi[id,2]\n",
    "c<-psi[id,3]\n",
    "\n",
    "<<ss_v>>\n",
    "\n",
    "a+<<spline_prod2>>}\n",
    "\n",
    "\n",
    "saemix.model<-saemixModel(model=model.sm,description='sitar_saem',\n", 
    "psi0=matrix(start,ncol=<<df+3>>,byrow=TRUE,dimnames=list(NULL,c(<<pars.saem>>))),\n", 
    "transform.par=transform.par,fixed.estim=fixed.estim, covariance.model=covariance.model, omega.init=omega.init,verbose = verbose)",
    .open = "<<",  .close = ">>"
  )
  
  saem.out <- eval(parse(text = saem.fitcode_p))
  
  attr(saem.out,"saemix.model") <- parse(text = saem.fitcode_p)
  attr(saem.out,"start") <- start
  
  return(saem.out)
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
# predict value function
#########################################


predict_saem<-function(newdata, xoffset,object,knots, bounds, df){
  
  idcoef<-object@results@map.psi[,1:3]
  idcoef<-data.frame(idcoef,id=as.factor(unique(object@data@data$id)))
  colnames(idcoef)<-c("a","b","c", "id")
  
  newdata<-merge(newdata,idcoef,by="id",all.x=TRUE)
  
  newdata$xnew=(newdata$x-xoffset-newdata$b)*(exp(newdata$c))
  newdata$xnew.pop=(newdata$x-xoffset-object@results@fixed.effects[2])*(exp(object@results@fixed.effects[3]))
  
  newns<-ns(newdata$xnew,knots = knots, Bound = bounds)
  newns.pop<-ns(newdata$xnew.pop,knots = knots, Bound = bounds)
  
  newdata$ipred<-newdata$a+newns%*%cbind(object@results@fixed.effects[4:(df+3)])
  newdata$ppred<-object@results@fixed.effects[1]+newns.pop%*%cbind(object@results@fixed.effects[4:(df+3)])
  newdata<- newdata[,c("id", "x", "ipred", "ppred")]
  
  return(newdata)
}


#########################################
# deriv_sitar_saemixModel function
#########################################


deriv_sitar_saemixModel  <-function(object, xoffset ,df, n.smooth=NULL,smooth.deriv=F,newdata=NULL){
  
  data<-object@data@data  
  x <- data$x
  y <- data$y
  id <- as.factor(data$id)
  # Define knots and bounds
  knots <-  quantile(x, (1:(df - 1)) / df)
  bounds <- range(x) + abs(0.04) * c(-1, 1) * diff(range(x))
  
  fulldata <- data.frame(x, y, id)
  
  ### calculate the a+a_i, b+b_i, c+c_i
  if(class(object)[1]!= "try-error"){
    idcoef<-object@results@map.psi[,1:3]
    idcoef<-data.frame(idcoef,id=as.factor(unique(data$id)))
    colnames(idcoef)<-c("a","b","c", "id")
    
    all_data<-merge(data, idcoef, by="id") %>% arrange(id, x)
    
    ################ deriv
    all_data<-all_data %>% mutate(z=(x-b)*(exp(c)), zpop=(x-object@results@fixed.effects[2])*(exp(object@results@fixed.effects[3])))
    
    dns_matrix<-JM::dns(all_data$z,knots = knots,Boundary.knots = bounds)
    dns_matrix.pop<-JM::dns(all_data$zpop,knots = knots,Boundary.knots = bounds)
    
    all_data$deriv=dns_matrix%*%cbind(object@results@fixed.effects[4:(df+3)])*(exp(all_data$c))
    all_data$derivpop=dns_matrix.pop%*%cbind(object@results@fixed.effects[4:(df+3)])*(exp(object@results@fixed.effects[3]))
    
    all_data<-all_data %>% mutate(x=x+xoffset)%>% dplyr::select(id,x,deriv,derivpop)
   
    
    if(smooth.deriv){
      
      if(is.null(n.smooth)) n.smooth=50
      if(is.null(newdata)){
        newdata1=data.frame(x=seq((min(data$x)+xoffset), (max(data$x)+xoffset), l=100), id=rep(unique(data$id), each=100))
      }else{ 
        newdata1=newdata
      }
      
      dat_pred=predict_saem(newdata=newdata1, xoffset,object,knots, bounds, df)
      dat_pred=dat_pred %>% mutate(x=x-xoffset) %>% rename(y=ipred)
      # dat_pred %>% ggplot(aes(x=x, y=y, group=id))+geom_line()
     
      
      idcoef<-object@results@map.psi[,1:3]
      idcoef<-data.frame(idcoef,id=as.factor(unique(data$id)))
      colnames(idcoef)<-c("a","b","c", "id")
      
      all_data_pred<-merge(dat_pred, idcoef, by="id") %>% arrange(id, x)
      
      ################ deriv
      all_data_pred<-all_data_pred %>% mutate(z=(x-b)*(exp(c)), zpop=(x-object@results@fixed.effects[2])*(exp(object@results@fixed.effects[3])))
      
      dns_matrix<-JM::dns(all_data_pred$z,knots = knots,Boundary.knots = bounds)
      dns_matrix.pop<-JM::dns(all_data_pred$zpop,knots = knots,Boundary.knots = bounds)
      
      all_data_pred$deriv=dns_matrix%*%cbind(object@results@fixed.effects[4:(df+3)])*(exp(all_data_pred$c))
      all_data_pred$derivpop=dns_matrix.pop%*%cbind(object@results@fixed.effects[4:(df+3)])*(exp(object@results@fixed.effects[3]))
      
      all_data_pred<-all_data_pred %>% mutate(age=x+xoffset)%>% dplyr::select(id,age,deriv,derivpop)
     
      

      deriv.smooth <- with(all_data_pred, 
                           by(all_data_pred, id, function(z) {
                             spl <- spline(z$age, y = z$deriv, n = n.smooth)
                             data.frame(x = spl$x, y = spl$y, id = z$id[1])  # Use z$id[1] instead of rep()
                           }, simplify = FALSE)  # Keep as a list
      )
      
   
      deriv.smooth <- do.call(rbind, deriv.smooth)
      

      # deriv.smooth_all %>% filter(id==83) %>% ggplot(aes(x=x, y=y, group=id)) + geom_line()
      
      return(list(all_data=all_data, deriv_smooth=deriv.smooth))
    }
  }
  else{
    all_data= data.frame( id="NA",age  ="NA", deriv="NA",derivpop="NA")
    
  }
  
  return(all_data)
  
}

