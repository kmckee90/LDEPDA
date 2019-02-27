
#--------------------------------------------------------------------------
#A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering
#Pt. 2a: Random-Effects Simulation, Diffusion Processes
#Author: Kevin McKee (mckeek@vcu.edu)
#Virginia Commonwealth University
#Virginia Institute of Psychiatric and Behavioral Genetics
#Last updated: 2/27/2019
#--------------------------------------------------------------------------

library(OpenMx)
library(signal)
library(deSolve)
library(DEoptim)
source("GLLAfunctions.R")

#Helper functions
smear<-function(x, t){
  x.s<-matrix(NA, nrow=length(x), ncol=t)
  x.s[,1]<-x
  for(i in 2:t){
    x.s[i:nrow(x.s),i]<-x[1:(length(x)-i+1)]
  }
  x.s.combine<-apply(x.s, 1, max, na.rm=T)
  return(x.s.combine)
}
genDifs<-function(x, kern){
  kern.mat<-matrix(kern, nrow=length(kern),ncol=1)
  a<-gllaEmbed(x, embed=length(kern), tau=1, idColumn = F)
  b<- a%*%kern.mat
  return(as.matrix(b))
}
#Outlier detection and removal algorithm
PDA<-function(model, data, thres.sd=1,AR.order, max.iter=10, plot=F){
  #Model variables
  tData<-as.matrix(data)
  D<-model$D$values[1]
  data.cut<-tData
  
  #include iteration count for max.iter stopping condition and minimum remaining rows of data.
  i<-1
  min.rows<-5
  
  #Initial values for max range stopping condition.
  iRange<-1
  iRange.prior<-0
  estEvent<-NULL
  par(mfrow=c(1,1), mai=c(.5,.5,.5,.5))
  
  while(iRange>iRange.prior & i<max.iter & nrow(tData)-length(estEvent) > min.rows) {
    
    model<-mxModel(model, mxData(data.cut, type="raw"))
    model<-mxRun(model, silent=T)
    if(model$output$status$code!=0) model<-mxTryHard(model, silent=T, extraTries = 15)
    
    m<-mahalanobis(tData[,1:D], center=model$fitfunction$info$expMean, cov=model$fitfunction$info$expCov)
    m.ar<-ar(m, order.max=AR.order,aic=FALSE, method="ols")$ar
    if(length(m.ar)==0) m.ar<-0
    m.d<- genDifs(m, c(-m.ar,1))
    
    m.d.thres<-as.numeric(mean(m.d)+thres.sd*sd(m.d))
    m.d.s<-matrix(c(smear(m.d, t=(D-1))))

    data.cut<-tData[m.d.s<m.d.thres, ]
    estEvent<- which(m.d.s>m.d.thres)
    iRange.prior<-iRange
    iRange<-range(m.d.s)[2]-range(m.d.s)[1]
    i<-i+1
    
    if(plot){
      plot(m.d.s, type="l")
      abline(h=m.d.thres)
      abline(v=estEvent, lty=3, col="blue")
    }
  }
  # Output ------------------------------------------------------------------
  return(list("fit"=model, "exclusions"=estEvent))
}

# Generate shot noise process ---------------------------------------------
DLOsim <- function(t, prevState, parms) {
  x <- prevState[1] 
  y <- prevState[2]

  with(as.list(parms), 
       {
         dx <- y
         dy <- parms[[1]]*x + parms[[2]]*y
         res<-c(dx, dy)
         list(res)
       }
  )
}

genSN<-function(pars, eventFreq=.1, eventType="level", iv=c(1,0)){
  # Generate events ---------------------------------------------------------
  eta<- pars[4]#-(2*pi/pars[[4]])^2
  zeta<-pars[5]#log(pars[[5]])*sqrt(-eta)/pi
  pars.conv<-c(eta, zeta)    
  t<-1:pars[[1]]
  
  shocks<-list()
  shocks<-rbinom(pars[[1]],1,eventFreq)*rnorm(pars[[1]]) #sample(c(-1, 1), 100, replace=T)#

  shocks.m<-NULL
  if(eventType=="slope"){
    shocks.m<-c(shocks.m, rep(0,pars[[1]]), shocks)
  }else{
    shocks.m<-c(shocks.m, shocks, rep(0,pars[[1]]))
  }
  eventdat <- data.frame(var = c(rep("x",pars[[1]]),rep("y",pars[[1]])),
                         time = rep(t, 2),
                         value = shocks.m,
                         method = rep("add", 2))
  eventdat<-eventdat[eventdat$value!=0,]
  names(iv)<-c("x","y")
  dat<-lsoda(y=iv, times=t, DLOsim, pars.conv, events = list(data=eventdat))
  x<-cbind(t, dat[,2]+rnorm(pars[[1]])*pars[[2]])
  
#  plot(x, type="l", col='darkgrey')
#  lines(dat[,1:2])
#  lines(eventdat[,2:3], type="h", col="red", lwd=2)

    return(x)
}

nIter<-200000

#---------------------------------------------------------------------------
# Simulation 1--------------------------------------------------------------
#---------------------------------------------------------------------------
  randEffects<-data.frame(
    "N"=sample(80:400, nIter, replace=T),
    "measErr"=runif(nIter, 1/32, 1/2),
    "dynErr"=.3,#runif(nIter, 1/32, 1/2),
    "Period"= -(2*pi/sample(8:40, nIter, replace=T))^2,
    "Damping"=0,
    "p"=runif(nIter, 0, 1),
    "type"=1#sample(c(1,2),nIter, replace=T)
  )
  randEffects$Damping<-log(runif(nIter, 0.05, .95))*sqrt(-randEffects$Period)/pi

cl<-makeCluster(4)

library(foreach)
library(doSNOW)

  registerDoSNOW(cl)
  
  #####
  continue<-F #If false, create new data set. If true, append existing set.
  #####
  
  foreach(i=1:nIter,
          .combine = rbind,
          .packages = c("OpenMx","deSolve"),
          .export=c("DLOsim","genSN","gllaEmbed","PDA","smear","randEffects", "genDifs"))%dopar%{
    
    parSet<-randEffects[i,]
  
    dat<-genSN(parSet[1:5], eventFreq = parSet$p, eventType = ifelse(parSet$type==1, "level","slope"), iv=c(1,0))[,2]
    
    
    #Embed data for the LDE
    D<-floor(round(2*pi/sqrt(-parSet$Period))/2)
    tEmbedded <- gllaEmbed(dat, embed=D, tau=1, idColumn=FALSE)
    colnames(tEmbedded)<-paste0("x1_",0:(D-1))
    tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))
    Lmat<-matrix((1:D)-mean(1:D),nrow=D, ncol=3)
    for(k in 1:3) Lmat[,k]<-Lmat[,k]^(k-1)/factorial(k-1)
    
    
    
  # LDE ---------------------------------------------------------------------
    LDE<-readRDS("LDE_model.RDS")
    LDE<-mxModel(LDE, 
                  mxMatrix("Full", nrow=D, 3, values=Lmat, name="L", dimnames=list(colnames(tEmbedded)[1:D], rownames(LDE$S$values))),
                  mxMatrix("Diag", nrow=D, ncol=D, values=.1, free=T, labels="e_x1", lbound=1e-10, name="U", dimnames=list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])),
                  mxMatrix("Full", nrow=1, ncol=D, values=.001, free=T, labels="mu_x1", name="M", dimnames=list(NULL,colnames(tEmbedded)[1:D])),
                  mxExpectationNormal("R","M"))
    dimnames(LDE$R)<-list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])
    LDE$D$values<-D
    LDEPDA<-LDE
    LDE<-mxModel(LDE, mxData(tEmbedded, type="raw"))
    LDE$compute<-omxDefaultComputePlan()
    LDE$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
    LDE<-mxRun(LDE)
    if(LDE$output$status$code!=0) LDE<-mxTryHard(LDE, silent=T, extraTries = 15)
    
  # LDE EA ------------------------------------------------------------------
    LDEPDA<-omxSetParameters(LDEPDA, "d2x1ts1", values=0,free=F,lbound=NA)
    LDEPDA$compute<-omxDefaultComputePlan()
    LDEPDA$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
    
    phi<- -0.152344    +0.004885*D    +0.002200*length(dat) #NEW: VARIANCE AND BIAS (QUICK CALIBRATION): Obj value 0.317456
    AR<-floor(D*sqrt(parSet$p))
    LDEPDA.fit<-PDA(LDEPDA, tEmbedded, thres.sd=phi,AR.order=AR, max.iter=25, plot=F)
  
  # SSM ---------------------------------------------------------------------
    ssm<-readRDS("SSM_model.RDS")
    ssm<-mxModel(ssm, mxData(data.frame("obs"=dat, "tim"=1:parSet$N), type="raw"))
    ssm$compute<-omxDefaultComputePlan()
    ssm$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
    ssm<-mxRun(ssm, silent=T)
    if(ssm$output$status$code!=0) ssm<-mxTryHard(ssm, silent=T, extraTries = 15)
   # Output ------------------------------------------------------------------
  out<-unlist(c(parSet, "LDE"=coef(LDE)[1:2], "LDEPDA"=coef(LDEPDA.fit$fit)[1:2], "SSM"=coef(ssm)[1:2],
  "SC.LDE"=LDE$output$status$code, "SC.LDE.OR"=LDEPDA.fit$fit$output$status$code, "SC.SSM"=ssm$output$status$code,
  "SE.LDE"=LDE$output$standardErrors[1:2,1], "SE.LDE.OR"=LDEPDA.fit$fit$output$standardErrors[1:2,1], "SE.SSM"=ssm$output$standardErrors[1:2,1])) 
  
  write.table(t(out), sep=",", file="outputShotNoise_Rand.csv",append=ifelse(i==1, continue, T), col.names=ifelse(i==1, !continue, F), row.names=F)
  }
stopCluster(cl)



