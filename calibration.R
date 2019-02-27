
#--------------------------------------------------------------------------
#A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering
#Pt. 1: Calibration of the Piecewise-Deterministic Algorithm (PDA)
#Author: Kevin McKee (mckeek@vcu.edu)
#Virginia Commonwealth University
#Virginia Institute of Psychiatric and Behavioral Genetics
#Last updated: 2/27/2019
#--------------------------------------------------------------------------

library(DEoptim)
library(OpenMx)
library(deSolve)
source("GLLAfunctions.R")


# Functions ---------------------------------------------------------------

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

#Generate successive differences using a kernel vector
genDifs<-function(x, kern){
  kern.mat<-matrix(kern, nrow=length(kern),ncol=1)
  a<-gllaEmbed(x, embed=length(kern), tau=1, idColumn = F)
  b<- a%*%kern.mat
  return(as.matrix(b))
}

#Outlier detection and removal algorithm
PDA<-function(model, data, thres.sd=1, AR.order=2,  max.iter=10, plot=F){
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
    if(model$output$status$code!=0) model<-mxTryHard(model, silent=T, extraTries = 5)
    
    m<-mahalanobis(tData[,1:D], center=model$fitfunction$info$expMean, cov=model$fitfunction$info$expCov)
    
    m.ar<-ar(m, order.max=AR.order,aic=FALSE, method="ols")$ar
    if(length(m.ar)==0) m.ar<-0
    m.d<- genDifs(m, c(-m.ar,1))
    
    m.d.thres<-mean(m.d)+thres.sd*sd(m.d)
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
  eta<- pars[3]
  # eta<- -(2*pi/pars[[3]])^2
  zeta<-pars[4]
  # zeta<- log(pars[[4]])*sqrt(-eta)/pi
  pars.conv<-c(eta, zeta)    
  t<-1:pars[[1]]
  
  shocks<-list()
  shocks<-rbinom(pars[[1]],1,eventFreq)*rnorm(pars[[1]])
  
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
  
  plot(x, type="l", col='darkgrey')
  lines(dat[,1:2])
  lines(eventdat[,2:3], type="h", col="red", lwd=2)
  
  return(x)
}



controlPars<-c("CONTROL.B1"= .15, "CONTROL.B2"= .001, "CONTROL.B3"= .01)

calSet<-list()
l<-0
calSet.par<-NULL
#Generate calibration data set
for(N in c(100, 200, 300))
  for(measErr in c(1/8))
    for(Period in -(2*pi/c(8, 20, 32))^2)
      for(Damping in log(c(.25, .75))*sqrt(-Period)/pi)
        for(p in c(.125, .25, .5))
          for(wType in 1){
            parSet<-data.frame(
              "N"=N,
              "measErr"=measErr,
              "Period"= Period,
              "Damping"=Damping,
              "p"=p,
              "type"=wType
            )
            for(z in 1:3){
              l<-l+1
              cat(l,"\r")
              calSet.par<-rbind(calSet.par, parSet)
              calSet[[l]]<-genSN(parSet[1:5], eventFreq = parSet$p, eventType = ifelse(parSet$type==1, "level","slope"), iv=c(1,0))[,2]
            }          
          }




# Model-fitting function: Take parameters values, fit to the data, output parameter error function --------
simIter<-function(par, trace=F){
  out<-matrix(0, nrow=length(calSet), ncol=1)
  for(i in 1:length(calSet)){
    parSet<-calSet.par[i,]
    dat<-calSet[[i]]
    
    #Embed data for the LDE and adjust the embedding dimensions
    D<-floor(round(2*pi/sqrt(-parSet$Period))/2)
    tEmbedded <- gllaEmbed(dat, embed=D, tau=1, idColumn=FALSE)
    colnames(tEmbedded)<-paste0("x1_",0:(D-1))
    tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))
    Lmat<-matrix((1:D)-mean(1:D),nrow=D, ncol=3)
    for(k in 1:3) Lmat[,k]<-Lmat[,k]^(k-1)/factorial(k-1)
    
    LDE<-readRDS("LDE_model.RDS")
    LDE<-mxModel(LDE, 
                 mxMatrix("Full", nrow=D, 3, values=Lmat, name="L", dimnames=list(colnames(tEmbedded)[1:D], rownames(LDE$S$values))),
                 mxMatrix("Diag", nrow=D, ncol=D, values=.1, free=T, labels="e_x1", lbound=1e-10, name="U", dimnames=list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])),
                 mxMatrix("Full", nrow=1, ncol=D, values=.001, free=T, labels="mu_x1", name="M", dimnames=list(NULL,colnames(tEmbedded)[1:D])),
                 mxExpectationNormal("R","M"),
                 mxData(tEmbedded, type="raw"))
    dimnames(LDE$R)<-list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])
    LDE$D$values<-D

    #Set up LDEPDA
    LDEPDA<-LDE
    LDEPDA<-omxSetParameters(LDEPDA, "d2x1ts1", values=0,free=F,lbound=NA)
    LDEPDA$compute<-mxComputeNelderMead(maxIter = 500000)
    phi<- par[1] + par[2]*D + par[3]*length(dat) 
    AR<-floor(D*sqrt(parSet$p))
    
    if(phi <= 0)return(Inf)
    LDEPDA.fit<-NULL
    try(LDEPDA.fit<-PDA(LDEPDA, tEmbedded, thres.sd=phi, AR.order = AR, max.iter=5, plot=F))
    if(is.null(LDEPDA.fit))return(Inf)
    
    estVal<-exp(pi*coef(LDEPDA.fit$fit)[2]/sqrt(-coef(LDEPDA.fit$fit)[1]))
    trueVal<-exp(pi*parSet[4]/sqrt(-parSet[3]))
    out[i, 1]<-estVal - as.numeric(trueVal)
    
    if(trace) cat(D, phi, AR, as.numeric(out[i,1]),"\t\t\t\t\r\n")
  }
  out[is.infinite(out)]<-NA
  return( sd(out) + mean(out[out^2 < 1], na.rm=T))
}

# simIter(c(.1, .1, .1), trace=T) #Test function

# Run global optimization on simIter --------------------------------------
# cl<-makeCluster(4)
DE<-DEoptim(simIter,
            lower= c(-100, -5, -.01),
            upper= c(100, 5, .01),
            control=list(
              itermax=5000,
              NP=20,
              strategy=4,
              # cluster=cl,
              parallelType=1,
              trace = T,
              CR=.8,
              F=.5,
              packages=list("OpenMx", "deSolve"),
              parVar=list("gllaEmbed","genSN","DLOsim","smear","PDA", "genDifs", "calSet", "calSet.par")
            )
)
# stopCluster(cl)

saveRDS(DE, "calibrationOutput.RDS")

