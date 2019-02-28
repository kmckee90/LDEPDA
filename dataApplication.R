
#--------------------------------------------------------------------------
#A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering
#Pt. 4: Application to postural control data
#Author: Kevin McKee (mckeek@vcu.edu)
#Virginia Commonwealth University
#Virginia Institute of Psychiatric and Behavioral Genetics
#Last updated: 2/27/2019
#--------------------------------------------------------------------------


# Functions ---------------------------------------------------------------
convPars<-function(pars){
  eta<-pars[1]; zeta<-pars[2]
  freq<- 1/(2*pi*sqrt(-1/(eta)+zeta^2/4))
  rat<-exp(pi*zeta/sqrt(-eta))
  return(c("Hz"=freq, "ratio"=rat))
}

stdz<-function(x) (x-mean(x))/sd(x,na.rm=T)



# Load and prepare data ---------------------------------------------------
dat<-read.table("sdata.txt",header=T)$COPy.cm.[2000:6000]
ds<-10# 100/ds Hz
dat<-dat[seq(1,length(dat),ds)]
dat<-stdz(dat)


# Spectral analysis -------------------------------------------------------
spec<-spectrum(dat, log="no", plot=T)
par(mfrow=c(2,1), mai=c(.75,.75,.1,.1))
plot(1:length(dat)*(ds/100), dat, type="l", xlab="", ylab="", cex.axis=1.2)
title(ylab="Center of Pressure (Lateral)", xlab="Time (Seconds)", cex.lab=1.5, line=2.2)
plot(spec$freq*(100/ds), spec$spec, type="l",xlim=c(0, 2), ylab="", xlab="", cex.axis=1.2, cex.lab=1.5)
title(ylab="Spectral Density", xlab="Frequency (Hz)", cex.lab=1.5, line=2.2)
#Oscillation period of 2s or so. 
#Sampling rate: 100Hz
#Half period: 100 obs

# Embed data for the LDE --------------------------------------------------
D<-6
tEmbedded <- gllaEmbed(dat, embed=D, tau=1, idColumn=FALSE)
colnames(tEmbedded)<-paste0("x1_",0:(D-1))
tEmbedded<-cbind(tEmbedded, "Occasion"=1:nrow(tEmbedded))
Lmat<-matrix((1:D)-mean(1:D) ,nrow=D, ncol=3)*(ds/100)
for(k in 1:3) Lmat[,k]<-Lmat[,k]^(k-1)/factorial(k-1)

# LDE ---------------------------------------------------------------------
LDE<-readRDS("LDE_model.RDS")
LDE<-mxModel(LDE, 
             mxMatrix("Full", nrow=D, 3, values=Lmat, name="L", dimnames=list(colnames(tEmbedded)[1:D], rownames(LDE$S$values))),
             mxMatrix("Diag", nrow=D, ncol=D, values=.1, free=T, labels="e_x1", lbound=1e-10, name="U", dimnames=list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])),
             mxMatrix("Full", nrow=1, ncol=D, values=.001, free=T, labels="mu_x1", name="M", dimnames=list(NULL,colnames(tEmbedded)[1:D])),
             mxExpectationNormal("R","M"))
LDE$A$lbound<-NA
LDE$A$ubound<-NA
LDE$etas$lbound<-NA
LDE$etas$ubound<-NA

dimnames(LDE$R)<-list(colnames(tEmbedded)[1:D],colnames(tEmbedded)[1:D])
LDE$D$values<-D
LDEPDA<-LDE
LDE<-mxModel(LDE, mxData(tEmbedded, type="raw"))
LDE$compute<-omxDefaultComputePlan()
LDE$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
LDE<-mxRun(LDE)
if(LDE$output$status$code!=0) LDE<-mxTryHard(LDE, silent=T, extraTries = 15)
coef(LDE)


# Boostrap ----------------------------------------------------------------
pars<-matrix(NA, 500, 4)
for(i in 1:nrow(pars)){
  rs<-tEmbedded[sample(1:nrow(tEmbedded), replace=T),]
  LDE<-mxModel(LDE, mxData(rs, type="raw"))
  LDE<-mxRun(LDE, silent=T)
  if(LDE$output$status$code!=0) LDE<-mxTryHard(LDE, silent=T, extraTries = 15)
  pars[i,]<-c(coef(LDE)[1:2], convPars(coef(LDE)[1:2]))
  cat("\r",i,"\t\t")
}
CI<-t(apply(pars, 2, function(x){ c(round( sort(x)[ round(length(x)*.025) ],2), round( sort(x)[ round(length(x)*.975) ],2) )} ))
apply(CI, 1, function(x) cat(paste0("(",x[1],", ",x[2],")"),"\n"))



# # LDEPDA ------------------------------------------------------------------
# LDEPDA<-readRDS("LDE_model.RDS")
LDEPDA<-omxSetParameters(LDEPDA, "d2x1ts1", values=0,free=F,lbound=NA)
LDEPDA$compute<-omxDefaultComputePlan()
LDEPDA$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
thres<- -0.152344    +0.004885*D    +0.002200*length(dat) #NEW: VARIANCE AND BIAS (QUICK CALIBRATION): Obj value 0.317456
AR<-floor(D*sqrt(.15)) #Guess: probability of occurrances is .15
LDEPDA.fit<-PDA(LDEPDA, tEmbedded, thres.sd=thres,AR.order=AR, max.iter=25, plot=T)

coef(LDEPDA.fit$fit)



# Boostrap ----------------------------------------------------------------
pars<-matrix(NA, 500, 4)
for(i in 1:nrow(pars)){
  rs<-tEmbedded[sample(1:nrow(tEmbedded), replace=T),]
  rs<-rs[order(rs[,7]),]
  LDEPDA.fit<-PDA(LDEPDA, rs, thres.sd=thres,AR.order=AR, max.iter=25, plot=F)
  pars[i,]<-c(coef(LDEPDA.fit$fit)[1:2], convPars(coef(LDEPDA.fit$fit)[1:2]))
  cat("\r",i,"\t\t")
}
CI<-t(apply(pars, 2, function(x){ c(round( sort(x)[ round(length(x)*.025) ],2), round( sort(x)[ round(length(x)*.975) ],2) )} ))
apply(CI, 1, function(x) cat(paste0("(",x[1],", ",x[2],")"),"\n"))



# SSM ---------------------------------------------------------------------
  ssm<-readRDS("SSM_model.RDS")
  ssm<-mxModel(ssm, mxData(data.frame("obs"=dat, "tim"=1:length(dat)*(ds/100)), type="raw"))
  ssm<-mxModel(ssm, 
               mxMatrix("Full", 1,1,values=0,free=T,labels="meanX", name="D"),
               mxMatrix("Full", 1,1,values=1,free=F, name="u"))
  ssm$A$values<-c(0,-5, 1, -1)
  ssm$R$lbound<-1e-10
  
  ssm$compute<-omxDefaultComputePlan()
  ssm$compute$steps$GD<-mxComputeNelderMead(maxIter = 500000)
  ssm<-mxTryHard(ssm, silent=F, extraTries = 15)
  
  coef(ssm)
  summary(ssm)
  
  
  # Boostrap ----------------------------------------------------------------
  blocksize<-80
  
  cl<-makeCluster(4)
  registerDoSNOW(cl)
  pars<-foreach(i = 1:500, .combine=rbind, .export = c("ssm","tEmbedded", "blocksize"),.packages = c("OpenMx"),.init=NULL )%dopar%{
    start<-sample(1:round(blocksize/2), 1)
    blockInds<-cbind(
      c(1,       seq(start,length(dat), blocksize)), 
      c(start-1, seq(start+blocksize-1,length(dat), blocksize), length(dat))  )
    seg<-list()
    for(j in 1:(nrow(blockInds))) seg[[j]]<-dat[blockInds[j,1]:blockInds[j,2]]
    rs<-c(seg[[1]], unlist( seg[ sample(3:(length(seg)-1))]), seg[[ length(seg) ]])
    
    ssm<-mxModel(ssm, mxData(data.frame("obs"=rs, "tim"=1:length(rs)*(ds/100)), type="raw"))
    ssm<-mxTryHard(ssm, silent=T, extraTries = 15)
    
    return(c(coef(ssm)[1:2], convPars(coef(ssm)[1:2])))
  }
  stopCluster(cl)
  
  CI<-t(apply(pars, 2, function(x){ c(round( sort(x)[ round(length(x)*.025) ],2), round( sort(x)[ round(length(x)*.975) ],2) )} ))
  apply(CI, 1, function(x) cat(paste0("(",x[1],", ",x[2],")"),"\n"))
  
  
  