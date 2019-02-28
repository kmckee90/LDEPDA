
#--------------------------------------------------------------------------
#A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering
#Pt. 3: Analysis of simulation results
#Author: Kevin McKee (mckeek@vcu.edu)
#Virginia Commonwealth University
#Virginia Institute of Psychiatric and Behavioral Genetics
#Last updated: 2/27/2019
#--------------------------------------------------------------------------


library(psych)
library(aplpack)
library(xtable)

# Functions ---------------------------------------------------------------
mround <- function(x,base){ 
  base*round(x/base) 
} 


# Load data ---------------------------------------------------------------
dp<-read.csv("outputDiffProc_Rand.csv", header=T)
sn<-read.csv("outputShotNoise_Rand.csv", header=T)

# Nonconvergence rates ----------------------------------------------------
ncr.dp<-matrix(unlist(apply(dp[,14:16], 2, function(x)table(x) )), 3,3, byrow=F)
round(ncr.dp/nrow(dp)*100,2)
ncr.dp

ncr.sn<-matrix(unlist(apply(sn[,14:16], 2, function(x)table(x) )), 3,3, byrow=F)
round(ncr.sn/nrow(sn)*100,2)
ncr.sn


# QC ----------------------------------------------------------------------
dp.qc<-dp[ apply(dp[,14:16], 1, sum)==0, ]
sn.qc<-sn[ apply(sn[,14:16], 1, sum)==0, ]

# qc.thres<-qchisq(.9, df=4)
# dp.qc<-dp.qc[outlier(dp.qc[,10:13], plot=F)<qc.thres,]
# sn.qc<-sn.qc[outlier(sn.qc[,10:13], plot=F)<qc.thres,]

# Convert parameters -------------------------------------------------------

dp.qc$Damping<-exp(dp.qc$Damping*pi/sqrt(-dp.qc$Period))
dp.qc$Period<-(2*pi/sqrt(-dp.qc$Period))
dp.qc$LDE.Zeta_x1ts1<-exp(dp.qc$LDE.Zeta_x1ts1*pi/sqrt(-dp.qc$LDE.Eta_x1ts1))
dp.qc$LDE.Eta_x1ts1<-(2*pi/sqrt(-dp.qc$LDE.Eta_x1ts1))
dp.qc$LDE.EA.Zeta_x1ts1<-exp(dp.qc$LDE.EA.Zeta_x1ts1*pi/sqrt(-dp.qc$LDE.EA.Eta_x1ts1))
dp.qc$LDE.EA.Eta_x1ts1<-(2*pi/sqrt(-dp.qc$LDE.EA.Eta_x1ts1))
dp.qc$SSM.ssm.A.2.2.<-exp(dp.qc$SSM.ssm.A.2.2.*pi/sqrt(-dp.qc$SSM.ssm.A.2.1.))
dp.qc$SSM.ssm.A.2.1.<-(2*pi/sqrt(-dp.qc$SSM.ssm.A.2.1.))

sn.qc$Damping<-exp(sn.qc$Damping*pi/sqrt(-sn.qc$Period))
sn.qc$Period<-(2*pi/sqrt(-sn.qc$Period))
sn.qc$LDE.Zeta_x1ts1<-exp(sn.qc$LDE.Zeta_x1ts1*pi/sqrt(-sn.qc$LDE.Eta_x1ts1))
sn.qc$LDE.Eta_x1ts1<-(2*pi/sqrt(-sn.qc$LDE.Eta_x1ts1))
sn.qc$LDE.EA.Zeta_x1ts1<-exp(sn.qc$LDE.EA.Zeta_x1ts1*pi/sqrt(-sn.qc$LDE.EA.Eta_x1ts1))
sn.qc$LDE.EA.Eta_x1ts1<-(2*pi/sqrt(-sn.qc$LDE.EA.Eta_x1ts1))
sn.qc$SSM.ssm.A.2.2.<-exp(sn.qc$SSM.ssm.A.2.2.*pi/sqrt(-sn.qc$SSM.ssm.A.2.1.))
sn.qc$SSM.ssm.A.2.1.<-(2*pi/sqrt(-sn.qc$SSM.ssm.A.2.1.))




# Clean up bad conversions ------------------------------------------------
sn.qc$SSM.ssm.A.2.1[is.nan(sn.qc$SSM.ssm.A.2.1.)]<-NA
sn.qc$SSM.ssm.A.2.2[is.nan(sn.qc$SSM.ssm.A.2.2.)]<-NA
sn.qc$LDE.EA.Eta_x1ts1[is.nan(sn.qc$LDE.EA.Eta_x1ts1)]<-NA
sn.qc$LDE.EA.Zeta_x1ts1[is.nan(sn.qc$LDE.EA.Zeta_x1ts1)]<-NA

sn.qc$SSM.ssm.A.2.1[is.infinite(sn.qc$SSM.ssm.A.2.1.)]<-NA
sn.qc$SSM.ssm.A.2.2[is.infinite(sn.qc$SSM.ssm.A.2.2.)]<-NA
sn.qc$LDE.EA.Eta_x1ts1[is.infinite(sn.qc$LDE.EA.Eta_x1ts1)]<-NA
sn.qc$LDE.EA.Zeta_x1ts1[is.infinite(sn.qc$LDE.EA.Zeta_x1ts1)]<-NA



dp.qc$SSM.ssm.A.2.1[is.nan(dp.qc$SSM.ssm.A.2.1.)]<-NA
dp.qc$SSM.ssm.A.2.2[is.nan(dp.qc$SSM.ssm.A.2.2.)]<-NA
dp.qc$LDE.EA.Eta_x1ts1[is.nan(dp.qc$LDE.EA.Eta_x1ts1)]<-NA
dp.qc$LDE.EA.Zeta_x1ts1[is.nan(dp.qc$LDE.EA.Zeta_x1ts1)]<-NA

dp.qc$SSM.ssm.A.2.1[is.infinite(dp.qc$SSM.ssm.A.2.1.)]<-NA
dp.qc$SSM.ssm.A.2.2[is.infinite(dp.qc$SSM.ssm.A.2.2.)]<-NA
dp.qc$LDE.EA.Eta_x1ts1[is.infinite(dp.qc$LDE.EA.Eta_x1ts1)]<-NA
dp.qc$LDE.EA.Zeta_x1ts1[is.infinite(dp.qc$LDE.EA.Zeta_x1ts1)]<-NA

sn.qc<-sn.qc[complete.cases(sn.qc[10:13]),]
dp.qc<-dp.qc[complete.cases(dp.qc[10:13]),]


# Overall performance comparison ------------------------------------------
#DP
LDE.better<-dp.qc[which((dp.qc$SSM.ssm.A.2.2.-dp.qc$Damping)^2 > (dp.qc$LDE.EA.Zeta_x1ts1-dp.qc$Damping)^2 ),]
SSM.better<-dp.qc[which((dp.qc$SSM.ssm.A.2.2.-dp.qc$Damping)^2 < (dp.qc$LDE.EA.Zeta_x1ts1-dp.qc$Damping)^2 ),]

nrow(LDE.better)/nrow(dp.qc)


LDE.better$measErr<-mround(LDE.better$measErr, .1)
SSM.better$measErr<-mround(SSM.better$measErr, .1)
LDE.better$Period<-mround(LDE.better$Period, 4)
SSM.better$Period<-mround(SSM.better$Period, 4)
LDE.better$Damping<-mround(LDE.better$Damping, .1)
SSM.better$Damping<-mround(SSM.better$Damping, .1)

pdf("GridsDP.pdf", width=9, height=3)

par(mfrow=c(1,3), mai=c(.7,.7,.2,.2))

compGrid<-table(LDE.better$measErr, LDE.better$Period)>table(SSM.better$measErr, SSM.better$Period)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Period",ylab="Measurement Error", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()
compGrid<-table(LDE.better$measErr, LDE.better$Damping)>table(SSM.better$measErr, SSM.better$Damping)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Amplitude Ratio",ylab="Measurement Error", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()
compGrid<-table(LDE.better$Period, LDE.better$Damping)>table(SSM.better$Period, SSM.better$Damping)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Amplitude Ratio",ylab="Period", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()
dev.off()


#SN
LDE.better<-sn.qc[which((sn.qc$SSM.ssm.A.2.2.-sn.qc$Damping)^2 > (sn.qc$LDE.EA.Zeta_x1ts1-sn.qc$Damping)^2 ),]
SSM.better<-sn.qc[which((sn.qc$SSM.ssm.A.2.2.-sn.qc$Damping)^2 < (sn.qc$LDE.EA.Zeta_x1ts1-sn.qc$Damping)^2 ),]

#Overall LDE performance
nrow(LDE.better)/nrow(sn.qc)

LDE.better$p<-mround(LDE.better$p, .1)
SSM.better$p<-mround(SSM.better$p, .1)
LDE.better$measErr<-mround(LDE.better$measErr, .1)
SSM.better$measErr<-mround(SSM.better$measErr, .1)
LDE.better$Period<-mround(LDE.better$Period, 4)
SSM.better$Period<-mround(SSM.better$Period, 4)
LDE.better$Damping<-mround(LDE.better$Damping, .1)
SSM.better$Damping<-mround(SSM.better$Damping, .1)


pdf("GridsSN.pdf", width=9, height=6)
par(mfrow=c(2,3), mai=c(.7,.7,.2,.2))

compGrid<-table(LDE.better$measErr, LDE.better$Period)>table(SSM.better$measErr, SSM.better$Period)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Period",ylab="Measurement Error", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()

compGrid<-table(LDE.better$measErr, LDE.better$Damping)>table(SSM.better$measErr, SSM.better$Damping)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Amplitude Ratio",ylab="Measurement Error", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()

compGrid<-table(LDE.better$Period, LDE.better$Damping)>table(SSM.better$Period, SSM.better$Damping)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Amplitude Ratio",ylab="Period", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()


compGrid<-table(LDE.better$p, LDE.better$Period)>table(SSM.better$p, SSM.better$Period)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Period",ylab="Outlier Probability", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()

compGrid<-table(LDE.better$p, LDE.better$measErr)>table(SSM.better$p, SSM.better$measErr)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Measurement Error",ylab="Outlier Probability", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()


compGrid<-table(LDE.better$p, LDE.better$Damping)>table(SSM.better$p, SSM.better$Damping)
image(t(compGrid), col=c("white","black"), axes=FALSE,xlab="Amplitude Ratio",ylab="Outlier Probability", useRaster=T, cex.lab=1.5)
axis(1, at = seq(0,1,length.out=ncol(compGrid)/1), labels=colnames(compGrid)[seq(1, ncol(compGrid), 1)],tick=TRUE, cex.axis=1.25)
axis(2, at = seq(0,1,length.out=nrow(compGrid)/1), labels=rownames(compGrid)[seq(1, nrow(compGrid), 1)],tick=TRUE, cex.axis=1.25)
box()

dev.off()


# Correlation of errors ---------------------------------------------------
cor((dp.qc$SSM.ssm.A.2.1.-dp.qc$Period), (dp.qc$LDE.EA.Eta_x1ts1-dp.qc$Period))
cor((sn.qc$SSM.ssm.A.2.1.-sn.qc$Period), (sn.qc$LDE.EA.Eta_x1ts1-sn.qc$Period))
cor((dp.qc$SSM.ssm.A.2.2.-dp.qc$Damping), (dp.qc$LDE.EA.Zeta_x1ts1-dp.qc$Damping))
cor((sn.qc$SSM.ssm.A.2.2.-sn.qc$Damping), (sn.qc$LDE.EA.Zeta_x1ts1-sn.qc$Damping))


# Multiple Regressions --------------------------------------------------------
  outlierThres<-2

  eLDE<- (dp.qc$LDE.Zeta_x1ts1-dp.qc$Damping)^2
  ePDA<- (dp.qc$LDE.EA.Zeta_x1ts1-dp.qc$Damping)^2
  eSSM<- (dp.qc$SSM.ssm.A.2.2.-dp.qc$Damping)^2
  
  lmLDE<-summary(lm(eLDE[eLDE<outlierThres]~ measErr +  N + Period + Damping, data=dp.qc[eLDE<outlierThres,] ))
  lmPDA<-summary(lm(ePDA[ePDA<outlierThres]~ measErr + N + Period + Damping, data=dp.qc[ePDA<outlierThres,] ), na.action=na.pass)
  lmSSM<-summary(lm(eSSM[eSSM<outlierThres]~ measErr +  N + Period + Damping, data=dp.qc[eSSM<outlierThres,] ))
  
  simEffects1a<-cbind(lmLDE$coefficients[,1], lmPDA$coefficients[,1], lmSSM$coefficients[,1])
  colnames(simEffects1a)<-c("LDE", "LDE+PDA","SSM")
  
  eLDE<- (sn.qc$LDE.Zeta_x1ts1-sn.qc$Damping)^2
  ePDA<- (sn.qc$LDE.EA.Zeta_x1ts1-sn.qc$Damping)^2
  eSSM<- (sn.qc$SSM.ssm.A.2.2.-sn.qc$Damping)^2
  
  lmLDE<-summary(lm(eLDE[eLDE<outlierThres]~ measErr +  N + Period + Damping + p, data=sn.qc[eLDE<outlierThres,] ))
  lmPDA<-summary(lm(ePDA[ePDA<outlierThres]~ measErr + N + Period + Damping+ p, data=sn.qc[ePDA<outlierThres,] ), na.action=na.pass)
  lmSSM<-summary(lm(eSSM[eSSM<outlierThres]~ measErr +  N + Period + Damping+ p, data=sn.qc[eSSM<outlierThres,] ))
  
  simEffects1b<-cbind(lmLDE$coefficients[,1], lmPDA$coefficients[,1], lmSSM$coefficients[,1])
  colnames(simEffects1b)<-c("LDE", "LDE+PDA","SSM")
  

  
  outlierThres<-20
  
  eLDE<- (dp.qc$LDE.Eta_x1ts1-dp.qc$Period)^2
  ePDA<- (dp.qc$LDE.EA.Eta_x1ts1-dp.qc$Period)^2
  eSSM<- (dp.qc$SSM.ssm.A.2.1.-dp.qc$Period)^2
  
  lmLDE<-summary(lm(eLDE[eLDE<outlierThres]~ measErr +  N + Period + Damping, data=dp.qc[eLDE<outlierThres,] ))
  lmPDA<-summary(lm(ePDA[ePDA<outlierThres]~ measErr + N + Period + Damping, data=dp.qc[ePDA<outlierThres,] ), na.action=na.pass)
  lmSSM<-summary(lm(eSSM[eSSM<outlierThres]~ measErr +  N + Period + Damping, data=dp.qc[eSSM<outlierThres,] ))
  
  simEffects2a<-cbind(lmLDE$coefficients[,1], lmPDA$coefficients[,1], lmSSM$coefficients[,1])
  colnames(simEffects2a)<-c("LDE", "LDE+PDA","SSM")

  eLDE<- (sn.qc$LDE.Eta_x1ts1-sn.qc$Period)^2
  ePDA<- (sn.qc$LDE.EA.Eta_x1ts1-sn.qc$Period)^2
  eSSM<- (sn.qc$SSM.ssm.A.2.1.-sn.qc$Period)^2

  lmLDE<-summary(lm(eLDE[eLDE<outlierThres]~ measErr +  N + Period + Damping + p, data=sn.qc[eLDE<outlierThres,] ))
  lmPDA<-summary(lm(ePDA[ePDA<outlierThres]~ measErr + N + Period + Damping+ p, data=sn.qc[ePDA<outlierThres,] ), na.action=na.pass)
  lmSSM<-summary(lm(eSSM[eSSM<outlierThres]~ measErr +  N + Period + Damping+ p, data=sn.qc[eSSM<outlierThres,] ))
  
  simEffects2b<-cbind(lmLDE$coefficients[,1], lmPDA$coefficients[,1], lmSSM$coefficients[,1])
  colnames(simEffects2b)<-c("LDE", "LDE+PDA","SSM")
  
  
  xtable(cbind(rbind(simEffects1a,simEffects1b),rbind(simEffects2a,simEffects2b)), digits=5)
  
    
      


  
  


