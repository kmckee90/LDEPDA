# LDEPDA
Scripts supporting manuscript "A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering" (2019)

These files were used to generate the results of the study "A Method of Correcting Estimation Failure in Latent Differential Equations with Comparisons to Kalman Filtering" (2019) by Kevin McKee, Michael Hunter, and Michael Neale.
A brief description of each file is given below.

GLLAfunctions.R - Functions for time-delay embedding data and generating GLLA projection matrices.

calibration.R - Script for calibrating the outlier classification threshold 'phi'

Simulation1_RandFX_DP.R - Diffusion process random effects simulations

Simulation1_RandFX_SN.R - Shot noise process random effects simulations

analysis.R - Analysis of the results of both simulations.

dataApplication.R - Analysis of a subset of lateral posture control data by Santos and Duarte (2016)

Files:

SSM_model.RDS - Pre-specified OpenMx state space model 

LDE_model.RDS - Pre-specified OpenMx latent differential equation model

sdata.txt - Lateral postural sway data in one individual

outputDiffProc_Rand.csv - Results of the random effects diffusion process simulation

outputShotNoise_Rand.csv - Results of the random effects shot noise process simulation

References:
Santos, D. A. and Duarte, M. (2016). A public data set of human balance evaluations.
PeerJ, 4. doi:10.7717/peerj.2648.


