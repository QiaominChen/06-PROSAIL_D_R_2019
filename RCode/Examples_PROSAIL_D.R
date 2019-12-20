# * Author:       Qiaomin Chen (qiaomin.chen@uq.net.au)
# * Created:      Sunday, 28 April 2019 (GMT+10)
# * Introduction: Run PROSAIL model
# * Description:  Rewritten in R on the base of the matlab version (PROSAIL_D_MATLAB_2017)

# ******************************************************************************************************
# Original info in matlab version
# This program allows modelling reflectance data from canopy.
# - modelling leaf optical properties with PROSPEcT-D (Feret et al., 2017)
# - modelling leaf inclination distribution function with the subroutine campbell(ALA) or dladgen(a,b)
# - modelling canopy reflectance with 4SAIL (verhoef et al., 2007)
#
# References:
# Verhoef et al. (2007) Unified Optical-Thermal Four-Stream Radiative Transfer Theory for Homogeneous 
# Vegetation Canopies, IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, 45(6), June 2007
# Féret et al. (2008) PROSPECT-4 and 5: Advances in the Leaf Optical Properties Model Separating
# Photosynthetic Pigments, REMOTE SENSING OF ENVIRONMENT
# 
# This specific absorption coefficient corresponding to brown pigment is provided by Frederic Baret
# (EMMAH, INRA Avignon, baret@avignon.inra.fr) and used with his autorization.
# The model PRO4SAIL is based on a version provided by 
# Wout Verhoef
# NLR
# April/May 2003
#
# The original 2-parameter LIDF model is developed by and described in:
# W. Verhoef, 1998, "Theory of radiative transfer models applied in optical remote sensing of vegetation canopies", 
# Wageningen Agricultural University, The Netherlands, 310 pp. (PhD thesis)
# The Ellipsoidal LIDF model is taken from:
# Campbell (1990) Derivtion of an angle density function for canopies with ellipsoidal leaf angle distribution, 
# Agricultural and Forest Meteorology, 49:173-176.
# ******************************************************************************************************


# Example One (main_PROSAIL) ----------------------------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)
source("RCode/model_PRO4SAIL.R")

# ******************************************************************************************************
# Leaf Chemical
# ******************************************************************************************************
Cab <- 40     # cholorophyll conent(ug/cm2)
Car <- 8      # carotenoid content (ug/cm2)
Ant <- 0      # Anthocyanins content (ug/cm2)
Cbrown <- 0.0 # brown pigment content (arbitrary units)
Cw <- 0.01    # EWT (cm)
Cm <- 0.009   # leaf matter weight (g/cm2)
N <- 1.5      # leaf structure coefficient

# ************************************************************************************
# 4SAIL canopy structure and geomery parameter
# ************************************************************************************
# LIDF model
# 2-parameter LIDF model: TypeLidf = 1
# LIDFa LIDF parameter a, which controls the average leaf slope
# LIDFb LIDF parameter b, which controls the distribution's bimodality
#	LIDF type 		LIDFa  LIDFb
#	Planophile 		    1		   0
#	Erectophile      -1	 	   0
#	Plagiophile 	    0		  -1
#	Extremophile 	    0		   1
#	Spherical     -0.35  -0.15
#	Uniform           0      0
# requirement: |LIDFa| + |LIDFb| < 1
#
# Ellipsoidal LIDF model: TypeLidf = 2
# LIDFa = average leaf angle (degree) [0 = planophile, 90 = erectophile]
# LIDFb = 0
# ************************************************************************************
TypeLidf=2;
if (TypeLidf == 1){
  LIDFa <-  -0.35
  LIDFb <-  -0.15
} else if(TypeLidf == 2){
  LIDFa <-  30
  LIDFb <-  0
} else{
  stop("Error: Nonexistent TypeLidf level!")
}
LAI <- 5       # leaf area index (m2/m2)
hspot <- 0.01  # hotspot size parameter
tts <- 30      # solar zenith angle (degree)
tto <- 10      # observer zenith angle (degree)
psi <- 90      # relative azimuth angle (degree)

# *************************************************************************************
# Soil reflectance properties
# *************************************************************************************
psoil <- 1     # soil factor (psoil=0: wet soil/ psoil=1: dry soil)

# ************************************************************************************
# Call PROSAIL model
# ************************************************************************************
df <- main_PROSAIL(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,psoil)

# resv: directional reflectance
ggplot(df) + 
  geom_line(aes(x = df$lambda, y = df$resv)) +
  labs(x = "wavelength (mm)", y = "canopy reflectance")

df <- df %>% 
  mutate_if(is.double,round,digits=6)

write.table(df, file = "Rcodes/Refl_CAN.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE)


# Example Two (pro4sail) --------------------------------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)
source("RCodes/function_PROSAIL.R")
source("RCodes/model_PRO4SAIL.R")

# ******************************************************************************************************
# LIDF model
# 2-parameter LIDF model: TypeLidf = 1
# LIDFa LIDF parameter a, which controls the average leaf slope
# LIDFb LIDF parameter b, which controls the distribution's bimodality
#	LIDF type 		LIDFa  LIDFb
#	Planophile 		    1		   0
#	Erectophile      -1	 	   0
#	Plagiophile 	    0		  -1
#	Extremophile 	    0		   1
#	Spherical     -0.35  -0.15
#	Uniform           0      0
# requirement: |LIDFa| + |LIDFb| < 1
#
# Ellipsoidal LIDF model: TypeLidf = 2
# LIDFa = average leaf angle (degree) [0 = planophile, 90 = erectophile]
# LIDFb = 0
# ******************************************************************************************************
TypeLidf=2;
if (TypeLidf == 1){
  LIDFa <-  -0.35
  LIDFb <-  -0.15
} else if(TypeLidf == 2){
  LIDFa <-  30
  LIDFb <-  0
} else{
  stop("Error: Nonexistent TypeLidf level!")
}

# ******************************************************************************************************
# Leaf Chemical and spetral properties
# ******************************************************************************************************
Cab <- 40     # cholorophyll conent(ug/cm2)
Car <- 8      # carotenoid content (ug/cm2)
Ant <- 0      # Anthocyanins content (ug/cm2)
Cbrown <- 0.0 # brown pigment content (arbitrary units)
Cw <- 0.01    # EWT (cm)
Cm <- 0.009   # leaf matter weight (g/cm2)
N <- 1.5      # leaf structure coefficient
data <- read_dataSpec_DB()

# ************************************************************************************
# 4SAIL canopy structure and geomery parameter
# ************************************************************************************
LAI <- 5       # leaf area index (m2/m2)
hspot <- 0.01  # hotspot size parameter
tts <- 30      # solar zenith angle (degree)
tto <- 10      # observer zenith angle (degree)
psi <- 90      # relative azimuth angle (degree)

# *************************************************************************************
# Soil reflectance properties
# *************************************************************************************
Rsoil1 <- data$Rsoil1   # dry soil reflectance property
Rsoil2 <- data$Rsoil2   # wet soil reflectance property
psoil <- 1              # soil factor (psoil=0: wet soil/ psoil=1: dry soil)
rsoil0 <- psoil*Rsoil1+(1-psoil)*Rsoil2


# ************************************************************************************
# Call PRO4SAIL model
# ************************************************************************************
df <- pro4sail(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0)
rdot <- df$rdot
rsot <- df$rsot
rddt <- df$rddt
rsdt <- df$rsdt

# ************************************************************************************
# direct / diffuse light
# The direct and diffuse light are taken into account as proposed by:
# Francois et al. (2002) Conversion of 400–1100 nm vegetation albedo measurements into
# total shortwave broadband albedo using a canopy radiative transfer model, Agronomie.
# ************************************************************************************
Es <- data$Es   # direct sun energy
Ed <- data$Ed   # diffuse sun energy
rd <- pi/180
# skyl: the proportion of diffuse energy
skyl <- 0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd) 
PARdiro <- (1-skyl)*Es
PARdifo <- skyl*Ed

# resv: directional reflectance
resv <- (rdot*PARdifo+rsot*PARdiro)/(PARdifo+PARdiro)
df <- tibble(lambda=data$lambda,resv=resv)

ggplot(df) + 
  geom_line(aes(x = lambda, y = resv)) +
  labs(x = "wavelength (mm)", y = "canopy reflectance")

df <- df %>% 
  mutate_if(is.double,round,digits=6)

write.table(df, file = "Rcodes/Refl_CAN.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE)