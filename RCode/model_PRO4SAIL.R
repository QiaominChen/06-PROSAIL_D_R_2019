# * Author:       Qiaomin Chen (qiaomin.chen@uq.net.au)
# * Created:      Sunday, 28 April 2019 (GMT+10)
# * Introduction: 4SAIL model
# * Description:  Rewritten in R on the base of the matlab version (PROSAIL_D_MATLAB_2017)

# ************************************************************************************************************
# Original info in matlab version:
# %	this model PRO4SAIL is based on a version provided by
# %	Wout Verhoef
# %	NLR
# %	April/May 2003,
# %	original version downloadable at http://teledetection.ipgp.jussieu.fr/prosail/
# %	Improved and extended version of SAILH model that avoids numerical singularities
# %	and works more efficiently if only few parameters change.
# % References:
# % 	Verhoef et al. (2007) Unified Optical-Thermal Four-Stream Radiative
# % 	Transfer Theory for Homogeneous Vegetation Canopies, IEEE TRANSACTIONS
# % 	ON GEOSCIENCE AND REMOTE SENSING, VOL. 45, NO. 6, JUNE 2007
# ************************************************************************************************************

library(ggplot2)
source("RCode/function_PROSAIL.R")
source("RCode/model_PROSPECT_DB.R")

#' PRO4SAIL model, coupling 4SAIL with PROsPECT-D
#' Calculate four canopy reflectance factors
#' @param N leaf structure paramter
#' @param Cab chlorophyll a+b content in µg/cm²
#' @param Car carotenoids content in µg/cm²
#' @param Anth anthocyanin content in µg/cm²
#' @param Cbrown brown pigments content in µg/cm²
#' @param Cw equivalent water thickness in g/cm² or cm
#' @param Cm dry matter content in g/cm²
#' @param LIDFa leaf inclination function parameter a or ALA
#' @param LIDFb leaf inclination function parameter b
#' @param TypeLidf type of leaf inclination function 
#' TypeLidf=1, leaf inclination function characteried by parameter a and b
#' TypeLidf=2, ellipsoidal distribution function characteried by average leaf angle in degree
#' @param lai leaf area index
#' @param q hotspot size parameter
#' @param tts solar zenith angle
#' @param tto viewing zenith angle
#' @param psi relative azimuth angle
#' @param rsoil soil reflectance
#' @return df a data frame with four columns: rdot,rsot,rddt,rsdt
#' rdot: hemispherical-directional reflectance factor for solar incident flux
#' rdot: bi-hemispherical reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rsot: bi-directional reflectance factor
pro4sail <- function(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,lai,q,tts,tto,psi,rsoil){
  # Leaf optical properties
  LRT <- prospect_DB(N,Cab,Car,Ant,Cbrown,Cw,Cm)
  rho <- LRT$refl
  tau <- LRT$tran
  
  # ggplot(LRT) +
  #   geom_line(aes(x=lambda,y=refl,color='refl')) +
  #   geom_line(aes(x=lambda,y=1-tran,color='1-tran')) + 
  #   labs(x = "wavelength (mm)", y = "leaf optical properties") +
  #   scale_color_manual("Legend", breaks = c("refl","1-tran"), values = c("blue","red"))
  
  # Geometric quantities
  rd <- pi/180
  cts <- cos(rd*tts)
  cto <- cos(rd*tto)
  ctscto <- cts*cto
  tants <- tan(rd*tts)
  tanto <- tan(rd*tto)
  cospsi <- cos(rd*psi)
  # dso <- sqrt(tants*tants+tanto*tanto-2*tants*tanto*cospsi)
  sqrtdso <- tants*tants+tanto*tanto-2*tants*tanto*cospsi
  sqrtdso <- replace(sqrtdso,sqrtdso<0,0)
  dso <- sqrt(sqrtdso)
  
  # Generate leaf angle distribution from average leaf angle (ellipsoidal) or (a,b) parameters
  if(TypeLidf == 1){
    df <- dladgen(LIDFa,LIDFb)
    lidf <- df$freq
    litab <- df$litab
  }
  if(TypeLidf == 2){
    df <- campbell(LIDFa)
    lidf <- df$freq
    litab <- df$litab
  }
  
  # angular distance, compensation of shadow length
  # Calculate geometric factors associated with extinction and scattering
  # Initialise sums
  ks <- 0
  ko <- 0
  bf <- 0
  sob <- 0
  sof <- 0
  
  # Weighted sums over LIDF
  na <- length(litab)
  for(i in 1:na){
    # leaf inclination discrete values
    ttl <- litab[i]
    ctl <- cos(rd*ttl)
    df <- volscatt(tts,tto,psi,ttl)
    chi_s <- df$chi_s
    chi_o <- df$chi_o
    frho <- df$frho
    ftau <- df$ftau
    # Extinction coefficients (ks,ko)
    ksli <- chi_s/cts
    koli <- chi_o/cto
    # Area scattering coefficient fractions
    sobli <- frho*pi/ctscto
    sofli <- ftau*pi/ctscto
    bfli <- ctl*ctl
    ks <- ks+ksli*lidf[i]
    ko <- ko+koli*lidf[i]
    bf <- bf+bfli*lidf[i]
    sob <- sob+sobli*lidf[i]
    sof <- sof+sofli*lidf[i]
  }
  # Geometric factors to be used later with rho and tau
  sdb <- 0.5*(ks+bf)
  sdf <- 0.5*(ks-bf)
  dob <- 0.5*(ko+bf)
  dof <- 0.5*(ko-bf)
  ddb <- 0.5*(1+bf)
  ddf <- 0.5*(1-bf)
  # Here rho and tau come in
  sigb <- ddb*rho+ddf*tau
  sigf <- ddf*rho+ddb*tau
  att <- 1-sigf
  m2 <- (att+sigb)*(att-sigb)
  m2 <- replace(m2,m2<=0,0)
  m <- sqrt(m2)
  sb <- sdb*rho+sdf*tau
  sf <- sdf*rho+sdb*tau
  vb <- dob*rho+dof*tau
  vf <- dof*rho+dob*tau
  w <- sob*rho+sof*tau
  # Here the LAI comes in
  # Outputs for the case LAI = 0
  if(lai<0){
    tss <- 1
    too <- 1
    tsstoo <- 1
    rdd <- 0
    tdd <- 1
    rsd <- 0
    tsd <- 0
    rdo <- 0
    tdo <- 0
    rso <- 0
    rsos <- 0
    rsod <- 0
    
    rddt <- rsoil
    rsdt <- rsoil
    rdot <- rsoil
    rsodt <- 0*rsoil
    rsost <- rsoil
    rsot <- rsoil
  } else{
    # Other case (LAI>0)
    e1 <- exp(-m*lai)
    e2 <- e1*e1
    rinf <- (att-m)/sigb
    rinf2 <- rinf*rinf
    re <- rinf*e1
    denom <- 1-rinf2*e2
    
    J1ks <- Jfunc1(ks,m,lai)
    J2ks <- Jfunc2(ks,m,lai)
    J1ko <- Jfunc1(ko,m,lai)
    J2ko <- Jfunc2(ko,m,lai)
    
    Ps <- (sf+sb*rinf)*J1ks
    Qs <- (sf*rinf+sb)*J2ks
    Pv <- (vf+vb*rinf)*J1ko
    Qv <- (vf*rinf+vb)*J2ko
    
    rdd <- rinf*(1-e2)/denom
    tdd <- (1-rinf2)*e1/denom
    tsd <- (Ps-re*Qs)/denom
    rsd <- (Qs-re*Ps)/denom
    tdo <- (Pv-re*Qv)/denom
    rdo <- (Qv-re*Pv)/denom
    
    tss <- exp(-ks*lai)
    too <- exp(-ko*lai)
    z <- Jfunc3(ks,ko,lai)
    g1 <- (z-J1ks*too)/(ko+m)
    g2 <- (z-J1ko*tss)/(ks+m)
    
    Tv1 <- (vf*rinf+vb)*g1
    Tv2 <- (vf+vb*rinf)*g2
    T1 <- Tv1*(sf+sb*rinf)
    T2 <- Tv2*(sf*rinf+sb)
    T3 <- (rdo*Qs+tdo*Ps)*rinf
    
    # Multiple scattering contribution to bidirectional canopy reflectance
    rsod <- (T1+T2-T3)/(1-rinf2)
    # Treatment of the hotspot-effect
    alf <- 1e6
    # Apply correction 2/(K+k) suggested by F.M. Bréon
    if(q>0){
      alf <- (dso/q)*2/(ks+ko)
    }
    if(alf>200){
      # inserted H. Bach 1/3/04
      alf <- 200
    }
    if(alf == 0){
      # The pure hotspot - no shadow
      tsstoo <- tss
      sumint <- (1-tss)/(ks*lai)
    } else{
      # Outside the hotspot
      # fhot <- lai*sqrt(ko*ks)
      sqrtfhot <- ko*ks
      sqrtfhot <- replace(sqrtfhot,sqrtfhot<0,0)
      fhot <- lai*sqrt(sqrtfhot)
      # Integrate by exponential Simpson method in 20 steps
      # the steps are arranged according to equal partitioning
      # of the slope of the joint probability function
      x1 <- 0
      y1 <- 0
      f1 <- 1
      fint <- (1-exp(-alf))*0.05
      sumint <- 0
      for(i in 1:20){
        if(i<20){
          x2 <- -log(1-i*fint)/alf
        } else{
          x2 <- 1
        }
        y2 <- -(ko+ks)*lai*x2+fhot*(1-exp(-alf*x2))/alf
        f2 <- exp(y2)
        sumint <- sumint+(f2-f1)*(x2-x1)/(y2-y1)
        x1 <- x2
        y1 <- y2
        f1 <- f2
      }
      tsstoo <- f1
    }
    # Bidirectional reflectance
    # Single scattering contribution
    rsos <- w*lai*sumint
    # Total canopy contribution
    rso <- rsos+rsod
    # Interaction with the soil
    dn <- 1-rsoil*rdd
    # rddt: bi-hemispherical reflectance factor
    rddt <- rdd+tdd*rsoil*tdd/dn
    # rsdt: directional-hemispherical reflectance factor for solar incident flux
    rsdt <- rsd+(tsd+tss)*rsoil*tdd/dn
    # rdot: hemispherical-directional reflectance factor for solar incident flux
    rdot <- rdo+tdd*rsoil*(tdo+too)/dn
    # rsot: bi-directional reflectance factor
    rsodt <- rsod+((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn
    rsost <- rsos+tsstoo*rsoil
    rsot <- rsost+rsodt
  }
  df <- data.frame(rdot=rdot,rsot=rsot,rddt=rddt,rsdt=rsdt)
  return(df)
}

#' Calculate canopy reflectance
#' @return df a data frame with two columns: lamda (wavelength), resv (canopy reflectance)
main_PROSAIL <- function(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,lai,q,tts,tto,psi,psoil){

  # Spetral properties
  data <- read_dataSpec_DB()

  # Soil reflectance properties
  Rsoil1 <- data$Rsoil1   # dry soil reflectance property
  Rsoil2 <- data$Rsoil2   # wet soil reflectance property
  rsoil0 <- psoil*Rsoil1+(1-psoil)*Rsoil2

  # Call 4SAIL model to calculate four reflectance factors
  df <- pro4sail(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,lai,q,tts,tto,psi,rsoil0)
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

  return(df)

}
