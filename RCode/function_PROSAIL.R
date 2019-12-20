# * Author:       Qiaomin Chen (qiaomin.chen@uq.net.au)
# * Created:      Sunday, 28 April 2019 (GMT+10)
# * Introduction: Provide general functions for calculation of PROSPECT-D and 4SAIL model
# * Description:  Rewritten in R on the base of the matlab version (PROSAIL_D_MATLAB_2017)

#' Read data from dataSpec_PDB.txt.
#' @return data the dataSpec_DB, which provides fixed input for PROSAIL model.
read_dataSpec_DB <- function(){
  path <- 'RCode/data_Spec_PDB.txt'
  data <- read.table(path, sep = "\t", skip = 26)
  data <- dplyr::rename(data,
                 lambda = V1,
                 nr = V2,
                 Kab = V3,
                 Kcar = V4,
                 Kant = V5,
                 KBrown = V6,
                 Kw = V7,
                 Km = V8,
                 Es = V9,
                 Ed = V10,
                 Rsoil1 = V11,
                 Rsoil2 = V12)
  return(data)
}

#' Compute transmission of isotropic radiation.
#' @param alfa 
#' @param nr 
#' @return tav transmission
calctav <- function(alfa,nr){
  rd <-  pi/180
  n2 <- nr^2
  np <- n2 + 1
  nm <- n2 - 1
  a <- (nr+1)*(nr+1)/2
  k <- -(n2-1)*(n2-1)/4
  sa <- sin(alfa*rd)
  
  # b1 <- (alfa!=90)*1*sqrt((sa^2-np/2)*(sa^2-np/2)+k)
  sqrtb1 <- (sa^2-np/2)*(sa^2-np/2)+k
  sqrtb1 <- replace(sqrtb1,sqrtb1<0,0)
  b1 <- (alfa!=90)*1*sqrt(sqrtb1)
  b2 <- sa^2-np/2
  b <- b1-b2
  b3 <- b^3
  a3 <- a^3
  ts <- (k^2/(6*b3)+k/b-b/2)-(k^2/(6*a3)+k/a-a/2)
  
  tp1 <- -2*n2*(b-a)/(np^2)
  tp2 <- -2*n2*np*log(b/a)/(nm^2)
  tp3 <- n2*(1/b-1/a)/2
  tp4 <- 16*n2^2*(n2^2+1)*log((2*np*b-nm^2)/(2*np*a-nm^2))/(np^3*nm^2)
  tp5 <- 16*n2^3*(1/(2*np*b-nm^2)-1/(2*np*a-nm^2))/(np^3)
  tp <- tp1+tp2+tp3+tp4+tp5
  
  tav <- (ts+tp)/(2*sa^2)
  return(tav)
}

#' A function is called in dladgen function.
dcum <- function(a,b,t){
  rd <- pi/180
  if(a >= 1){
    f <- 1-cos(rd*t)
  } else{
    eps <- 1e-8
    delx <- 1
    x <- 2*rd*t
    p <- x
    while (delx >= eps) {
      y <- a*sin(x)+0.5*b*sin(2.0*x)
      dx <- 0.5*(y-x+p)
      x <- x+dx
      delx <- abs(dx)
    }
    f <- (2.0*y+p)/pi
  }
  return(f)
}

#' LIDF function with parameter a and b
#' @param a 
#' @param b 
#' @return df a data frame, the first column is freq and the second is litab.
dladgen <- function(a,b){
  litab <- c(5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 81.0, 83.0, 85.0, 87.0, 89.0)
  freq <- vector(length = 13)
  for (i1 in 1:8) {
    t <- i1*10
    freq[i1] <- dcum(a,b,t)
  }
  for (i2 in 9:12) {
    t <- 80.0+(i2-8)*2
    freq[i2] <- dcum(a,b,t)
  }
  freq[13] <- 1
  for (i in 13:2) {
    freq[i] <- freq[i]-freq[i-1]
  }
  df <- data.frame(freq=freq,litab=litab)
  return(df)
}

#' LIDF function with parameter ALA
#' Computation of the leaf angle distribution function value (freq) 
#' Ellipsoidal distribution function caracterised by the average leaf 
#' inclination angle in degree (ala), Campbell 1986 
#' @param ala average leaf inclination angle
#' @return df a data frame, the first column is freq0 and the second is litab
campbell <- function(ala){
  tx1 <- c(10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,82.0,84.0,86.0,88.0,90.0)
  tx2 <- c(0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,82.0,84.0,86.0,88.0)
  
  litab <- (tx1+tx2)/2
  n <- length(litab)
  t11 <- tx1*(pi/180)
  t12 <- tx2*(pi/180)
  excent <- exp(-1.6184e-5*ala^3+2.1145e-3*ala^2-1.2390e-1*ala+3.2491)
  sum0 <- 0
  
  freq <- vector(length = 13)
  for(i in 1:n){
    # x1 <- excent/sqrt(1+excent^2*tan(t11[i])^2)
    # x2 <- excent/sqrt(1+excent^2*tan(t12[i])^2)
    sqrtx1 <- 1+excent^2*tan(t11[i])^2
    sqrtx1 <- replace(sqrtx1,sqrtx1<0,0)
    x1 <- excent/sqrt(sqrtx1)
    sqrtx2 <- 1+excent^2*tan(t12[i])^2
    sqrtx2 <- replace(sqrtx2,sqrtx2<0,0)
    x2 <- excent/sqrt(sqrtx2)
    if(excent == 1){
      freq[i] <- abs(cos(t11[i])-cos(t12[i]))
    } else{
      alpha <- excent/sqrt(abs(1-excent^2))
      alpha2 <- alpha^2
      x12 <- x1^2
      x22 <- x2^2
      if(excent >1 ){
        # alpx1 <- sqrt(alpha2+x12)
        # alpx2 <- sqrt(alpha2+x22)
        sqrtalpx1 <- alpha2+x12
        sqrtalpx1 <- replace(sqrtalpx1,sqrtalpx1<0,0)
        alpx1 <- sqrt(sqrtalpx1)
        sqrtalpx2 <- alpha2+x22
        sqrtalpx2 <- replace(sqrtalpx2,sqrtalpx2<0,0)
        alpx2 <- sqrt(sqrtalpx2)
        dum <- x1*alpx1+alpha2*log(x1+alpx1)
        freq[i] <- abs(dum-(x2*alpx2+alpha2*log(x2+alpx2)))
      } else{
        # almx1 <- sqrt(alpha2-x12)
        # almx2 <- sqrt(alpha2-x22)
        sqrtalmx1 <- alpha2-x12
        sqrtalmx1 <- replace(sqrtalmx1,sqrtalmx1<0,0)
        almx1 <- sqrt(sqrtalmx1)
        sqrtalmx2 <- alpha2-x22
        sqrtalmx2 <- replace(sqrtalmx2,sqrtalmx2<0,0)
        almx2 <- sqrt(sqrtalmx2)
        dum <- x1*almx1+alpha2*asin(x1/alpha)
        freq[i] <- abs(dum-(x2*almx2+alpha2*asin(x2/alpha)))
      }
    }
  }
  sum0 <- sum(freq)
  freq0 <- freq/sum0
  df <- data.frame(freq=freq0,litab=litab)
  return(df)
}

#' Compute coefficients of interception function and volumn scattering function
#' @param tts solar zenith angle
#' @param tto veiwing zenith angle
#' @param psi relative azimuth angle
#' @param ttl leaf inclination angle
#' @return df a data fram with four columns: chi_s, chi_o, frho, ftau
#' @return chi_s interception functions
#' @return chi_o interception functions
#' @return frho function to be multiplied by leaf reflectance rho
#' @return ftau function to be multiplied by leaf transmittance tau
#' @description 
#' % Compute volume scattering functions and interception coefficients
#' % for given solar zenith, viewing zenith, azimuth and leaf inclination angle.
#' % chi_s and chi_o are the interception functions.
#' % frho and ftau are the functions to be multiplied by leaf reflectance rho and
#' % leaf transmittance tau, respectively, in order to obtain the volume scattering function.
#' % Wout Verhoef, april 2001, for CROMA
#' % REAL(KIND=8),INTENT(in) :: tts
#' % REAL(KIND=8),INTENT(in) :: tto
#' % REAL(KIND=8),INTENT(in) :: psi
#' % REAL(KIND=8),INTENT(in) :: ttl
#' % REAL(KIND=8),INTENT(inout) :: chi_s
#' % REAL(KIND=8),INTENT(inout) :: chi_o
#' % REAL(KIND=8),INTENT(inout) :: frho
#' % REAL(KIND=8),INTENT(inout) :: ftau
#' % 
#' % REAL(KIND=8) costs,costo,sints,sinto,cospsi
#' % REAL(KIND=8) psir
#' % REAL(KIND=8) costl,sintl,cs,co,ss,so,ds
#' % REAL(KIND=8) cosbts,cosbto,bts,bto
#' % REAL(KIND=8) btran1,btran2,bt1,bt2,bt3,t1,t2
#' % REAL(KIND=8) doo
#' % REAL(KIND=8) denom
volscatt <- function(tts, tto, psi, ttl){
  rd <- pi/180
  costs <- cos(rd*tts)
  costo <- cos(rd*tto)
  sints <- sin(rd*tts)
  sinto <- sin(rd*tto)
  cospsi <- cos(rd*psi)
  psir <- rd*psi
  costl <- cos(rd*ttl)
  sintl <- sin(rd*ttl)
  cs <- costl*costs
  co <- costl*costo
  ss <- sintl*sints
  so <- sintl*sinto
  # ..................................................................................
  # betas -bts- and betao -bto- computation
  # Transition angles (beta) for solar (betas) and view (betao) directions
  # if thetav+thetal>pi/2, bottom side of the leaves is observed for leaf azimuth
  # interval betao+phi<leaf azimut<2pi-betao+phi.
  # if thetav+thetal<pi/2, top side of the leaves is always observed, betao=pi
  # same consideration for solar direction to compute betas
  # ...................................................................................
  
  cosbts <- 5
  if(abs(ss)>1e-6){
    cosbts <- -cs/ss
  }
  cosbto <- 5
  if(abs(so)>1e-6){
    cosbto <- -co/so
  }
  
  if(abs(cosbts)<1){
    bts <- acos(cosbts)
    ds <- ss
  } else{
    bts <- pi
    ds <- cs
  }
  chi_s <- 2/pi*((bts-pi*0.5)*cs+sin(bts)*ss)
  
  if(abs(cosbto)<1){
    bto <- acos(cosbto)
    doo <- so
  } else if(tto<90){
    bto <- pi
    doo <- co
  } else{
    bto <- 0
    doo <- -co
  }
  chi_o <- 2/pi*((bto-pi*0.5)*co+sin(bto)*so)
  
  # ................................................................................
  # Computation of auxiliary azimuth angles bt1, bt2, bt3 used for
  # the computation of the bidirectional scattering coefficient w.
  # ................................................................................
  
  btran1 <- abs(bts-bto)
  btran2 <- pi-abs(bts+bto-pi)
  
  if(psir<=btran1){
    bt1 <- psir
    bt2 <- btran1
    bt3 <- btran2
  } else{
    bt1 <- btran1
    if(psir<=btran2){
      bt2 <- psir
      bt3 <- btran2
    } else{
      bt2 <- btran2
      bt3 <- psir
    }
  }
  
  t1 <- 2*cs*co+ss*so*cospsi
  t2 <- 0
  if(bt2>0){
    t2 <- sin(bt2)*(2*ds*doo+ss*so*cos(bt1)*cos(bt3))
  }
  denom <- 2*pi*pi
  frho <- ((pi-bt2)*t1+t2)/denom
  ftau <- (-bt2*t1+t2)/denom
  
  if(frho<0){
    frho <- 0
  }
  if(ftau<0){
    ftau <- 0
  }
  df <- data.frame(chi_s=chi_s, chi_o=chi_o, frho=frho, ftau=ftau)
  return(df)
}

#' These three functions (Jfun1, Jfun2, Jfun3) are used to avoid singularities in 4SAIL model.
Jfunc1 <- function(k,l,t){
  del <- (k-l)*t
  Jout <- vector(length = length(del))
  for(i in seq(del)){
    if(abs(del[i])>1e-3){
      Jout[i] <- (exp(-l[i]*t)-exp(-k*t))/(k-l[i])
    } else{
      Jout[i] <- 0.5*t*(exp(-k*t)+exp(-l[i]*t)*(1-del[i]*del[i])/12)
    }
  }
  return(Jout)
}
# Jfunc2
Jfunc2 <- function(k,l,t){
  Jout <- (1-exp(-(k+1)*t))/(k+l)
  return(Jout)
}
# Jfunc3
Jfunc3 <- function(k,l,t){
  Jout <- (1-exp(-(k+l)*t))/(k+l)
  return(Jout)
}