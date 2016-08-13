###functions for the CO2 model**********************************
#********v2016-08-13********************************************
#***************************************************************

#===============================================================
#Differentials
#alkalinity exchange
falk <- function(rr,rp)
{
    vAlk <- -0.15*(rr-rp)
    return (vAlk)
}

#CO2 exchange
fCO2 <- function(patm,kh,zCO2,k600,rr,rp)
{
  vCO2 <- 1000*(patm*kh-zCO2)*k600+rr-rp
  return (vCO2)
}

#DIC exchange
fDIC <- function (patm,kh,zCO2,rr,rp)
{
  vDIC <- 1000*(patm*kh-zCO2)*k600+rr-rp
  return (vDIC)
}

#isotope exchange for DIC
fRdic <- function(zdic,patm,kh,Ratm,zCO2,zRdic,k600,rr,rp)
{
  t1 <-rr
  t2 <- rp
  vRdic <- 1/(zdic)*((1000*(patm*kh*Ratm-zCO2*zRdic/adic)*akin*aatm*k600+rr*Rrr-rp*zRdic*ap)-zRdic*(1000*k600*(patm*kh-zCO2)+rr-rp))
  return (vRdic)
}

#isotope exchange for DIC
fRdic1 <- function(zdic,patm,kh,Ratm,zCO2,zRdic,k600,rr,rp)
{
  vRdic1 <- 1/(zdic)*((1000*(patm*kh*Ratm-zCO2*zRdic/adic)*aatm*k600)-zRdic*1000*k600*(patm*kh-zCO2))
  return (vRdic1)
}

#-----------------------------------------------------------------

#=================================================================
#speciation carbonate system
#calculate H from dic and co2
fH <- function(zdic,zCO2,k1,k2)
{
  #x(n+1)=x(n)-f(xn)/f'(xn)
  #whole equation
  #h^2*(zdic-zco2)-zco2*frmain.k1*h-zco2*frmain.k1*frmain.k2=0
  a1 <- (zdic-zCO2)
  a2 <- zCO2*k1
  a3 <- zCO2*k1*k2
  
  b1 <- 2*(zdic-zCO2)
  b2 <- zCO2*k1
  
  #precision of calculation
  prec <- 10^-10
  
  #newton rhapson solver here
  xn <- 1E-7 #guess pH start
  xn1 <- 0 #calc pH
  Da <- abs(xn-xn1)
  while(Da>prec)
  {
    z1 <- a1*xn^2
    z2 <- a2*xn
    z3 <- a3
    z4 <- b1*xn
    z5 <- b2
    
    xn1 <- xn-(a1*xn^2-a2*xn-a3)/(b1*xn-b2)
    Da <- abs(xn-xn1)
    xn <- xn1
  } 
  return (xn)
}

#calculate HCO3 from CO2 and H+
fHCO3 <- function(k1,zCO2,zh)
{
  #calculate HCO3
  #fhco3=(zDIC-(zDIC^2-zAlk*(2*zDIC-zAlk)*(1-4*frmain.k2/frmain.k1))^0.5)/(1-4*frmain.k2/frmain.k1)
  vhco3 <- k1*zCO2/zh
  return (vhco3)
}

#calulate CO3 from HCO3 and H
fCO3 <- function(k2,zhco3,zh)
  {
  co3 <- k2*zhco3/zh
  return (co3)
  }

#-------------------------------------------------------------------

#===================================================================
#CONSTANTS
#K1
fk1 <- function(t,sal)
  {
  #absolute temp
  t <- t+273.15
  
  #calculate K1
  a0 <- 13.4191
  a1 <- 0.0331
  a2 <- -0.0000533
  a3 <- -530.1228
  a4 <- -6.103
  a5 <- -2.0695
  a6 <- 0
  A <- a0*sal^0.5+a1*sal+a2*sal^2
  B <- a3*sal^0.5+a4*sal
  C <- a5*sal^0.5+a6*sal
  pk10 <- -126.34048+6320.813/t+19.568224*log(t)
  pk1_millero <- pk10+A+B/t+C*log(t)
  vk1 <<- 10^-pk1_millero
  return (vk1)
} 

#K2
fk2 <- function(t,sal)
  {
  #absolute temp
  t <- t+273.15
  
  #calculate K2
  a0 <-  21.0894
  a1 <-  0.1248
  a2 <- -0.0003687
  a3 <- -772.483
  a4 <- -20.051
  a5 <- -3.32254
  a6 <- 0
  A <- a0 * sal ^ 0.5 + a1 * sal + a2 * sal ^ 2
  B <- a3 * sal ^ 0.5 + a4 * sal
  C <- a5 * sal ^ 0.5 + a6 * sal
  pK20 <- -90.18333 + 5143.692 / t + 14.613358 * log(t)
  pK2_Millero <- pK20 + A + B / t + C * log(t)
  vk2 <<- 10 ^ -pK2_Millero
  return (vk2)
}

#Kh
fkh <- function(tc)
{
  #calculate KH
  pKH <- -0.0000000024372 * tc ^ 5 + 0.00000034041 * tc ^ 4 - 0.000015787 * tc ^ 3 + 0.00017829 * tc ^ 2 + 0.015355 * tc + 1.1102
  khv <- 10 ^ -pKH
  return (khv)
}

fifracdic <- function(fco2,fhco3,fco3)
  {
  
  #alternative way to calculate DIC form the distribution
  #of each species in DIC, CO2, HCO3 and CO3
  edic <<- fhco3*ehco3+fco3*eco3+eatm 
  adic <<-(edic/1000+1) 
  
}

#give all isotope fractionation factors
fifrac <- function(tc)
  {
  Rstd <- 0.0112372
  
  #fractionation factors
  #-0.0049*temp-1.31
  #for ratios this means (-0.0049*temp/1000+1)*Rst
  eatm <<- -0.0049*tc-1.31
  aatm <<- (eatm/1000+1)
  aco2 <<- aatm
  
  #kinetic isotope effect
  #Zhang et al (1995) -0.81 at 21 degrees and -0.95 at 5 degrees
  ekin <<- 0 #should be -0.9
  akin <<- (ekin/1000+1)
  
  #read isotope fractionation during production
  apc <<- (-13/1000+1)*Rstd
  ap <<- apc/Rstd
  
  #isotope fractionation co2(g) to HCO3
  #from Zhang (1995)
  #-0.141*temp+10.78, adding co2(g) to co2(aq)
  ehco3 <<- -0.141*tc+10.78-eatm
  ahco3 <<- (ehco3/1000+1)
  
  #isotope fractionation co2(g) to CO3
  #from Zhang (1995)
  #-0..052*temp+7.22 adding co2(g) to co2(aq)
  eco3 <<- -0.052*tc+7.22-eatm
  aco3 <<- (eco3/1000+1)
  
  #isotope fractionation CO2(g) to DIC
  #from Zhang (1995)
  #-0.014*tc*fco3-0.105*T+10.53
  #since fco3 is small compared to the effect of temp fco3 is omitted
  edic <<- 10.53-0.105*tc
  adic <<- (edic/1000+1)
  
}

#-------------------------------------------------------------

#=============================================================
#HELP FUNCTIONS
ival <- function(Rs)
  {
  Rstd <- 0.0112372
  iv <- (Rs/Rstd-1)*1000
  return (iv)
}

#Calculate d13C of respired DOC
fDOCin <- function(sal)
{
  #from Barbaras paper, end members, marine -23.1/terr -28.1 and mixing
  #calculation of d13C respired from salinity relation
  Rstd <- 0.0112372
  fracT <- 1.0688*sal^-0.25246
  fracM <- 1-fracT
  if (fracT>1)
    {
    fracT <- 1
    fracM <- 1-fracT
    }
  calc_d_resp <- fracT*-28.1+fracM*-23.1
  calcR <- (calc_d_resp/1000+1)*Rstd
  return (calcR)
}

#calculate d13C in CO2
f13CO2 <- function(iDIC,fCO2,fhco3,fco3)
{
  ivCO2 <- iDIC-fhco3*ehco3-fco3*eco3
  return (ivCO2)
} 

#--------------------------------------------------------------
#==============================================================&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
fCO2modrp <- function(x0)
{
  
  rp.conc <- x0[2]
  rr.conc <- x0[1]
  
  #respiration per day/m2 and volume and volume of water
  #rr.conc <- 0.7 #umol/L correspond to 184 gr C per m2 and year in BP E Gustafsson
  #average depth in BP is 61.7 m
  rr.mass <- rr.conc*1e-6*V #respiration rate per day mol/m2
  rrc.conc <<- rr.conc*1e-6 #in mol/L
  #===============================
  
  #CONSTANT#
  #production per day mol/m2 and volume of water
  rp.mass <- rp.conc*1e-6*V #production rate per day/m2
  rpc.conc <<- rp.conc*1e-6 #in mol/L
  #===============================
  
  ssum <- 0 #for sum of squares of d13C in DIC
  t <- 0 #time
  
  #start loop for differentials
  while(no<=mxt)
  {
    
    #eq 1 CO2 exchange
    #fCO2 <- function(patm,kh,zCO2,k600,rr,rp)
    odCO2.mass <- dCO2.mass
    dCO2.conc <- odCO2.mass/V
    dCO2.mass <- odCO2.mass+dt*fCO2(patm,kh,dCO2.conc,k600,rr.mass,rp.mass)
    #===============================
    
    #eq 2 Alk 
    #falk <- function(rr,rp)
    oalk.mass <- dalk.mass
    dalk.conc <- oalk.mass/V
    dalk.mass <- oalk.mass+dt*falk(rr.mass,rp.mass)
    #===============================
    
    #eq 3 DIC change with time
    #fDIC <- function (patm,kh,zCO2,rr,rp)
    oddic.mass <- ddic.mass
    ddic.conc <- oddic.mass/V
    ddic.mass <- oddic.mass+dt*fDIC(patm,kh,dCO2.conc,rr.mass,rp.mass)
    #===============================
    
    #eq 4 calculate pH as a function of Alk and DIC
    #fH <- function(zdic,zCO2,k1,k2)
    dh <- fH(ddic.conc,dCO2.conc,k1,k2)
    dph <- -log10(dh)
    #===============================
    
    #eq 5
    #calculate hco3
    #fHCO3 <- function(k1,zCO2,zh)
    dhco3.conc <- fHCO3(k1,dCO2.conc,dh)
    #===============================
    
    #eq 6
    #calculate CO3
    #fCO3 <- function(k2,zhco3,zh)
    dco3.conc <- fCO3(k2,dhco3.conc,dh.conc)
    #===============================
    
    #calculate fractions
    fco2 <- dCO2.conc/ddic.conc #fraction CO2
    fhco3 <- dhco3.conc/ddic.conc #fraction HCO3
    fco3 <- dco3.conc/ddic.conc #fraction CO3
    fdic <- fco2+fhco3+fco3 #sum of above for check
    #===============================
    
    #eq 7
    #rate of change for R in DIC
    #(zdic,patm,kh,Ratm,zCO2,zRdic,k600,rr,rp)
    odRdic <- dRdic
    dRdic <- odRdic+dt*fRdic(ddic.mass,patm,kh,Ratm,dCO2.conc,odRdic,k600,rr.mass,rp.mass)
    #===============================
    
    #calculate d13C in DIC, CO2, HCO3, CO3
    #ival <- function(Rs)
    Rs <- dRdic
    ivdic <- ival(Rs)
    #===============================
    
    #d13C in CO2, HCO3 and CO3
    ivCO2 <- f13CO2(ivdic,fco2,fhco3,fco3)
    ivhco3 <- ivCO2+ehco3
    ivco3 <- ivCO2+eco3
    #===============================
    
    #store data in array for latest run
    dact[no,1] <<- t #time in days
    dact[no,2] <<- ivdic.start #start for d13C in DIC
    dact[no,3] <<- ivdic #isotope value for DIC
    dact[no,4] <<- ivco2.start #start for d13C in DIC
    dact[no,5] <<- ivCO2 #isotope value for DIC
    dact[no,6] <<- dCO2.conc*1e6 # to get umol/L
    dact[no,7] <<- ddic.conc*1e6 # to get umol/L
    dact[no,8] <<- dhco3.conc*1e6 # to get umol/L
    dact[no,9] <<- dco3.conc*1e6 # to get umol/L
    dact[no,10] <<- fco2 #fractions for each ion
    dact[no,11] <<- fhco3 #fractions for each ion
    dact[no,12] <<- fco3 #fractions for each ion
    dact[no,13] <<- dph #pH
    dact[no,14] <<- rr.mass #to get umol/m2 and day
    dact[no,15] <<- rp.mass # to get umol/m2 and day
    dact[no,16] <<- eatm #to get umol/m2 and day
    dact[no,17] <<- ehco3 # to get umol/m2 and day
    dact[no,18] <<- eco3 #to get umol/m2 and day
    dact[no,19] <<- ekin # to get umol/m2 and day
    dact[no,20] <<- fco2 #to get umol/m2 and day
    dact[no,21] <<- fhco3 # to get umol/m2 and day
    dact[no,22] <<- fco3 #to get umol/m2 and day
    dact[no,23] <<- fdic # to get umol/m2 and day
    #===============================
    
    #takes sum of isotope values and CO2 conc
    ivss <- ivdic.start
    ivcs <- ivdic
    
    #normalize to start values
    ivss <- ivss/ivss
    ivcs <- ivcs/ivss
    
    #take sum of squares
    ssum1 <<- (ivcs-ivss)^2
    #ssum1 <- 0 #only conc is used for this run
    
    #get data for sums
    co2ss <- co2.start*1e6
    co2cs <- dCO2.conc*1e6
    
    #normalize to start values
    co2ss <- co2ss/co2ss
    co2cs <- co2cs/co2ss
    
    #take sum of squares
    #ssum2 <<- (co2cs-co2ss)^2
    ssum2 <- 0
    
    #take final sum
    ssum <- ssum+ssum1+ssum2
    
    t <- t+dt
    no <- no+1
    #===============================
  }
  
  return(ssum)
    
}