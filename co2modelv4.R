#==========CO2 model for sea water================
#==========v2016-08-13================================
#*****************************************************
#*****************************************************

#important settings
#Line 64 number of days
#line 73 depth of water
#line 78 k600
#line 161 run sample no

#set working directory
#setwd("/Users/CMM/Dropbox/Rcode/CO2Code") #macbook pro 15 and iMac
setwd("/Users/CMM/Dropbox/Dropbox/Rcode/CO2Code") #macbook pro 13
#setwd("/Volumes/Duo/Dropbox/Rcode/CO2Code") #mac pro
#===============================

#call CO2f functions needed to run the model
source ("co2fv4.R")

#library for optimisation routine
library(nloptr)

#load se carb library
library(seacarb)
#===============================

#read in data frame
da <- read.csv("bbsw3.csv", sep=",")

#get number of rows in data frame da
nrmax <- nrow(da)
#input format
#1  id
#2	Temp
#3	sal
#4	wind
#5	co2 air
#6	co2 flux
#7	TA
#8	DIC
#9	pH
#10	d13C CO2 water
#11	d13C CO2 gas
#12	d13C Air
#13	d13C hco3
#14	d13C co3
#15	f co2
#16	f hco3
#17	f co3
#18	pco2
#19	hco3
#20	co3
#21	co2
#22	d13C.calc.DIC
#23	Basin
#24	BasinID
#===============================

#iteration
dt <- 0.1 #time step
h <- dt #traditional Euler step size
#===============================

#Run from day 0 to xx
FromDay <- 0 #day
ToDay <- 50 #days
mxt <- (ToDay-FromDay)/dt #number of days to run divided by time step gives array
#===============================

#isotope ratio of 13C/12C in pdb
Rstd <- 0.011237
#===============================

#volume of water and depth
d <- 20 #depth in meters
V <- d*1000 #depth multiplied by 1000 to give volume in liters
#===============================

#piston velocity for CO2
k600 <- 2 #m/day
#===============================

#PCO2 in air
co2_atm <- 10^-3.41
patm <- co2_atm
#===============================

#atmospheric d13C CO2 and ratio 13C/12C 
ivatm <- -8.5
Ratm=(ivatm/1000+1)*Rstd #atmospheric isotope ratio and value
#===============================

#estimated production umol/day and m2
#estrp <- 2 #umol production per L and day
#this number is for umol/L and day to get production per day
#multipli with volume (V) and 1E-6 to get mol and 12 to get gr
#rp.C.mass <- estrp*1e-6*12*V
#rp.conc <- estrp
#===============================

#estimated respiration umol/day and m2
#estrr <- 1 #0.75 umol /day per L
#this number is for umol/L and day to get production per day
#multiply with volume (V) and 1E-6 to get mol and 12 to get gr
#rr.C.mass <- estrr*1e-6*12*V
#rr.conc <- estrr
#===============================

#create array for run and copy with results for final values
dact <- data.frame(matrix(0, ncol = 23, nrow = mxt)) #latest run
#===============================

#set colnames
colnames(dact)[1] <- "t"
colnames(dact)[2] <- "d13C.dic.in"
colnames(dact)[3] <- "d13C.dic"
colnames(dact)[4] <- "d13C.co2.in"
colnames(dact)[5] <- "d13C.co2"
colnames(dact)[6] <- "CO2"
colnames(dact)[7] <- "dic"
colnames(dact)[8] <- "hco3"
colnames(dact)[9] <- "co3"
colnames(dact)[10] <- "fco2"
colnames(dact)[11] <- "fhco3"
colnames(dact)[12] <- "fco3"
colnames(dact)[13] <- "pH"
colnames(dact)[14] <- "rr"
colnames(dact)[15] <- "rp"
colnames(dact)[16] <- "eatm"
colnames(dact)[17] <- "ehco3"
colnames(dact)[18] <- "eco3"
colnames(dact)[19] <- "ekin"
colnames(dact)[20] <- "fco2"
colnames(dact)[21] <- "fhco3"
colnames(dact)[22] <- "fco3"
colnames(dact)[23] <- "fdic"
#===============================

#add extra column names 25 to 34
da[25] <- 0
da[26] <- 0
da[27] <- 0
da[28] <- 0
da[29] <- 0
da[30] <- 0
da[31] <- 0
da[32] <- 0
da[33] <- 0
da[34] <- 0
colnames(da)[25] <- "cpH.start"
colnames(da)[26] <- "cdic.conc.start"
colnames(da)[27] <- "cd13C.dic.start"
colnames(da)[28] <- "cd13C.dic.end"
colnames(da)[29] <- "cd13C.CO2.start"
colnames(da)[30] <- "cd13C.CO2.end"
colnames(da)[31] <- "cCO2.conc.end"
colnames(da)[32] <- "cdic.conc.end"
colnames(da)[33] <- "crr"
colnames(da)[34] <- "crp"
#===============================

#run sample no
nr <- 1
#===============================

#temperature in celsius=tc
tc <- da[nr,2] #read from frame
#===============================
    
#read salinity from table
sal <- da[nr,3] #read form frame
#===============================
  
#given ratio of 13C in DOC (respired material)
Rrr <- fDOCin(sal)
#===============================
  
#calc carbonate system
kh <- fkh(tc)
k1 <- fk1(tc,sal)
k2 <- fk2(tc,sal)
#===============================
  
#set isotope fractionation factors
fifrac(tc)
#===============================
  
#read from input sheet
dCO2in <- da[nr,18]*1e-6
dCO2.conc <- dCO2in*kh #CO2 dissolved in water start value
dCO2.mass <- dCO2.conc*V #mass CO2 in volume of water
#===============================
  
#read TA
dalk.conc <- da[nr,7]*1e-6  #alk in water
dalk.mass <- dalk.conc*V #mass of alkality in water
#===============================

#read DIC
ddic.conc <- da[nr,8]*1e-6 #conc of DIC in water
ddic.mass <- ddic.conc*V #mass of DIC in water
#===============================

#pH values and H+
#calculate pH from DIC and CO2
dh.conc <- fH(ddic.conc,dCO2.conc,k1,k2)
dph <- -log10(dh.conc)
da[nr,25] <- dph #store calculated pH
#===============================

#calculate conc of HCO3 and CO3
dhco3.conc <- fHCO3(k1,dCO2.conc,dh.conc)
dco3.conc <- fCO3(k2,dhco3.conc,dh.conc)
#===============================

#calculate DIC and Alk
dalk.conc <- dhco3.conc+2*dco3.conc-dh.conc
ddic.conc <- dCO2.conc+dhco3.conc+dco3.conc
da[nr,26] <- ddic.conc*1000000 #store dic
#===============================

#calculate fraction of each species
fco2 <- dCO2.conc/ddic.conc #fraction CO2
fhco3 <- dhco3.conc/ddic.conc #fraction HCO3
fco3 <- dco3.conc/ddic.conc #fraction CO3
fdic <- fco2+fhco3+fco3 #check that sum becomes 1
fifracdic(fco2,fhco3,fco3)
#===============================

#respiration per day/m2 and volume and volume of water
#rr.mass <- rr.conc*1e-6*V #respiration rate per day mol/m2
#rr.conc <- rr.conc*1e-6 #in mol/L
#===============================

#production per day mol/m2 and volume of water
#rp.mass <- rp.conc*1e-6*V #producwarningtion rate per day/m2
#rp.conc <- rp.conc*1e-6 #in mol/L
#===============================

#read isotope start values from sheet
ivCO2 <- da[nr,10] #read d13C CO2 in water
dRCO2 <- (ivCO2/1000+1)*Rstd #start value isotope ratio CO2 in water
ivCO2g <- ivCO2-eatm #(-0.0049*tc-1.31) start value gas CO2
#===============================

#calculate fractions for carbonate system
ivdic <- fco2*ivCO2+fhco3*(ivCO2+ehco3)+fco3*(ivCO2+eco3)
dRdic <- (ivdic/1000+1)*Rstd #start value DIC
#===============================

#target
#no change in dic data, i.e. d13C in dic should be the same as
#the start value
ivdic.start <<- ivdic
ivco2.start <<- ivCO2
dic.start <<- ddic.conc
co2.start <<- dCO2.conc
#===============================
  
no <- 1
#check that it is overpressure
#if finding correct rr here is is where the loop should start
#*********************************************************************************************
#*********************************************************************************************
#*********************************************************************************************

  
      #call model
      #fCO2modrr(rr.conc)
      #x0 <- c(0.2,0.4)
      #sres <- cobyla(x0, fCO2modrr, 0.2, 5, hin=NULL, nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000)) 
  
    #call model
    #fCO2modrp(rp.conc)
    #x0 <- c(0.2,0.4)
    
    # Solve using NLOPT_LN_COBYLA without gradient information
    #x0, x0[1]=rp, x0[2]=rr

    #start values and lower and upper bound
    #x0 <- c(0.25,0.1)
    #lb <- c(0.1,0.1)
    #ub <- c(1,0.1)

    res1 <- nloptr(x0=c(0.25,0.1),
                    eval_f=fCO2modrp,
                    lb=c(0.1,0.1),
                    ub=c(1,0.1),
                    opts = list("algorithm"="NLOPT_LN_COBYLA",
                                "xtol_rel" = 1e-8,
                                "maxeval" = 1000,
                                "print_level" = 1)
                   )
    print(res1)
    
    #sres <- cobyla(x0, fCO2modrp, 0.2, 5, hin=NULL, nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000)) 

#############END INTEGRATION LOOP#########################################
  
  #add from model run
  #write data to array final numbers
  da[nr,27] <- ivdic.start
  da[nr,28] <- ivdic
  da[nr,29] <- ivco2.start
  da[nr,30] <- ivCO2
  da[nr,31] <- dCO2.conc*1000000
  da[nr,32] <- ddic.conc*1000000
  da[nr,33] <- rrc.conc*1000000
  da[nr,34] <- rpc.conc*1000000
  #===============================

  #write da and dact data to file
  write.csv(da,file="da.csv")
  write.csv(dact,file="dact.csv")

  quartz("CO2")
  plot(dact$t,dact$CO2)
  
  quartz("DIC")
  plot(dact$t,dact$dic)
  
  quartz("d13C CO2")
  plot(dact$t,dact$d13C.co2)
  
  quartz("d13C DIC")
  plot(dact$t,dact$d13C.dic)
  
  quartz("pH")
  plot(dact$t,dact$pH)
