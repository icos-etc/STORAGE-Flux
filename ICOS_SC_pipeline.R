# *****************************************************************************
# 
# ICOS_SC_pipeline
# Processing pipeline code, settings and calls to task-specific code
# 
# author: Giacomo Nicolini
# contact: g.nicolini@unitus.it
# date: 2020-01-19
# For license information see LICENSE file
# 
# *****************************************************************************



# *****************************************************************************
## Libraries & Functions 

# Libraries
kpacks <- c(
  'xts',
  'gplots',
  'data.table',
  'birk',
  'googledrive',
  'data.table',
  'tidyverse',
  'crayon'
)   
new.packs <- kpacks[!(kpacks %in% installed.packages()[,"Package"])]
if(length(new.packs)){install.packages(new.packs)}
lapply(kpacks, require, character.only=T)
remove(kpacks, new.packs)

# Environment options
options(scipen = 999) 


# Functions
source(paste0(getwd(), '/ICOS_SC_Functions.R'))



# *****************************************************************************

## Read the configuration file
source(paste0(getwd(), '/ICOS_SC_config.R'))



# Data and metadata collection ================================================

## Meteo data .................................................................

if('ICOS_SC_meteo_data.R' %in% list.files(getwd())){
  source('ICOS_SC_meteo_data.R')
  Meteo_dwnld(site.ID, t0 = t0, tX = tX, out_dir = BM.data.dir, remove_raw = TRUE, remove_intermediate = TRUE)
}

# Check and list the meteo file
meteo.files <- list.files(BM.data.dir, pattern = paste0(site.ID, '_meteo_sc.csv'), full.names = TRUE)

if(length(meteo.files) == 0){
  cat(red('ERROR: meteo data must be provided, covering the processing period.\n'))
  cat(red('The file must be a comma separeted ASCII with the following variables/columns:\n'))
  cat(silver('TIMESTAMP_START, TIMESTAMP_END, PA_H_V_R, RH_H_V_R, TA_H_V_R\n'))
  cat(red('where "_H_V_R" are the position indexes.\n'))
}


## Storage (profile) data .....................................................

if('ICOS_SC_storage_data.R' %in% list.files(getwd())){
  source('ICOS_SC_storage_data.R')
  # download and manipulate profile data
  Storage_dwnld(site.ID, t0 = t0, tX = tX, out_dir = SC.data.dir)
}
# check and list the storage file
SC.files <- list.files(SC.data.dir, full.names = TRUE)

if(length(SC.files) == 0){
  cat(red('ERROR: storage profile data must be provided, covering the processing period.\n'))
  cat(red('The file must be a comma separeted ASCII with - at least - the following variables/columns:\n'))
  cat(silver('TIMESTAMP, CO2, H2O, LEVEL, T_CELL, PRESS_CELL\n'))
  cat(red('Please see the storage data instruction document for detail.\n'))
}


# BADM ........................................................................

if('ICOS_SC_badm.R' %in% list.files(getwd())){
  # download and manipulate storage, meteo and EC BADM
  source('ICOS_SC_badm.R')
  BADM_dwnld(site.ID)
}

# check and list the storage file
BADM.files <- list.files(MD.dir, full.names = FALSE)

if(!all(grepl('CO2|H2O|ECMD|BMMD', BADM.files))){
  cat(red('ERROR: at least the BADMs for the storage system of CO2 and H2O, and the BADMs for the EC and meteo (TA, RH and PA) systems must be provided.\n'))
  # NOTE: EXAMPLE MUST BE PROVIDED 
}else{
  BADM <- list()
  BADM <- lapply(list.files(MD.dir, full.names = TRUE), function(x) fread(x, header = T, sep = ',', data.table = F))
  names(BADM) <- sub('|_', '', sub(site.ID, '', sub('.csv', '', BADM.files)))
}



# Meteo profile characteristics ...............................................

## Air temperature
TA.indx <- grep('TA', BADM$'BMMD'$'BM_VARIABLE_H_V_R')
# Number of sampling levels
N.LEV.TA <- length(BADM$'BMMD'$'BM_HEIGHT'[TA.indx])
# Heights of sampling levels
z.i.TA <- BADM$'BMMD'$'BM_HEIGHT'[TA.indx]
# Heights of the level-specifc boundary layers
zl.i.TA <- c(0, c(z.i.TA[-N.LEV.TA] + z.i.TA[-1])/2, max(z.i.TA))
# Boundary air layer depth
Dz.i.TA <- diff(zl.i.TA)

## Air relative humidity
RH.indx <- grep('RH', BADM$'BMMD'$'BM_VARIABLE_H_V_R')
# Number of sampling levels
N.LEV.RH <- length(BADM$'BMMD'$'BM_HEIGHT'[RH.indx])
# Heights of sampling levels
z.i.RH <- BADM$'BMMD'$'BM_HEIGHT'[RH.indx]
# Heights of the level-specifc boundary layers
zl.i.RH <- c(0, c(z.i.RH[-N.LEV.RH] + z.i.RH[-1])/2, max(z.i.RH))
# Boundary air layer depth
Dz.i.RH <- diff(zl.i.RH)

## Air pressure
PA.indx <- grep('PA', BADM$'BMMD'$'BM_VARIABLE_H_V_R')
# Number of sampling levels
N.LEV.PA <- length(BADM$'BMMD'$'BM_HEIGHT'[PA.indx])
# Heights of sampling levels
z.i.PA <- BADM$'BMMD'$'BM_HEIGHT'[PA.indx]
# Heights of the level-specifc boundary layers
zl.i.PA <- c(0, c(z.i.PA[-N.LEV.PA] + z.i.PA[-1])/2, max(z.i.PA))
# Boundary air layer depth
Dz.i.PA <- diff(zl.i.PA)


# Constants ...................................................................

AP.0 <- 101.325 * exp(-BADM$'ECMD'$'ALTITUDE'[1]/8200) # # Air pressure at the tower base [kPa]
R. <- 8.3144598        # Universal gas constant [J K-1 mol-1 or m3 Pa K-1 mol-1]
M.air.d <- 0.0289647   # Molar mass of dry air [kg mol-1] (28.9647 g mol-1)
M.H2O <- 0.0180153     # Molar mass of H2O [kg mol-1] (18.01528 g mol-1 )
R.H2O <- R. / M.H2O    # Gas constant for water vapour [J kg-1 K-1]
R.d <- R. / M.air.d    # Gas constant for dry air [J kg-1 K-1]
g <- 9.8               # Gravitational acceleration [m s-2]
Dt <- 1800             # EC averaging time [s]




# ******************************************************************************
# METEO PROCESSING =============================================================
# . compute half-hourly aggregated values used for SC computation
# . calculate HEAT storage
# ******************************************************************************

# Load data
BM.data <- fread(meteo.files, header = TRUE, na.strings = '-9999', fill = TRUE, blank.lines.skip = TRUE, data.table = FALSE)


# Air temperature profile [°C]
TA <- BM.data[, grep('TA_', names(BM.data), fixed = TRUE)]
# cut-out data over the EC system
TA <- TA[, names(TA) %in% BADM$'BMMD'$'BM_VARIABLE_H_V_R'[TA.indx]]
# order the variable according to the (ascending) sampling height (from N = ground to 1 = tower-top)
TA <- TA[, BADM$'BMMD'$'BM_VARIABLE_H_V_R'[TA.indx]]
# check wheter is measured along a profile or not. * NOTE: it must be a profile *
TA.on.profile <- ifelse(ncol(TA) == 1, FALSE, TRUE)
if(!TA.on.profile){cat(red('Error: air temperature must be measured along a profile.\n')); stop()}

# Air relative humidity (can be a profile) [%]
RH <- BM.data[, grep('RH_', names(BM.data), fixed = TRUE)]
# if measured along a profile
if(is.data.frame(RH)) {
  # cut-out data over the EC system
  RH <- RH[, names(RH) %in% BADM$'BMMD'$'BM_VARIABLE_H_V_R'[RH.indx]]
  # order the variable according to the (ascending) sampling height (from N = ground to 1 = tower-top)
  RH <- RH[, BADM$'BMMD'$'BM_VARIABLE_H_V_R'[RH.indx]]
  RH.on.profile <- TRUE
}else{
  RH.on.profile <- FALSE
}

# Air pressure (can be a profile) [kPa]
PA <- BM.data[, grep('PA_', names(BM.data), fixed = TRUE)]
# if measured along a profile
if(is.data.frame(PA)) {
  # cut-out data over the EC system
  PA <- PA[, names(PA) %in% BADM$'BMMD'$'BM_VARIABLE_H_V_R'[PA.indx]]
  # order the variable according to the (ascending) sampling height (from N = ground to 1 = tower-top)
  PA <- PA[, BADM$'BMMD'$'BM_VARIABLE_H_V_R'[PA.indx]]
  PA.on.profile <- TRUE
}else{
  PA.on.profile <- FALSE
}
# check if PA is actually reported in kPa
# 1 atm = 98.0665 kPa, 980.665 hPa, 98066.5 Pa
if(abs(1 - mean(PA, na.rm=T)/100000) < 0.01){ # pascal
  PA <- PA * 0.001
}else if((abs(1 - mean(PA, na.rm=T)/1000) < 0.01)){  # hectopascal
  PA <- PA * 0.1
}


# Meteorological data aggregation .............................................
# NOTE: meteo data are averaged over the full half-hour (referring to the FOLLOWING one)

# Air temperature
# Average air temperature profile - weighted by the depths of the air layers representative of the individual measurement points
TA.avg <- apply(TA, 1, function(x) weighted.mean(x, w = Dz.i.TA, na.rm = TRUE))


# Air relative humidity 
if(!is.null(RH)){ # 1. if air pressure is reported somehow
  if(RH.on.profile){ # 1.1 if measured along the profile
    # Average air pressure along the profile - weighted by the depths of the air layers representative of the individual measurement points
    RH.avg <- apply(RH, 1, function(x) weighted.mean(x, w = Dz.i.RH, na.rm = TRUE))
  }else{ # 1.2 if measured at a single point
    # Average air pressure along the profile estimated at the RH measuring point
    RH.avg <- RH
  }
}else{           # 2. if air humidity is not reported
  cat(red('Error: at least one air relative humidity must be reported'))
  stop()
}


# Air pressure
if(!is.null(PA)){ # 1. if air pressure is reported somehow
  if(PA.on.profile) { # 1.1 if measured along the profile
    # Average air pressure along the profile - weighted by the depths of the air layers representative of the individual measurement points
    PA.avg <- apply(PA, 1, function(x) weighted.mean(x, w = Dz.i.PA, na.rm = TRUE))
  }else{ # 1.2 if measured at a single point
    # Average air pressure along the profile estimated at the profile MIDDLE point
    if(max(z.i.TA) < z.i.RH){
      D.H.PA <- max(z.i.TA)/2 + z.i.PA # height difference between the PA sampling level and the profile middle
      PA.avg <- PA * (1 + D.H.PA / (R.d * c(TA.avg + 273.15) / g))  # note that the formula starts from AP at the top and the gradient is therefore added (downward increment)
    }else{
      D.H.PA <- max(z.i.TA)/2 - z.i.PA # height difference between the PA sampling level and the profile middle
      PA.avg <- PA * (1 - D.H.PA / (R.d * c(TA.avg + 273.15) / g))  # note that the formula starts from AP at the top and the gradient is therefore added (downward increment)
    }
  }
}else{           # 2. if air pressure is not reported
  # 2.1 # Average air pressure along the profile estimated at the profile MIDDLE point
  # estimated based on site altitude (Campbell and Norman, 1998)
  D.H.PA <- max(z.i.TA)/2 # profile middle point
  PA.avg <- AP.0 * (1 - D.H.PA / (R.d * c(TA.avg + 273.15) / g))
}


# H2O molar fraction - resulting from hygrometers data - Profile weighted average
# NOTE: if measured along the profile, must be used to compute the dry air molar density and correct the storage flux
#       (it can be used to cross-check the irga data anyway)
H2O.MF.avg.meteo <- (6088.484 * exp(TA.avg/(TA.avg+237.6429)*17.31303) / 18.0153) * RH.avg / (PA.avg/R.*(TA.avg+273.15))


# Save halfhourly and vertically integaretd meteo data for storage flux compuation
# it can be either mantained or removed after storage processing, according to user's choice
BM.data.avg <- data.frame(NULL)
BM.data.avg <- data.frame(
  TIMESTAMP = BM.data[, 'TIMESTAMP_END'],
  TA = unname(TA.avg),
  PA = unname(PA.avg),
  RH = unname(RH.avg),
  H2O_MF_meteo = H2O.MF.avg.meteo
)
write.csv(BM.data.avg, paste0(BM.data.dir, site.ID, '_meteo_sc_AVG.csv', sep=''), quote = FALSE, row.names = FALSE)





# *****************************************************************************
# STORAGE PROCESSING ==========================================================
# Data are processed according to the actual sampling scheme
# *****************************************************************************


# . SENSIBLE HEAT storage  ----------------------------------------------------

# source('ICOS_SCH.R') # UNDER DEVELOPMENT




# . GHG storage  ---------------------------------------------------------------
# Verify how many gases are actually measured basing on provided SC BADM
sc.gases <- grep('CO2|H2O|CH4|N2O', sub('_', '', sub('...D', '', sub('|_', '', sub(site.ID, '', sub('.csv', '', BADM.files))))), value = T)

# create output data structures
GAS.avg.pro.c <- list(NULL) # AVERAGE PROFILE CONCENTRATIONS
SC.GAS <- list(NULL) # STORAGE FLUXES (GHGs, H and LE)

i.gas <- 0
for(i.gas in 1:length(sc.gases)){
  
  # set current gas
  cur.GAS <- sc.gases[i.gas]
  
  # set current gas/profile BADM 
  cur.BADM.SCMD <- BADM[[grep(sc.gases[i.gas], names(BADM))]]
  
  # storage sampling scheme
  STO_CONFIG <- unique(cur.BADM.SCMD$'STO_CONFIG')[1]
  
  ## Profile characteristics ...................................................
  # Number of sampling levels
  N.LEV <- length(na.omit(cur.BADM.SCMD$'STO_PROF_LEVEL'))
  # Levels ID
  ID.LEV <- c(na.omit(cur.BADM.SCMD$'STO_PROF_LEVEL'))
  # Heights of sampling levels
  z.i <- c(na.omit(cur.BADM.SCMD$'STO_PROF_HEIGHT'))
  # Heights of the level-specifc boundary layers
  zl.i <- c(0, c(z.i[-N.LEV] + z.i[-1])/2, max(z.i))
  # Boundary air layer depth
  Dz.i <- diff(zl.i)
  names(Dz.i) <- c(na.omit(cur.BADM.SCMD$'STO_PROF_LEVEL'))
  
  # Compute line volumes from the lines switch to the GA (closed path only) [m3]
  if(!is.null(cur.BADM.SCMD$'STO_GA_TUBE_LENGTH')){ 
    line.vol.dw <- cur.BADM.SCMD$'STO_GA_TUBE_LENGTH' * (cur.BADM.SCMD$'STO_GA_TUBE_DIAM' * 10^-3)^2 * pi
  }
  
  # Set colors for profile levels
  prof.col <- rich.colors(N.LEV)
  
  
  # list the storage files according to the sampling scheme
  if(!grepl('sep', STO_CONFIG, ignore.case = TRUE)){ # 'SEQUENTIAL' and 'SIMULTANEOUS' setup
    
    # Already done
  
    }else{ # 'Separate' setup
    
    # WIP 
    
  }
  
  # create receiving data structures
  PRO.GAS.avg.ls <- list(NULL) # level-averaged [c] profile at half-hour
  PRO.GAS.SD.ls <- list(NULL)  # level standard deviation [c] profile at half-hour
  SC.GAS.un <- c(NULL)         # UNCORRECTED STORAGE FLUX 

  # Data import and manipulation **********************************************
  
  if(!grepl('separate', STO_CONFIG, ignore.case = TRUE)){ # 'Sequential' and 'Simultaneous' setup
    
    SC.data <- data.frame(NULL)
    SC.data <- as.data.frame(do.call('rbind', lapply(SC.files, function(x) fread(x, header = T, sep = ',', na.strings = '-9999', data.table = FALSE))))
    
  }else{ # 'Separate' setup
    
    # WIP
  }
  
  # Control on actual levels according to the EC measurement height and possibly subset 
  levels.in.data <- c(na.omit(unique(SC.data$'LEVEL')))
  levels.in.BADM <- c(na.omit(cur.BADM.SCMD$'STO_PROF_LEVEL'))
  if(length(levels.in.data) > length(levels.in.BADM)){
    more.level.on.data <- TRUE
    # subset the dataset
    SC.data[!SC.data$'LEVEL' %in% levels.in.BADM, 2:ncol(SC.data)] <- NA
  }else{
    more.level.on.data <- FALSE
  }
  
  ## Time specification
  # time object in UTC (POSIXlt)
  PRO.time <- strptime(SC.data$'TIMESTAMP', '%Y%m%d%H%M%S', tz='GMT')
  # refer the timestamps to the nearest half-hour
  PRO.time.hh <- floor.POSIXct.near(PRO.time, "30 mins")
  # Make half-hourly timestamp for storage flux
  SC.time <- seq(first(PRO.time.hh), last(PRO.time.hh), 1800)

  
  
  # Auxiliary data handling ***************************************************
  # 1 Hz mandatory auxiliary variables
  
  ## Flow rate
  # note that the check allows for user choice in case of missing data: continue using nominal flow rate or stop processing
  if(any(grepl('FLOW_VOLRATE', names(SC.data))) & !all(is.na(SC.data$'FLOW_VOLRATE'))){
    FLOW_VOLRATE <- SC.data[, grep('FLOW_VOLRATE', names(SC.data))] # [L min-1]
  }else{
    FR.prompt <- readline(prompt = cat(yellow("WARNING: In-line flow rate is not reported in the data. Doyou want to use nominal flow rate instead (Y/N)? ")))
    if(FR.prompt == 'Y'){
      FLOW_VOLRATE <- rep(na.omit(unique(cur.BADM.SCMD$'STO_PROF_BUFFER_FLOWRATE')), nrow(SC.data))
      cat(green('OK, processing will continue...'))
    }else{
      cat(red('OK, then processing will terminate here'))
      stop()
    }
  }
  
  ## Levels' number
  LEVEL <- SC.data[, grep('LEVEL', names(SC.data))] # [#] 
  
  ## GA pressure and temperature 
  if(any(grepl('GA_CP', cur.BADM.SCMD$'STO_GA_MODEL'))){ # of the GA cell in case of closed path analyzer
    PRESS_CELL <- SC.data[, grep('PRESS_CELL', names(SC.data))] # [kPa]
    T_CELL <- SC.data[, grep('T_CELL', names(SC.data))] # [°C]
  }else{ # open path analyzer
    stop('ETC DEBUG: PRESS and TEMP in case of OP are not yet handled')
  }
  
  ## Level's sampling time
  if(STO_CONFIG == 'Sequential'){
    # Actual level-specific average sampling time [s]
    # Note that is computed by median to allow for some flexibility and because of some site which repeat the sampling of the first level
    tau.LEVEL <- round(median(diff(which(diff(LEVEL) != 0)), na.rm = T), 0)
    # Time to complete a full profile sampling (Instruction: 60 - 300 s)
    tau.prof.integ <- tau.LEVEL * length(na.omit(cur.BADM.SCMD$'STO_PROF_LEVEL'))
  }
  
  
  # Concentrations or densities coversions ************************************
  
  # H2O .......................................................................
  if(i.gas == 1){ # do it once
    if(any(grepl('H2O_CONC', names(SC.data)))){ # gas data are expressed in concentration density [mmol m-3]
      H2O.CONC <- SC.data[, grep('H2O_CONC', names(SC.data), ignore.case = TRUE, value = TRUE)]
      H2O.MF <- H2O.CONC * R. * (T_CELL + 273.15) / (PRESS_CELL * 1000)
      H2O.MR <- H2O.MF / (1 - H2O.MF * 10^-3)
    }else if(any(grepl('H2O_DRY', names(SC.data)))){ # gas data are expressed in molar fraction in dry air [μmol mol-1]
      H2O.MR <- SC.data[, grep('H2O_DRY', names(SC.data), ignore.case = TRUE, value = TRUE)]
      H2O.MF <- H2O.MR * PRESS_CELL * 1000 / (R. * (T_CELL + 273.15))
    }else{ # gas data are expressed in molar fraction in humid air [μmol mol-1]
      H2O.MF <- SC.data[, grep('H2O', names(SC.data), ignore.case = TRUE, value = TRUE)]
    } 
  }
  
  # GAS ........................................................................
  if(cur.GAS != 'H2O'){
    cur.GAS.var <- grep(cur.GAS, names(SC.data), ignore.case = TRUE, value = TRUE)
    if(any(grepl('CONC', cur.GAS.var))){ # gas data are expressed in concentration density [mmol m-3]
      GAS.CONC <- SC.data[, cur.GAS.var]
      GAS.MF <- GAS.CONC * R. * (T_CELL + 273.15) / (PRESS_CELL * 1000)
      GAS.MR <- GAS.MF / (1 - H2O.MF * 10^-3)
    }else if(any(grepl('DRY', cur.GAS.var))){   # gas data are expressed in molar fraction in dry air [μmol mol-1]
      GAS.MR <- SC.data[, cur.GAS.var]
    }else{                                  # gas data are expressed in molar fraction in humid air [μmol mol-1]
      GAS.MF <- SC.data[, cur.GAS.var]
      GAS.MR <- GAS.MF / (1 - H2O.MF * 10^-3)
    } 
  }
  
  
  # Profile weighted average H2O molar fraction - from GA data - **************
  # NOTE: if RH is NOT measured along the profile, it MUST be used to compute the dry air molar density and correct the storage flux
  
  if(i.gas == 1){ # do it once
    
    if(STO_CONFIG == 'Separate'){
      H2O.MF.avg.i <- apply(H2O.MF, 2, tapply(x, PRO.time.hh, mean, na.rm = TRUE))
    }else{
      H2O.MF.avg.i <- aggregate(H2O.MF, by = list(LEVEL, as.character(PRO.time.hh)), FUN = mean, na.rm = TRUE, na.action = na.pass)
      names(H2O.MF.avg.i) <- c('level','timestamp','H2O')
    }
    
    # Profile average H2O weighted by the depths of the air layers representative of the individual measurement points
    if(STO_CONFIG == 'Separate'){
      H2O.MF.avg <- apply(H2O.MF.avg.i, 1, function(x) weighted.mean(x, w = rev(Dz.i), na.rm = TRUE))
    }else if(STO_CONFIG == 'Sequential'){
      H2O.MF.avg <- tapply(H2O.MF.avg.i$'H2O', H2O.MF.avg.i$'timestamp', function(x) weighted.mean(x, w = rev(Dz.i)[1:length(x)], na.rm = TRUE))
      # NOTE: it ahs been added "[1:length(x)]" to comply with possible differences in time blocks length (it is not correct but it can be accpeted as is referred to 1 or two data)
    }else{
      H2O.MF.avg <- H2O.MF.avg.i
    }
    
    # # crosscheck with 'meteo' H2O.MF
    # windows()
    # plot(BM.time, H2O.MF.avg.meteo, type='o', ylim=c(0,20))
    # points(SC.time, H2O.MF.avg, type='o', col=4, pch=19)
    # legend('topright', c('MF[H2O]_meteo','MF[H2O]_GA'), col=c(1,4), lty=1, pch=c(1,19))
    
  }
  
  
  # Some controls and checks ............................
  
  # |> Uppermost storage profile level
  # It should be at the same level of the EC system (tolerance: +- 5% of the EC height)
  if(any(max(z.i) < unique(BADM$'ECMD'$'SA_HEIGHT') - (unique(BADM$'ECMD'$'SA_HEIGHT') * 0.05) | 
         max(z.i) > unique(BADM$'ECMD'$'SA_HEIGHT') + (unique(BADM$'ECMD'$'SA_HEIGHT') * 0.05))){
    ST.H.prompt <- readline(prompt = cat(yellow("\nWARNING: The storage profile top level is not at the same height of the EC system. Do you want to continue processing anyway? (Y/N) ... ")))
    if(ST.H.prompt == 'Y'){
      cat(green('OK, processing will continue...'))
    }else{
      stop(red('Ok then processing terminates here because the storage profile top level is not at the same height of the EC system'))
    }
  }
  
  # Specifc tests to the sequential sampling scheme
  if(STO_CONFIG == 'Sequential'){
    
    # - Level-specific average sampling time [s]
    if(tau.LEVEL != median(cur.BADM.SCMD$'STO_PROF_SAMPLING_TIME', na.rm = T)){
      cat(yellow('Warning: actual sampling period is different from the nominal one\n'))
    }else{
      cat(green('OK, actual sampling period is on average coherent with the nominal one\n'))
    }
    
    # - Profile integration period (or return time) [s]
    if(tau.prof.integ > 300){
      cat(yellow('Warning: profile integration period is longer than those suggested in ICOS\n'))
    }else{
      cat(green('OK, profile integration period is in line with ICOS requirements\n'))
    }
    
    
  }
  
  # Here some other checks on submitted data: WIP
  
  
  
  # Profile averaged dry air density ******************************************
  # the average dry air density for the air column is used as scaling variable
  # to convert the kinematic storage term into a mass-based quantity
  # NOTE: the time-averaging is performed at half-hourly scale
  
  if(i.gas == 1){
    BM.data.avg <- fread(paste0(BM.data.dir, site.ID, '_meteo_sc_AVG.csv', sep=''), sep = ',', header = TRUE, fill = TRUE, blank.lines.skip = TRUE, data.table = FALSE)
    BM.time <- strptime(BM.data.avg$'TIMESTAMP', '%Y%m%d%H%M', tz='GMT')
    BM.data.avg$'TIMESTAMP' <- as.double(BM.data.avg$'TIMESTAMP')
    
    # align H2O MF and master (SC) timestamps and in case of missing values, fill with NA
    BM.data.avg.xts <- xts(BM.data.avg, order.by = BM.time) # covert to xts object
    BM.data.avg.xts <- merge(BM.data.avg.xts, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -(ncol(BM.data.avg.xts)+1)] # align
    BM.data.avg <- as.data.frame(coredata(BM.data.avg.xts)) # and get back to vector
    rm(BM.data.avg.xts)
    
    # Air molar volume [m3 mol-1] (should be around 0.025 m3 mol-1 at 25 °C)
    amv <- R. * (BM.data.avg$'TA' + 273.15) / (BM.data.avg$'PA'*10^3) # PA from KPa to Pa
    
    # Water vapor mass density [kg m-3] (e.g. range: 0.005 kgm-3 at 0 °C  - 0.051 kg m-3 at 40 °C)
    # NOTE: if RH is measured along the profile use the resulting H2O molar fraction, otherwise use the H2O measured by the GA
    if(RH.on.profile){
      rho.H2O <- c(BM.data.avg$'H2O_MF_meteo'*10^-3) * M.H2O / amv  # [H2O] from mmol to mol
    }else{
      rho.H2O <- c(H2O.MF.avg*10^-3) * M.H2O / amv  # [H2O] from mmol to mol
    }
    
    # Water vapor partial pressure [kPa] (should be around 2.338 kPa at 20 °C)
    e.avg <- rho.H2O * R.H2O*10^-3 * (BM.data.avg$'TA' + 273.15) # R.H2O from J kg-1 K-1 (i.e. Pa m3 kg-1 K-1) to kPa m3 kg-1 K-1
    # alternative formula : 0.61094*exp(17.625 * TA.avg / (TA.avg + 243.04)) # August-Roche-Magnus equation
    
    # Dry air partial pressure [kPa]
    PA.avg.d <- BM.data.avg$'PA' - e.avg
    
    ## Average dry air molar density [mol m-3]
    # NOTE: the dry air density is expressed in MOLAR density instead of MASS density
    rho.d <- PA.avg.d*10^3 / (R. * (BM.data.avg$'TA' + 273.15)) # PA form kPa to Pa
    
  }
  
  
  # Specific evaporation heat (lambda) [J kg-1] *******************************
  # calculated as a function of ambient air temperature
  if(i.gas == 1){
    seh <- 10^3 * (3147.5 - 2.37 * (BM.data.avg$'TA' + 273.15))
  }
  
  
  
  # ***************************************************************************
  # Storage flux conditional processing
  
  # set the current scalar data to be processed
  
  if(cur.GAS == 'H2O'){
    cur.scalar <- H2O.MF
  }else{
    cur.scalar <- GAS.MR
  }
  
  
  # .. SEPARATE sampling scheme -----------------------------------------------
  # Only 10 s of C values, crossing each half-hour, will be used for averaging
  # Time reference example:
  # Assuming that turbulent flux and storage flux estimates are to be computed for timestamps
  # 08:00:00, then by convention the 08:00:00 timestamp represents the flux corresponding to
  # observations between 07:30:00 and 08:00:00. SC will then be computed from the time 
  # average of measurements from 07:59:55 to 08:00:05 minus the time average
  # of measurements from 07:29:55 to 07:30:05
  
  if(STO_CONFIG == 'Separate'){
    
    # UNDER DEVELOPMNET
  }
  
  
  # .. SIMULTANEOUS sampling scheme -------------------------------------------
  # Only 10 s of C values, crossing each half-hour, will be used for averaging
  # Time reference example:
  # Assuming that turbulent flux and storage flux estimates are to be computed for timestamps
  # 08:00:00, then by convention the 08:00:00 timestamp represents the flux corresponding to
  # observations between 07:30:00 and 08:00:00. SC will then be computed from the time 
  # average of measurements from 07:59:55 to 08:00:05 minus the time average
  # of measurements from 07:29:55 to 07:30:05
  
  if(STO_CONFIG == 'Simultaneous'){
    
    # UNDER DEVELOPMNET
  }
  
  
  # .. SEQUENTIAL sampling scheme ---------------------------------------------
  # 
  # Time reference example: TODO
  # 
  
  if(STO_CONFIG == 'Sequential'){
    
    # Average profile concentrations ********************************************
    
    cur.GAS.avg.pro.c <- aggregate(cur.scalar, list(LEVEL, PRO.time.hh), median, na.rm = TRUE)
    cur.GAS.avg.pro.c.list <- list(NULL)
    cur.GAS.avg.pro.c.list <- split(cur.GAS.avg.pro.c, cur.GAS.avg.pro.c$Group.1)
    cur.GAS.avg.pro.c.list.r <- Reduce(function(x, y) merge(x, y, by = "Group.2", all = TRUE), cur.GAS.avg.pro.c.list)
    cur.GAS.avg.pro.c.list.r[, grep('Group.1', names(cur.GAS.avg.pro.c.list.r))] <- NULL
    names(cur.GAS.avg.pro.c.list.r) <- c('TIMESTAMP', paste0(cur.GAS, '_', rev(ID.LEV)))
    GAS.avg.pro.c[[i.gas]] <- cur.GAS.avg.pro.c.list.r
    GAS.avg.pro.c[[i.gas]]$'TIMESTAMP' <- format(GAS.avg.pro.c[[i.gas]]$'TIMESTAMP', '%Y%m%d%H%M%S', tz='GMT')
    rm(cur.GAS.avg.pro.c, cur.GAS.avg.pro.c.list, cur.GAS.avg.pro.c.list.r)
    
    
    # Storage flux computation ************************************************
    
    PRO.GAS.avg.ls <- list(NULL)    # level-averaged [c] profile at half-hour
    PRO.GAS.SD.ls <- list(NULL)     # level standard deviation [c] profile at half-hour
    SC.GAS.un <- c(NULL)            # UNCORRECTED STORAGE FLUX (no air density correction)
    
    j <- 1  # loop on each half-hour
    # 1 half-hour step back at the end to avoid out of bounds errors at the last half-hour
    for(j in 1:(length(SC.time)-1)){ 
      
      # find the current series indexes, respective of the WHOLE current half-hour (1800 records)
      hh.indxs <- which(PRO.time.hh %in% SC.time[j])
      
      if(length(hh.indxs) != 0 & !all(is.na(LEVEL[hh.indxs])) & j != 1){
        
        # Check on lines flow-rates (requirements: 10% maximum difference among lines), [L min-1]
        lines.FR <- tapply(FLOW_VOLRATE[hh.indxs], LEVEL[hh.indxs], median, na.rm = TRUE)
        # A. between current levels (may be less than the wahole profile levels)
        if(any(abs(diff(c(lines.FR, lines.FR[1])))/median(lines.FR) > 0.1*median(lines.FR))){
          cat(lines.FR[abs(diff(c(lines.FR,lines.FR[1])))/median(lines.FR) > 0.1])
          cat(yellow('WARNING: There is more than 10% difference among lines flow rates'))
        }
        # B. nominal flow-rate, maximum allowed 20%
        cur.lev.indx <- which(as.character(cur.BADM.SCMD$'STO_PROF_LEVEL') %in% unlist(dimnames(lines.FR)))
        if(any(abs(lines.FR - cur.BADM.SCMD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx])/lines.FR > 0.2*median(cur.BADM.SCMD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx]))){
          cat(lines.FR[abs(lines.FR - cur.BADM.SCMD$'STO_BUFFER_FLOW_RATE')/lines.FR > 0.2])
          cat(yellow('WARNING: There is more than 20% difference among actual and theoretical lines flow rates'))
        }
        
        # Flushing time: according to line volume and flow rate, [s] for each level
        tau.flush <- line.vol.dw[cur.lev.indx] / (lines.FR/60/1000)
        
        # Averaging time: seconds used for [c] averaging at each level (Protocol: 5 - 15 s) for each level
        tau.avging <- round(tau.LEVEL - tau.flush)
        # check:
        if(!all(is.na(tau.avging)) & any(tau.avging < 5)){
          cat(tau.avging[tau.avging < 5])
          cat(red('Levels averaging time is too short..'))
          stop()
        }
        # if averaging time is long enough, extend the flushing time and limit the averaging to 25 seconds
        if(!all(is.na(tau.avging)) & any(tau.avging > 25)){
          tau.avging[tau.avging > 25] <- 25
        }
        
        # Profile of average [c] for each half-hour
        GAS.c <- list(NULL) # list with the last [c] (N=tau.avging) of each level (used for averaging)
        
        # Find which level crosses the current half-hour
        # NOTE: the case where at the current half-hour there are not any level sampled, it is not currently handled
        crossing.level.inx <- NULL
        crossing.level.inx <- which(PRO.time == SC.time[j])
        # level crossing the current half-hour
        crossing.level <- LEVEL[crossing.level.inx]
        
        # if it is found, proceed with SC computation
        if(any(crossing.level.inx) & !is.na(crossing.level)){
          
          # if(see.data){ # visualize the current [CO2] frame
          #   # windows()
          #   plot.win <- c((crossing.level.inx-200) : (crossing.level.inx+200)) # to fit to the current half-hour (200 points windows)
          #   plot(as.double(PRO.time[plot.win]), unlist(cur.scalar)[plot.win],
          #        col=prof.col[LEVEL[plot.win]-min(LEVEL, na.rm = TRUE)+1], pch=19, cex=1.5)
          #   points(as.double(PRO.time[hh.indxs]), unlist(cur.scalar)[hh.indxs])
          #   abline(v=as.double(SC.time[j]), lwd=9, col=adjustcolor('grey', 0.6))
          #   legend('topleft', legend=c(rev(ID.LEV), 'ref. h/2'), col=c(prof.col, adjustcolor('grey', 0.4)),
          #          pch=c(rep(19,N.LEV), NA), pt.cex=c(rep(1.5,N.LEV), NA), lwd=c(rep(NA,N.LEV), 9), bty='n')
          # }
          
          jj <- 1 # loop on each level 
          for(jj in 1:N.LEV){
            
            # 1. find the timestamp of the record of this level which is closest to the current half-hour
            cur.lev.ts <- PRO.time[LEVEL == ID.LEV[jj]]
            level.ts <- cur.lev.ts[which.closest(cur.lev.ts, SC.time[j])]
            # and respective index
            level.indx <- which(PRO.time == level.ts)
            
            # if(see.data){abline(v=as.double(level.ts), col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=1, lty=2)}
            
            if(ID.LEV[jj] != crossing.level){ # IF THE CURRENT LEVEL IS NOT THE CROSSER
              
              # 2. extract the level [C] series (at 1 Hz) to average
              
              # 2.1. IF the TS is BEFORE the current half-hour:
              if(level.ts <  SC.time[j]){ 
                GAS.c[[jj]] <-  cur.scalar[(level.indx - (tau.avging[jj]-1)) : level.indx] # consider the PREVIOUS N seconds for averaging
                # if(see.data){
                #   abline(v=as.double(PRO.time[(level.indx - (tau.avging[jj]-1))]), lty=1, col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=2)
                #   abline(v=as.double(PRO.time[level.indx]), lty=1, col=prof.col[ID.LEV[jj]], lwd=2)
                #   rect(xleft = as.double(PRO.time[(level.indx - (tau.avging[jj]-1))]), ybottom = 300,
                #        xright = as.double(PRO.time[level.indx]), ytop = 800, col=adjustcolor(prof.col[ID.LEV[jj]-min(ID.LEV)+1], 0.3))
                # }
                
                # 2.2.IF the TS is AFTER the current half-hour:
              }else{ # if(level.ts >=  TS.ST.master.double.OP[j]) 
                GAS.c[[jj]] <- cur.scalar[(level.indx + tau.LEVEL - tau.avging[jj]) : (level.indx + (tau.LEVEL-1))] # consider the LAST N seconds (of that level) for averaging
                # if(see.data){
                #   abline(v=as.double(PRO.time[(level.indx + tau.LEVEL - tau.avging[jj])]), lty=1, col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=2)
                #   abline(v=as.double(PRO.time[(level.indx + (tau.LEVEL-1))]), lty=1, col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=2)
                #   rect(xleft = as.double(PRO.time[(level.indx + tau.LEVEL - tau.avging[jj])]), ybottom = 300,
                #        xright = as.double(PRO.time[(level.indx + (tau.LEVEL-1))]), ytop = 800, col=adjustcolor(prof.col[ID.LEV[jj]-min(ID.LEV)+1], 0.3), border = NA)
                # }
              }
              
            }else{  # IF THE CURRENT LEVEL IS THE CROSSER
              
              # find the swithching point to the NEXT level
              sw.pts <- which(diff(LEVEL[level.indx:(level.indx+tau.LEVEL*3/2)]) != 0) - 1 # -1 because of diff
              
              if(length(sw.pts)>1){ # this loop is to skip situations in which a level is sampled twice
                level.indx.c <- level.indx + sw.pts[2] - tau.LEVEL
              }else{
                level.indx.c <- level.indx + sw.pts
              }
              
              GAS.c[[jj]] <- cur.scalar[(level.indx.c - (tau.avging[jj]-1)) : level.indx.c] # consider the PREVIOUS N seconds for averaging
              # if(see.data){
              #   abline(v=as.double(PRO.time[level.indx.c]), lty=1, col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=2)
              #   abline(v=as.double(PRO.time[(level.indx.c - (tau.avging[jj]-1))]), lty=1, col=prof.col[ID.LEV[jj]-min(ID.LEV)+1], lwd=2)
              #   rect(xleft = as.double(PRO.time[(level.indx.c - (tau.avging[jj]-1))]), ybottom = 300,
              #        xright = as.double(PRO.time[level.indx.c]), ytop = 800, col=adjustcolor(prof.col[ID.LEV[jj]-min(ID.LEV)+1], 0.3), border = NA)
              # }
            }
            
          }
          
          # 3. current [c] average profile
          names(GAS.c) <- ID.LEV
          PRO.GAS.avg.ls[[j]] <- unlist(lapply(GAS.c, function(x) median(unlist(x), na.rm = TRUE)))
          # and respective standard deviation (2 sigma)
          PRO.GAS.SD.ls[[j]] <- 2 * unlist(lapply(GAS.c, function(x) sd(unlist(x), na.rm = TRUE)))
          
          
        }else{ 
          
          PRO.GAS.avg.ls[[j]] <- rep(NA, N.LEV); names(PRO.GAS.avg.ls[[1]]) <- ID.LEV
          PRO.GAS.SD.ls[[j]] <- rep(NA, N.LEV); names(PRO.GAS.SD.ls[[1]]) <- ID.LEV
          
        }
        
      }else{
        
        PRO.GAS.avg.ls[[j]] <- rep(NA, N.LEV); names(PRO.GAS.avg.ls[[1]]) <- ID.LEV
        PRO.GAS.SD.ls[[j]] <- rep(NA, N.LEV); names(PRO.GAS.SD.ls[[1]]) <- ID.LEV
        
      }
      
    } # HALF-HOUR LOOP END
    
    # assign time reference
    names(PRO.GAS.avg.ls) <- as.character(SC.time[-c(length(SC.time))])
    
    # UNCORRECTED STORAGE FLUX ************************************************
    D.GAS <- apply(do.call(rbind, PRO.GAS.avg.ls), 2, diff)
    SC.GAS.un <- apply(D.GAS, 1, function(x) sum((x/Dt) * Dz.i))

    # Time referencing the time series
    # convert to xts object
    SC.GAS.un <- xts(SC.GAS.un, order.by = as.POSIXct(names(SC.GAS.un), '%Y-%m-%d %H:%M:%S', tz='GMT'), tzone = 'GMT') # covert to xts object
    # and order according to the SC timestamp, possibly adding NAs
    SC.GAS.un <- merge(SC.GAS.un, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -2] # align
    SC.GAS.un <- coredata(SC.GAS.un)

  }
  
  
  
  # CORRECTED STORAGE FLUX (mass-based quantity) ****************************
  
  # CO2, CH4 or N2O storege flux
  if(any(grepl('CO2|CH4|N2O', cur.GAS))){
    SC.GAS[[i.gas]] <- SC.GAS.un * rho.d
  }
  
  # Latent heat storage flux [Wm-2]
  if(cur.GAS == 'H2O'){
    SC.GAS[[i.gas]] <- SC.GAS.un * rho.d * seh * M.H2O *10^-3
  }
  
  # Sensible heat storage flux [Wm-2]
  # WIP 
  
  
  windows()
  plot(SC.time, SC.GAS[[i.gas]], type='o', pch =19)
  
  
}






# Output file -----------------------------------------------------------------

# renaming SC list 
names(SC.GAS) <- paste0('SC', sc.gases)
# initialize the full output dataframe
SC.output <- data.frame(NULL)
SC.output <- data.frame(
  TIMESTAMP = format(SC.time, '%Y%m%d%H%M%S', tz='GMT')
)
# add storage fluxes
SC.output <- as.data.frame(cbind(SC.output, do.call(cbind, SC.GAS)))
names(SC.output)[2:ncol(SC.output)] <- paste0('SC', sc.gases)
# add average concentration profiles
for(i in 1:length(GAS.avg.pro.c)){
  SC.output <- merge(SC.output, GAS.avg.pro.c[[i]], by='TIMESTAMP')
}
# cut the dataset to the desired period
if(is.null(t0) | is.null(tX)){ # last day only
  yesterday <- paste0(format(Sys.Date()-1, '%Y%m%d'), '030000')
  SC.output <- SC.output[which(as.character(SC.output$'TIMESTAMP') == yesterday) : nrow(SC.output), ]
}else{ # defined period
  SC.output <- SC.output[which(SC.output$'TIMESTAMP' == as.character(t0)) : which(SC.output$'TIMESTAMP' == as.character(tX)), ]
}
# save it
fwrite(SC.output, paste0(OUT.dir, 'ICOS_SC_fulloutput.csv'), sep = ',')

