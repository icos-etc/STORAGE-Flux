
#'' ***************************************************************************
# STORAGE FLUX COMPUTATION
# + SEQUENTIAL SAMPLING SCHEME  +
# 
# author: Giacomo Nicolini
# contact: g.nicolini@unitus.it
# date: 2020-06-02
# For license information see LICENSE file
# 
#'' ***************************************************************************

stfcomp_SEQ <- function (
  site.ID, # station ID
  gas, # current gas name (e.g. 'CO2')
  Dt = 1800, # time interval between the two profile samplings (generally 1800 s, the default)
  PRO.data.path, # path to profile data folder
  SC.MD,  # path to profile meta-data folder
  EC.height, # height of the eddy covariance system
  BM.data.path,  # path to raw meteo data folder
  BM.data.AVG.path,  # path to profile averaged meteo data folder
  RH.profile = FALSE, # set TRUE if air relative humidity is sampled along a profile
  unc = 'SD', # storage flux uncertainty. Can be quantified either from the standard deviation ('SD', default) or standard error ('SE') of the concentrations.
  plot.profile.data = FALSE, # set TRUE to save an image if the averaged concentrations profile. Warning one figure is created for each averaging period, use it for debugging purposes only.
  plot.conc.data = FALSE, # set TRUE to save an image if the actually averaged concentrations. Warning one figure is created for each averaging period, use it for debugging purposes only.
  plot.SC.flux = TRUE, # set TRUE to save an image of the storage flux time series and daily cycle.
  warns = FALSE, # set TRUE to display warning messages related to controls on profile data.
  write.output = TRUE, # set TRUE to save numerical output files.
  output.dir = NULL  # output folder path
){
  
  # Install required packages if they are not already installed ..............
  if (!require("data.table")) {install.packages("data.table"); library(data.table)}
  if (!require("lubridate")) {install.packages("lubridate"); library(lubridate)}
  if (!require("bigleaf")) {install.packages("bigleaf"); library(bigleaf)}
  if (!require("birk")) {install.packages("birk"); library(birk)}
  if (!require("xts")) {install.packages("xts"); library(xts)}
  if (!require("gplots")) {install.packages("gplots"); library(gplots)}
  if (!require("crayon")) {install.packages("crayon"); library(crayon)}
  
  options(scipen=9999)
  
  # Useful constants
  Rgas <- bigleaf.constants()$Rgas  # Universal gas constant [J K-1 mol-1 or m3 Pa K-1 mol-1]
  Rd <- bigleaf.constants()$Rd      # Gas constant for dry air (Rgas/Md) [J kg-1 K-1]
  Rv <- bigleaf.constants()$Rv      # Gas constant for water vapour [J kg-1 K-1]
  Md <- bigleaf.constants()$Md      # Molar mass of dry air [kg mol-1] (28.9647 g mol-1)
  Mw <- bigleaf.constants()$Mw      # Molar mass of H2O [kg mol-1] (18.01528 g mol-1 )
  
  cat(cyan('Calculating SC flux for ' %+% gas %+% ' ...\n'))
  
  
  # **************************************************************************
  # Profile data import and manipulation
  
  cat(cyan('Getting and manipulating profile data ...\n'))
  
  SC.files <- list.files(PRO.data.path, pattern = '_processed', full.names = TRUE)
  if (length(SC.files) == 0) {
    cat(red('ERROR: 1 Hz profile data must be provided, covering the processing period.\n'))
    cat(yellow('The files must be daily, comma separeted ASCII with - at least - the following variables/columns:\n'))
    cat(silver('TIMESTAMP, CO2, H2O, LEVEL, FLOW_VOLRATE, T_CELL, PRESS_CELL\n'))
    cat(silver('Please see the storage data instruction document for details.\n'))
    stop()
  }
  
  ## Whole data period dataset (about 2.5 GB prof-1 Y-1)
  PRO.data <- data.frame(NULL)
  PRO.data <- as.data.frame(data.table::rbindlist(lapply(SC.files, function(x) 
    suppressWarnings(fread(x, na.strings = '-9999', header = TRUE, data.table = FALSE))
  )))
  PRO.time <- as.POSIXct(as.character(PRO.data$'TIMESTAMP'), '%Y%m%d%H%M%S', tz='GMT') # Profile data original timestamp, time object in UTC (POSIXct)
  PRO.data.f <- median(diff(PRO.time))     # native data frequency [s], it must be 1 second
  stopifnot(PRO.data.f == 1L)
  rm(SC.files)
  
  # Check and possibly correct for time-series jumps
  system.time({
    if (any(diff(PRO.time) != PRO.data.f)) {
      
      PRO.time.cont <- c(NULL)
      PRO.time.cont <- data.frame(TIMESTAMP_tmp = seq(PRO.time[1], PRO.time[length(PRO.time)], by = PRO.data.f)) # Profile data desired timestamp (pure 1 Hz frequency)
      PRO.data$'TIMESTAMP_tmp' <- PRO.time
      PRO.data <- merge(PRO.data, PRO.time.cont, by = 'TIMESTAMP_tmp', all = TRUE) # Align
      PRO.time <- c(NULL)  # Reset the profile data original timestamp
      PRO.time <- PRO.time.cont$'TIMESTAMP_tmp'
      PRO.data$'TIMESTAMP' <- format(PRO.time, '%Y%m%d%H%M%S')
      PRO.data$'TIMESTAMP_tmp' <- NULL
      rm(PRO.time.cont)

      # Find and drop possible duplicated timestamps values (WARNING: this may remove good - although badly acquired - data)
      if (any(duplicated(PRO.data$'TIMESTAMP'))){ PRO.data <- PRO.data[!duplicated(PRO.data$'TIMESTAMP'), ] }

    }
  })
  
  # Refer the timestamps to the *next* period (generally 30 minutes)
  PRO.time.Dt <- lubridate::ceiling_date(PRO.time, paste(Dt/60, 'minutes'))

  # Benchmark timestamp for the storage flux (generally 30 minuts resolution)
  SC.time <- seq(PRO.time.Dt[1], PRO.time.Dt[length(PRO.time.Dt)], Dt)

  cat(cyan('Done.\n'))
  
  
  # **************************************************************************
  # Check & set profile characteristics
  
  cat(cyan('Get/set profile characteristics ...\n'))
  
  ## Possibly limit the profile to the EC system height (up to the closest level)
  profile.top <- SC.MD$'STO_PROF_HEIGHT'[birk::which.closest(SC.MD$'STO_PROF_HEIGHT', EC.height)] 
  if (any(na.omit(SC.MD$'STO_PROF_HEIGHT') > profile.top)) {
    SC.MD <- SC.MD[- which(SC.MD$'STO_PROF_HEIGHT' > profile.top), ]
  }
  
  ## Levels ID (possibly recoding values in data)
  if (!all(is.na(SC.MD$'STO_PROF_LEVEL_ICOS'))) { # if level's ID were corrected in the BADM check
    # Recode data values
    level_key <- c(SC.MD$'STO_PROF_LEVEL_ICOS')[!is.na(SC.MD$'STO_PROF_LEVEL')]; names(level_key) <- SC.MD$'STO_PROF_LEVEL'[!is.na(SC.MD$'STO_PROF_LEVEL')]
    PRO.data[, grep('LEVEL', names(PRO.data))] <- dplyr::recode(PRO.data[, grep('LEVEL', names(PRO.data))], !!!level_key)
    # log the change
    cat(yellow('WARNING: levels ID are not consistent with levels heights and have been recoded as:\n' ))
    cat(silver(paste0(SC.MD$'STO_PROF_LEVEL'[!is.na(SC.MD$'STO_PROF_LEVEL')], ' >> ', SC.MD$'STO_PROF_LEVEL_ICOS'[!is.na(SC.MD$'STO_PROF_LEVEL')], "\n")))
    cat(yellow('Levels ID on data have been recoded accordingly.\n' ))
    fwrite(data.frame(ID_orig = SC.MD$'STO_PROF_LEVEL', ID_ICOS = SC.MD$'STO_PROF_LEVEL_ICOS'), paste0(output.dir,'Profile_recoding.txt'), sep='\t', quote = F, row.names = F)
    
    levels.in.BADM <- c(na.omit(SC.MD$'STO_PROF_LEVEL_ICOS'))    # store levels badm ID
    N.LEV <- length(na.omit(SC.MD$'STO_PROF_LEVEL_ICOS')) # Number of sampling levels
    ID.LEV <- c(na.omit(SC.MD$'STO_PROF_LEVEL_ICOS')) # Levels ID
    
  } else {
    
    levels.in.BADM <- c(na.omit(SC.MD$'STO_PROF_LEVEL'))    # store levels badm ID
    N.LEV <- length(na.omit(SC.MD$'STO_PROF_LEVEL')) # Number of sampling levels
    ID.LEV <- c(na.omit(SC.MD$'STO_PROF_LEVEL')) # Levels ID
    
  }
  
  # set colors for profile levels
  prof.col <- gplots::rich.colors(N.LEV) 
  
  # Calculate referece heights
  z.i <- c(na.omit(SC.MD$'STO_PROF_HEIGHT')); names(z.i) <- ID.LEV # Heights of sampling levels
  zl.i <- c(0, c(z.i[-N.LEV] + z.i[-1])/2, max(z.i)) # Heights of the level-specifc boundary layers
  Dz.i <- diff(zl.i); names(Dz.i) <- ID.LEV # Boundary air layer depth

  # Levels ID data series
  LEVEL <- PRO.data[, grep('LEVEL', names(PRO.data))] # [#] 
  
  
  ## Control on actual levels number 
  # according those reported in the badm, and possibly trim the dataset accordingly
  levels.in.data <- c(na.omit(unique(PRO.data$'LEVEL')))
  if (!all(levels.in.data %in% levels.in.BADM)) {
    more.level.on.data <- TRUE
    PRO.data[!PRO.data$'LEVEL' %in% levels.in.BADM, 2:ncol(PRO.data)] <- NA
    cat(silver('NOTE: the levels reported in the data files are more than the required ones according the the EC height.\n'))
    cat(silver('The levels above the EC system will not be used for storage flux computation.\n'))
  } else {
    more.level.on.data <- FALSE
  }
  
  ## Flow rate (mandatory auxiliary variables) series
  if (!all(grepl('_OP', SC.MD$'STO_GA_MODEL'))) { 
    if(any(grepl('FLOW_VOLRATE', names(PRO.data))) & !all(is.na(PRO.data$'FLOW_VOLRATE'))){
      FLOW_VOLRATE <- PRO.data[, grep('FLOW_VOLRATE', names(PRO.data))] # [L min-1]
    }else{
      FLOW_VOLRATE <- rep(na.omit(unique(SC.MD$'STO_PROF_BUFFER_FLOWRATE')), nrow(PRO.data))
      if(warns) cat(yellow("WARNING: In-line flow rate is not reported in the data. Nominal flow rate was used instead\n"))
    }
  }

  ## Line volumes from the lines switch to the GA (closed path only) [m3]
  if (!is.null(SC.MD$'STO_GA_TUBE_LENGTH')) { 
    line.vol.dw <- na.omit(unique(SC.MD$'STO_GA_TUBE_LENGTH' * (SC.MD$'STO_GA_TUBE_DIAM' * 10^-3)^2 * pi))[1]
    line.vol.dw <- rep(line.vol.dw, N.LEV)
  }
  
  ## Level's sampling time [s], calculated as median of the level's frequency
  tau.LEVEL <- round(median(diff(which(diff(LEVEL) != 0)), na.rm = T), 0)
  if (tau.LEVEL != median(SC.MD$'STO_PROF_SAMPLING_TIME', na.rm = T) & warns) {
    cat(yellow(paste('WARNING: the actual sampling period (', tau.LEVEL, ' s) is different from the nominal one (', median(SC.MD$'STO_PROF_SAMPLING_TIME', na.rm = T),' s)\n', sep='')))
  }
  
  # Profile integration period (or return time). Average time to complete a full profile sampling (instruction: 60 - 300 s)
  tau.PROF <- tau.LEVEL * N.LEV
  if(tau.PROF > 330 & warns){ # if more then 10 % of ICOS recommendations
    cat(yellow(paste('WARNING: the current profile integration period (', tau.PROF, ' s) is more then 10% longer than the one suggested in ICOS (300 s)\n', sep=''))) # currently disabled as it will go in a log-file
  }
  
  
  # TODO: Further controls to implement
  # * buffer volume tau >= tau.PROF
  # * x < tau.LEVEL IQr < y : plot(table(diff(which(diff(LEVEL) != 0)))); abline(v=quantile(diff(which(diff(LEVEL) != 0)))[c(2,4)], col=2)
  # * ...   
  
  
  ## GA pressure and temperature 
  if (!all(grepl('_OP', SC.MD$'STO_GA_MODEL'))) { # of the GA cell in case of closed path analyzer
    PRESS_CELL <- PRO.data[, grep('PRESS_CELL', names(PRO.data))] # [kPa]
    T_CELL <- PRO.data[, grep('T_CELL', names(PRO.data))] # [°C]
  }else{ # open path analyzer
    stop('ETC WORK IN PROGRESS: PRESS and TEMP in case of OP are not handled yet')
  }
  
  
  cat(cyan('Done.\n'))
  
  
  
  # **************************************************************************
  # Concentrations or densities coversions
  cat(cyan('Concentrations data and auxiliary variables calculation ...\n'))
  
  # H2O ......................................................................
  if (any(grepl('H2O_CONC', names(PRO.data)))) { # gas data are expressed in concentration density [mmol m-3]
    H2O.CONC <- PRO.data[, grep('H2O_CONC', names(PRO.data), ignore.case = TRUE, value = TRUE)]
    H2O.MF <- H2O.CONC * Rgas * (T_CELL + 273.15) / (PRESS_CELL * 1000)
    H2O.MR <- H2O.MF / (1 - H2O.MF * 10^-3)
  } else if (any(grepl('H2O_DRY', names(PRO.data)))) { # gas data are expressed in molar fraction in dry air [μmol mol-1]
    H2O.MR <- PRO.data[, grep('H2O_DRY', names(PRO.data), ignore.case = TRUE, value = TRUE)]
    H2O.MF <- H2O.MR * PRESS_CELL * 1000 / (Rgas * (T_CELL + 273.15))
  } else { # gas data are expressed in molar fraction in humid air [μmol mol-1]
    H2O.MF <- PRO.data[, grep('H2O', names(PRO.data), ignore.case = TRUE, value = TRUE)]
  }
  
  
  # GAS ......................................................................
  if(gas != 'H2O'){
    gas.var <- grep(gas, names(PRO.data), ignore.case = TRUE, value = TRUE)
    if(any(grepl('CONC', gas.var))){ # gas data are expressed in concentration density [mmol m-3]
      GAS.CONC <- PRO.data[, gas.var]
      GAS.MF <- GAS.CONC * Rgas * (T_CELL + 273.15) / (PRESS_CELL * 1000)
      GAS.MR <- GAS.MF / (1 - H2O.MF * 10^-3)
    }else if(any(grepl('DRY', gas.var))){   # gas data are expressed in molar fraction in dry air [μmol mol-1]
      GAS.MR <- PRO.data[, gas.var]
    }else{                                  # gas data are expressed in molar fraction in humid air [μmol mol-1]
      GAS.MF <- PRO.data[, gas.var]
      GAS.MR <- GAS.MF / (1 - H2O.MF * 10^-3)
    } 
  }
  
  # free some memory
  rm(PRO.data)
  
  
  
  # **************************************************************************
  # Profile weighted average H2O molar fraction - from GA data - 
  # NOTE: if RH is NOT measured along the profile, this MUST be used to compute the dry air molar density and correct the storage flux
  
  ## Drop concentration data during line flushing 
  # set to NA concentration values during flush (considering the last x seconds of measurements in each level)
  lev.switch.indx <- which(diff(data.table::nafill(LEVEL, 'nocb')) != 0) + 1
  # nominal flushing time [s] in case of closed-path GA
  if (!all(grepl('_OP', SC.MD$'STO_GA_MODEL'))) {
    tau.flush.nom <- ceiling(median(line.vol.dw, na.rm = TRUE) / (median(FLOW_VOLRATE, na.rm = TRUE)/60/1000))
  } else {
    tau.flush.nom <- 0
  }
  # guessed time to stabilize the concentration at each level sampling, 25% of level's sampling time
  tau.stab <- ceiling(tau.LEVEL * 0.25) 
  if (tau.LEVEL - (tau.flush.nom + tau.stab) > 15) { # if at least 15 seconds remains to average 
    flushing.indx <- c(sapply(lev.switch.indx, function(x) seq(x, x + tau.flush.nom + tau.stab))) # consider both times
  } else { # otherwise
    flushing.indx <- c(sapply(lev.switch.indx, function(x) seq(x, x + tau.flush.nom))) # only flushing time
  }
  # Drop flushing values
  H2O.MF[flushing.indx] <- NA
  
  # Level specific averages 
  # NOTE: it takes time
  H2O.MF.avg.i <- data.frame(NULL)
  H2O.MF.avg.i <- aggregate(H2O.MF, list(LEVEL, PRO.time.Dt), median, na.rm = TRUE)
  names(H2O.MF.avg.i) <- c('level','timestamp','H2O')
  
  ## Profile average H2O weighted by the depths of the air layers representative of the individual measurement points
  # H2O.MF.avg.0 <- tapply(H2O.MF.avg.i$'H2O', H2O.MF.avg.i$'timestamp', function(x) weighted.mean(x, w = rev(Dz.i), na.rm = TRUE))
  H2O.MF.avg.0 <- tapply(H2O.MF.avg.i$'H2O', H2O.MF.avg.i$'timestamp', function(x) weighted.mean(x, w = rev(Dz.i)[1:length(x)], na.rm = TRUE))
  # NOTE: the "[1:length(x)]" constraint has been added avoid possible crashes given by differences in time blocks length (it is not correct but it happen very rarely and affect a few data points)
  
  # Possibly align the timeline to SC timestamp (fill passible NAs)
  H2O.MF.avg.xts <- xts(H2O.MF.avg.0, order.by = as.POSIXct(names(H2O.MF.avg.0), '%Y-%m-%d %H:%M:%S', tz='GMT')) 
  H2O.MF.avg <- merge(H2O.MF.avg.xts, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -(ncol(H2O.MF.avg.xts)+1)] # align
  # row.names(H2O.MF.avg) <- format(time(H2O.MF.avg.xts), '%Y%m%d%H%M%S')
  colnames(H2O.MF.avg) <- 'H2O.MF.avg'
  rm(H2O.MF.avg.xts, H2O.MF.avg.0, H2O.MF.avg.i)
  # object.size(H2O.MF.avg)/1024^2
  
  
  # **************************************************************************
  ## Compute auxiliary micrometeorological variables 
  # Based profile-wheighted half-hourly averages as computed according to the "process_BM" function
  
  BM.data.avg <- fread(paste0(BM.data.AVG.path, site.ID, '_BM_SC_AVG.csv', sep=''), sep = ',', header = TRUE, fill = TRUE, blank.lines.skip = TRUE, data.table = FALSE)
  
  # Possibly align H2O MF and master (SC) timestamps and in case of missing values, fill with NA
  BM.time <- as.POSIXct(as.character(BM.data.avg$'TIMESTAMP'), '%Y%m%d%H%M', tz='GMT')
  # BM.data.avg$'TIMESTAMP' <- as.double(BM.data.avg$'TIMESTAMP')
  BM.data.avg.xts <- xts(BM.data.avg, order.by = BM.time) # covert to xts object
  BM.data.avg <- merge(BM.data.avg.xts, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -(ncol(BM.data.avg.xts)+1)] # align
  rm(BM.data.avg.xts)
  
  # crosscheck with H2O.MF from GA
  # windows()
  # plot(SC.time, H2O.MF.avg, type='o', col=4, pch=19, ylim=c(mean(H2O.MF.avg, na.rm = T)-15,mean(H2O.MF.avg, na.rm = T)+15))
  # points(BM.time, BM.data.avg$'H2O_MF_meteo', type = 'o')
  # legend('bottomleft', c('MF[H2O]_meteo','MF[H2O]_GA'), col=c(1,4), lty=1, pch=c(1,19))
  # plot(coredata(H2O.MF.avg), BM.data.avg$'H2O_MF_meteo', xlim=c(0,20), ylim=c(0,20)); abline(0,1)
  # plot(coredata(H2O.MF.avg) - coredata(BM.data.avg$'H2O_MF_meteo'), type='s'); abline(h=mean(H2O.MF.avg - coredata(BM.data.avg$'H2O_MF_meteo'), na.rm=T), col=2, lty=32)
  
  
  # Dry air density ..........................................................
  # the average dry air density for the air column is used as scaling variable to convert the kinematic storage term into a mass-based quantity
  # NOTE: the time-averaging is performed at half-hourly scale
  
  # Air molar volume [m3 mol-1] (should be around 0.025 m3 mol-1 at 25 °C)
  amv <- Rgas * (BM.data.avg$'TA' + 273.15) / (BM.data.avg$'PA'*10^3) # PA from KPa to Pa
  
  # Water vapor mass density [kg m-3] (e.g. range: 0.005 kgm-3 at 0 °C  - 0.051 kg m-3 at 40 °C)
  # NOTE: if RH is measured along the profile use the resulting H2O molar fraction, otherwise use the H2O measured by the GA
  if(RH.profile){
    rho.H2O <- c(BM.data.avg$'H2O_MF_meteo'*10^-3) * Mw / amv  # [H2O] from mmol to mol
  }else{
    rho.H2O <- c(H2O.MF.avg*10^-3) * Mw / amv  # [H2O] from mmol to mol
  }
  
  # Water vapor partial pressure [kPa] (should be around 2.338 kPa at 20 °C)
  e.avg <- rho.H2O * Rv*10^-3 * (BM.data.avg$'TA' + 273.15) # Rv from J kg-1 K-1 (i.e. Pa m3 kg-1 K-1) to kPa m3 kg-1 K-1
  
  # Dry air partial pressure [kPa]
  PA.avg.d <- BM.data.avg$'PA' - e.avg
  
  ## Average *dry* air *molar* density [mol m-3]
  rho.d <- PA.avg.d*10^3 / (Rd * (BM.data.avg$'TA' + 273.15) * Md) 
  names(rho.d) <- 'rho'
  
  # Free some memory
  rm(PA.avg.d, e.avg, rho.H2O, amv)
  
  
  # Specific evaporation heat (lambda) [J kg-1] ..............................
  seh <- 10^3 * (3147.5 - 2.37 * (BM.data.avg$'TA' + 273.15))
  names(seh) <- 'seh'
  
  
  # Dry air heat capacity at constant pressure (cp) [J kg-1 K-1] .............
  cp <- 1005 + ((BM.data.avg$'TA' + 23.12)^2) / 3364
  names(cp) <- 'cp'
  
  
  # Set the actual scalar data to process
  if(gas == 'H2O'){
    cur.scalar <- H2O.MF
  }else{
    cur.scalar <- GAS.MR
  }
  
  
  # ***************************************************************************
  ## Average profile concentrations
  # NOTE: this are not the **average** values used for the flux calculation, as are averaged over the full 30 half-hours.
  #       They are meant to be used for profile analysis.
  
  gas.pro <- data.frame(cur.scalar.or = cur.scalar, 
                        cur.scalar, 
                        LEVEL, 
                        PRO.time.Dt)
  # head(cbind(gas.pro, unname(replace(rep(0,nrow(gas.pro)), lev.switch.indx, 1))), 100)
  
  # Drop concentration data during line flushing (already calculated for H2O), set flushing values to NA
  gas.pro[flushing.indx, 'cur.scalar'] <- NA
  
  # save illustrative figure for the gas averaging strategy 
  jpeg(paste0(output.dir, site.ID, '_[', gas, ']_averaging_ex.jpeg'), quality = 100, res = 300, width=480*5, height=480*3)
  cpstart <- sample(1:nrow(gas.pro), 1)
  cpframe <- c(cpstart:(cpstart+500))
  plot(PRO.time[cpframe], gas.pro$cur.scalar.or[cpframe], las=1, xlab='Time [minutes]', ylab=paste0('[', gas,']'))
  abline(v=PRO.time[lev.switch.indx], lty=32)
  points(PRO.time[cpframe], gas.pro$cur.scalar[cpframe], col=2, pch=19)
  dev.off()
  
  # Concentration average by level and time (computed over the whole averaging period and referred to the end of the averaging period)
  gas.pro.avg.0 <- data.frame(NULL)
  gas.pro.avg.0 <- aggregate(cur.scalar ~ LEVEL + PRO.time.Dt, data = gas.pro, median, na.rm = TRUE)
  
  # Concentrations uncertainty expressed as standard deviation (or standard error) by level and time
  # (computed over the whole averaging period and referred to the end of the averaging period) 
  cur.scalar.unc <- numeric(nrow(gas.pro.avg.0))
  if (unc == 'SD') {
      cur.scalar.unc <- aggregate(cur.scalar ~ LEVEL + PRO.time.Dt, data = gas.pro, sd, na.rm = TRUE)[, 3]
  } else {
    cur.scalar.unc <- aggregate(cur.scalar ~ LEVEL + PRO.time.Dt, data = gas.pro, function(x) sd(x, na.rm = T)/sqrt(length(x)))$cur.scalar
  }
  gas.pro.avg.0$cur.scalar.unc <- cur.scalar.unc
  
  # rearrange columns
  gas.pro.avg.list <- split(gas.pro.avg.0, gas.pro.avg.0$'LEVEL')    # split by levels
  gas.pro.avg <- suppressWarnings(Reduce(function(x, y) merge(x, y, by = "PRO.time.Dt", all = TRUE), gas.pro.avg.list))    # merge back columnwise by 
  gas.pro.avg[, grep('LEVEL', names(gas.pro.avg))] <- NULL    # remove level's columns and rename
  names(gas.pro.avg) <- c('TIMESTAMP', sort(c(cbind(paste0(gas, '_', rev(ID.LEV)), paste0(gas, '_', rev(ID.LEV), '_sd.bt')))))
  gas.pro.avg$'TIMESTAMP' <- format(gas.pro.avg$'TIMESTAMP', '%Y%m%d%H%M%S', tz='GMT')
  
  
  # Concentrations uncertainty summed in quadrature by time: 
  # ** this will used as the "between" component of the storage flux uncertainty **
  SC.UNC.bt <- numeric(nrow(gas.pro.avg))
  SC.UNC.bt <- tapply(gas.pro.avg.0$cur.scalar.unc,  gas.pro.avg.0$PRO.time.Dt, function(x) sqrt(sum(x^2, na.rm = TRUE)))
 
  SC.UNC.bt <- xts(SC.UNC.bt, order.by = strptime(names(SC.UNC.bt), '%Y-%m-%d %H:%M:%S', tz='GMT')) # covert to xts object
  SC.UNC.bt <- merge(SC.UNC.bt, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -(ncol(SC.UNC.bt)+1)] # align
  
  # Free some memory
  rm(gas.pro, gas.pro.avg.list, gas.pro.avg.0, cur.scalar.unc) 
  
  
  # ***************************************************************************
  # STORAGE FLUX COMPUTATION
  cat(cyan('Computing storage flux ...\n'))
  
  
  # Create receiving data structures
  PRO.GAS.avg.ls <- list(NULL) # level-averaged [c] profile at half-hour
  PRO.GAS.unc.ls <- list(NULL)  # level standard deviation [c] profile at half-hour
  SC.GAS.k <- numeric(nrow(gas.pro.avg))         # UNCORRECTED STORAGE FLUX 
  SC.GAS <- numeric(nrow(gas.pro.avg))  # STORAGE FLUXES (GHGs and LE)
  
  # find the indexes of the profile timeline with respect to the SC timeline
  hh.indx.series <- c(NULL)
  hh.indx.series <- which(PRO.time %in% SC.time)
  
  # system.time({
    
    i.hh <- 12  # loop on each half-hour
    pb <- txtProgressBar(min = 0, max = length(SC.time), style = 3, char = '+', width = 50)
    for(i.hh in 1:length(SC.time)){ 
      
      setTxtProgressBar(pb, i.hh)
      
      # Find the current series indexes, respective of the WHOLE current period (1800 records in case of 30 minutes)
      # hh.indxs <- which(PRO.time.Dt %in% SC.time[i.hh]) # if PRO.time.Dt is rounded to the *closest* period
      cur.hh.indx <- integer()
      cur.hh.indx <- hh.indx.series[i.hh]
      hh.indxs <- integer(Dt)
      if (i.hh != length(SC.time)) {
          hh.indxs <- c((cur.hh.indx - Dt/2) : (cur.hh.indx + Dt/2))
      } else {
        hh.indxs <- (length(PRO.time)-Dt) : length(PRO.time)
      }
      
      if(length(hh.indxs) != 0 & !all(is.na(LEVEL[hh.indxs])) & i.hh != 1){
        
        # Check on lines flow-rates 
        # (requirements: 10% maximum difference among lines), [L min-1]
        lines.FR <- c(NULL)
        lines.FR <- tapply(FLOW_VOLRATE[hh.indxs], LEVEL[hh.indxs], median, na.rm = TRUE)
        # a. between current levels (may be less than the wahole profile levels) # TODO: probably send to a log-file
        if (any(abs(diff(c(lines.FR, lines.FR[1])))/median(lines.FR) > 0.1) & warns) {
          cat(lines.FR[abs(diff(c(lines.FR,lines.FR[1])))/median(lines.FR) > 0.1])
          cat(yellow(paste('WARNING: In the', SC.time[i.hh],'half-hour there is more than 10% difference among lines flow rates\n')))
        }
        # b. nominal flow-rate, maximum allowed difference 20%   # TODO: probably send to a log-file
        cur.lev.indx <- c(NULL)
        cur.lev.indx <- which(as.character(SC.MD$'STO_PROF_LEVEL') %in% unlist(dimnames(lines.FR)))
        if (any(abs(lines.FR - SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx]) / SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx] > 0.2) & warns) {
          levFR.out <- which(abs(lines.FR - SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx]) / SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx] > 0.2)
          cat(yellow(paste('WARNING: in the sampling line ', levFR.out , ' there is more than 20% difference (', 
                           round(abs(lines.FR - SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx]) / SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx], 2)*100,
                           '%) among actual and expected flow rate\n', sep='')))
        }
		# possibly fill missing values in flow rates by using the nominal ones
        if (any(is.na(lines.FR))) {
          lines.FR[is.na(lines.FR)] <- SC.MD$'STO_PROF_BUFFER_FLOWRATE'[cur.lev.indx][is.na(lines.FR)]
        }
        
        # Actual flushing time: according to line volume and flow rate, [s] for each level
        tau.flush <- c(NULL)
        tau.flush <- line.vol.dw[cur.lev.indx] / (lines.FR/60/1000)
        
        ### Level's averaging time
        # seconds used for [c] averaging at each level (Protocol: 5 - 15 s) for each level
        tau.avging <- c(NULL)
        tau.avging <- floor(tau.LEVEL - tau.flush)
        # check it:
        if (!all(is.na(tau.avging)) & any(tau.avging < 5) & warns) {
          cat(red(paste('ERROR: The resulting level averaging time of', tau.avging[tau.avging < 5], 's in line #', unlist(dimnames(lines.FR)), 'is too short!\n')))
          stop()
        }
        # if averaging time is long enough, extend the flushing time and limit the averaging to 15 seconds
        if(!all(is.na(tau.avging)) & any(tau.avging > 15)){
          tau.avging[tau.avging > 15] <- 15
        }
        
        # Profile of average [c] for each half-hour ............................
        GAS.c <- vector(mode = "list", length = N.LEV) # list with the last [c] (N = tau.avging) of each level (this is used for flux computation)
        
        # Level crossing the current half-hour
        # NOTE: the case where at the current half-hour there are not any level sampled, it is not currently handled
        crossing.level <- LEVEL[cur.hh.indx]
        
        if (any(cur.hh.indx) & !is.na(crossing.level)) {
          
          if(plot.conc.data){ # visualize the current [CO2] frame
            # windows()
            jpeg(paste0(output.dir, site.ID, '_', gas, '_conc_used.jpeg'), quality = 100, res = 300, width=480*5, height=480*5)
            plot.win <- c((cur.hh.indx-200) : (cur.hh.indx+200)) # to fit to the current half-hour (200 points windows)
            plot(PRO.time[plot.win], cur.scalar[plot.win], 
                 # col=prof.col[LEVEL[plot.win]-min(LEVEL, na.rm = TRUE)+1], pch=19, cex=1.5, xlab='Time [s]')
                 col=prof.col[LEVEL[plot.win]], pch=19, cex=1.5, xlab='Time [s]')
            points(PRO.time[hh.indxs], cur.scalar[hh.indxs])
            abline(v=SC.time[i.hh], lwd=9, col=adjustcolor('grey', 0.6))
            legend('topleft', legend=c(rev(ID.LEV), 'ref. h/2'), col=c(prof.col, adjustcolor('grey', 0.4)),
                   pch=c(rep(19,N.LEV), NA), pt.cex=c(rep(1.5,N.LEV), NA), lwd=c(rep(NA,N.LEV), 9), bty='n')
          }
          
          # Extract the level-specific concentrations (at 1 Hz) to average
          
          # system.time({
            i.lev <- 1
            for (i.lev in 1:N.LEV) {
              
              # 1. find the timestamp of the current level record which is closest to the current half-hour
              cur.lev.ts <- c(NULL)
              cur.lev.ts <- PRO.time[LEVEL == ID.LEV[i.lev]]
              level.ts <- numeric(1L)
              level.ts <- cur.lev.ts[birk::which.closest(cur.lev.ts, SC.time[i.hh])]
              # and respective index
              level.indx <- numeric(1L)
              level.indx <- which(PRO.time == level.ts)
              
              GAS.c[[i.lev]] <- c(NULL)
              cT0 <- cTX <- integer(1L)
              
              if (ID.LEV[i.lev] != crossing.level) { # If the current level IS NOT the crosser
                # 2.1. IF the TS is BEFORE the current half-hour, consider the PREVIOUS N seconds for averaging:
                if (level.ts <  SC.time[i.hh]) { 
                  cT0 <- level.indx - (tau.avging[i.lev]-1)
                  cTX <- level.indx
                } else { 
                # 2.2. IF the TS is AFTER the current half-hour, consider the LAST N seconds (of that level) for averaging:
                  cT0 <- level.indx + tau.LEVEL - tau.avging[i.lev]
                  cTX <- level.indx + (tau.LEVEL-1)
                }
              } else {  # If the current level IS the crosser
                # Find the switching point to the NEXT level
                sw.pts <- integer(1L)
                sw.pts <- which(diff(LEVEL[level.indx:(level.indx+tau.LEVEL*3/2)]) != 0) - 1 # -1 because of diff
                if (length(sw.pts) == 0) {sw.pts <- 0} # last half-hour
                if (length(sw.pts) > 1) { # this loop is to skip situations in which a level is sampled twice
                  level.indx.c <- level.indx + sw.pts[2] - tau.LEVEL
                } else {
                  level.indx.c <- level.indx + sw.pts
                }
                cT0 <- level.indx.c - (tau.avging[i.lev]-1)
                cTX <- level.indx.c
              }
              
              # level-specific concentrations 
              if (!any(is.na(c(cT0, cTX)))) {
                GAS.c[[i.lev]] <- cur.scalar[cT0 : cTX]
              } else {
                GAS.c[[i.lev]] <- rep(NA, median(tau.avging, na.rm=T))
              }
              
              if (plot.conc.data) {
                abline(v=level.ts, col=prof.col[ID.LEV[i.lev]-min(ID.LEV)+1], lwd=1, lty=2)
                abline(v=as.double(PRO.time[cTX]), lty=1, col=prof.col[ID.LEV[i.lev]-min(ID.LEV)+1], lwd=2)
                abline(v=as.double(PRO.time[cT0]), lty=1, col=prof.col[ID.LEV[i.lev]-min(ID.LEV)+1], lwd=2)
                rect(xleft = as.double(PRO.time[cT0]), ybottom = 0,
                     xright = as.double(PRO.time[cTX]), ytop = 800, col=adjustcolor(prof.col[ID.LEV[i.lev]-min(ID.LEV)+1], 0.3), border = NA)
              }
              
             # cat(paste('level', ID.LEV[i.lev], 'processed at', Sys.time(), '\n')) 
            }
          # }) # system.time close
            if (plot.conc.data) {dev.off()}
            
          ## Current average [c] profile, used for SC flux computation
          names(GAS.c) <- ID.LEV
          PRO.GAS.avg.ls[[i.hh]] <- c(NULL)
          PRO.GAS.avg.ls[[i.hh]] <- unlist(lapply(GAS.c, mean, na.rm = TRUE))
          
          # Current profile [c] standard deviation
          PRO.GAS.unc.ls[[i.hh]] <- c(NULL)
          if (unc == 'SD') {
            PRO.GAS.unc.ls[[i.hh]] <- unlist(lapply(GAS.c, sd, na.rm = TRUE))
          } else {
            PRO.GAS.unc.ls[[i.hh]] <- unlist(lapply(GAS.c, function(x) sd(x, na.rm = TRUE)/sqrt(length(x))))
          }
      
          
          if (plot.profile.data) {
            # windows()
            jpeg(paste0(output.dir, site.ID, '_', gas, '_profile.jpeg'), quality = 100, res = 300, width=480*5, height=480*5)
            plot(PRO.GAS.avg.ls[[i.hh]], SC.MD$'STO_PROF_HEIGHT'[1:N.LEV], 
                 xlim = mean(PRO.GAS.avg.ls[[i.hh]], na.rm = T) + c(-10,10),
                 type = 'o', pch=19, col = rev(prof.col),
                 xlab='Concentrations', ylab='Height (m)', las=1)
            polygon(c(c(PRO.GAS.avg.ls[[i.hh]]-PRO.GAS.unc.ls[[i.hh]]), rev(PRO.GAS.avg.ls[[i.hh]]+PRO.GAS.unc.ls[[i.hh]])),
                    c(SC.MD$'STO_PROF_HEIGHT'[1:N.LEV], SC.MD$'STO_PROF_HEIGHT'[N.LEV:1]),
                    col=adjustcolor(8, 0.3), border = NA)
            mtext(SC.time[i.hh], 3, 0)
            dev.off()
          }
          
          
        }else{ 
          
          PRO.GAS.avg.ls[[i.hh]] <- rep(NA, N.LEV); names(PRO.GAS.avg.ls[[1]]) <- ID.LEV
          PRO.GAS.unc.ls[[i.hh]] <- rep(NA, N.LEV); names(PRO.GAS.unc.ls[[1]]) <- ID.LEV
          
        }
        
      }else{
        
        PRO.GAS.avg.ls[[i.hh]] <- rep(NA, N.LEV); names(PRO.GAS.avg.ls[[1]]) <- ID.LEV
        PRO.GAS.unc.ls[[i.hh]] <- rep(NA, N.LEV); names(PRO.GAS.unc.ls[[1]]) <- ID.LEV
        
      }
      
    } 
    
    close(pb)
    
  # })# system.time closing
  
  # Assign the time reference
  names(PRO.GAS.avg.ls) <- as.character(SC.time)
  names(PRO.GAS.unc.ls) <- as.character(SC.time)
  
  
  ### Concentrations average profiles used for SC flux computation
  gas.pro.avg.SC <- do.call(rbind, PRO.GAS.avg.ls)
  
  ### Concentrations uncertanty profiles
  gas.pro.unc.SC <- do.call(rbind, PRO.GAS.unc.ls)
  
  ### Concentrations profiles uncertanty (standard deviations or standard error) summed in quadrature by time (in vertical) and then by two consecutive profile: 
  # ** this will used as the "within" component of the storage flux uncertainty **
  SC.UNC.wt.i <- apply(gas.pro.unc.SC, 1, function(x) sqrt(sum(x^2, na.rm = TRUE))) # profile-specific
  SC.UNC.wt <- data.table::frollapply(SC.UNC.wt.i, n = 2, function(x) sqrt(sum(x^2, na.rm = TRUE))) # quadrature for differentials
  # summary(SC.UNC.wt)
  
  # Store if and which level is missing whithin each profile
  missing.lev.N <- apply(gas.pro.avg.SC, 1, function(x) length(which(is.na(x))))
  missing.lev.ID <- apply(gas.pro.avg.SC, 1, function(x) paste(ID.LEV[which(is.na(x))], collapse = ';'))
  
  
  
  # KINEMATIC STORAGE TERM **************************************************
  
  # Differential concentration profiles
  D.GAS <- apply(gas.pro.avg.SC, 2, diff)
  
  # Kinematic SC flux
  SC.GAS.k <- apply(D.GAS, 1, function(x) sum((x/Dt) * Dz.i))
  
  
  # Time referencing the time series
  SC.GAS.k <- xts(SC.GAS.k, order.by = as.POSIXct(names(SC.GAS.k), '%Y-%m-%d %H:%M:%S', tz='GMT'), tzone = 'GMT') # covert to xts object
  SC.GAS.k <- merge(SC.GAS.k, xts(rep(1, length(SC.time)), order.by = SC.time, tzone = 'GMT'))[, -2] # align
  
  
  # STORAGE FLUX (mass-based quantity) ***************************************
  
  # CO2, CH4 or N2O storege flux
  if (any(grepl('CO2|CH4|N2O', gas))) {
    SC.GAS <- rho.d * SC.GAS.k  
  }
  
  # Latent heat storage flux [Wm-2]
  if(gas == 'H2O'){
    SC.GAS <- rho.d * seh * Mw *10^-3 * SC.GAS.k
  }
  
  names(SC.GAS) <- paste0('SC_', gas)
  
  # STORAGE FLUX COMBINED UNCERTAINTY ****************************************
  SC.GAS.UNC <- sqrt(SC.UNC.wt^2 + SC.UNC.bt^2)
  names(SC.GAS.UNC) <- 'SC_CombUNC'
  # summary(SC.GAS.UNC)
  
  # Relative Flux Frror
  RFE <- SC.GAS.UNC/abs(SC.GAS)
  # summary(RFE)
  # windows()
  # hist(na.omit(RFE))
  # plot(density(na.omit(RFE)), xlim=c(0,5))
  
  
  
  # Storage flux daily cycle ************************************************* 
  
  SC.GAS.dc <- tapply(SC.GAS, format(SC.time, '%H:%M'), mean, na.rm=T)
  # Standard error of the mean
  SC.GAS.sem <- tapply(SC.GAS, format(SC.time, '%H:%M'), function(x) sd(x, na.rm = T)/sqrt(length(x))) 
  # 95% confidence intervals of the mean
  SC.GAS.ci <- data.frame(ll = SC.GAS.dc - 2 * SC.GAS.sem, ul = SC.GAS.dc + 2 * SC.GAS.sem)
  SC.GAS.dc.ts <- strptime(names(SC.GAS.dc), '%H:%M', tz='GMT')
  
  
  cat(cyan('Storage flux successfully computed.\n'))
  
  # Possibly plot the SC time series
  if (plot.SC.flux) {
    
    if(!is.null(output.dir)){
      jpeg(paste0(output.dir, site.ID, '_SC-', ifelse(gas == 'H2O', 'LE', gas), '.jpeg'), quality = 100, res = 300, width=480*5, height=480*5)
    } else {
      windows()
    }
    
    nf <- layout(matrix(c(0,1,0,0,2,0), ncol=2), widths = c(0.7,0.3))
    # layout.show(nf)
    par(mar=c(.2,.1,.1,.2), oma=c(4.5,4,4.5,0))
    sc.ylim <- range(SC.GAS, na.rm = T) * 1.5
    # Storage flux time series
    plot(SC.time, SC.GAS, type='l', xlab = 'Time', las=1, ylim = sc.ylim); abline(h=0, col=8)
    mtext('Time', 1, 2.5)
    if(gas == 'CO2'){mtext(bquote(CO[2] ~ 'SC Flux [' ~ mu ~ 'mol' ~ m^-2 ~ s^-1 ~ ']'), 2, 2.5); shdw.col <- 1}
    if(gas == 'H2O'){mtext(bquote('Latent heat SC Flux [W' ~ m^-2 ~ ']'), 2, 2.5); shdw.col <- 4}
    if(gas == 'CH4'){mtext(bquote(CH[4] ~ 'SC Flux [nmol' ~ m^-2 ~ s^-1 ~ ']'), 2, 2.5); shdw.col <- 2}
    if(gas == 'N2O'){mtext(bquote(N[2]*'O' ~ 'SC Flux [nmol' ~ m^-2 ~ s^-1 ~ ']'), 2, 2.5); shdw.col <- 3}
    if (unc == 'SD'){
      polygon(c(SC.time, rev(SC.time)),
              replace(c(c(coredata(SC.GAS-SC.GAS.UNC)), rev(c(coredata(SC.GAS+SC.GAS.UNC)))), is.na(c(c(coredata(SC.GAS-SC.GAS.UNC)), rev(c(coredata(SC.GAS+SC.GAS.UNC))))), 0),
              col=adjustcolor(shdw.col, 0.2), border=NA)
    } else {
      polygon(c(SC.time, rev(SC.time)),
              replace(c(c(coredata(SC.GAS-2*SC.GAS.UNC)), rev(c(coredata(SC.GAS+2*SC.GAS.UNC)))), is.na(c(c(coredata(SC.GAS-SC.GAS.UNC)), rev(c(coredata(SC.GAS+SC.GAS.UNC))))), 0),
              col=adjustcolor(shdw.col, 0.2), border=NA)
    }
    
    points(SC.time, SC.GAS, type='l')
    mtext(paste(site.ID, ' ', paste(range(index(SC.GAS)), collapse = ' - ')), 3, 0)
    
    # Storage flux daily cycle
    plot(SC.GAS.dc.ts, SC.GAS.dc, ylim = sc.ylim, type='l', xaxt ='n', yaxt ='n'); abline(h=0, col=8)
    polygon(c(SC.GAS.dc.ts, rev(SC.GAS.dc.ts)),
            c(SC.GAS.ci$ll, rev(SC.GAS.ci$ul)),
            col=adjustcolor(shdw.col, 0.2), border=NA)
    points(SC.GAS.dc.ts, SC.GAS.dc, type='l')
    axis.POSIXct(1, at=SC.GAS.dc.ts[c(T,F,F)], labels=format(SC.GAS.dc.ts,"%H:%M")[c(T,F,F)], las=1)
    mtext('Time [HH:MM]', 1, 2.5)
    
    if(!is.null(output.dir)){dev.off()}
    
  }
  
  
  if (write.output) { # write results on file
    # SC flux
    SC.df.out <- data.frame(TIMESTAMP_END = format(index(SC.GAS), '%Y%m%d%H%M%S'), 
                            SC = c(unname(coredata(SC.GAS))),
                            SC_UNC = c(unname(coredata(SC.GAS.UNC))),
                            SC_UNC_b = c(unname(coredata(SC.UNC.bt))),
                            SC_UNC_w = SC.UNC.wt,
                            RFE = c(unname(coredata(RFE))),
                            missing_lev_N = unname(missing.lev.N),
                            missing_lev_ID = unname(missing.lev.ID))
    fwrite(SC.df.out, paste0(output.dir, site.ID, '_SC-', ifelse(gas == 'H2O', 'LE', gas), '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    
    # SC flux daily cycle
    fwrite(data.frame(TIMESTAMP_START = format(SC.GAS.dc.ts,"%H:%M"), SC = unname(SC.GAS.dc), SC_cill = SC.GAS.ci$ll, SC_ciul = SC.GAS.ci$ul), 
           paste0(output.dir, site.ID, '_SCdc-', ifelse(gas == 'H2O', 'LE', gas), '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    
    # Concentration profile
    fwrite(gas.pro.avg, 
           paste0(output.dir, site.ID, '_AvgProfile-', ifelse(gas == 'H2O', 'LE', gas), '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    
  } 
  
  
  # output list
  SC.list <- list(SC.GAS, SC.GAS.dc, gas.pro.avg)
  names(SC.list) <- c(paste0('SC_', ifelse(gas == 'H2O', 'LE', gas)), paste0('SC_', ifelse(gas == 'H2O', 'LE', gas), '_dc'), paste0('PRO_', ifelse(gas == 'H2O', 'LE', gas), '_AvgC'))
  return(SC.list)
  
  
  
  
}

