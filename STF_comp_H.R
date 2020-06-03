

#'' ***************************************************************************
# STORAGE FLUX COMPUTATION
# + SEQUENTIAL SAMPLING SCHEME  +
# 
# author: Giacomo Nicolini
# contact: g.nicolini@unitus.it
# date: 2020-04-10
# For license information see LICENSE file
# 
#'' ***************************************************************************

stfcomp_H <- function (
  site.ID,
  BM.data.path,
  BM.data.AVG.path,
  dfn.pattern = 'profile', # data file name pattern to match
  BM.MD,
  TA.tau = 300, # Air temperautre data averaging time [seconds]. Default is 150 seconds before and 150 seconds after the halfhour
  EC.height,
  # RH.profile,
  # plot.profile.data = FALSE,
  # plot.conc.data = FALSE,
  plot.SCH.flux = TRUE,
  warns = FALSE,
  write.output = FALSE,
  output.dir = NULL
){
  
  # Install required packages if they are not already installed ..............
  if (!require("data.table")) {install.packages("data.table"); library(data.table)}
  if (!require("lubridate")) {install.packages("lubridate"); library(lubridate)}
  if (!require("bigleaf")) {install.packages("bigleaf"); library(bigleaf)}
  if (!require("xts")) {install.packages("xts"); library(xts)}
  
  options(scipen=9999)
  
  cat(cyan('Calculating the latent heat storage flux ...\n'))
  
  
  # Meteo profile characteristics ********************************************
  # From BADM system
  # NOTE: the variables ID are renamed according the level's respective heights
  cat(cyan('Get and manipulate meteo profile metadata ...\n'))
  
  ## Air temperature .........................................................
  BM.MD.TA <- BM.MD[grep('TA', BM.MD$'BM_VARIABLE_H_V_R'), ]
  BM.MD.TA <- BM.MD.TA[order(BM.MD.TA$BM_HEIGHT), ]
  # Scale TA variables vertical index according the the effective heights and rename
  TA.new.V.indx <- order(BM.MD$'BM_HEIGHT'[grep('TA', BM.MD$'BM_VARIABLE_H_V_R')], decreasing = TRUE)
  BM.MD.TA$'BM_VARIABLE_H_V_R' <- paste0('TA_', TA.new.V.indx)
  # Limit the profile to the EC system height (for meteo is taken the closest level)
  BM.MD.TA <- BM.MD.TA[1: which.min(abs(BM.MD.TA$'BM_HEIGHT' - EC.height)), ]
  # TA.indx <- grep('TA', BM.MD.TA$'BM_VARIABLE_H_V_R')
  # Number of sampling levels
  # N.LEV.TA <- length(BM.MD$'BM_HEIGHT'[TA.indx])
  N.LEV.TA <- length(BM.MD.TA$'BM_HEIGHT')
  # Heights of sampling levels
  # z.i.TA <- BM.MD$'BM_HEIGHT'[TA.indx]
  z.i.TA <- BM.MD.TA$'BM_HEIGHT'
  # names(z.i.TA) <- paste0('TA_', N.LEV.TA:1) # renamed as in the profile data
  names(z.i.TA) <- BM.MD.TA$'BM_VARIABLE_H_V_R' # renamed as in the profile data
  # Heights of the level-specifc boundary layers
  zl.i.TA <- c(0, c(z.i.TA[-N.LEV.TA] + z.i.TA[-1])/2, max(z.i.TA))
  # Boundary air layer depth
  Dz.i.TA <- diff(zl.i.TA)
  # names(Dz.i.TA) <- paste0('TA_', N.LEV.TA:1) # renamed as in the profile data
  names(Dz.i.TA) <- BM.MD.TA$'BM_VARIABLE_H_V_R'  # renamed as in the profile data
  
  ## Air relative humidity
  BM.MD.RH <- BM.MD[grep('RH', BM.MD$'BM_VARIABLE_H_V_R'), ]
  BM.MD.RH <- BM.MD.RH[order(BM.MD.RH$BM_HEIGHT), ]
  # Scale RH variable vertical index according the the effective heights
  RH.new.V.indx <- order(BM.MD$'BM_HEIGHT'[grep('RH', BM.MD$'BM_VARIABLE_H_V_R')], decreasing = TRUE)
  BM.MD.RH$'BM_VARIABLE_H_V_R' <- paste0('RH_', RH.new.V.indx)
  # Limit the profile to the EC system height
  BM.MD.RH <- BM.MD.RH[1: which.min(abs(BM.MD.RH$'BM_HEIGHT' - EC.height)), ]
  # RH.indx <- grep('RH', BM.MD.RH$'BM_VARIABLE_H_V_R')
  # Number of sampling levels
  # N.LEV.RH <- length(BM.MD$'BM_HEIGHT'[RH.indx])
  N.LEV.RH <- length(BM.MD.RH$'BM_HEIGHT')
  # Heights of sampling levels
  # z.i.RH <- BM.MD$'BM_HEIGHT'[RH.indx]
  z.i.RH <- BM.MD.RH$'BM_HEIGHT'
  # names(z.i.RH) <- paste0('RH_', N.LEV.RH:1) # renamed as in the profile data
  names(z.i.RH) <- BM.MD.RH$'BM_VARIABLE_H_V_R' # renamed as in the profile data
  # Heights of the level-specifc boundary layers
  zl.i.RH <- c(0, c(z.i.RH[-N.LEV.RH] + z.i.RH[-1])/2, max(z.i.RH))
  # Boundary air layer depth
  Dz.i.RH <- diff(zl.i.RH)
  # names(Dz.i.RH) <- paste0('RH_', N.LEV.RH:1) # renamed as in the profile data
  names(Dz.i.RH) <- BM.MD.RH$'BM_VARIABLE_H_V_R'  # renamed as in the profile data
  
  ## Air pressure
  BM.MD.PA <- BM.MD[grep('PA', BM.MD$'BM_VARIABLE_H_V_R'), ]
  BM.MD.PA <- BM.MD.PA[order(BM.MD.PA$BM_HEIGHT), ]
  # Scale PA variable vertical index according the the effective heights
  PA.new.V.indx <- order(BM.MD$'BM_HEIGHT'[grep('PA', BM.MD$'BM_VARIABLE_H_V_R')], decreasing = TRUE)
  BM.MD.PA$'BM_VARIABLE_H_V_R' <- paste0('PA_', PA.new.V.indx)
  # Limit the profile to the EC system height
  BM.MD.PA <- BM.MD.PA[1: which.min(abs(BM.MD.PA$'BM_HEIGHT' - EC.height)), ]
  # PA.indx <- grep('PA', BM.MD.PA$'BM_VARIABLE_H_V_R')
  # Number of sampling levels
  # N.LEV.PA <- length(BM.MD$'BM_HEIGHT'[PA.indx])
  N.LEV.PA <- length(BM.MD.PA$'BM_HEIGHT')
  # Heights of sampling levels
  # z.i.PA <- BM.MD$'BM_HEIGHT'[PA.indx]
  z.i.PA <- BM.MD.PA$'BM_HEIGHT'
  # names(z.i.PA) <- paste0('PA_', N.LEV.PA:1) # renamed as in the profile data
  names(z.i.PA) <- BM.MD.PA$'BM_VARIABLE_H_V_R' # renamed as in the profile data
  # Heights of the level-specifc boundary layers
  zl.i.PA <- c(0, c(z.i.PA[-N.LEV.PA] + z.i.PA[-1])/2, max(z.i.PA))
  # Boundary air layer depth
  Dz.i.PA <- diff(zl.i.PA)
  # names(Dz.i.PA) <- paste0('PA_', N.LEV.PA:1) # renamed as in the profile data
  names(Dz.i.PA) <- BM.MD.PA$'BM_VARIABLE_H_V_R'  # renamed as in the profile data
  

  # **************************************************************************
  # Meteo profile data import and manipulation
  
  BMSC.files <- list.files(BM.data.path, pattern = dfn.pattern, full.names = TRUE)
  if (length(BMSC.files) == 0) {
    cat(red('ERROR: Meteorological data (temperature rofile, air relative humididty and pressure) data must be provided, covering the processing period.\n'))
    cat(yellow('The files must be daily, comma separeted ASCII with - at least - the following variables/columns:\n'))
    cat(silver('TIMESTAMP_START,TIMESTAMP_END,TA_1,...,TA_n,RH_1,...,RH_n,PA_1,...,PAH_n\n'))
    cat(silver('Please see the storage data instruction document for details.\n'))
    stop()
  }
  
  # Create a unique dataframe for the whole data period (about 1.5-2 GB prof-1 y-1)
  BM.data.or <- data.frame(NULL)
  BM.data.or <- as.data.frame(data.table::rbindlist(lapply(BMSC.files, function(x) fread(x, na.strings = '-9999', header = TRUE, data.table = FALSE)), fill = TRUE))
  
  # Drop the variables (columns) sampled above the EC system
  # Timestamps
  BM.data <- BM.data.or[, 1:2]
  # TA 
  BM.data <- cbind(BM.data, BM.data.or[, names(BM.data.or) %in% names(Dz.i.TA)])
  # RH
  if (any(names(BM.data.or) %in% names(Dz.i.RH))){ # profile
    BM.data <- cbind(BM.data, BM.data.or[, names(BM.data.or) %in% names(Dz.i.RH)])
  } else { # single sample
    BM.data <- cbind(BM.data, data.frame('RH' = BM.data.or[, grep('RH', names(BM.data.or))]))
  }
  # PA
  if (any(names(BM.data.or) %in% names(Dz.i.PA))){ # profile
    BM.data <- cbind(BM.data, BM.data.or[, names(BM.data.or) %in% names(Dz.i.PA)])
  } else { # single sample
    BM.data <- cbind(BM.data, data.frame('PA' = BM.data.or[, grep('PA', names(BM.data.or))]))
  }
  rm(BM.data.or)
  
  
  ## Air temperature profile [°C] ............................................
  TA <- BM.data[, grep('TA_', names(BM.data), fixed = TRUE)]
  # ordering TA bottom-up (from N = ground to 1 = tower-top)
  TA <- TA[, order(as.numeric(substr(names(TA), 4, 6)), decreasing = TRUE)] 
  
  # check wheter is measured along a profile or not. * NOTE: it must be a profile so in case negative, stop processing *
  if (!ifelse(ncol(TA) == 1, FALSE, TRUE) & warns) {cat(red('ERROR: air temperature must be measured along a profile.\n')); stop()}
  
  ## Time specification
  # Profile data original timestamp
  BM.time <- strptime(BM.data$'TIMESTAMP_END', '%Y%m%d%H%M', tz='GMT') # time object in UTC (POSIXlt)
  # native data frequency [s]
  BM.data.f <- as.numeric(names(which.max(table(diff(as.numeric(BM.time)))))) 
  
  ## Check and possibly correct for time-series jumps
  if (any(diff(as.numeric(BM.time)) != BM.data.f)) {
    
    # Profile data desired timestamp (pure 1 Hz frequency)
    BM.time.cont <- data.frame(TIMESTAMP_tmp = seq(first(BM.time), last(BM.time), by = BM.data.f))
    
    # Alining
    BM.data$'TIMESTAMP_tmp' <- as.POSIXct(BM.time)
    BM.data <- merge(BM.data, 
                      BM.time.cont, 
                      by = 'TIMESTAMP_tmp', 
                      all = TRUE)
    
    # Reset the profile data original timestamp
    BM.time <- c(NULL)
    BM.time <- BM.data$'TIMESTAMP_tmp'
    BM.data$'TIMESTAMP_END' <- format(BM.time, '%Y%m%d%H%M%S')
    BM.data$'TIMESTAMP_tmp' <- NULL
    
    # Find and drop possible duplicated timestamps values (WARNING: this may remove good - although badly acquired - data)
    if (any(duplicated(BM.data$'TIMESTAMP_END'))){
      BM.data <- BM.data[!duplicated(BM.data$'TIMESTAMP_END'), ]
    }
    
  }
  
  # Refer the timestamps to the nearest half-hour
  BM.time.hh <- as.POSIXct(lubridate::round_date(strptime(BM.data$'TIMESTAMP_END', '%Y%m%d%H%M', tz='GMT'), "30 minutes"))
  # BM.time.hh <- as.POSIXct(lubridate::ceiling_date(strptime(BM.data$'TIMESTAMP_END', '%Y%m%d%H%M', tz='GMT'), "30 minutes"))
  
  # Benchmark timestamp for the storage flux (30 minuts resolution)
  SCH.time <- seq(xts::first(BM.time.hh)+1800, xts::last(BM.time.hh), 1800)
  
  
  # Useful constants
  Rd <- bigleaf.constants()$Rd      # Gas constant for dry air (Rgas/Md) [J kg-1 K-1]
  Rgas <- bigleaf.constants()$Rgas  # Universal gas constant [J K-1 mol-1 or m3 Pa K-1 mol-1]
  Rv <- bigleaf.constants()$Rv      # Gas constant for water vapour [J kg-1 K-1]
  Md <- bigleaf.constants()$Md      # Molar mass of dry air [kg mol-1] (28.9647 g mol-1)
  Mw <- bigleaf.constants()$Mw      # Molar mass of H2O [kg mol-1] (18.01528 g mol-1 )
  if(!exists('Dt')){Dt <- 1800}
  
  
  # **************************************************************************
  ## Compute auxiliary micrometeorological variables 
  # Based profile-wheighted half-hourly averages as computed according to the "process_BM" function
  
  # Load profile wheithed half-hourly averages
  BM.data.avg <- fread(paste0(BM.data.AVG.path, site.ID, '_BM_SC_AVG.csv', sep=''), sep = ',', header = TRUE, fill = TRUE, blank.lines.skip = TRUE, data.table = FALSE)
  
  # Possibly align H2O MF and master (SCH) timestamps and in case of missing values, fill with NA
  BM.data.avg.ts <- strptime(BM.data.avg$'TIMESTAMP', '%Y%m%d%H%M', tz='GMT')
  BM.data.avg$'TIMESTAMP' <- as.double(BM.data.avg$'TIMESTAMP')
  BM.data.avg.xts <- xts(BM.data.avg, order.by = BM.data.avg.ts) # covert to xts object
  BM.data.avg <- merge(BM.data.avg.xts, xts(rep(1, length(SCH.time)), order.by = SCH.time, tzone = 'GMT'))[, -(ncol(BM.data.avg.xts)+1)] # align
  # BM.data.avg <- as.data.frame(coredata(BM.data.avg.xts)) # and get back to vector
  rm(BM.data.avg.xts)
  
  
  # Dry air density ..........................................................
  # the average dry air density for the air column is used as scaling variable
  # to convert the kinematic storage term into a mass-based quantity
  # NOTE: the time-averaging is performed at half-hourly scale
  
  # Air molar volume [m3 mol-1] (should be around 0.025 m3 mol-1 at 25 °C)
  amv <- Rgas * (BM.data.avg$'TA' + 273.15) / (BM.data.avg$'PA'*10^3) # PA from KPa to Pa
  
  # Water vapor mass density [kg m-3] (e.g. range: 0.005 kgm-3 at 0 °C  - 0.051 kg m-3 at 40 °C)
  # NOTE: if RH is measured along the profile use the resulting H2O molar fraction as calculated by ,
  #       otherwise try to compute it with the data provided here
  if ('H2O_MF_meteo' %in% names(BM.data.avg)) {
    rho.H2O <- c(BM.data.avg$'H2O_MF_meteo'*10^-3) * Mw / amv  # [H2O] from mmol to mol
  } else {
    # Average air temperature profile - weighted by the depths of the air layers representative of the individual measurement points
    TA.avg <- apply(TA, 1, function(x) stats::weighted.mean(x, w = Dz.i.TA, na.rm = TRUE))
    TA.avg.h <- tapply(TA.avg, BM.time.hh, mean, na.rm=T)
    # Average relative humidity along the profile - weighted by the depths of the air layers representative of the individual measurement points
    RH.avg <- apply(BM.data[, grep('RH_', names(BM.data), fixed = TRUE)], 1, function(x) stats::weighted.mean(x, w = Dz.i.RH, na.rm = TRUE))
    RH.avg.h <- tapply(RH.avg, BM.time.hh, mean, na.rm=T)
    # Average relative humidity along the profile - weighted by the depths of the air layers representative of the individual measurement points
    if (length(grep('PA', names(BM.data), fixed = TRUE))> 1) {
      PA.avg <- apply(BM.data[, grep('PA', names(BM.data), fixed = TRUE)], 1, function(x) stats::weighted.mean(x, w = Dz.i.PA, na.rm = TRUE))
    } else {
      PA.avg <- BM.data[, grep('PA', names(BM.data), fixed = TRUE)]
    }
    PA.avg.h <- tapply(PA.avg, BM.time.hh, mean, na.rm=T)
    # H2O molar fraction - resulting from hygrometers data 
    H2O.MF.avg.meteo <- (6088.484 * exp(TA.avg.h/(TA.avg.h+237.6429)*17.31303) / 18.0153) * RH.avg.h / (PA.avg.h/Rgas*(TA.avg.h+273.15))
    # Water vapor mass density
    rho.H2O <- c(H2O.MF.avg.meteo*10^-3) * Mw / amv  # [H2O] from mmol to mol
  }
  
  # Water vapor partial pressure [kPa] (should be around 2.338 kPa at 20 °C)
  e.avg <- rho.H2O * Rv*10^-3 * (BM.data.avg$'TA' + 273.15) # Rv from J kg-1 K-1 (i.e. Pa m3 kg-1 K-1) to kPa m3 kg-1 K-1
  
  # Dry air partial pressure [kPa]
  PA.avg.d <- BM.data.avg$'PA' - e.avg
  
  ## Average *moist* air *mass* density [Kg m-3]
  rho.a <- PA.avg.d*10^3 / (Rd * (BM.data.avg$'TA' + 273.15)) 
  names(rho.a) <- 'rho'
  
  
  # Dry air heat capacity at constant pressure (cp) [J kg-1 K-1] .............
  cp <- 1005 + ((BM.data.avg$'TA' + 23.12)^2) / 3364
  names(cp) <- 'cp'
  
  cat(cyan('Done.\n'))
  
  
  # ***************************************************************************
  # STORAGE FLUX COMPUTATION
  # NOTE: The averaging is based on user defined seconds of TA values, crossing each half-hour
  cat(cyan('Computing storage flux ...\n'))
  
  # Time reference example:
  # Assuming that turbulent flux and storage flux estimates are to be computed for timestamps
  # 08:00:00, then by convention the 08:00:00 timestamp represents the flux corresponding to
  # observations between 07:30:00 and 08:00:00. DC will then be computed from the time 
  # average of measurements from 07:59:55 to 08:00:05 minus the time average
  # of measurements from 07:29:55 to 07:30:05
  
  # convert TA to K
  TA <- TA + 273.15
  
  # Find the profile specific timpestamps indexes matching the half-hours
  TA.win.s <- which(as.double(BM.time) %in% as.double(SCH.time)) # steps
  
  # Set the averaging time (HALF) windows according to the original data resolution (as number of rows) and the specified averaging time
  TA.win.w <- length(which(BM.time >= SCH.time[1] - lubridate::seconds(trunc(TA.tau/2)) & BM.time <= SCH.time[1] + lubridate::seconds(trunc(TA.tau/2)))) # window width
  
  # Create an index list over which calculating the averages
  avg.indx.df <- data.frame(start = TA.win.s - TA.win.w, stop = TA.win.s + TA.win.w)
  avg.indx.ls <- lapply(seq(nrow(avg.indx.df)), function(x) avg.indx.df[x,1]:avg.indx.df[x,2])
  
  # Averaging of [c] at each half-hour (average profiles) - data by row
  TA.i.avg <- simplify2array(
    do.call(cbind, apply(TA, 2, function(xx) lapply(avg.indx.ls, function(x) mean(xx[x], na.rm = TRUE)))
    )
  )
  TA.i.avg <- apply(TA.i.avg, 2, function(x) unlist(x))
  row.names(TA.i.avg) <- as.character(SCH.time)
  
  # Calculate DTA for each level (differential profiles)
  D.TA.i <- data.frame(NULL)
  D.TA.i <- as.data.frame(apply(TA.i.avg, 2, function(x) diff(x)))
  D.TA.i[D.TA.i == "NaN"] <- NA
  
  # interpolate possible NAs between levels
  if(any(is.na(D.TA.i))){
    D.TA.i <- t(apply(D.TA.i, 1, function(x) data.table::nafill(x, 'nocb')))
  }
  
  # Align data with the flux timestamp (add the first row as NA removed by diff)
  D.TA.i.xts <- xts(D.TA.i, order.by = as.POSIXct(rownames(D.TA.i), '%Y-%m-%d %H:%M:%S', tz='GMT'), tzone = 'GMT') # covert to xts object
  D.TA.i.xts <- merge(D.TA.i.xts, xts(rep(1, length(SCH.time)), order.by = SCH.time, tzone = 'GMT'))[, -c(ncol(D.TA.i.xts)+1)] # align
  

  # KINEMATIC STORAGE TERM ***************************************************
  
  SCH.k <- apply((D.TA.i.xts), 1, function(x) sum((x/Dt) * Dz.i.TA))
  
  
  # STORAGE TERM (mass-based quantity) ***************************************
  
  SCH <- SCH.k * rho.a * cp

  cat(cyan('Sensible heat storage flux succesfully computed.\n'))
  
  
  # Storage flux daily cycle ************************************************* 
  
  SCH.dc <- tapply(SCH, format(SCH.time, '%H:%M'), mean, na.rm=T)
  # Standard error of the mean
  SCH.sem <- tapply(SCH, format(SCH.time, '%H:%M'), function(x) sd(x, na.rm = T)/sqrt(length(x))) 
  # 95% confidence intervals of the mean
  SCH.ci <- data.frame(ll = SCH.dc - 2 * SCH.sem, ul = SCH.dc + 2 * SCH.sem)
  SCH.dc.ts <- strptime(names(SCH.dc), '%H:%M', tz='GMT')
  
  
  
  # Possibly plot the SC time series
  if (plot.SCH.flux) {
    
    if(!is.null(output.dir)){
      jpeg(paste0(output.dir, site.ID, '_SC-H.jpeg'), quality = 100, res = 300, width=480*5, height=480*5)
    } else {
      windows()
    }
    
    nf <- layout(matrix(c(0,1,0,0,2,0), ncol=2), widths = c(0.7,0.3))
    # layout.show(nf)
    par(mar=c(.2,.1,.1,.2), oma=c(4.5,4,4.5,0))
    sc.ylim <- range(SCH, na.rm = T)
    # Storage flux time series
    plot(SCH.time, SCH, type='l', xlab = 'Time', las=1, ylim = sc.ylim); abline(h=0, col=8)
    mtext('Time', 1, 2.5)
    mtext(bquote('Sensible heat SC Flux [W' ~ m^-2 ~ ']'), 2, 2.5)
    mtext(paste(site.ID, ' ', paste(range(index(SCH)), collapse = ' - ')), 3, 0)  
    # Storage flux daily cycle
    plot(SCH.dc.ts, SCH.dc, ylim = sc.ylim, type='l', xaxt ='n', yaxt ='n'); abline(h=0, col=8)
    polygon(c(SCH.dc.ts, rev(SCH.dc.ts)), 
            c(SCH.ci$ll, rev(SCH.ci$ul)), 
            col=adjustcolor('orange', 0.5), border=NA)
    points(SCH.dc.ts, SCH.dc, type='l')
    axis.POSIXct(1, at=SCH.dc.ts[c(T,F,F)], labels=format(SCH.dc.ts,"%H:%M")[c(T,F,F)], las=1)
    mtext('Time [HH:MM]', 1, 2.5)
    if(!is.null(output.dir)){dev.off()}
    
  }
  
  cat(cyan('Done.\n'))
  
  
  if (write.output) { # write results on file
    # SC flux
    fwrite(
	data.frame(TIMESTAMP_END = format(index(SCH), '%Y%m%d%H%M%S'), SC = coredata(unname(SCH))), 
	paste0(output.dir, site.ID, '_SC-H', '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    # SC flux daily cycle
    fwrite(data.frame(TIMESTAMP_START = format(SCH.dc.ts,"%H:%M"), SC = unname(SCH.dc), SC_cill = SCH.ci$ll, SC_ciul = SCH.ci$ul), 
           paste0(output.dir, site.ID, '_SCdc-H', '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    # # Concentration profile
    # fwrite(GAS.avg.pro.c, 
    #        paste0(output.dir, site.ID, '_AvgProfile-', gas, '.csv'), sep=',', row.names = F, quote = F, na = '-9999')
    
  } 
  
  # output list
  SCH.list <- list(
    SCH = SCH, 
    SCH_dc = SCH.dc)
  return(SCH.list)
  
  
  
  
}
