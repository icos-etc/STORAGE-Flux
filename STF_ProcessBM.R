


# ******************************************************************************

proc_profilebm <- function(
  site.ID,
  BM.data.path,
  BM.MD,
  EC.height,
  Dt = 1800
){
  
  # install required packages if they are not already installed
  if (!require("bigleaf")) {install.packages("bigleaf"); library(bigleaf)}
  if (!require("data.table")) {install.packages("data.table"); library(data.table)}
  
  
  # ******************************************************************************
  # METEO PROCESSING =============================================================
  # Compute half-hourly aggregated values used for SC computation
  # For profile data, profile averages are weighted by the depths of the air layers 
  # representative of the individual measurement points
  # ******************************************************************************
  cat(cyan('Meteo data processing ...\n'))
  
  # useful constants
  Rd <- bigleaf.constants()$Rd      # Gas constant for dry air (Rgas/Md) [J kg-1 K-1]
  g <- bigleaf.constants()$g        # Gravitational acceleration [m s-2]
  Rgas <- bigleaf.constants()$Rgas  # Universal gas constant [J K-1 mol-1 or m3 Pa K-1 mol-1]
  
  # Check and list the meteo file
  BM.files <- list.files(BM.data.path, pattern = '_profile', full.names = TRUE)
  
  if (length(BM.files) == 0) {
    cat(red('ERROR: 1 minute agregated meteo data must be provided, covering the processing period.\n'))
    cat(yellow('The files must be daily, comma separeted ASCII with the following variables/columns:\n'))
    cat(silver('TIMESTAMP_START, TIMESTAMP_END, PA_H_V_R, RH_H_V_R, TA_H_V_R\n'))
    cat(silver('where "_H_V_R" are the position indexes.\n'))
  }
  
  
  
  # Meteo profile characteristics ********************************************
  # From BADM system
  # NOTE: the variables ID are possibly renamed according the level's respective heights
  
  
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
  
  
  
  # Meteo data ***************************************************************
  
  # Create a unique dataframe for the whole data period (about 30-45 MB y-1)
  BM.data.or <- data.frame(NULL)
  BM.data.or <- as.data.frame(data.table::rbindlist(lapply(BM.files, function(x) 
    fread(x, na.strings = '-9999', header = TRUE, data.table = F)), fill = TRUE))
  
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
  
  
  ### Time specs
  BM.time.nat <- strptime(BM.data$'TIMESTAMP_END', '%Y%m%d%H%M', tz='GMT') # native timestamp
  # native data frequency [s]
  BM.data.f <- as.numeric(names(which.max(table(diff(as.numeric(BM.time.nat)))))) # it must be 60 seconds
  
  # Check and possibly correct for time jumps
  if (any(as.numeric(names(table(diff(as.numeric(BM.time.nat))))) != BM.data.f)) {
    
    # creaate a temporary timeline
    BM.data$'TIMESTAMP_tmp' <- as.POSIXct(BM.time.nat)
    # align
    BM.data <- merge(BM.data, 
                     data.frame(TIMESTAMP_tmp = seq(first(BM.time.nat), last(BM.time.nat), by = BM.data.f)), 
                     by = 'TIMESTAMP_tmp', all = TRUE)
    # adjust the meteo reference times
    BM.time.nat <- BM.data$'TIMESTAMP_tmp'
    BM.data$'TIMESTAMP_END' <- format(BM.time.nat, '%Y%m%d%H%M')
    BM.data$'TIMESTAMP_START' <- format(BM.time.nat-BM.data.f, '%Y%m%d%H%M')
    BM.data$'TIMESTAMP_tmp' <- NULL
    
  }
  
  # Timestam rounding to 30 minutes
  BM.time.nat_30 <- lubridate::ceiling_date(BM.time.nat, paste(Dt/60, 'minutes'))
  
  
  # Air temperature profile [°C]
  TA <- BM.data[, grep('TA_', names(BM.data), fixed = TRUE)]
  # ordering TA bottom-up (from N = ground to 1 = tower-top)
  TA <- TA[, order(as.numeric(substr(names(TA), 4, 6)), decreasing = TRUE)] 
  # check wheter is measured along a profile or not. NOTE: it must be a profile
  TA.on.profile <- ifelse(ncol(TA) == 1, FALSE, TRUE)
  if (!TA.on.profile) {
    cat(red('Error: air temperature must be measured along a profile.\n')); stop()
    }
  
  # Air relative humidity (can be a profile) [%]
  RH <- BM.data[, grep('RH', names(BM.data), fixed = TRUE)]
  # if measured along a profile
  if(is.data.frame(RH)) {
    # ordering TA bottom-up (from N = ground to 1 = tower-top)
    RH <- RH[, order(as.numeric(substr(names(RH), 4, 6)), decreasing = TRUE)] 
    RH.on.profile <- TRUE
  }else{
    RH.on.profile <- FALSE
  }
  
  # Air pressure (can be a profile) [kPa]
  PA <- BM.data[, grep('PA', names(BM.data), fixed = TRUE)]
  # if measured along a profile
  if(is.data.frame(PA)){
    # ordering TA bottom-up (from N = ground to 1 = tower-top)
    PA <- PA[, order(as.numeric(substr(names(PA), 4, 6)), decreasing = TRUE)] 
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
  

  
  
  # Meteo data aggregation ......................................................
  # Meteo data are averaged over the full half-hour (referring to the FOLLOWING one)
  
  ### Native resolution data
  
  ## Air temperature
  # Average air temperature profile - weighted by the depths of the air layers representative of the individual measurement points
  TA.avg <- apply(TA, 1, function(x) stats::weighted.mean(x, w = Dz.i.TA, na.rm = TRUE))
  
  ## Air relative humidity 
  if(!is.null(RH)){ # 1. if air pressure is reported somehow
    if(RH.on.profile){ # 1.1 if measured along the profile
      # Average relative humidity along the profile - weighted by the depths of the air layers representative of the individual measurement points
      RH.avg <- apply(RH, 1, function(x) stats::weighted.mean(x, w = Dz.i.RH, na.rm = TRUE))
    }else{ # 1.2 if measured at a single point
      # Average air pressure along the profile estimated at the RH measuring point
      RH.avg <- RH
    }
  }else{           # 2. if air humidity is not reported
    cat(red('Error: at least one air relative humidity must be reported'))
    stop()
  }
  
  ## Air pressure
  if(!is.null(PA)){ # 1. if air pressure is reported somehow
    if(PA.on.profile) { # 1.1 if measured along the profile
      # Average air pressure along the profile - weighted by the depths of the air layers representative of the individual measurement points
      PA.avg <- apply(PA, 1, function(x) stats::weighted.mean(x, w = Dz.i.PA, na.rm = TRUE))
    }else{ # 1.2 if measured at a single point
      # Average air pressure along the profile estimated at the profile MIDDLE point
      if(max(z.i.TA) < z.i.PA){
        D.H.PA <- max(z.i.TA)/2 + z.i.PA # height difference between the PA sampling level and the profile middle
        PA.avg <- PA * (1 + D.H.PA / (Rd * c(TA.avg + 273.15) / g))  # note that the formula starts from AP at the top and the gradient is therefore added (downward increment)
      }else{
        D.H.PA <- max(z.i.TA)/2 - z.i.PA # height difference between the PA sampling level and the profile middle
        PA.avg <- PA * (1 - D.H.PA / (Rd * c(TA.avg + 273.15) / g))  # note that the formula starts from AP at the top and the gradient is therefore added (downward increment)
      }
    }
  }else{           # 2. if air pressure is not reported
    # 2.1 # Average air pressure along the profile estimated at the profile MIDDLE point
    # estimated based on site altitude (Campbell and Norman, 1998)
    D.H.PA <- max(z.i.TA)/2 # profile middle point
    PA.avg <- AP.0 * (1 - D.H.PA / (Rd * c(TA.avg + 273.15) / g))
  }
  
  
  # H2O molar fraction - resulting from hygrometers data 
  # Profile weighted average
  # NOTE: if measured along the profile, must be used to compute the dry air molar density and correct the storage flux
  #       (it can be used to cross-check the irga data anyway)
  H2O.MF.avg.meteo <- (6088.484 * exp(TA.avg/(TA.avg+237.6429)*17.31303) / 18.0153) * RH.avg / (PA.avg/Rgas*(TA.avg+273.15))
  
  
  ### Aggregate to 30 minutes (rounds time up to the nearest half-hour)
  #   and save halfhourly and vertically integaretd meteo data used in the storage flux compuation
  BM.data.avg <- data.frame(NULL)
  BM.data.avg <- data.frame(
    TIMESTAMP = format(unique(BM.time.nat_30), '%Y%m%d%H%M%S'),
    TA = tapply(TA.avg, as.double(BM.time.nat_30), mean, na.rm = TRUE),
    RH = tapply(RH.avg, as.double(BM.time.nat_30), mean, na.rm = TRUE),
    PA = tapply(PA.avg, as.double(BM.time.nat_30), mean, na.rm = TRUE),
    H2O_MF_meteo = tapply(H2O.MF.avg.meteo, as.double(BM.time.nat_30), mean, na.rm = TRUE)
  )
  
  data.table::fwrite(BM.data.avg, paste0(BM.data.path, site.ID, '_BM_SC_AVG.csv', sep=''), quote = FALSE, row.names = FALSE, sep=',')
  
  rm(BM.data.avg) # TODO: add other variables to remove
  
  cat(cyan('Done!\n'))
  
  
}