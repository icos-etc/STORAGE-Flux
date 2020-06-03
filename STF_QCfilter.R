

#'' ***************************************************************************
# STORAGE FLUX Quality Control FILTER
# 
# author: Giacomo Nicolini
# contact: g.nicolini@unitus.it
# date: 2020-06-03
# For license information see LICENSE file
# 
#'' ***************************************************************************

stfQC <- function(
  SC.data.path,
  abs.lim.thrs. = NULL,
  unc.Q.thrs. = NULL,
  weekly = FALSE,
  plot.res.ts = FALSE,
  plot.res.dc = FALSE
){
  
  # Install required packages if they are not already installed ..............
  if (!require("data.table")) {install.packages("data.table"); library(data.table)}
  
  SC.df <- fread(SC.data.path, header = T, data.table = F, na.strings = '-9999' )
 
  # Local data structures
  SC.or <- SC.qc <- SC.df$'SC'
  SC.unc <- SC.df$'SC_UNC'
  SC.time <- strptime(SC.df$'TIMESTAMP_END', '%Y%m%d%H%M%S', tz = 'GMT')
  
  
  # Absolute limit filter ****************************************************
  
  if (is.null(abs.lim.thrs.)) {
    
    if (grepl('-CO2', SC.data.path)) {sc.abslim <- c(-100,100)}
    if (grepl('-LE', SC.data.path)) {sc.abslim <- c(-200,200)}
    # if (grepl('CH4', SC.data.path)) {sc.ablim <- c(-Y,Y)}
    # if (grepl('N2O', SC.data.path)) {sc.ablim <- c(-Y,Y)}
      
  }
  
  # Absolute limit flag
  abslim.flag <- c(SC.qc < sc.abslim[1] | SC.qc > sc.abslim[2])
  # Absolute limit filter
  SC.qc[abslim.flag] <- NA
  
  
  # Uncertainty filter *******************************************************
  # filtering values above the 90° % of UNC 
  
  if (weekly) {
    
    # Set uncertainty quantiles at weekly scale
    SC.unc.Q <- tapply(SC.unc, format(SC.time, '%Y-%V'), 
                       function(x) quantile(x, seq(0,1,.02), na.rm = TRUE))
    
    # Split uncertainty at weekly scale
    SC.unc.l <- tapply(SC.unc, format(SC.time, '%Y-%V'), function(x) x)
    
    # Uncertainty threshold, if not specifed, set to 95° quantile
    if (!is.null(unc.Q.thrs.)) {
      unc.Q.thrs <- lapply(SC.unc.Q, function(x) x[grep(unc.Q.thrs., names(x))])
    } else {
      unc.Q.thrs <- lapply(SC.unc.Q, function(x) x['98%'])
    }
    
    # Uncertainty flag
    i <- 1
    unc.flag <- numeric()
    for (i in 1:length(SC.unc.l)) {
      unc.flag <- c(unc.flag, SC.unc.l[[i]] > unc.Q.thrs[[i]])
    }
    unc.flag <- unc.flag == 1
    
    # Uncertainty filter
    SC.qc[unc.flag] <- NA
    
  } else {
    
    # Set uncertainty quantiles at weekly scale
    SC.unc.Q <- numeric()
    SC.unc.Q <- quantile(SC.unc, seq(0, 1,.02), na.rm = T)
    
    # Uncertainty threshold, if not specifed, set to 95° quantile
    if (!is.null(unc.Q.thrs.)) {
      unc.Q.thrs <- SC.unc.Q[grep(unc.Q.thrs., names(SC.unc.Q))]
    } else {
      unc.Q.thrs <- SC.unc.Q['98%']
    }
    
    # Uncertainty flag
    unc.flag <- SC.unc > unc.Q.thrs
    
    # Uncertainty filter
    SC.qc[unc.flag] <- NA
  }
  
  
  if (plot.res.ts) {
    
    uncColors <- adjustcolor(2, 0.5)
    
    windows()
    plot(SC.time, SC.or, type='l', col=8, ylim=c(-50,50))
    # absolute limit falgged data
    points(SC.time[abslim.flag], SC.or[abslim.flag], pch=19, col=2)
    # uncertainty flagged data
    points(SC.time[unc.flag], SC.or[unc.flag], pch=19, col=uncColors)
    points(SC.time, SC.qc, type='l', col=1)
    
  }
  
  if (plot.res.dc) {
    SC.qcd.df <- data.frame(SC = SC.qc, 
                            SC.or = SC.or,
                            SC.unc = SC.unc,
                            months.N = as.integer(format(SC.time, '%m')))
    SC.qcd.dc.avg.seas <- aggregate(SC.qcd.df, list(format(SC.time, '%B'), format(SC.time, '%R')), mean, na.rm=T)
    SC.qcd.dc.avg.seas <- SC.qcd.dc.avg.seas[order(SC.qcd.dc.avg.seas$months.N), ]
    
    time.groups <- unique(SC.qcd.dc.avg.seas$Group.1)
    if (grepl('-CO2', SC.data.path)) {ylim.sc <- c(-5,5); ylim.unc <- c(-15,15)}
    if (grepl('-LE', SC.data.path)) {ylim.sc <- c(-10,10); ylim.unc <- c(-5,5)}
    windows()
    par(mfrow=c(4,3), mar=c(.1,.1,.1,.1), oma=c(3,3,1,3))
    i <- 1
    for(i in 1:length(time.groups)) {
      
      cur.indx <- which(SC.qcd.dc.avg.seas$Group.1 == time.groups[i])
      
      # SC fluxes
      plot(SC.qcd.dc.avg.seas$SC[cur.indx], type='l', col=1, lwd=1, ylim=ylim.sc, las=1, axes = F, frame.plot = F)
      points(SC.qcd.dc.avg.seas$SC.or[cur.indx], type='l', col=8, lwd=2)
      points(SC.qcd.dc.avg.seas$SC[cur.indx], type='l', col=1, lwd=1)
      
      abline(h=0, col=8)
      mtext(time.groups[i], side = 3, -1.5)
      if(i == 1){legend('bottomleft', c('SC_or', 'SC', 'SC_unc'), lty=1, col=c(8,1,2), bty='n')}
      if(i %in% c(10,11,12)){axis(1)}
      if(i %in% c(1,4,7,10)){axis(2, las=1)}
      
      # Uncertainty
      par(new=T)
      plot(SC.qcd.dc.avg.seas$SC.unc[cur.indx], type='l', col=2, lwd=1, ylim=ylim.unc, las=1, axes = F, frame.plot = F)
      
      if(i %in% c(3,6,9,12)){axis(4, at=pretty(0:max(ylim.unc)), labels=pretty(0:max(ylim.unc)), las=1, col=2)}
      
      
    }
  }
  
  return(SC.qc)
  
}






