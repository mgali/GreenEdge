# Calculate relative weight of MIZ EDMS compared to open ocean EDMS, in 5 degrees latitude bands

datapath <- '~/Desktop/Artic_DOSES/Arctic_DOSES_timeseries/time_series_3Dmatrix/'

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_MIZ_EDMS/"
plotres <- 600

# Read data
ooedms <- list(year = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_year_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F),
               mjja = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_MJJA_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F),
               jj = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_JJ_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F))

mizedms <- list(A = read.csv(paste0(datapath,'MIZ_INTMIZ_caseB_5_x_years_65N85N_20032014_8d.csv'), header = F), # NOTE HERE I'M SWITCHING A/B NAMES
                B = read.csv(paste0(datapath,'MIZ_INTMIZ_caseA_5_x_years_65N85N_20032014.csv'), header = F))

# Scaling factors to convert MIZ area to fluxes
fixedflux <- 1 # µmol m-2 s-1 = mol km-2 d-1
ndays <- 8 # days
cf <- 32e-9 # Gg S / mol S

# ----------------------------------------------------------------------
# # Quick calculation of mean and std emission from comparable domains (done prior to conversion)
# selmiz <- 'A'
# selperiod <- 'mjja'
# 
# for (j in 1:4) {
#   mmiz <- mean(as.numeric(mizedms[[selmiz]][j,]), na.rm = T) * fixedflux * ndays * cf
#   moo <- mean(as.numeric(ooedms[[selperiod]][j+4,]), na.rm = T)
#   
#   print(paste0('Mean MIZ EDMS = ',mmiz))
#   print(paste0('Mean OO EDMS = ',moo))
#   print(paste0('Mean MIZ/(OO+MIZ) % = ',100*mmiz/(moo+mmiz)))
# }

# Scale mizedms by fixedflux, ndays and units conversion factor (INTMIZ contains only the sum of MIZ area in a given period in km2)
mizedms <- lapply(mizedms, function(x) x * fixedflux * ndays * cf)

# Name columns (latbands) and rows (years)
ooedms <- lapply(ooedms, function(x) {rownames(x) <- seq(45,85,5); colnames(x) <- seq(2003,2016,1); return(x)})
mizedms <- lapply(mizedms, function(x) {rownames(x) <- seq(65,80,5); colnames(x) <- seq(2003,2014,1); return(x)})

# Restrict ooedms to same latbands (65-70, 70-75,75-80,80,85) and years
ooedms <- lapply(ooedms, function(x) {x <- x[5:8,1:12]; return(x)})

# ----------------------------------------------------------------------
# Plot
pletS <- letters[1:dim(mizedms$A)[1]]
names(pletS) <- rownames(mizedms$A)
xp <- colnames(mizedms$A)       # x axis of bar plot
xl <- seq(1,length(xp)) - 0.5   # x axis of lines plot
colA <- "cadetblue3"
colB <- "cadetblue2"

for (iseason in names(ooedms)) {
  
  if (exportimg) {png(filename = paste0(opath,"Fig8_MIZ_EDMS_v0_",iseason,".png"), width = 8, height = 12, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 5)
  m <- rbind(m0+1,m0+2,m0+3,m0+4)
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  for (ilat in rownames(mizedms$A)) {
    
    yA <- mizedms$A[ilat,]
    yB <- mizedms$B[ilat,]
    yBmax <- max(yB, na.rm = T)
    ybarA <- t(yA * 100 / ooedms[[iseason]][ilat,])
    ybarB <- t(yB * 100 / ooedms[[iseason]][ilat,])    
    ybarBmax <- max(ybarB, na.rm = T)
    sfactor <- ybarBmax / yBmax
    
    plet <- pletS[[ilat]]
    
    par(mar = c(3,5,2,5))
    barplot(height = ybarB,
            names.arg = xp,
            ylim = c(0, ybarBmax * 1.1),
            col = colB, border = "black", lwd = 0.5,
            width = 1, space = 0, beside = T,
            legend.text = NULL,
            axes = F)
    barplot(height = ybarA,
            names.arg = xp, 
            add = T,
            col = colA, border = "black", lwd = 0.5,
            width = 1, space = 0, beside = T,
            legend.text = NULL,
            axes = F)
    ytickmax <- round(ybarlim[2], -1)
    lines(xl, yA * sfactor, col = "black", lwd = 1, lty = 1)
    lines(xl, yB * sfactor, col = "black", lwd = 1, lty = 2)
    points(xl, yA * sfactor, pch = 19, col = colA, cex = 1.5)
    points(xl, yB * sfactor, pch = 19, col = colB, cex = 1.5)
    points(xl, yA * sfactor, pch = 1, col = "black", lwd = 1, cex = 1.5)
    points(xl, yB * sfactor, pch = 1, col = "black", lwd = 1, cex = 1.5)
    # Axes
    ybarticks <- seq(0, ytickmax, length.out = 5)
    axis(side = 2, at = ybarticks, las = 1)
    yscale <- max(ybarticks / sfactor, na.rm = T)
    axis(side = 4, at = ybarticks, labels = round(ybarticks / sfactor, -round(log10(yscale) - 1)), las = 1)
    
    if (plet == "c") {
      mtext(side = 2, "                           100 * EDMS_MIZ / EDMS_OW", line = 3)
      mtext(side = 4, "                          MIZ DMS emission, Gg S / year", line = 3)
    }
  } # end loop on latitude bands
  
  if (exportimg) {dev.off()}
  
} # end loop on periods (annual, MJJA, JJ)

