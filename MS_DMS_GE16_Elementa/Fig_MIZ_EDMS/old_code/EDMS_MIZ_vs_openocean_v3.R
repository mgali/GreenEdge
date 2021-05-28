# Calculate relative weight of MIZ EDMS compared to open ocean EDMS, in 5 degrees latitude bands

datapath <- '~/Desktop/DMS-SAT/Artic_DOSES/Arctic_DOSES_timeseries/time_series_3Dmatrix/'

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
fixedflux <- 5 # µmol m-2 d-1 = mol km-2 d-1
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
names(pletS) <- rev(rownames(mizedms$A))
xp <- colnames(mizedms$A)       # x axis of bar plot
xl <- seq(1,length(xp)) - 0.5   # x axis of lines plot
colA <- "cadetblue3"
# colB <- "cadetblue2"
colB <- "azure1"
m0 <- matrix(data = 0, nrow = 4, ncol = 5)
m <- rbind(m0+1,m0+2,m0+3,m0+4)

for (iseason in names(ooedms)[2]) {
  
  if (exportimg) {png(filename = paste0(opath,"Fig9_MIZ_EDMS_v3_",iseason,".png"), width = 8, height = 10, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

  # Multipanel setup
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  # Common left y axes for all panels
  ybarBmax <- max(mizedms$B, na.rm = T) * 1.4
  
  for (ilat in rev(rownames(mizedms$A))) {
    
    ybarA <- mizedms$A[ilat,]
    ybarB <- mizedms$B[ilat,]  
    yfracA <- t(ybarA * 100 / ooedms[[iseason]][ilat,])
    yfracB <- t(ybarB * 100 / ooedms[[iseason]][ilat,])
    
    # Different right y axes
    yfracBmax <- max(yfracB, na.rm = T) * 1.4
    sfactor <- ybarBmax / yfracBmax
    
    # Output
    print(mean(yfracA))
    print(sd(yfracA))
    print(mean(yfracB))
    print(sd(yfracB))
    
    par(mar = c(2,5,1,5))
    barplot(height = unlist(ybarB),
            names.arg = xp,
            ylim = c(0, ybarBmax),
            col = colB, border = "black", lwd = 0.5,
            width = 1, space = 0, beside = T,
            legend.text = NULL,
            axes = F)
    barplot(height = unlist(ybarA),
            names.arg = xp, 
            add = T,
            col = colA, border = "black", lwd = 0.5,
            width = 1, space = 0, beside = T,
            legend.text = NULL,
            axes = F)
    lines(xl, yfracA * sfactor, col = "black", lwd = 1, lty = 1)
    lines(xl, yfracB * sfactor, col = "black", lwd = 1, lty = 2)
    points(xl, yfracA * sfactor, pch = 19, col = colA, cex = 1.2)
    points(xl, yfracB * sfactor, pch = 19, col = colB, cex = 1.2)
    points(xl, yfracA * sfactor, pch = 1, col = "black", lwd = 1, cex = 1.2)
    points(xl, yfracB * sfactor, pch = 1, col = "black", lwd = 1, cex = 1.2)
    
    # Axes
    nticks <- 4
    yint <- round(ybarBmax / nticks, -round(log10(ybarBmax) - 1))
    ybarBticks <- seq(0, ybarBmax, yint)
    yintB <- round(yfracBmax / nticks, -round(log10(yfracBmax) - 1))
    yfracBticks <- seq(0, yfracBmax, yintB)
    axis(side = 2, at = ybarBticks, las = 1, cex.axis = 1.1)
    axis(side = 4, at = yfracBticks * sfactor, labels = yfracBticks, las = 1, cex.axis = 1.1)
    
    # Axis and panel labels
    plet <- pletS[[ilat]]
    t1 <- paste0(plet, ") ") 
    t2 <- paste0(ilat, "-", as.character(as.numeric(ilat)+5))
    t3 <- expression('             '*degree*'N')
    text(0.1, max(ybarBticks) * 0.96, t1, cex = 1.2)
    mtext(side = 3, t2, line = -1.5, cex = 0.8, font.axis=2)
    mtext(side = 3, t3, line = -1.5, cex = 0.8, font.axis=2)
    if (plet == "c") {
      mtext(side = 2, expression("                          E"[DMS_MIZ]*", Gg S yr"^-1), line = 3, cex = 1.1)
      mtext(side = 4, expression("                          -o-   100 x E"[DMS_MIZ]*" / E"[DMS_OW]), line = 3.5, cex = 1.1)
    }
  } # end loop on latitude bands
  
  if (exportimg) {dev.off()}
  
} # end loop on periods (annual, MJJA, JJ)


# ------------------------------------------
# # Old plot settings
# ybarBmaxS <- c(10, 10, 50, 100)
# names(ybarBmaxS) <- rownames(mizedms$A)
# yfracBmaxS <- c(1.25, 1.25, 1.25, 0.5)
# names(yfracBmaxS) <- rownames(mizedms$A)
# maxmax <- max(unlist(mizedms), na.rm = T)


# Rate of increase north of 80N
x <- seq(2003,2014)/10
y <- unlist(mizedms$B[4,])
fit <- lm(y ~ x)
fit$coefficients[2] * 100 / mean(y)


# Pan-Arctic sums
colsums <- lapply(mizedms, colSums)




