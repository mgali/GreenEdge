# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(tidyverse)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_temporal_OWDaxis/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
df <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# ---- Remove rows with no DMS data
df <- df[!is.na(df$dms),]

# Remove columns with no data
df <- as.matrix(df) # Convert to matrix for col selection
df <- df[,colSums(is.nan(df)) < dim(df)[1]-2]
df <- as.data.frame(df)

# Add TChla and cpsmooth by merging with prof.all
prof.all <- prof.all[prof.all$depth<2,c('cast','tchla','cpsmooth1')]
# df <- select(df, -c(cpsmooth1,TChla))
df <- merge(x = df, y = prof.all, all.y = T, all.x = F, by = 'cast')

# Scale fdms by sic
df$fdms <- df$fdmsW97c24*(1-df$SICday)

# Calculate mean PAR in hBD layer
df$kd <- log(df$PAR_at_3m_E_m2_d / 0.415) / df$isolume_m_at_0415
df$par_1m <- df$PAR_at_3m_E_m2_d * exp((3-1)*0.1)
df$par_srd <- df$par_1m * (1 / (df$kd*df$hBD_m)) * (1 - exp(-df$kd*df$hBD_m))

# Sort by OWD
df <- df[order(df$OWD),]

# Interpolate by daily OWD
dxi <- 5 # important: smaller interp interval dxi makes line more spiky, and viceversa
owd <- df$OWD
owdi <- seq(min(owd, na.rm = T), max(owd, na.rm = T), dxi)
dfi <- lapply(df, function(x) {y <- approx(owd, x, owdi); return(y$y)})
dfi <- as.data.frame(dfi)

# Merge dfi with df, output is dfi with some gaps left by interpolation filled
# dfim <- merge(x = df, y = dfi, by = 'OWD', all.x = F, all.y = T)

# Set xlim and ylim for DMS and other variables (take into account scaling factor)
xl <- c(min(df$OWD, na.rm = T), max(df$OWD, na.rm = T))
ylfdms <- c(0,35)
yldms <- c(0,30)
yldmspt <- c(0,300)
yltchla <- c(0,4)
ylsst <- c(-2,4)
ylpar <- c(0,25)
ylhBD <- c(40,0)

# Set scaling factors
f.dmspt <- 0.1
f.cp <- 4
f.par <- 0.1

# Scale salinity
linmod <- lm(c(-2,4) ~ c(31,34))
print(linmod$coefficients)

# X ticks. If specified, override prevously defined ticks
xticks <- seq(-20, 35, 5)
xl <- c(min(xticks, na.rm = T), max(xticks, na.rm = T))

# Set colors
col.fdms <- rgb(red = 1, green = .6, blue = 0, alpha = 1)
col.fdms.line <- rgb(red = 1, green = .6, blue = 0, alpha = 0.3)
col.dms <- rgb(red = 1, green = 0, blue = 0, alpha = 1)
col.dms.line <- rgb(red = 1, green = 0, blue = 0, alpha = 0.3)
col.dmspt <- rgb(red = 0, green = 0, blue = 1, alpha = 1)
col.dmspt.line <- rgb(red = 0, green = 0, blue = 1, alpha = 0.3)
col.tchla <- rgb(red = 76/255, green = 153/255, blue = 0, alpha = 1)
col.tchla.line <- rgb(red = 76/255, green = 153/255, blue = 0, alpha = 0.3)
col.cp <- rgb(red = 204/255, green = 204/255, blue = 0, alpha = 1)
col.cp.line <- rgb(red = 204/255, green = 204/255, blue = 0, alpha = 0.3)
col.sst <- rgb(red = 1, green = 0, blue = .5, alpha = 1)
col.sst.line <- rgb(red = 1, green = 0, blue = .5, alpha = 0.3)
col.sss <- rgb(red = .6, green = 0, blue = .6, alpha = 1)
col.sss.line <- rgb(red = .6, green = 0, blue = .6, alpha = 0.3)
col.par <- rgb(red = 1, green = .9, blue = 0, alpha = 1)
col.par.line <- rgb(red = 1, green = .9, blue = 0, alpha = 0.3)
col.hBD <- rgb(red = 0, green = 0, blue = 0, alpha = 1)
col.hBD.line <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
col.mld <- rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 1)
col.mld.line <- rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5)
col.dbm <- col.cp
col.dbm.line <- col.cp.line
# col.dbm <- rgb(red = 204/255, green = 204/255, blue = 0, alpha = 1)
# col.dbm.line <- rgb(red = 204/255, green = 204/255, blue = 0, alpha = 0.3)
# col.ncline <- rgb(red = 102/255, green = 102/255, blue = 0, alpha = 1)
# col.ncline.line <- rgb(red = 102/255, green = 102/255, blue = 0, alpha = 0.3)

# ------------------------- PLOT ------------------------- 
# png(filename = paste0(genpath,pdir,'Fig4','.png'),
#     width = 8, height = 16, units = 'cm', pointsize = 8, # height 10 for 3 panels, 13 for 4 panels
#     bg = 'white', res = 300, type = 'cairo')

# Layout
m <- rbind(matrix(data=1,nrow = 4, ncol = 6), matrix(data=2,nrow = 4, ncol = 6), matrix(data=3,nrow = 4, ncol = 6))
m <- rbind(m,matrix(data=4,nrow = 4, ncol = 6)) # Add a fourth panel
m <- rbind(m,matrix(data=5,nrow = 5, ncol = 6)) # Add a fifth panel
layout(m) #, widths = lcm(12), heights = lcm(4))
# layout.show(5)
par(oma = c(1,1,0.5,0.5))

# Panel A
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$fdms, axes = F, xlim = xl, ylim = ylfdms, col = col.fdms, pch = 19, ylab = 'DMS, DMSPt/10 (nM)', cex.lab = 1.1)
lines(spline(dfi$OWD, dfi$fdms, n = 201, method = 'natural'), col = col.fdms.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(0,30,10))
points(df$OWD, df$dms, col = col.dms, pch = 19)
lines(spline(dfi$OWD, dfi$dms, n = 201, method = 'natural'), col = col.dms.line, lwd = 3)
points(df$OWD, df$dmspt * f.dmspt, ylim = yldmspt, col = col.dmspt, pch = 19)
lines(spline(dfi$OWD, dfi$dmspt * f.dmspt, n = 201, method = 'natural'), col = col.dmspt.line, lwd = 3)
axis(4, labels = T, tcl = -0.3, at = seq(0,30,10))
mtext('FDMS (µmol/m2/d)', side = 4, cex = .8, line = 2.5, srt = 180)

# Panel B
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$tchla, axes = F, xlim = xl, ylim = yltchla, col = col.tchla, pch = 19, ylab = 'TChla (µg/L)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$tchla, n = 201, method = 'natural'), col = col.tchla.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3)
points(df$OWD, df$cpsmooth1 * f.cp, col = col.cp, pch = 19)
lines(spline(dfi$OWD, dfi$cpsmooth1 * f.cp, n = 201, method = 'natural'), col = col.cp.line, lwd = 3)
axis(4, labels = T, tcl = -0.3, at = seq(0,4,1))
mtext('Cp (1/m)', side = 4, cex = .8, line = 2.5, srt = 180)

# Panel C
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$sst, axes = F, xlim = xl, ylim = ylsst, col = col.sst, pch = 19, ylab = 'SST (C)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$sst, n = 201, method = 'natural'), col = col.sst.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(-2,4,2))
points(df$OWD, df$sal * 2 - 64, ylim = yldmspt, col = col.sss, pch = 19)
lines(spline(dfi$OWD, dfi$sal  * 2 - 64, n = 201, method = 'natural'), col = col.sss.line, lwd = 3)
axis(4, tcl = -0.3, at = seq(-2,4,2), labels = seq(31,34,1))
mtext('SSS', side = 4, cex = .8, line = 2.5, srt = 180)

# Panel D
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$par_srd, axes = F, xlim = xl, ylim = ylpar, col = col.par, pch = 19, ylab = 'PAR (mol/m2/d)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$par_srd, n = 201, method = 'natural'), col = col.par.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(0,25,5))

# Panel E
par(mar = c(4,4,0.5,4))
plot(df$OWD, df$hBD_m, axes = F, xlim = xl, ylim = ylhBD, col = col.hBD, pch = 19, xlab = 'Open Water Days', ylab = 'Depth (m)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$hBD_m, n = 201, method = 'natural'), col = col.hBD.line, lwd = 3)
axis(1, labels = T, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3)
points(df$OWD, df$mld03, col = col.mld, pch = 19)
lines(spline(dfi$OWD, dfi$mld03, n = 201, method = 'natural'), col = col.mld.line, lwd = 3)
points(df$OWD, df$dbm, col = col.dbm, pch = 19)
lines(spline(dfi$OWD, dfi$dbm, n = 201, method = 'natural'), col = col.dbm.line, lwd = 3)
# points(df$OWD, df$Nitracline_m, col = col.ncline, pch = 19)
# lines(spline(dfi$OWD, dfi$Nitracline_m, n = 201, method = 'natural'), col = col.ncline.line, lwd = 3)


# dev.off()
