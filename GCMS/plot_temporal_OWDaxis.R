# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(tidyverse)
library(zoo)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_temporal_OWDaxis/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
df <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Binning window OWD width
bb <- 5
smooth_bin <- 'b'

if (smooth_bin=='b') {
  fname <- paste0(genpath,pdir,'Fig4_bin',bb,'owd.png')
} else if (smooth_bin=='s') {
  fname <- paste0(genpath,pdir,'Fig4_smooth',bb,'owd.png')
}

# Exporting image?
exportimg <- T

if (exportimg) {
  png(filename = fname,
      width = 8, height = 16, units = 'cm', pointsize = 8, # for 4 panels, paper
      # width = 8, height = 8, units = 'cm', pointsize = 8, # for 2 panels, ppt
      bg = 'white', res = 600, type = 'cairo')
}

# ---- Remove rows with no DMS data
df <- df[!is.na(df$dms),]

# Remove columns with no data
df <- as.matrix(df) # Convert to matrix for col selection
df <- df[,colSums(is.nan(df)) < dim(df)[1]-2]
df <- as.data.frame(df)
print(summary(df$fdmsW97c24))

# Add TChla by merging with prof.all (formerly used to add cpsmooth1, now corrected)
# IMPROVE TO AVOID DUPLICATIONS OF DATA!
prof.all <- prof.all[prof.all$depth<4,c('cast','tchla')]
# df <- select(df, -c(cpsmooth1,TChla))
df <- merge(x = df, y = prof.all, all.y = T, all.x = F, by = 'cast')

# Scale fdms by sic
df$fdms <- df$fdmsW97c24*(1-df$SICday)
print(summary(df$fdms))

# Calculate mean PAR in hBD layer
df$kd <- log(df$PAR_at_3m_E_m2_d / 0.415) / df$isolume_m_at_0415
df$par_1m <- df$PAR_at_3m_E_m2_d * exp((3-1)*0.1)
df$par_srd <- df$par_1m * (1 / (df$kd*df$hBD_m)) * (1 - exp(-df$kd*df$hBD_m))

# Sort by OWD
df <- df[order(df$OWD),]

# Remove data from stn 400 which has nothing to do with other tr.400 stations
df <- df[df$stn > 400,]

# Interpolate by daily OWD, then bin or smooth (bb width)
dxi <- 0.5 # important: smaller interp interval dxi makes line more spiky, and viceversa. Use 3 to 6
owd <- df$OWD
owdi <- seq(min(owd, na.rm = T), max(owd, na.rm = T), dxi)
# owdi <- seq(-17, 33, dxi)
dfi <- lapply(df, function(x) {y <- approx(owd, x, owdi); return(y$y)})
dfi <- as.data.frame(dfi)

# Bin or smooth to selected OWD period
if (smooth_bin=='b') {
  b4vec <- floor((dfi$OWD + bb/2) / bb) * bb
  dfi <- aggregate(dfi, by = list(OWD = b4vec), function(x) mean(x, na.rm = T), simplify = TRUE) # This can be done with rollapply with argument by = bb
} else if (smooth_bin=='s') {
  dfi <- as.data.frame(rollapply(data = dfi, width = round(bb/dxi), mean))
}

# Merge dfi with df, output is dfi with some gaps left by interpolation filled
# dfim <- merge(x = df, y = dfi, by = 'OWD', all.x = F, all.y = T)

# Set xlim and ylim for DMS and other variables (take into account scaling factor)
xl <- c(min(df$OWD, na.rm = T), max(df$OWD, na.rm = T))
ylfdms <- c(0,30)
yldms <- c(0,30)
yldmspt <- c(0,300)
yltchla <- c(0,5)
ylsst <- c(-2,4)
ylpar <- c(0,25)
ylhBD <- c(50,0)

# Set scaling factors
f.dmspt <- 0.1
f.cp <- 10
f.par <- 1/5
f.no3 <- 5

# Scale salinity
linmod <- lm(c(-2,4) ~ c(31,34))
# print(linmod$coefficients)

# X ticks. If specified, override prevously defined ticks
xticks <- seq(-20, 35, 5)
xl <- c(min(xticks, na.rm = T), max(xticks, na.rm = T))

# Set colors
col.mld <- rgb(0.5, 0.5, 0.5, alpha = 1)
col.mld.line <- rgb(0.5, 0.5, 0.5, alpha = 0.5)
col.fdms <- rgb(255/255, 165/255, 0, alpha = 1)
col.fdms.line <- rgb(255/255, 165/255, 0, alpha = 0.3)
# col.fdms <- rgb(1, .63, .48, alpha = 1)
# col.fdms.line <- rgb(1, .63, .48, alpha = 0.3)
# col.fdms <- rgb(1, .6, 0, alpha = 1)
# col.fdms.line <- rgb(1, .6, 0, alpha = 0.3)
col.dms <- rgb(1, 0, 0, alpha = 1)
col.dms.line <- rgb(1, 0, 0, alpha = 0.3)
# col.dmspt <- rgb(0.54, 0.17, .89, alpha = 1)
# col.dmspt.line <- rgb(0.54, 0.17, .89, alpha = 0.3)
col.dmspt <- rgb(0, 0, .85, alpha = 1)
col.dmspt.line <- rgb(0, 0, .85, alpha = 0.3)
col.tchla <- rgb(30/255, 150/255, 0/255, alpha = 1)
col.tchla.line <- rgb(30/255, 150/255, 0/255, alpha = 0.3)
col.cp <- rgb(204/255, 204/255, 0, alpha = 1)
col.cp.line <- rgb(204/255, 204/255, 0, alpha = 0.3)
# col.cp <- rgb(65/255, 105/255, 225/255, alpha = 1)
# col.cp.line <- rgb(65/255, 105/255, 225/255, alpha = 0.3)
col.sst <- rgb(1, 0, .5, alpha = 1)
col.sst.line <- rgb(1, 0, .5, alpha = 0.3)
col.sss <- rgb(0, 0, 0, alpha = 1) #rgb(.6, 0, .6, alpha = 1)
col.sss.line <- rgb(0, 0, 0, alpha = 0.5) #rgb(.6, 0, .6, alpha = 0.3)
col.par <- col.mld
col.par.line <- col.mld.line
# col.par <- rgb(0, .75, 1, alpha = 1)
# col.par.line <- rgb(0, .75, 1, alpha = 0.3)
col.hBD <- rgb(0, 0, 0, alpha = 1)
col.hBD.line <- rgb(0, 0, 0, alpha = 0.5)
col.dbm <- col.cp
col.dbm.line <- col.cp.line
col.ncline <- rgb(255/255, 165/255, 0, alpha = 1)
col.ncline.line <- rgb(255/255, 165/255, 0, alpha = 0.3)
# col.ncline <- rgb(102/255, 102/255, 0, alpha = 1)
# col.ncline.line <- rgb(102/255, 102/255, 0, alpha = 0.3)
# col.ncline <- col.sst
# col.ncline.line <- col.sst.line

# ------------------------- PLOT ------------------------- 

# Layout for paper
m <- rbind(matrix(data=1,nrow = 4, ncol = 6), matrix(data=2,nrow = 4, ncol = 6), matrix(data=3,nrow = 4, ncol = 6))
m <- rbind(m,matrix(data=4,nrow = 5, ncol = 6)) # Add a fourth panel
# m <- rbind(m,matrix(data=5,nrow = 5, ncol = 6)) # Add a fifth panel

# # Layout for ppt
# m <- rbind(matrix(data=1,nrow = 4, ncol = 6), matrix(data=2,nrow = 5, ncol = 6))

layout(m) #, widths = lcm(12), heights = lcm(4))
# layout.show(5)
par(oma = c(1,1,0.5,0.5))

# ---- Panel A
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$fdms, axes = F, xlim = xl, ylim = ylfdms, col = col.fdms, pch = 20, ylab = "")
lines(spline(dfi$OWD, dfi$fdms, n = 201, method = 'natural'), col = col.fdms.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(0,30,10))
mtext(side = 2, expression(paste('DMS (nM), FDMS (µmol ',m^-2,' ',d^-1,')')), cex = 0.8, line = 3)
points(df$OWD, df$dms, col = col.dms, pch = 20)
lines(spline(dfi$OWD, dfi$dms, n = 201, method = 'natural'), col = col.dms.line, lwd = 3)
points(df$OWD, df$dmspt * f.dmspt, ylim = yldmspt, col = col.dmspt, pch = 20)
lines(spline(dfi$OWD, dfi$dmspt * f.dmspt, n = 201, method = 'natural'), col = col.dmspt.line, lwd = 3)
axis(4, tcl = -0.3, at = seq(0,30,10), labels = seq(0,300,100))
mtext('DMSPt (nM)', side = 4, cex = .8, line = 3, srt = 180)
grid(nx = NULL, ny = NA)

# Legend: DMS, FDMS, DMSPt
legend(x = 23, y = 30, pch = 20,
       cex = 0.9,
       lwd = rep(2,3),
       legend = c('DMS','FDMS','DMSPt'),
       col = c(col.dms,col.fdms,col.dmspt),
       bg= "gray97", box.col = "gray97")

# ---- Panel B
# For paper
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$tchla, axes = F, xlim = xl, ylim = yltchla, col = col.tchla, pch = 20, ylab = "")
# # For ppt
# par(mar = c(4,4,0.5,4))
# plot(df$OWD, df$tchla, axes = F, xlim = xl, ylim = yltchla, col = col.tchla, pch = 20, ylab = "", xlab = "Open Water Days")
# mtext(expression(paste('TChla (µg ',L^-1,')')), side = 2, cex = 0.8, line = 3) # Put ylabel out of plot to control position (line)
mtext(expression(paste('TChla (µg ',L^-1,'), Cp*10 (',m^-1,')')), side = 2, cex = 0.8, line = 3) # Put ylabel out of plot to control position (line)
lines(spline(dfi$OWD, dfi$tchla, n = 201, method = 'natural'), col = col.tchla.line, lwd = 3, cex.lab = 1.2)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(0,5,1))
points(df$OWD, df$par_srd * f.par, col = col.par, pch = 20)
lines(spline(dfi$OWD, dfi$par_srd * f.par, n = 201, method = 'natural'), col = col.par.line, lwd = 3)
axis(4, tcl = -0.3, at = seq(0,5,1), labels = seq(0,25,5))
mtext(expression(paste(PAR[hBD],' (mol ',m^-2,' ',d^-1,')')), side = 4, cex = 0.8, line = 3)
points(df$OWD, df$cpsmooth1 * f.cp, col = col.cp, pch = 20)
lines(spline(dfi$OWD, dfi$cpsmooth1 * f.cp, n = 201, method = 'natural'), col = col.cp.line, lwd = 3)
# axis(4, tcl = -0.3, at = seq(0,4,1), labels = seq(0,0.8,0.2))
# mtext(expression(paste('Cp (',m^-1,')')), side = 4, cex = .8, line = 3, srt = 180)
grid(nx = NULL, ny = NA)

# Legend: TChla, Cp
legend(x = -20, y = 5, pch = 20,
       cex = 0.9,
       lwd = rep(2,2),
       legend = c('TChla','Cp*10',expression(paste(PAR[hBD]))),
       col = c(col.tchla,col.cp,col.par),
       bg= "gray97", box.col = "gray97")

# ---- Panel C
par(mar = c(1,4,0.5,4))
plot(df$OWD, df$sst, axes = F, xlim = xl, ylim = ylsst, col = col.sst, pch = 20, ylab = 'SST (C)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$sst, n = 201, method = 'natural'), col = col.sst.line, lwd = 3)
axis(1, labels = F, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3, at = seq(-2,4,2))
points(df$OWD, df$sal * 2 - 64, ylim = yldmspt, col = col.sss, pch = 20)
lines(spline(dfi$OWD, dfi$sal  * 2 - 64, n = 201, method = 'natural'), col = col.sss.line, lwd = 3)
axis(4, tcl = -0.3, at = seq(-2,4,2), labels = seq(31,34,1))
mtext('SSS', side = 4, cex = .8, line = 2.5, srt = 180)
grid(nx = NULL, ny = NA)

# Legend: SST, SSS
legend(x = -20, y = 4, pch = 20,
       cex = 0.9,
       lwd = rep(2,2),
       legend = c('SST','SSS'),
       col = c(col.sst,col.sss),
       bg= "gray97", box.col = "gray97")

# # ---- Panel D
# par(mar = c(1,4,0.5,4))
# plot(df$OWD, df$par_srd, axes = F, xlim = xl, ylim = ylpar, col = col.par, pch = 20, ylab = "")
# mtext(expression(paste(PAR[hBD],' (mol ',m^-2,' ',d^-1,')')), side = 2, cex = 0.8, line = 3)
# lines(spline(dfi$OWD, dfi$par_srd, n = 201, method = 'natural'), col = col.par.line, lwd = 3)
# axis(1, labels = F, tcl = -0.3, at = xticks)
# axis(2, labels = T, tcl = -0.3, at = seq(0,25,5))
# points(df$OWD, df$NO3 * f.no3, col = col.sss, pch = 20)
# lines(spline(dfi$OWD, dfi$NO3 * f.no3, n = 201, method = 'natural'), col = col.sss.line, lwd = 3)
# axis(4, tcl = -0.3, at = seq(0,25,5), labels = seq(0,5,1))
# mtext('Nitrate (µM)', side = 4, cex = .8, line = 2.5, srt = 180)
# grid(nx = NULL, ny = NA)

# ---- Panel bottom
par(mar = c(4,4,0.5,4))
plot(df$OWD, df$hBD_m, axes = F, xlim = xl, ylim = ylhBD, col = col.hBD, pch = 20, xlab = 'Open Water Days', ylab = 'Depth (m)', cex.lab = 1.2)
lines(spline(dfi$OWD, dfi$hBD_m, n = 201, method = 'natural'), col = col.hBD.line, lwd = 3)
axis(1, labels = T, tcl = -0.3, at = xticks)
axis(2, labels = T, tcl = -0.3)
points(df$OWD, df$mld03, col = col.mld, pch = 20)
lines(spline(dfi$OWD, dfi$mld03, n = 201, method = 'natural'), col = col.mld.line, lwd = 3)
points(df$OWD, df$dbm, col = col.dbm, pch = 20)
lines(spline(dfi$OWD, dfi$dbm, n = 201, method = 'natural'), col = col.dbm.line, lwd = 3)
points(df$OWD, df$Nitracline_m, col = col.ncline, pch = 20)
lines(spline(dfi$OWD, dfi$Nitracline_m, n = 201, method = 'natural'), col = col.ncline.line, lwd = 3)
grid(nx = NULL, ny = NA)

# Legend: depth of hBD, mld_0.03, dbm
legend(x = -20, y = 35, pch = 20,
       cex = 0.9,
       lwd = rep(2,3),
       legend = c('hBD','MLD0.03','DBM','NO3-cline'),
       col = c(col.hBD,col.mld,col.dbm,col.ncline),
       bg= "gray97", box.col = "gray97")

if (exportimg) {dev.off()}
