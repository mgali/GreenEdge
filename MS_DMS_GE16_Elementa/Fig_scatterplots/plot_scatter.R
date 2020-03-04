# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(tidyverse)
library(classInt) # for function classIntervals

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)


# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_scatterplots/"

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
# pal <- colorRampPalette(colors = c("red","white","blue"))       # Color palette for ANP
# col <- pal(n = 101)
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

# ---------------------
# DATA WRANGLING

# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[(!is.na(prof.all$dms) | !is.na(prof.all$dmspt)) & !is.na(prof.all$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(3.93,79.9)
prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)

# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Change units of N2 from s-2 to h-1
pplot$N2 <- sqrt(pplot$N2) * 3600

# Calculate photosynthetic and photoprotective carotenoids (Bricaud 2004)
pplot$psc <- rowSums(pplot[,c("fuco","peri","but19_like","hex","hex19_likeSUM")], na.rm = T)
pplot$ppc <- rowSums(pplot[,c("zea","anthera","viola","diadino","diato","allo","tcar")], na.rm = T)
pplot$dd <- rowSums(pplot[,c("diadino","diato")], na.rm = T)
pplot$vaz <- rowSums(pplot[,c("zea","anthera","viola")], na.rm = T)
pplot$tpig <- rowSums(pplot[,seq(59,88,1)], na.rm = T)

# ---------------------
# Bin profiles by station categories
df2bin <- pplot
pplot$z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3)) # labels = c("0_9","10_20","21_40","41_80")
pplot$owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW"))

class <- classIntervals(var = pplot$OWD, n = 3, style="fixed", fixedBreaks = c(-35,-3.5,3.5,35))
ccol <- findColours(class, col)

# ---------------------
# 6 PANELS

xvarS <- list(a = "cpsmooth1", d = "cpsmooth1", b = "tchla", e = "tchla", c = "dmspt", f = "dmspt")
yvarS <- list(a = "dmspt", d = "dms", b = "dmspt", e = "dms", c = "dms", f = "dms")
xlabS <- list(a = "", d = expression('Cp (m'^-1*')'), b = "", e = expression('TChla (µg L'^-1*')'), c = "", f = "DMSPt (nM)")
ylabS <- list(a = "DMSPt (nM)", d = "DMS (nM)", b = "", e = "", c = "DMS (nM)", f = "")
refline <- c(1e-4, 1e4)

if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels",".png"), width = 17, height = 12, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

# Multipanel setup
m0 <- matrix(data = 0, nrow = 4, ncol = 4)
mc1 <-  rbind(m0+1,m0+2)
m <- cbind(mc1,mc1+2,mc1+4)
layout(m)
par(oma = c(1,1,0.5,0.5))

for (p in names(xvarS)) {
  
  par(mar = c(4,3,1,0.5))
  tp <- data.frame(x = pplot[,xvarS[[p]]],
                   y = pplot[,yvarS[[p]]],
                   z_class = pplot$z_class,
                   owd_class = pplot$owd_class)
  
  # Remove NA
  nona <- !is.na(tp$x) & !is.na(tp$y)
  tp <- tp[nona, ]
  pcol <- ccol[nona]
  
  # Remove duplicated values variable by variable
  dx <- as.logical(c(1,diff(tp$x)))
  dy <- as.logical(c(1,diff(tp$y)))
  tp <- tp[dx&dy,]
  pcol <- pcol[dx&dy]
  
  # Limits
  xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
  yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
  
  plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
  box()
  axis(side = 1, cex.axis = 1.1)
  axis(side = 2, cex.axis = 1.1)
  mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
  mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
  
  for (j in 10^seq(-4,4)) {
    lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
    lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
  }
  
  points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.2, axes = F)
  points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.2, axes = F)
  points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.2, axes = F)
  points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.2, axes = F)
  
  text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
  text(0.3*xl[2], 4.5*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.2)
  text(0.3*xl[2], 3*yl[1], paste0("r_s = ", round(cor(tp$x, tp$y, "pairwise", "spear"), 2)), cex = 1.2)
  text(0.3*xl[2], 2*yl[1], paste0("r_p = ", round(cor(tp$x, tp$y, "pairwise", "pears"), 2)), cex = 1.2)
  text(0.3*xl[2], (4/3)*yl[1], paste0("r2 = ", round((cor(tp$x, tp$y, "pairwise", "pears"))^2, 2)), cex = 1.2)
  
}
if (exportimg) {dev.off()}


# ---------------------
# 2 PANELS

xvarS <- list(a = "cpsmooth1", b = "dmspt")
yvarS <- list(a = "dmspt", b = "dms")
xlabS <- list(a = expression('Cp (m'^-1*')'), b = "DMSPt (nM)")
ylabS <- list(a = "DMSPt (nM)", b = "DMS (nM)")
refline <- c(1e-4, 1e4)

if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels",".png"), width = 17, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

# Multipanel setup
m <-  cbind(matrix(data = 1, nrow = 4, ncol = 4), matrix(data = 2, nrow = 4, ncol = 4))
layout(m)
par(oma = c(1,1,0.5,0.5))

for (p in names(xvarS)) {
  
  par(mar = c(4,4,1,1))
  tp <- data.frame(x = pplot[,xvarS[[p]]],
                   y = pplot[,yvarS[[p]]],
                   z_class = pplot$z_class,
                   owd_class = pplot$owd_class)
  
  # Remove NA
  nona <- !is.na(tp$x) & !is.na(tp$y)
  tp <- tp[nona, ]
  pcol <- ccol[nona]
  
  # Remove duplicated values variable by variable
  dx <- as.logical(c(1,diff(tp$x)))
  dy <- as.logical(c(1,diff(tp$y)))
  tp <- tp[dx&dy,]
  pcol <- pcol[dx&dy]
  
  # Limits
  xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
  yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
  
  plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
  box()
  axis(side = 1, cex.axis = 1.1)
  axis(side = 2, cex.axis = 1.1)
  mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
  mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
  
  for (j in 10^seq(-4,4)) {
    lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
    lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
  }
  med <- median(tp$y/tp$x, na.rm = T)
  lines(refline, refline*med, lty = 1, lwd = 0.5, col = "black")
  
  points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.5, axes = F)
  points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.5, axes = F)
  points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.5, axes = F)
  points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.5, axes = F)
  
  text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
  text(0.3*xl[2], 4.5*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.2)
  text(0.3*xl[2], 3*yl[1], paste0("r_s = ", round(cor(tp$x, tp$y, "pairwise", "spear"), 2)), cex = 1.2)
  text(0.3*xl[2], 2*yl[1], paste0("r_p = ", round(cor(tp$x, tp$y, "pairwise", "pears"), 2)), cex = 1.2)
  text(0.3*xl[2], (4/3)*yl[1], paste0("r2 = ", round((cor(tp$x, tp$y, "pairwise", "pears"))^2, 2)), cex = 1.2)
  
}
if (exportimg) {dev.off()}
