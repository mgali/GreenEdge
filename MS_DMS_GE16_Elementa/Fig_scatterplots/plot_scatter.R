# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(RColorBrewer)
library(tidyverse)
library(classInt) # for function classIntervals

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)


# Exporting image?
exportimg <- F
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

# Calculate POC
cf <- 500 # [mg C/m2], Cetinic et al. 2012
pplot$poc <- pplot$cpsmooth1 * cf /12 # [mmol C/m3] = µM

# ---------------------
# Bin profiles by station categories
df2bin <- pplot
pplot$z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3)) # labels = c("0_9","10_20","21_40","41_80")
pplot$owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW"))

class <- classIntervals(var = pplot$OWD, n = 3, style="fixed", fixedBreaks = c(-35,-3.5,3.5,35))
ccol <- findColours(class, col)


# ---------------------
# FOR WIDER VIEW, NOT PAPER
# ---------------------
# 6 PANELS

# xvarS <- list(a = "cpsmooth1", d = "cpsmooth1", b = "tchla", e = "tchla", c = "dmspt", f = "dmspt")
# yvarS <- list(a = "dmspt", d = "dms", b = "dmspt", e = "dms", c = "dms", f = "dms")
# xlabS <- list(a = "", d = expression('Cp (m'^-1*')'), b = "", e = expression('TChla (µg L'^-1*')'), c = "", f = "DMSPt (nM)")
# ylabS <- list(a = "DMSPt (nM)", d = "DMS (nM)", b = "", e = "", c = "DMS (nM)", f = "")
# refline <- c(1e-4, 1e4)
# 
# if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels",".png"), width = 17, height = 12, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
# 
# # Multipanel setup
# m0 <- matrix(data = 0, nrow = 4, ncol = 4)
# mc1 <-  rbind(m0+1,m0+2)
# m <- cbind(mc1,mc1+2,mc1+4)
# layout(m)
# par(oma = c(1,1,0.5,0.5))
# 
# for (p in names(xvarS)) {
# 
#   par(mar = c(4,3,1,0.5))
#   tp <- data.frame(x = pplot[,xvarS[[p]]],
#                    y = pplot[,yvarS[[p]]],
#                    z_class = pplot$z_class,
#                    owd_class = pplot$owd_class)
# 
#   # Remove NA
#   nona <- !is.na(tp$x) & !is.na(tp$y)
#   tp <- tp[nona, ]
#   pcol <- ccol[nona]
# 
#   # Remove duplicated values variable by variable
#   dx <- as.logical(c(1,diff(tp$x)))
#   dy <- as.logical(c(1,diff(tp$y)))
#   tp <- tp[dx&dy,]
#   pcol <- pcol[dx&dy]
# 
#   # Limits
#   xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
#   yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
# 
#   plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
#   box()
#   axis(side = 1, cex.axis = 1.1)
#   axis(side = 2, cex.axis = 1.1)
#   mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
#   mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
# 
#   for (j in 10^seq(-4,4)) {
#     lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
#     lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
#   }
# 
#   points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.2)
#   points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.2)
#   points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.2)
#   points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.2)
# 
#   text(0.3*xl[2], 4.5*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.2)
#   text(0.3*xl[2], 3*yl[1], paste0("r_s = ", round(cor(tp$x, tp$y, "pairwise", "spear"), 2)), cex = 1.2)
#   text(0.3*xl[2], 2*yl[1], paste0("r_p = ", round(cor(tp$x, tp$y, "pairwise", "pears"), 2)), cex = 1.2)
#   text(0.3*xl[2], (4/3)*yl[1], paste0("r2 = ", round((cor(tp$x, tp$y, "pairwise", "pears"))^2, 2)), cex = 1.2)
# 
#   text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
# 
# }
# if (exportimg) {dev.off()}


# # ---------------------
# # FOR PAPER
# # ---------------------
# # 2 PANELS WITH DMSPt vs. Cp and DMS vs. DMSPt
# 
# xvarS <- list(a = "cpsmooth1", b = "dmspt")
# yvarS <- list(a = "dmspt", b = "dms")
# xlabS <- list(a = expression('Cp (m'^-1*')'), b = "DMSPt (nM)")
# ylabS <- list(a = "DMSPt (nM)", b = "DMS (nM)")
# xtickS <- list(a = c(.005,.01,.02,.05,.1,.2,.5,1,2),
#                b = c(2,5,10,20,50,100,200,500))
# ytickS <- list(a = c(2,5,10,20,50,100,200,500),
#                b = c(.05,.1,.2,.5,1,2,5,10,20,50,100))
# refline <- c(1e-4, 1e4)
# reflabS <- list(a = data.frame(x = c(.15,.45,1.5,3),
#                                y = c(550,550,550,350),
#                                t = as.character(c(3000,1000,300,100))),
#                 b = data.frame(x = c(60,200,400,500,500),
#                                y = c(75,75,50,18,6),
#                                t = as.character(c(1,.3,.1,.03,.01))))
# 
# if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels",".png"), width = 17, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
# 
# # Multipanel setup
# m <-  cbind(matrix(data = 1, nrow = 4, ncol = 4), matrix(data = 2, nrow = 4, ncol = 4))
# layout(m)
# par(oma = c(1,1,0.5,0.5))
# 
# for (p in names(xvarS)) {
#   
#   par(mar = c(4,4,1,1))
#   tp <- data.frame(stn = pplot$stn,
#                    cast = pplot$cast,
#                    depth = pplot$depth,
#                    x = pplot[,xvarS[[p]]],
#                    y = pplot[,yvarS[[p]]],
#                    z_class = pplot$z_class,
#                    owd_class = pplot$owd_class)
#   
#   # Remove NA
#   nona <- !is.na(tp$x) & !is.na(tp$y)
#   tp <- tp[nona, ]
#   pcol <- ccol[nona]
#   
#   # Remove duplicated values variable by variable
#   dx <- as.logical(c(1,diff(tp$x)))
#   dy <- as.logical(c(1,diff(tp$y)))
#   tp <- tp[dx&dy,]
#   pcol <- pcol[dx&dy]
#   
#   # Calculate ratio
#   tp$y2x <- tp$y/tp$x
#   
#   # Limits
#   xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
#   yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
#   
#   # Background plot with all points
#   plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
#   box()
#   axis(side = 1, cex.axis = 1.1, at = xtickS[[p]], labels = xtickS[[p]])
#   axis(side = 2, cex.axis = 1.1, at = ytickS[[p]], labels = ytickS[[p]])
#   mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
#   mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
#   
#   # Ref lines and their labels
#   for (j in 10^seq(-4,4)) {
#     lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
#     lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
#   }
#   med <- median(tp$y2x, na.rm = T)
#   lines(refline, refline*med, lty = 1, lwd = 0.5, col = "black")
#   rl <- reflabS[[p]]
#   # Ref line labels
#   for (j in length(rl$x)) {
#    text(rl$x, rl$y, rl$t, srt = 42, col = "darkgray") 
#   }
#   print(paste(p,xlabS[[p]],round(med,2),sep = " "))
#   
#   # Colored scatterplot
#   points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.5)
#   points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.5)
#   points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.5)
#   points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.5)
#   
# text(0.5*xl[2], 3*yl[1], bquote( r[S] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "spear"), 2)) ), cex = 1.4)
# text(0.5*xl[2], 2*yl[1], bquote( r[P] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "pears"), 2)) ), cex = 1.4)
# text(0.45*xl[2], (4/3)*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.4)
#   
#   # Labels and legends
#   text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
#   if (p  == "a") {
#     legend(x = xl[1],
#            y = 0.7*yl[2],
#            pch = 16,
#            legend = c("ICE","MIZ","OW"),
#            cex = 1.3,
#            col = col,
#            bg = "white", box.lwd = 0)
#     legend(x = xl[1],
#            y = 0.2*yl[2],
#            pch = c(16,17,15,18),
#            legend = c("0_9","10_20","21_40","41_80"),
#            cex = 1.3,
#            col = "black",
#            bg = "white", box.lwd = 0)
#   }
#   
#   # Outlier labelling with 1.5*IQR criterion (in log10 space if lognormal distribution of the ratio: checked beforehand)
#   if (p%in%c("a","b")) {tp$y2x <- log10(tp$y2x)}
#   qq <- quantile(tp$y2x, c(0.25, 0.75), na.rm = T)
#   iqr <- diff(qq)
#   io <- tp$y2x < qq[1]-1.5*iqr | tp$y2x > qq[2]+1.5*iqr
#   if (p%in%c("a","b")) {tp$y2x <- 10^(tp$y2x)}
#   # View(tp[io,])
#   
#   # Highlight outliers
#   points(tp$x[tp$z_class==0 & io], tp$y[tp$z_class==0 & io], pch = 1, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==1 & io], tp$y[tp$z_class==1 & io], pch = 2, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==2 & io], tp$y[tp$z_class==2 & io], pch = 0, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==3 & io], tp$y[tp$z_class==3 & io], pch = 5, col = "black", cex = 1.5)
#   
#   # Label outliers
#   for (j in which(io)) {
#     text(tp$x[j], tp$y[j], tp$stn[j], pos = 4)
#   }
#   
# }
# if (exportimg) {dev.off()}


# ---------------------
# FOR PAPER
# ---------------------
# 2 PANELS WITH DMSPt vs. Cp and DMS vs. DMSPt

xvarS <- list(a = "poc", b = "dmspt")
yvarS <- list(a = "dmspt", b = "dms")
xlabS <- list(a = "POC (µM)", b = "DMSPt (nM)")
ylabS <- list(a = "DMSPt (nM)", b = "DMS (nM)")
xtickS <- list(a = c(.02,.05,.1,.2,.5,1,2,5,10,20,50),
               b = c(2,5,10,20,50,100,200,500))
ytickS <- list(a = c(2,5,10,20,50,100,200,500),
               b = c(.05,.1,.2,.5,1,2,5,10,20,50,100))
refline <- c(1e-4, 1e4)
reflabS <- list(a = data.frame(x = c(.15,.45,1.5,3)*30,
                               y = c(550,550,550,350),
                               t = as.character(c(3000,1000,300,100)/30)),
                b = data.frame(x = c(60,200,400,500,500),
                               y = c(75,75,50,18,6),
                               t = as.character(c(1,.3,.1,.03,.01))))

if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels",".png"), width = 17, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

# Multipanel setup
m <-  cbind(matrix(data = 1, nrow = 4, ncol = 4), matrix(data = 2, nrow = 4, ncol = 4))
layout(m)
par(oma = c(1,1,0.5,0.5))

for (p in names(xvarS)) {
  
  par(mar = c(4,4,1,1))
  tp <- data.frame(stn = pplot$stn,
                   cast = pplot$cast,
                   depth = pplot$depth,
                   x = pplot[,xvarS[[p]]],
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
  
  # Calculate ratio
  tp$y2x <- tp$y/tp$x
  
  # Limits
  xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
  yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
  
  # Background plot with all points
  plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
  box()
  axis(side = 1, cex.axis = 1.1, at = xtickS[[p]], labels = xtickS[[p]])
  axis(side = 2, cex.axis = 1.1, at = ytickS[[p]], labels = ytickS[[p]])
  mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
  mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
  
  # Ref lines and their labels
  for (j in 10^seq(-4,4)) {
    lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
    lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
  }
  med <- median(tp$y2x, na.rm = T)
  lines(refline, refline*med, lty = 1, lwd = 0.5, col = "black")
  rl <- reflabS[[p]]
  # Ref line labels
  for (j in length(rl$x)) {
    text(rl$x, rl$y, rl$t, srt = 42, col = "darkgray") 
  }
  print(paste(p,xlabS[[p]],round(med*100*5*0.9/1000,2),sep = " "))
  print(paste(p,xlabS[[p]],round(max(tp$y2x, na.rm = T)*100*5*0.9/1000,2),sep = " "))
  
  # Colored scatterplot
  points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.5)
  points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.5)
  points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.5)
  points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.5)
  
  text(0.5*xl[2], 3*yl[1], bquote( r[S] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "spear"), 2)) ), cex = 1.4)
  text(0.5*xl[2], 2*yl[1], bquote( r[P] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "pears"), 2)) ), cex = 1.4)
  text(0.45*xl[2], (4/3)*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.4)
  
  # Labels and legends
  text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
  if (p  == "a") {
    legend(x = xl[1],
           y = 0.7*yl[2],
           pch = 16,
           legend = c("ICE","MIZ","OW"),
           cex = 1.3,
           col = col,
           bg = "white", box.lwd = 0)
    legend(x = xl[1],
           y = 0.2*yl[2],
           pch = c(16,17,15,18),
           legend = c("0_9","10_20","21_40","41_80"),
           cex = 1.3,
           col = "black",
           bg = "white", box.lwd = 0)
  }
  
  # Outlier labelling with 1.5*IQR criterion (in log10 space if lognormal distribution of the ratio: checked beforehand)
  if (p%in%c("a","b")) {tp$y2x <- log10(tp$y2x)}
  qq <- quantile(tp$y2x, c(0.25, 0.75), na.rm = T)
  iqr <- diff(qq)
  io <- tp$y2x < qq[1]-1.5*iqr | tp$y2x > qq[2]+1.5*iqr
  if (p%in%c("a","b")) {tp$y2x <- 10^(tp$y2x)}
  # View(tp[io,])
  print(round((10^qq)*100*5*0.9/1000,2))
  print(round((range(tp$y2x[tp$y2x<28], na.rm = T))*100*5*0.9/1000,2))
  
  # Highlight outliers
  points(tp$x[tp$z_class==0 & io], tp$y[tp$z_class==0 & io], pch = 1, col = "black", cex = 1.5)
  points(tp$x[tp$z_class==1 & io], tp$y[tp$z_class==1 & io], pch = 2, col = "black", cex = 1.5)
  points(tp$x[tp$z_class==2 & io], tp$y[tp$z_class==2 & io], pch = 0, col = "black", cex = 1.5)
  points(tp$x[tp$z_class==3 & io], tp$y[tp$z_class==3 & io], pch = 5, col = "black", cex = 1.5)
  
  # Label outliers
  for (j in which(io)) {
    text(tp$x[j], tp$y[j], tp$stn[j], pos = 4)
  }
  
}
if (exportimg) {dev.off()}


# ---------------------
# FOR SUPPORTING INFO
# ---------------------
# 2 PANELS WITH DMSPt or DMS vs. TChla

# xvarS <- list(a = "tchla", b = "tchla")
# yvarS <- list(a = "dmspt", b = "dms")
# xlabS <- list(a = expression('TChla (µg L'^-1*')'), b = expression('TChla (µg L'^-1*')'))
# ylabS <- list(a = "DMSPt (nM)", b = "DMS (nM)")
# xtickS <- list(a = c(.05,.01,.02,.05,.1,.2,.5,1,2,5,10),
#                b = c(.05,.01,.02,.05,.1,.2,.5,1,2,5,10))
# ytickS <- list(a = c(2,5,10,20,50,100,200,500),
#                b = c(.05,.1,.2,.5,1,2,5,10,20,50,100))
# refline <- c(1e-4, 1e4)
# reflabS <- list(a = data.frame(x = c(.15,.45,1.5,3),
#                                y = c(550,550,550,350),
#                                t = as.character(c(3000,1000,300,100))),
#                 b = data.frame(x = c(60,200,400,500,500),
#                                y = c(75,75,50,18,6),
#                                t = as.character(c(1,.3,.1,.03,.01))))
# 
# if (exportimg) {png(filename = paste0(opath,"Fig4_",length(xvarS),"panels_tchla",".png"), width = 17, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
# 
# # Multipanel setup
# m <-  cbind(matrix(data = 1, nrow = 4, ncol = 4), matrix(data = 2, nrow = 4, ncol = 4))
# layout(m)
# par(oma = c(1,1,0.5,0.5))
# 
# for (p in names(xvarS)) {
#   
#   par(mar = c(4,4,1,1))
#   tp <- data.frame(stn = pplot$stn,
#                    cast = pplot$cast,
#                    depth = pplot$depth,
#                    x = pplot[,xvarS[[p]]],
#                    y = pplot[,yvarS[[p]]],
#                    z_class = pplot$z_class,
#                    owd_class = pplot$owd_class)
#   
#   # Remove NA
#   nona <- !is.na(tp$x) & !is.na(tp$y)
#   tp <- tp[nona, ]
#   pcol <- ccol[nona]
#   
#   # Remove duplicated values variable by variable
#   dx <- as.logical(c(1,diff(tp$x)))
#   dy <- as.logical(c(1,diff(tp$y)))
#   tp <- tp[dx&dy,]
#   pcol <- pcol[dx&dy]
#   
#   # Calculate ratio
#   tp$y2x <- tp$y/tp$x
#   
#   # Limits
#   xl <- c(0.9*min(tp$x, na.rm = T), 1.1*max(tp$x, na.rm = T))
#   yl <- c(0.9*min(tp$y, na.rm = T), 1.1*max(tp$y, na.rm = T))
#   
#   # Background plot with all points
#   plot(tp$x, tp$y, pch = 20, col = pcol, cex = 0.1, axes = F, xlab = "", ylab = "", log = "xy", xlim = xl, ylim = yl)
#   box()
#   axis(side = 1, cex.axis = 1.1, at = xtickS[[p]], labels = xtickS[[p]])
#   axis(side = 2, cex.axis = 1.1, at = ytickS[[p]], labels = ytickS[[p]])
#   mtext(side = 1, xlabS[[p]], cex = 0.9, line = 3)
#   mtext(side = 2, ylabS[[p]], cex = 0.9, line = 2.8)
#   
#   # Ref lines and their labels
#   for (j in 10^seq(-4,4)) {
#     lines(refline, refline*j, lty = 1, lwd = 0.5, col = "gray")
#     lines(refline, refline*j*3, lty = 3, lwd = 0.5, col = "gray")
#   }
#   med <- median(tp$y2x, na.rm = T)
#   lines(refline, refline*med, lty = 1, lwd = 0.5, col = "black")
#   rl <- reflabS[[p]]
#   # Ref line labels
#   for (j in length(rl$x)) {
#     text(rl$x, rl$y, rl$t, srt = 42, col = "darkgray") 
#   }
#   print(paste(p,xlabS[[p]],round(med,2),sep = " "))
#   
#   # Colored scatterplot
#   points(tp$x[tp$z_class==0], tp$y[tp$z_class==0], pch = 16, col = pcol[tp$z_class==0], cex = 1.5)
#   points(tp$x[tp$z_class==1], tp$y[tp$z_class==1], pch = 17, col = pcol[tp$z_class==1], cex = 1.5)
#   points(tp$x[tp$z_class==2], tp$y[tp$z_class==2], pch = 15, col = pcol[tp$z_class==2], cex = 1.5)
#   points(tp$x[tp$z_class==3], tp$y[tp$z_class==3], pch = 18, col = pcol[tp$z_class==3], cex = 1.5)
#   
# text(0.5*xl[2], 3*yl[1], bquote( r[S] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "spear"), 2)) ), cex = 1.4)
# text(0.5*xl[2], 2*yl[1], bquote( r[P] ~ ' = ' ~ .(round(cor(tp$x, tp$y, "pairwise", "pears"), 2)) ), cex = 1.4)
# text(0.45*xl[2], (4/3)*yl[1], paste0("n = ", round(length(tp$x), 0)), cex = 1.4)
#   
#   # Labels and legends
#   text(1.1*xl[1], 0.9*yl[2], paste0(p,")"), cex = 1.5)
#   if (p  == "a") {
#     legend(x = xl[1],
#            y = 0.7*yl[2],
#            pch = 16,
#            legend = c("ICE","MIZ","OW"),
#            cex = 1.3,
#            col = col,
#            bg = "white", box.lwd = 0)
#     legend(x = xl[1],
#            y = 0.2*yl[2],
#            pch = c(16,17,15,18),
#            legend = c("0_9","10_20","21_40","41_80"),
#            cex = 1.3,
#            col = "black",
#            bg = "white", box.lwd = 0)
#   }
#   
#   # Outlier labelling with 1.5*IQR criterion (in log10 space if lognormal distribution of the ratio: checked beforehand)
#   if (p%in%c("a","b")) {tp$y2x <- log10(tp$y2x)}
#   qq <- quantile(tp$y2x, c(0.25, 0.75), na.rm = T)
#   iqr <- diff(qq)
#   io <- tp$y2x < qq[1]-1.5*iqr | tp$y2x > qq[2]+1.5*iqr
#   if (p%in%c("a","b")) {tp$y2x <- 10^(tp$y2x)}
#   View(tp[io,])
#   
#   # Highlight outliers
#   points(tp$x[tp$z_class==0 & io], tp$y[tp$z_class==0 & io], pch = 1, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==1 & io], tp$y[tp$z_class==1 & io], pch = 2, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==2 & io], tp$y[tp$z_class==2 & io], pch = 0, col = "black", cex = 1.5)
#   points(tp$x[tp$z_class==3 & io], tp$y[tp$z_class==3 & io], pch = 5, col = "black", cex = 1.5)
#   
#   # Label outliers
#   for (j in which(io)) {
#     text(tp$x[j], tp$y[j], tp$stn[j], pos = 4)
#   }
#   
# }
# if (exportimg) {dev.off()}

