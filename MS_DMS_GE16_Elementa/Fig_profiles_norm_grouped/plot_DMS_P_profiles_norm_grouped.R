# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(classInt) # for function classIntervals
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_profiles_norm_grouped/"

# ---------------------

# Color palette (Consistent with profile and surface plots made with Matlab)
pal <- colorRampPalette(brewer.pal(9, "Spectral"))
col <- pal(n = 21)[c(21,18,6)]

# Rename DMS variable
prof.all$dms <- prof.all$dms_consens_cf68

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[!is.na(prof.all$dms) | !is.na(prof.all$dmspt),]

# Add MIZ classification
icecrit1 <- 0.15
icecrit2 <- 0.85
ICE <- surf.all[,c("SICm2d","SICm1d","SICday")]
icemin <- apply(ICE, 1, min, na.rm = T)
icemax <- apply(ICE, 1, max, na.rm = T)
surf.all$sic_class <- "MIZ"
surf.all$sic_class[icemax<icecrit1] <- "OW"
surf.all$sic_class[icemin>icecrit2] <- "ICE"
surf.all$sic_class[is.na(rowSums(ICE))] <- NA

# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Hide data from transect 400
pplot[pplot$stn<=400,] <- NA

# Add chlc3 to tchla variable
pplot$chlc3_2_tchla <- pplot$chlc3/pplot$tchla

# Add ratios
pplot$xcp <- (pplot$diadino+pplot$diato)/pplot$tchla  # xantophyll cycle pigments to tchla
pplot$dms2dmspt <- pplot$dms/pplot$dmspt              # dms/dmspt ratio
pplot$dmspt2tchla <- pplot$dmspt/pplot$tchla          # dmspt/tchla ratio
pplot$dmspt2cp <- pplot$dmspt/pplot$cpsmooth1         # dmspt/cp ratio

# Bin profiles by station categories
df2bin <- pplot
z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3))

pplot.bin <- list(mean = aggregate.data.frame(df2bin,
                                              by = list(
                                                Z_CLASS=z_class, SIC_CLASS=pplot$sic_class),
                                              FUN = mean,
                                              na.rm = T),
                  median = aggregate.data.frame(df2bin,
                                                by = list(
                                                  Z_CLASS=z_class, SIC_CLASS=pplot$sic_class),
                                                FUN = median,
                                                na.rm = T),
                  count = aggregate.data.frame(df2bin,
                                               by = list(
                                                 Z_CLASS=z_class, SIC_CLASS=pplot$sic_class),
                                               FUN = function(x) sum(!is.na(x), na.rm = T))
)

# ---------------------
# Plot settings

xvar <- "dms"; xlabel <- "DMS (nM)"
# xvar <- "dmspt"; xlabel <- "DMSPt (nM)"
# xvar <- "tchla"; xlabel <- "TChla (µg/L)"
# xvar <- "cpsmooth1"; xlabel <- "Cp (1/m)"
# xvar <- "chlc3_2_tchla"; xlabel <- "Chlc3/TChla"
# xvar <- "xcp"; xlabel <- "(Dt+Dd)/TChla"
# xvar <- "dms2dmspt"; xlabel <- "DMS/DMSPt"
# xvar <- "dmspt2tchla"; xlabel <- "DMSPt/TChla"
xvar <- "dmspt2cp"; xlabel <- "DMSPt/Cp"

xl <- c(0, 1.1*max(pplot.bin$mean[,xvar], na.rm = T))
yvar <- "depth"

# ---------------------
if (exportimg) {png(filename = paste0(opath,xvar,".png"), width = 6, height = 6, units = 'cm', pointsize = 6, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
     y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
     ylim = c(70,0), xlim = xl,
     pch = 19, col = col[1], cex = 1.9,
     xlab = xlabel, ylab = "Depth", cex.lab = 1.2)
points(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="MIZ",xvar],
      y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="MIZ",yvar],
      pch = 19,  col = col[2], cex = 1.9)
points(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="OW",xvar],
       y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="OW",yvar],
       pch = 19, col = col[3], cex = 1.9)
lines(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
     y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
     col = col[1], lwd = 2)
lines(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="MIZ",xvar],
       y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="MIZ",yvar],
       col = col[2], lwd = 2)
lines(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="OW",xvar],
       y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="OW",yvar],
       col = col[3], lwd = 2)
lines(x = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",xvar],
      y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
      col = col[1], lwd = 1.5, lty = 3)
lines(x = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="MIZ",xvar],
      y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="MIZ",yvar],
      col = col[2], lwd = 1.5, lty = 3)
lines(x = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="OW",xvar],
      y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="OW",yvar],
      col = col[3], lwd = 1.5, lty = 3)

if (exportimg) {dev.off()}

