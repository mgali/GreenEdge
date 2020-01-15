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

# Add MIZ classification
icecrit1 <- 0.10
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
pplot[pplot$stn<500,] <- NA

# if (exportimg) {png(filename = paste0(opath,"XXX.png"), width = 6, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

xvar <- "dms_consens_cf68"
# xvar <- "dmspt"
xl <- c(0, quantile(pplot[,xvar],0.95,na.rm = T))
yvar <- "depth"
  
plot(x = pplot[pplot$sic_class=="ICE",xvar],
     y = pplot[pplot$sic_class=="ICE",yvar],
     ylim = c(80,0),
     xlim = xl,
     # log = "x",
     pch = 19,
     col = "blue")
points(x = pplot[pplot$sic_class=="MIZ",xvar],
     y = pplot[pplot$sic_class=="MIZ",yvar],
     pch = 21,
     col = "green")
points(x = pplot[pplot$sic_class=="OW",xvar],
       y = pplot[pplot$sic_class=="OW",yvar],
       pch = 20,
       col = "orange")


plot(x = pplot[pplot$anp >= 0.2,xvar],
     y = pplot[pplot$anp >= 0.2,yvar],
     ylim = c(80,0),
     xlim = xl,
     # log = "x",
     pch = 19,
     col = "blue")
points(x = pplot[pplot$anp < 0.2,xvar],
       y = pplot[pplot$anp < 0.2,yvar],
       pch = 21,
       col = "green")

plot(x = pplot[pplot$AW_ArW_clustering_coefficient >= 0.5,xvar],
     y = pplot[pplot$AW_ArW_clustering_coefficient >= 0.5,yvar],
     ylim = c(80,0),
     xlim = xl,
     # log = "x",
     pch = 19,
     col = "blue")
points(x = pplot[pplot$AW_ArW_clustering_coefficient < 0.5,xvar],
       y = pplot[pplot$AW_ArW_clustering_coefficient < 0.5,yvar],
       pch = 21,
       col = "green")


# if (exportimg) {dev.off()}

       