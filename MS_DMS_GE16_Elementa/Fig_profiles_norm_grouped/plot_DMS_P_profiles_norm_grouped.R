# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)
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
col <- pal(n = 21)[c(21,18,5)]

# Rename DMS variable
prof.all$dms <- prof.all$dms_consens_cf68

# Remove unnecessary DMS variables
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[(!is.na(prof.all$dms) | !is.na(prof.all$dmspt)) & !is.na(prof.all$depth),]

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

# Remove duplicated variables
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Hide data from transect 400
pplot[pplot$stn<=400,] <- NA

# Calculate photosynthetic and photoprotective carotenoids (Bricaud 2004)
pplot$psc <- rowSums(pplot[,c("fuco","peri","but19_like","hex","hex19_likeSUM")], na.rm = T)
pplot$ppc <- rowSums(pplot[,c("zea","anthera","viola","diadino","diato","allo","tcar")], na.rm = T)
pplot$dd <- rowSums(pplot[,c("diadino","diato")], na.rm = T)
pplot$vaz <- rowSums(pplot[,c("zea","anthera","viola")], na.rm = T)
pplot$tpig <- rowSums(pplot[,seq(59,88,1)], na.rm = T)

# Add ratios
pplot$cp2tchla <- pplot$cpsmooth1/pplot$tchla                           # dms/dmspt ratio
pplot$dms2dmspt <- pplot$dms/pplot$dmspt                                # dms/dmspt ratio
pplot$dmspt2tchla <- pplot$dmspt/pplot$tchla                            # dmspt/tchla ratio
pplot$dmspt2cp <- pplot$dmspt/pplot$cpsmooth1                           # dmspt/cp ratio
pplot$ppc2psc <- pplot$ppc/pplot$psc                                    # PPC to PSC
pplot$ppc2tchla <- pplot$ppc/pplot$tchla                                # PPC to TChla
pplot$psc2tchla <- pplot$psc/pplot$tchla                                # PSC to tchla
pplot$npp <- pplot$ppc/pplot$tpig                                       # PPC to TPig
pplot$dd2tchla <- pplot$dd/pplot$tchla                                  # D+D xantophyll cycle pigments to tchla
pplot$vaz2tchla <- pplot$vaz/pplot$tchla                                # VAZ xantophyll cycle pigments to tchla
pplot$chlc3_2_tchla <- pplot$chlc3/pplot$tchla                          # chlc3 to tchla (Phaeocystis proxy?)
pplot$chlc3_2_psc <- pplot$chlc3/pplot$psc                              # chlc3 to tchla (Phaeocystis proxy?)
pplot$phaeo2chl <- pplot$phaeo_Tu_ugL/pplot$chla_Tu_ugL                 # Phaeopigments to Chl (Turner)

# Remove phaeopigments outlier
pplot[pplot$phaeo2chl > 3 & !is.na(pplot$phaeo2chl),c("phaeo_Tu_ugL","chla_Tu_ugL","phaeo2chl")] <- NA

# Bin profiles by station categories
df2bin <- pplot
z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3))
# st_class <- pplot$sic_class
st_class <- list(sic_class = pplot$sic_class,
                 owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,4.5,35), labels = c("ICE","MIZ","OW")))

# ---------------------
# Plot settings
# xvarS <- list(dms = "DMS (nM)",
# #               dmspt = "DMSPt (nM)",
# #               tchla = "TChla (µg/L)",
#               # cpsmooth1 = "Cp (1/m)",
#               cp2tchla = "Cp/TChla (m2/mg)",
# #               dms2dmspt = "DMS/DMSPt",
# #               dmspt2tchla = "DMSPt/TChla",
# #               dmspt2cp = "DMSPt/Cp",
# #               temp = "Temperature (C)",
# #               sal = "Salinity",
# #               sigt = "sigma-t (kg/m3)",
# #               anp = "ANP",
#               par_d_p24h_ein_m_2_day_1 = "PAR (µE/m2/d)")
# xvarS <- list(diat_pelagic_mg_L = "Diatoms (mg C/L)",
#               melo_mg_L = "Melosira (mg C/L)",
#               phaeo_mg_L = "Phaeocystis (mg C/L)",
#               prym_clumped_mg_L = "Prym (mg C/L)",
#               detritus_mg_L = "Detritus (mg C/L)",
#               dino_mg_L = "Dinoflagellates (mg C/L)",
#               dino_athec = "Dinoflagellates, athecate (cells/mL)",
#               dino_thec = "Dinoflagellates, thecate (cells/mL)",
#               Phaeo = "Phaeocystis (cells/mL)",
#               flag = "Flagellates (cells/mL)",
#               crypt = "Cryptophytes (cells/mL)",
#               hetero = "HNF (cells/mL)",
#               choano = "Choanoflagellates (cells/mL)",
#               cilli = "Cilliates (cells/mL)")
xvarS <- list(ppc2psc = "PPC/PSC",
              # npp = "PPC/TPig",
              # ppc2tchla = "PPC/TChla",
              psc2tchla = "PSC/TChla",
              # chlc3_2_tchla = "Chlc3/TChla",
              # chlc3_2_psc = "Chlc3/PSC",
              # phaeo2chl = "Phaeopigments/Chla (Turner)",
              dd = "(Dd+Dt)/TChla",
              vaz = "(Vi+Anth+Zea)/TChla")
yvar <- "depth"

# ---------------------
# Loop on different station classifications

for (sc in names(st_class)) {
  
  rm(pplot.bin)
  pplot.bin <- list(mean = aggregate.data.frame(df2bin,
                                                by = list(
                                                  Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                                FUN = mean,
                                                na.rm = T),
                    sd = aggregate.data.frame(df2bin,
                                                by = list(
                                                  Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                                FUN = sd,
                                                na.rm = T),
                    median = aggregate.data.frame(df2bin,
                                                  by = list(
                                                    Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                                  FUN = median,
                                                  na.rm = T),
                    min = aggregate.data.frame(df2bin,
                                                  by = list(
                                                    Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                                  FUN = min,
                                                  na.rm = T),
                    max = aggregate.data.frame(df2bin,
                                               by = list(
                                                 Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                               FUN = max,
                                               na.rm = T),
                    count = aggregate.data.frame(df2bin,
                                                 by = list(
                                                   Z_CLASS=z_class, SIC_CLASS=st_class[[sc]]),
                                                 FUN = function(x) sum(!is.na(x), na.rm = T))
  )
  
  for (xvar in names(xvarS)) {
    
    if (exportimg) {png(filename = paste0(opath,paste(sc,xvar,sep = "_"),".png"), width = 6, height = 6, units = 'cm', pointsize = 6, bg = 'white', res = 600, type = 'cairo')}
    
    print(xvar)
    xl <- c(min(c(0,1.1*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))),
            1.1*max(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))
    if (xvar  %in% c("sal","sigt")) {xl[1] <- 0.9*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T)}
    if (xvar  == "anp") {xl <- rev(xl)}
    print(xl)
    
    plot(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
         y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
         ylim = c(70,0), xlim = xl,
         pch = 19, col = col[1], cex = 1.9,
         xlab = xvarS[[xvar]], ylab = "Depth", cex.lab = 1.2)
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
  }
}

# --------------------------------------------------
# View some tables, compute some summary stats

# View(pplot[,c("stn","OWD","dms","dmspt")])
# View(pplot.bin$count[,c("stn","OWD","dms","dmspt")])
# View(pplot.bin$mean[,grep("SIC",names(pplot.bin$mean))]) # equivalent to: View(pplot.bin$mean[,c("SIC_CLASS","SICm2d","SICm1d","SICday")])
# View(pplot.bin$mea[,c("stn","OWD","dms","dmspt")])
# View(pplot.bin$mean[,c("SIC_CLASS","OWD")])
# 
# a <- as.matrix(pplot.bin$mean[,c("SICm2d","SICm1d","SICday")])
# mean(a[seq(5,6),1])
# 
# b <- as.matrix(pplot.bin$mean[,"OWD"])
# mean(b[seq(1,2),1])

dmean <- pplot.bin$mean[pplot.bin$mean$Z_CLASS==0,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm")]
dmin <- pplot.bin$min[pplot.bin$min$Z_CLASS==0,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm")]
dmax <- pplot.bin$max[pplot.bin$max$Z_CLASS==0,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm")]
View(cbind(dmean,dmin,dmax))

