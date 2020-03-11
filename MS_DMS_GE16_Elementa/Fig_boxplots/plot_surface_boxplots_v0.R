# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
surf <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- F
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_boxplots/"

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

# ---------------------
# Replace NaN by NA
surf[is.nan(as.matrix(surf))] <- NA

# Remove data where no DMS or DMSPt are available
surf <- surf[!is.na(surf$dms) | !is.na(surf$dmspt),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# surf[surf$stn==519,c("dms","dmspt")] <- c(3.93,79.9)
surf[surf$stn==519,c("dms","dmspt")] <- c(11.42,79.9)

# Add MIZ classification by SIC
icecrit1 <- 0.15
icecrit2 <- 0.70
ICE <- surf[,c("SICm2d","SICm1d","SICday")]
icemin <- apply(ICE, 1, min, na.rm = T) # Min-max SIC criterion
icemax <- apply(ICE, 1, max, na.rm = T) # Min-max SIC criterion
icemean <- apply(ICE, 1, min, na.rm = T) # Mean concentration criterion
surf$sic_class <- NA
surf$sic_class[icemax<icecrit1] <- "OW"
surf$sic_class[icemin>icecrit2] <- "ICE"
surf$sic_class[icemax>=icecrit1 & icemin<=icecrit2] <- "MIZ"

# Add MIZ classification by OWD
surf$owd_class = cut(surf$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW"))

# # Hide data from stn 400
# surf[surf$stn<=400,] <- NA

# Decide if using 3-day or day SIC and wind speed (NOW USING DAILY MEAN DATA)
surf$wsp <- surf$wsp24
surf$SIC <- surf$SICday

# CORRECT FDMS FOR SIC
surf$fdms <- surf$fdmsW97c24 * (1 - surf$SIC)

# Compute additional variables
surf$dms2dmspt <- surf$dms / surf$dmspt
surf$kvent <- surf$fdmsW97c24 / surf$hBD_m / surf$dms

# --------------------------------------------------
# FINAL PLOTTING FOR PAPER

# Plot settings
xvarS <- list(fdms = expression('FDMS (µmol m'^-2*' d'^-1*')'),
              SIC = "SIC",
              wsp = expression('Wind speed (m s'^-1*')'),
              sst = expression('SST ('*degree*'C)'),
              dms = "DMS (nM)",
              kvent = expression('k'[vent]*' (d'^-1*')'),
              dms2dmspt = "DMS:DMSPt",
              dmspt = "DMSPt (nM)")

if (exportimg) {png(filename = paste0(opath,"Fig7_boxplots_v0.png"), width = 16, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}

# Multipanel setup
m0 <- matrix(data = 0, nrow = 4, ncol = 5)
mtop <- cbind(m0+1,m0+2,m0+3,m0+4)
m <- rbind(mtop, mtop+4)
layout(m)
par(oma = c(1,1,0.5,0.5))

for (xv in names(xvarS)) {
  
  par(mar = c(3,5,1,1))
  toplot <- list(ICE = surf[surf$owd_class=="ICE",xv],
                 MIZ = surf[surf$owd_class=="MIZ",xv],
                 OW = surf[surf$owd_class=="OW",xv])
  # boxplot(toplot,
  boxplot(surf[,xv] ~ surf$owd_class, # Equivalent to using "toplot" list created above
          col = col, # color-filled
          # border = col, # colored lines
          notch = F,
          lwd = 0.5, # 0.5 for filled, 1 for lines
          ylab = xvarS[[xv]],
          las = 1)
  points(c(1,2,3), unlist(lapply(toplot, mean, na.rm = T)), col = "white", pch = 15, cex = 0.9)
  points(c(1,2,3), unlist(lapply(toplot, mean, na.rm = T)), col = "black", pch = 0, cex = 1)
  if (xv == "fdms") {print(unlist(lapply(toplot, mean, na.rm = T)))}
  
} # end loop on variables

if (exportimg) {dev.off()}

# ----------------------------------------------------------------
# # Linear regression model between FDMS and its controlling factors
# surf$IF <- 1 - surf$SIC
# yreg <- surf$fdms 
# Xreg <- surf[,c("dms","wsp","sst","IF")]
# Xreg <- as.data.frame(scale(Xreg, center = T, scale = T))
# full.model <- lm(yreg ~., data = Xreg)
# print(summary(full.model))

