# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
surf <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_boxplots/"

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

# ===============================================================
# Correction of wind speed height, from 16 m measurement to 10 m standard

# Drag coefficients for neutral atmospheric stability (cdn) and 0 degress air temp taken from Smith 1988 JGR Table 1
# Meteo data (GreenEdge/meteo/GreenEdge_MET_2016_MGT.eps) shows air temperature was typically 0 +/- 2
cdn <- data.frame(u10 = c(2,5,10,15,20,25),
                  cdn10 = c(0.98,1.03,1.30,1.56,1.80,2.04))
# Find empirical relationship between cdn and u10
plot(cdn)
# Since we are almost always in the 5-11 m/s wind regime, I use constant cd (Large & Pond 1981 JGR)
cd <- 0.0012 # neutral drag coefficient for 5-11 m/s wind regime (Large & Pond 1981 JGR)
k <- 0.41 # von Karman constant
# So we can calculate u10 from u16 as:
# u16 = (u_star/k) * ln(z/z0)
# u_star = sqrt(cd) * u10
# u16 = ( (sqrt(cd) * u10) / k ) * ln(z/z0)
# u10 = ( u16 * k / sqrt(cd) ) / ln(z/z0)
f_u10 <- function(z, uz) {
  cd <- 0.0012 # neutral drag coefficient for 5-11 m/s wind regime (Large & Pond 1981 JGR)
  k <- 0.41 # von Karman constant
  z0 <- 0.0002 # Stull 1988 quoted in Soren Ejling Larsen PhD courses DTU  (document on MBP BSC)
  u10 <- ( uz * k / sqrt(cd) ) / log(z/z0)
  return(u10)
}
windprofile <- data.frame(z = seq(0,20,0.1))
windprofile$u10 <- f_u10(windprofile$z, 5)
plot(windprofile$u10, windprofile$z)
# STILL INCORRECT, but getting there
# Should we perhaps get rid of proportionality constants k and cd? Then, with simple logarihmic scaling
d <- 10
z <- seq(10,30,0.1) - d
z0 <- 0.0002
# cd <- 0.0012 # neutral drag coefficient for 5-11 m/s wind regime (Large & Pond 1981 JGR)
# k <- 0.41 # von Karman constant
zz <- log(z/z0)
# zz <- zz*k/sqrt(cd)
plot(zz, z)
print(zz[which(z==10)]/zz[which(z==16)])
cf <- zz[which(z==10)]/zz[which(z==16)]
# It seems that multiplying u16 by cf, or by cf^2 in quadratic model (???) should be enough to correct fluxes approximately

# ===============================================================

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

# Decide if using 3-day or day SIC and wind speed (NOW USING DAILY MEAN DATA). Select flux variable
surf$wsp <- surf$wsp24
surf$SIC <- surf$SICday
fvar <- "fdmsW97c24"

# Print FDMS prior to correction
toprint <- list(ICE = surf[surf$owd_class=="ICE",fvar],
               MIZ = surf[surf$owd_class=="MIZ",fvar],
               OW = surf[surf$owd_class=="OW",fvar])
print(lapply(toprint, summary))
print(lapply(toprint, function(x) {return(sum(!is.na(x)))}))

# CORRECT FDMS FOR SIC
surf$fdms <- surf[,fvar] * (1 - surf$SIC)

# Compute additional variables
surf$dms2dmspt <- surf$dms / surf$dmspt
# surf$kvent <- surf$fdmsW97c24 / surf$hBD_m / surf$dms
surf$kvent <- surf$fdmsW97c24 / surf$mld03 / surf$dms

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
  if (xv == "kvent") {print(lapply(toplot, summary))}
  
} # end loop on variables

if (exportimg) {dev.off()}

# ----------------------------------------------------------------
# # Linear regression model between FDMS and its controlling factors
# surf$IF <- 1 - surf$SIC
# yreg <- surf$fdms
# Xreg <- surf[,c("dms","wsp","sst","IF")]
# # Xreg <- surf[,c("wsp","dms")]
# Xreg <- as.data.frame(scale(Xreg, center = T, scale = T))
# full.model <- lm(yreg ~., data = Xreg)
# print(summary(full.model))

