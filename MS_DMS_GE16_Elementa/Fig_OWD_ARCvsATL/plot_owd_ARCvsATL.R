# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)
library(tidyverse)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
surf <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_OWD_ARCvsATL/"

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
surf$OWD_zone = cut(surf$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW"))

# Add horizontal domain classification by clustering coefficient
surf$Domain = cut(surf$AW_ArW_clustering_coefficient, breaks = c(0,0.4,1), labels = c("Arctic","Atlantic"))

# Remove data with unknown Domain
surf <- surf[!is.na(surf$Domain),]

# # Hide data from stn 400
# surf[surf$stn<=400,] <- NA

# Decide if using 3-day or day SIC and wind speed (NOW USING DAILY MEAN DATA). Select flux variable
surf$wsp <- surf$wsp24
surf$SIC <- surf$SICday
fvar <- "fdmsW97c24"

# CORRECT FDMS FOR SIC
surf$fdms <- surf[,fvar] * (1 - surf$SIC)

# Compute additional variables
surf$dms2dmspt <- surf$dms / surf$dmspt
# surf$kvent <- surf$fdmsW97c24 / surf$hBD_m / surf$dms
surf$kvent <- surf$fdmsW97c24 / surf$mld03 / surf$dms

# --------------------------------------------------
# FINAL PLOTTING FOR PAPER

# Plot settings
xvarS <- list(idmspt_z60 = expression('DMSPt (µmol m'^-2*')'),
              idms_z60 = expression('DMS (µmol m'^-2*')'),
              dms = "DMS (nM)",
              fdms = expression('FDMS (µmol m'^-2*' d'^-1*')'))
pchS <- list(Arctic = 15,
             Atlantic = 16)


# prof.plots <- data.frame(xvar = surf[,"OWD"],
#                          yvar = surf[,"idmspt_z60"])
# PUT ALL Y VARIABLES AS FACTOR COLUMN AND SPECIFY VARIABLE UNITS IN LABELLER (NOT SURE AXIS LIMITS CAN DIFFER AMONG PANELS...)

xl <- "OWD (d)"
yl <- expression('DMSPt (µmol m'^-2*')')

p <- ggplot(surf, aes(x = OWD, y = idmspt_z60, shape = Domain, colour = OWD_zone)) + geom_point(size = 3)
p + scale_color_manual(values = col) + facet_wrap(vars(yvariable), labeller = labellet(yvariable = mylabels)) + theme()

# if (exportimg) {png(filename = paste0(opath,"Figowd_ARCvsATL.png"), width = 16, height = 5, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
# 
# # Multipanel setup
# m0 <- matrix(data = 0, nrow = 4, ncol = 5)
# m <- cbind(m0+1,m0+2,m0+3,m0+4)
# layout(m)
# par(oma = c(1,1,0.5,0.5))

# for (xv in names(xvarS)) {
#   
#   for (clco in levels(surf$Domain)) {
#     
#     xy <- surf[ surf$Domain == clco , ]
#     
#     par(mar = c(3,5,1,1))
#     
#     xy1 = xy[xy$OWD_zone=="MIZ",c("OWD",xv)]
#     plot(xy1[,], xy1[,xv],
#          pch = pchS[[clco]],
#          col = col[2], # color-filled
#          ylab = xvarS[[xv]])
#     xy2 = xy[xy$OWD_zone=="ICE",c("OWD",xv)]
#     if (sum(!is.na(xy2[,xv])) != 0) {points(xy2[,"OWD"], xy2[,xv], col = col[1], pch = pchS[[clco]])}
#     xy3 = xy[xy$OWD_zone=="OW",c("OWD",xv)]
#     if (sum(!is.na(xy3[,xv])) != 0) {points(xy3[,"OWD"], xy3[,xv], col = col[3], pch = pchS[[clco]])}
#     
#   } # end loop on horizontal domains
# } # end loop on variables
# 
# if (exportimg) {dev.off()}
