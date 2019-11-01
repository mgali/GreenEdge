# Martí Galí Tàpias, May 2019

# Install packages if needed

# # If using https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
# install.packages("tidyverse")                # tidyR, dplyr, ggplot2
# install.packages(c("devtools", "stringr"))   # basic packages
# install.packages(c("maps", "mapdata"))       # some standard map packages.
# devtools::install_github("dkahle/ggmap")     # the github version of ggmap

# # If using https://hansenjohnson.org/post/bathymetric-maps-in-r/
# install.packages("marmap")
# install.pacakges("oce")

# Load packages
library(marmap)
library(oce)
library(ocedata)
library(RNetCDF)

# NOTE: oce compilation fails at first
# "make: gfortran-4.8: No such file or directory"
# Fix proposed in https://stackoverflow.com/questions/23916219/os-x-package-installation-depends-on-gfortran-4-8
# But adding the -L flag to curl as proposed in https://github.com/arq5x/bedtools2/issues/189
# So the fix involves typng in the terminal
# "curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2"
# "sudo tar fvxj gfortran-4.8.2-darwin13.tar.bz2 -C /"
# And then intall.packages(oce) works

# Load data
data("coastlineWorldFine")
fpath_ge_stn <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig1_map/"
ge_stn_all <- read.csv(file=paste0(fpath_ge_stn,"GE_Amundsen_Station_Coordinates.csv"), header=T)
ge_stn_1b <- read.csv(file=paste0(fpath_ge_stn,"GE_Amundsen_Station_Coordinates_leg1b.csv"), header=T)
ge_fdms <- read.csv(file=paste0(fpath_ge_stn,"GE_Amundsen_fdms_sicday.csv"), header=T)
sat1file <- open.nc(paste0(fpath_ge_stn,"GE_SatObs_2016176_2016180.nc"))
sat2file <- open.nc(paste0(fpath_ge_stn,"GE_SatObs_2016186_2016190.nc"))
sat1 <- read.nc(sat1file)
sat2 <- read.nc(sat2file)

# Remove some ice data in southwest corner which do not display nicely (they overlap land)
mask <- matrix(data = F, nrow = length(sat1$lon), ncol = length(sat1$lat))
latmat <- t(matrix(data = sat1$lat, ncol = length(sat1$lon), nrow = length(sat1$lat)))
lonmat  <- matrix(data = sat1$lon, nrow = length(sat1$lon), ncol = length(sat1$lat))
mask[(lonmat < -64.9 & latmat < 68.0) | (lonmat > -56 & latmat > 69.0)] <- T
tmp1 <- sat1$`sea-ice`
tmp1[mask] <- NA
sat1$`sea-ice` <- tmp1
tmp2 <- sat2$`sea-ice`
tmp2[mask] <- NA
sat2$`sea-ice` <- tmp2

# Scale fdms by sic
ge_fdms$fdms <- ge_fdms$fdms*(1-ge_fdms$sicday)

# Remove repeated rows
d1 <- duplicated(ge_fdms$dms)
d2 <- duplicated(ge_fdms$fdms)
ge_fdms <- ge_fdms[!(d1 & d2),]

# Define plotting region and plot coastline (no projection)
mlon = mean(ge_stn_all$lon)
mlat = mean(ge_stn_all$lat)
span = 565 # 400 to 600 works fine


# Exporting image? --------------
exportimg <- T

if (exportimg) {
  png(filename = '~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig1_map/Fig1_1panel.png',
      width = 15, height = 12, units = 'cm', pointsize = 8, # height 10 for 3 panels, 13 for 4 panels
      bg = 'white', res = 600, type = 'cairo')
}

## ---- Plot ----
plot(coastlineWorldFine, clon = mlon, clat = mlat, span = span, # span = c(2000,1000)
     col = 'tan', bg = rgb(red = 1, green = 1, blue = 1, alpha = 0),
     axes = F) # xlab="Longitude", ylab="Latitude")
axis(side = 1, at = c(-67,-53))
axis(side = 1, at = c(-65,-60,-55), labels = c(expression(paste('',65*degree,'W')),expression(paste('',60*degree,'W')),expression(paste('',55*degree,'W'))))
axis(side = 2, at = c(67,72))
axis(side = 2, at = seq(68,71,1), labels = c(expression(paste('',68*degree,'N')),expression(paste('',69*degree,'N')),expression(paste('',70*degree,'N')),expression(paste('',71*degree,'N'))))
box()

contour(sat1$lon,sat1$lat,sat1$`sea-ice`,
        levels = c(0.10,0.80, 0.85, 0.90),
        lwd = c(2,3,4,5),
        lty = c(3,3,3),
        drawlabels = F, add = TRUE, col = 'skyblue1') # turquoise2
contour(sat2$lon,sat2$lat,sat2$`sea-ice`,
        levels = c(0.10,0.80, 0.85, 0.90),
        lwd = c(2,3,4,5),
        lty = c(1,1,1),
        drawlabels = F, add = TRUE, col = 'skyblue3') #

# Plots station points
points(x=ge_stn_all$lon, y=ge_stn_all$lat, col='red', cex=0.4, pch=4) # unsampled
points(x=ge_fdms$lon, y=ge_fdms$lat, col='black', cex=0.5, pch=4) # sampled
# Plots station fdms
ff <- 1.7
col.fdms <- rgb(red = .6, green = .6, blue = .6, alpha = .7)
# col.fdms <- rgb(red = 1, green = .8, blue = .5, alpha = .7)
# col.fdms <- rgb(red = 1, green = 193/255, blue = 7/255, alpha = .5)
# col.dms <- rgb(red = .8, green = .6, blue = .7, alpha = .5)
col.dms <- rgb(red = 0, green = 0, blue = 0, alpha = 1)
# col.dms <- rgb(red = 216/255, green = 27/255, blue = 96/255, alpha = .5)
points(x=ge_fdms$lon, y=ge_fdms$lat, col=col.fdms, cex=ff*sqrt(ge_fdms$fdms), pch=19)
points(x=ge_fdms$lon, y=ge_fdms$lat, col=col.dms, cex=ff*sqrt(ge_fdms$dms), pch=1, lwd=3)

# Add legends -----------------

# Ice contours
legend("bottomleft", seg.len = 3, cex = 1,
       # legend(x = -66, y = 69, seg.len = 3, cex = 1,
       lwd = rep(c(2,3,4,5),2),
       lty = c(3,3,3,3,1,1,1,1),
       legend = c("start 0.10","        0.80", "        0.85", "        0.90", "end  0.10", "        0.80", "        0.85", "        0.90"),
       col = c(rep('skyblue1',4),rep('skyblue3',4)),
       title = "Sea ice concentration",
       bg= "gray97", box.col = "gray95")

# # FDMS and DMS overlaid
# legend("topright", pch = 19,
#        pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
#        legend = c("    0.2","    2","","    20","      "),
#        col = rep(col.fdms,5),
#        title = "FDMS or DMS",
#        ncol = 1,
#        bg= "gray97", box.lwd = 0, box.col = "gray95")
# legend("topright", pch = 1, lty = 0, lwd = 3,
#        pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
#        legend = c("       ","     ","","      ","      "),
#        col = rep(col.dms,5),
#        title = "FDMS or DMS",
#        text.col = rgb(red = 1, green = 1, blue = 1, alpha = 0) ,
#        ncol = 1,
#        bg=rgb(red = 1, green = 1, blue = 1, alpha = 0), box.lwd = 0, box.col = "white")

# FDMS and DMS, separate
xpos <- -56.15
legend(x = xpos, y = 71.07, pch = 19,
       pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
       legend = c("    0.2","    2","","    20","      "),
       col = rep(col.fdms,5),
       title = expression(paste('FDMS (µmol ',m^-2,' ',d^-1,')')),
       ncol = 1,
       bg= "gray97", box.lwd = 0, box.col = "gray95")
legend(x = xpos, y = 70.27, pch = 1, lty = 0, lwd = 3,
       pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
       legend = c("    0.2","    2","","    20","      "),
       col = rep(col.dms,5),
       title = "DMS (nM)                   ",
       ncol = 1,
       bg= "gray97", box.lwd = 0, box.col = "gray95")

if (exportimg) {dev.off()}
