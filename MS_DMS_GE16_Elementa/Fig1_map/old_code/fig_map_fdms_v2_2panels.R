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
sat2file <- open.nc(paste0(fpath_ge_stn,"GE_SatObs_2016181_2016185.nc"))
sat3file <- open.nc(paste0(fpath_ge_stn,"GE_SatObs_2016186_2016190.nc"))
sat1 <- read.nc(sat1file)
sat2 <- read.nc(sat2file)
sat3 <- read.nc(sat3file)
gebcofile <- open.nc("~/Desktop/GreenEdge/GEBCO_2019_01_Nov_2019_89f1161f595f/gebco_2019_n76.0_s59.0_w-80.0_e-45.0.nc")
gebco <- read.nc(gebcofile)
gebco$elevation[gebco$elevation > -50] <- NA

# # Correct error in cordinates: already corrected in csv file
# ge_fdms$lon[ge_fdms$stn==507] <- -59.123

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
tmp3 <- sat3$`sea-ice`
tmp3[mask] <- NA
sat3$`sea-ice` <- tmp3

# Scale fdms by sic
ge_fdms$fdms <- ge_fdms$fdms*(1-ge_fdms$sicday)

# # Mean SST
# sst <- array(c(sat1$`sst`, sat2$`sst`, sat3$`sst`), dim = c(dim(sat3$`sst`)[1],dim(sat3$`sst`)[2],3))
# sstmean <- apply(sst, c(1,2), mean, na.rm = T)

# Remove repeated rows
d1 <- duplicated(ge_fdms$dms)
d2 <- duplicated(ge_fdms$fdms)
ge_fdms <- ge_fdms[!(d1 & d2),]

# Define plotting region and plot coastline (no projection)
mlon = mean(ge_stn_all$lon)
mlat = mean(ge_stn_all$lat)
span = 565 # 400 to 600 works fine

# Colors
col.land <- 'gray95' # 'tan'
col.fdms <- rgb(.6, .6, .6, alpha = .7)
col.dms <- rgb(0, 0, 0, alpha = 1)
col.ice1 <- 'skyblue'
col.ice2 <- 'skyblue2'
col.ice3 <- 'skyblue3'

# Exporting image? --------------
exportimg <- T

if (exportimg) {
  png(filename = '~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig1_map/Fig1_2panels_v3.png',
      width = 15, height = 7, units = 'cm', pointsize = 8, # height 10 for 3 panels, 13 for 4 panels
      bg = 'white', res = 600, type = 'cairo')
}

# Layout
m <- cbind(matrix(data = 1, nrow = 7, ncol = 6), matrix(data = 2, nrow = 7, ncol = 9))
layout(m)

## ---- Plot ----

# Panel 1
# ===========================================================
plot(coastlineWorldFine, clon = mlon-1, clat = mlat, span = c(2950,1850), #c(2000,2000)
     col = col.land, bg = rgb(red = 1, green = 1, blue = 1, alpha = 0), lwd = 0.5,
     axes = F) # xlab="Longitude", ylab="Latitude")
# image(x = sat1$lon, y = sat1$lat, z = sstmean, add = T, col = rev(heat.colors(100)))

# Bathy
image(gebco$lon,gebco$lat,gebco$elevation, add = T, col = oce.colorsGebco(11))
contour(gebco$lon,gebco$lat,gebco$elevation,
        levels = c(-1000),
        lwd = 0.5, lty = 1, cex = 0.4,
        drawlabels = T, add = TRUE, col = "#1E8EB3")

# Stations and close-up rectangle
points(x=ge_stn_all$lon, y=ge_stn_all$lat, col='red', cex=0.3, pch=20) # unsampled
points(x=ge_fdms$lon, y=ge_fdms$lat, col='black', cex=0.4, pch=20) # sampled
rect(-66.6, 67.4, -53.5, 71.1, col = NULL, border = 'black', lty = 1, lwd = 1)
       
# Annotations
text(x = -77, y = 74.5, labels = 'a)', offset = 0.5, cex = 1.3, srt = 0, col = 'black')
text(x = -74, y = 70, labels = 'Baffin Island', offset = 0.5, cex = 1, srt = -45, col = 'black')
text(x = -50, y = 73, labels = 'Greenland', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -56.7, y = 61, labels = 'Labrador Sea', offset = 0.5, cex = 1, srt = 0, col = 'white')
text(x = -57.5, y = 66, labels = 'Davis Strait', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -67, y = 74, labels = 'Baffin Bay', offset = 0.5, cex = 1, srt = 0, col = 'black')

axis(side = 1, at = c(-80,-40))
axis(side = 1, at = c(-70,-60,-50), labels = c(expression(paste('',70*degree,'W')),expression(paste('',60*degree,'W')),expression(paste('',50*degree,'W'))))
axis(side = 2, at = c(55,85))
axis(side = 2, at = seq(60,75,5), labels = c(expression(paste('',60*degree,'N')),expression(paste('',65*degree,'N')),expression(paste('',70*degree,'N')),expression(paste('',75*degree,'N'))))
box()

# Panel 2
# ===========================================================
plot(coastlineWorldFine, clon = mlon, clat = mlat, span = 565,
     col = col.land, bg = rgb(red = 1, green = 1, blue = 1, alpha = 0), lwd = 0.5,
     axes = F) # xlab="Longitude", ylab="Latitude")
# image(x = sat3$lon, y = sat3$lat, z = sat3$`sea-ice`, add = T, col = oce.colors(10))

axis(side = 1, at = c(-67,-53))
axis(side = 1, at = c(-65,-60,-55), labels = c(expression(paste('',65*degree,'W')),expression(paste('',60*degree,'W')),expression(paste('',55*degree,'W'))))
axis(side = 2, at = c(67,72))
axis(side = 2, at = seq(68,71,1), labels = c(expression(paste('',68*degree,'N')),expression(paste('',69*degree,'N')),expression(paste('',70*degree,'N')),expression(paste('',71*degree,'N'))))
box()
# arrows(x0 = -63, y0 = 71, x1 = -65, y1 = 70.8, col = col.ice2, lwd = 1)
text(x = -66.2, y = 70.9, labels = 'b)', offset = 0.5, cex = 1.3, srt = 0, col = 'black')

contour(sat1$lon,sat1$lat,sat1$`sea-ice`,
        levels = c(0.10),
        lwd = c(1,2),
        lty = 3,
        drawlabels = F, add = TRUE, col = col.ice1)
contour(sat2$lon,sat2$lat,sat2$`sea-ice`,
        levels = c(0.10),
        lwd = c(1,2),
        lty = 2,
        drawlabels = F, add = TRUE, col = col.ice2)
contour(sat3$lon,sat3$lat,sat3$`sea-ice`,
        levels = c(0.10,0.50,0.80,0.85),
        lwd = c(1,2,3,4),
        lty = 1,
        drawlabels = F, add = TRUE, col = col.ice3)

# Plots station points
points(x=ge_stn_all$lon, y=ge_stn_all$lat, col='red', cex=0.4, pch=4) # unsampled
points(x=ge_fdms$lon, y=ge_fdms$lat, col='black', cex=0.5, pch=4) # sampled

# Plots station fdms
ff <- 1.7
points(x=ge_fdms$lon, y=ge_fdms$lat, col=col.fdms, cex=ff*sqrt(ge_fdms$fdms), pch=19)
points(x=ge_fdms$lon, y=ge_fdms$lat, col=col.dms, cex=ff*sqrt(ge_fdms$dms), pch=1, lwd=1)

# Add legends -----------------

# Ice contours
legend("bottomleft", seg.len = 2, cex = 1,
       # legend(x = -66, y = 69, seg.len = 3, cex = 1,
       lwd = c(1,1,1,2,3,4),
       lty = c(3,2,1,1,1,1),
       legend = c("d176, 0.10","d181, 0.10","d186, 0.10", "          0.50", "          0.80", "          0.85"),
       col = c(col.ice1,col.ice2,rep(col.ice3,4)),
       title = "Sea ice concentration",
       bg= "gray97", box.col = "gray95")

# FDMS and DMS, separate
xpos <- -56.25
legend(x = xpos, y = 71.00, pch = 19,
       pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
       legend = c("    0.2","    2","","    20","      "),
       col = rep(col.fdms,5),
       title = expression(paste('FDMS (µmol ',m^-2,' ',d^-1,')')),
       ncol = 1,
       bg= "gray97", box.lwd = 0, box.col = "gray95")
legend(x = xpos, y = 70.10, pch = 1, lty = 0, lwd = 1,
       pt.cex = ff*sqrt(c(0.2,2,NA,20,NA)),
       legend = c("    0.2","    2","","    20","      "),
       col = rep(col.dms,5),
       title = "DMS (nM)                   ",
       ncol = 1,
       bg= "gray97", box.lwd = 0, box.col = "gray95")

# Labels for transects and dates
text(x = -62.8, y = 68.40, labels = 'T4', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -62.7, y = 68.27, labels = 'd176 ->', offset = 0.5, cex = 1, srt = 0, col = 'red')
text(x = -56.4, y = 68.13, labels = '-> d181', offset = 0.5, cex = 1, srt = 0, col = 'red')

text(x = -57, y = 70.23, labels = 'T5', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -57, y = 70.1, labels = '<- d181', offset = 0.5, cex = 1, srt = 0, col = 'red')
text(x = -63.7, y = 70.0, labels = 'd184 <-', offset = 0.5, cex = 1, srt = 0, col = 'red')

text(x = -65, y = 70.63, labels = 'T6', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -64.9, y = 70.5, labels = 'd185 ->', offset = 0.5, cex = 1, srt = 0, col = 'red')
text(x = -58.2, y = 70.5, labels = '-> d188', offset = 0.5, cex = 1, srt = 0, col = 'red')

text(x = -57, y = 69.63, labels = 'T7', offset = 0.5, cex = 1, srt = 0, col = 'black')
text(x = -57, y = 69.5, labels = '<- d189', offset = 0.5, cex = 1, srt = 0, col = 'red')
text(x = -64.1, y = 69.5, labels = 'd191 <-', offset = 0.5, cex = 1, srt = 0, col = 'red')

if (exportimg) {dev.off()}
