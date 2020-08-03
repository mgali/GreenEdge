# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggplot2)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
surf <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_OWD_ARCvsATL/"

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi

# ---------------------
# Replace NaN by NA
surf[is.nan(as.matrix(surf))] <- NA

# Remove data where no DMS or DMSPt are available
surf <- surf[!is.na(surf$dms) | !is.na(surf$dmspt),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# surf[surf$stn==519,c("dms","dmspt")] <- c(3.93,79.9)
surf[surf$stn==519,c("dms","dmspt")] <- c(11.42,79.9)

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

# Remove rows where surface sulfur variables are not available
surf <- surf[!is.na(surf$dms) & !is.na(surf$fdms),] # !is.na(surf$idms_z60) & !is.na(surf$idmspt_z60) & 

# Change units of vertical integrals
surf$idms_z60 <- surf$idms_z60 / 1000
surf$idmspt_z60 <- surf$idmspt_z60 / 1000


# --------------------------------------------------
# Calculate vertical integrals for putative Phaeocystis diagnostic pigments

# Exception for diagnostic pigments and ratios: replace NAs with zeros (= take into account, and plot, values below DL)
diapigs <- c("chlc3","but","hex","but19_like")
set2zero <- is.na(prof[,diapigs])
prof[,diapigs][set2zero] <- 0

# Repeat TChla integral from same station and different cast if not available for sulfur cast
for (cc in surf$stn) {
  tmp <- surf[surf$stn==cc,"iTchla_z60"]
  tmp[is.na(tmp)] <- mean(tmp, na.rm = T)
  surf[surf$stn==cc,"iTchla_z60"] <- tmp
}
# View(surf[,c("stn","cast","iTchla_z60")]) # debug

for (vv in diapigs) {
  newVar <- paste0("i",vv,"_z60")
  surf[,newVar] <- NA
  for (cc in surf$stn) {
    xy <- prof[prof$stn==cc, c("depth",vv)]
    xy <- xy[!duplicated(xy[,"depth"]),]
    # if (cc == 512 & vv == "but19_like") { # debug
    #   print(xy)
    #   View(prof[,c("stn","cast","depth",vv)])
    # }
    
    if ( sum(!is.na(xy[ xy$depth<=60 , vv ])) >= 4 ) {
      xy <- rbind(xy,c(100,0))
      xyout <- approx(x = xy$depth, y = xy[,vv], xout = seq(0.5,99.5,1),
                     method = "linear", rule = 2, ties = mean)
      surf[surf$stn==cc, newVar] <- sum(xyout$y[1:60], na.rm = T)
    }
    # if (cc == 512 & vv == "but19_like") { # debug
    #   print(sum(xyout$y[1:60], na.rm = T))
    # }
  }
  # Add ratios: Phaeocystis proxies?
  surf[,paste0("i",vv,"_2_tchla_z60")] <- surf[,newVar] / surf$iTchla_z60
}

# --------------------------------------------------
# FINAL PLOTTING FOR PAPER
# NOTE: I tried to put 4 variables into a single column and converting variable names to factors to be able to
# use facet_wrap to make the 4x4 plot matrix. It didn't work because all variables share the same y axis limits
# toplot <- list()
# for (svar in names(svarS)) {
#   toplot[[svar]] <- data.frame(station=surf$station, OWD_zone=surf$OWD_zone, Domain=surf$Domain, OWD=surf$OWD, yvar=surf[[svar]], groupvar=svar)
# }
# toplot <- data.table::rbindlist(toplot, use.names = F, fill = F) # Nearly equivalent: toplot <- do.call("rbind", toplot)
# ...facet_wrap(vars(groupvar), labeller = labeller(yvariable = svarS))

yvarS <- list(icp_z60 = expression(paste(sum(Cp[0-60]),' (-)')),
              iTchla_z60 = expression(paste(sum(TChl_a[0-60]),' (mg ',m^-2,')')),
              idmspt_z60 = expression(paste(sum(DMSPt[0-60]),' (mmol ',m^-2,')')),
              idms_z60 = expression(paste(sum(DMS[0-60]),' (mmol ',m^-2,')')),
              dms = expression('<DMS>'[0-5]*' (nM)'),
              fdms = expression('FDMS (µmol m'^-2*' d'^-1*')'),
              # ichlc3_z60 = expression(paste(sum(Chlc3[0-60]),' (mg ',m^-2,')')),
              # ibut_z60 = expression(paste(sum(But[0-60]),' (mg ',m^-2,')')),
              # ihex_z60 = expression(paste(sum(Hex[0-60]),' (mg ',m^-2,')')),
              # ibut19_like_z60 = expression(paste(sum(But-like[0-60]),' (mg ',m^-2,')')),
              # ichlc3_2_tchla_z60 = expression(paste(sum(Chl_c3/TChl_a[0-60]),' (-)')),
              # ibut_2_tchla_z60 = expression(paste(sum(But/TChl_a[0-60]),' (-)')),
              # ihex_2_tchla_z60 = expression(paste(sum(Hex/TChl_a[0-60]),' (-)')),
              ibut19_like_2_tchla_z60 = expression(paste(sum(But-like/TChl_a[0-60]),' (-)')))
surf$station <- as.character(surf$station)

for (yvar in names(yvarS)) {
  
  par(mar = c(3,5,1,1))
  toplot <- data.frame(station=surf$station, OWD_zone=surf$OWD_zone, Domain=surf$Domain, OWD=surf$OWD, yvar=surf[[yvar]])
  
  # Remove duplicated rows
  toplot <- toplot[!duplicated(toplot$station) & !duplicated(toplot$yvar)  & !is.na(toplot$yvar),]
  
  # Remove labels for selected y variables and conditions
  if (yvar %in% c("icp_z60","iTchla_z60","idmspt_z60","idms_z60","dms","fdms")) {
    toplot$station <- ifelse(toplot$yvar > quantile(toplot$yvar, 0.35, na.rm = T), # Max labels is quantile 0.2 just for FDMS
                             as.character(toplot$station), "")
  }
  
  p <- ggplot(toplot, aes(x=OWD, y=yvar, shape=Domain, colour=OWD_zone))
  p + geom_point(size = 3) +
    geom_text(aes(label=station), hjust=-0.3, vjust=0.2, show.legend = F, size = 11/4, check_overlap = T, color = "gray60") + # Setting size to x/4 is to maintain proportion with default ggplot of 15/4
    scale_color_manual(values = col) +
    scale_shape_manual(values = c(16,17)) +
    xlim(c(-23,37)) +
    xlab("Open water days") +
    ylab(yvarS[[yvar]]) +
    ylim(c(0, 1.05*max(toplot$yvar, na.rm = T))) + 
    theme_bw()
  
  if (exportimg) {
    ggsave(
      filename = paste0(yvar,".png"),
      plot = last_plot(),
      device = NULL,
      path = opath,
      scale = 1.6,
      width = 6.5,
      height = 4,
      units = "cm",
      dpi = plotres
    )
  }
}
