# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggplot2)

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

# Change units of vertical integrals
surf$idms_z60 <- surf$idms_z60 / 1000
surf$idmspt_z60 <- surf$idmspt_z60 / 1000

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

yvarS <- list(idmspt_z60 = expression(paste(sum(DMSPt[0-60]),' (mmol ',m^-2,')')),
              idms_z60 = expression(paste(sum(DMS[0-60]),' (mmol ',m^-2,')')),
              dms = expression('<DMS>'[0-5]*' (nM)'),
              fdms = expression('FDMS (µmol m'^-2*' d'^-1*')'))
surf$station <- as.character(surf$station)
xl <- "OWD (d)"
yl <- expression('DMSPt (µmol m'^-2*')')

for (yvar in names(yvarS)) {
  
  par(mar = c(3,5,1,1))
  toplot <- data.frame(station=surf$station, OWD_zone=surf$OWD_zone, Domain=surf$Domain, OWD=surf$OWD, yvar=surf[[yvar]])
  
  # Remove labels for selected y variables and conditions
  if (yvar %in% c("dms","fdms")) {
    toplot$station <- ifelse(toplot$yvar > quantile(toplot$yvar, 0.4, na.rm = T), as.character(toplot$station), "")
  } else if (yvar == "idms_z60") { # if (yvar == "idms_z60")
    toplot$station <- ifelse(toplot$yvar > quantile(toplot$yvar, 0.2, na.rm = T), as.character(toplot$station), "")
  }
  
  p <- ggplot(toplot, aes(x=OWD, y=yvar, shape=Domain, colour=OWD_zone))
  p + geom_point(size = 3) +
    geom_text(aes(label=station), hjust=-0.3, vjust=0.2, show.legend = F, size = 11/4, check_overlap = T, color = "gray60") + # Setting size to x/4 is to maintain proportion with default ggplot of 15/4
    scale_color_manual(values = col) +
    scale_shape_manual(values = c(16,17)) +
    ylab(yvarS[[yvar]]) +
    ylim(c(0, 1.05*max(toplot$yvar, na.rm = T))) + 
    theme_bw()

  if (exportimg) {
    ggsave(
      filename = paste0(yvar,".png"),
      plot = last_plot(),
      device = NULL,
      path = opath,
      scale = 2,
      width = 6,
      height = 4,
      units = "cm",
      dpi = 600)
  }
}
