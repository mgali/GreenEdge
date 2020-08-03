# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library("readxl")
library(RColorBrewer)
library(dplyr)

# Load data
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
prof.all <- read.csv(file = "~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_GreenEdge_ANP_v0.1.csv", header = T)
surf.all <- read.csv(file = "~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_Ice-edge-paper_per-station_v0.1.csv", header = T)
pre.qpcr <- read_excel(path = "~/Desktop/GreenEdge/Results_methylotrophy story_ODV2-1_mg.xlsx", sheet = "Results_methylotrophy story_ODV")

# Exporting image?
exportimg <- T # Figure with binned profiles by OWD categories
exportfig <- T # Figure with profiles in T5 where all data available
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_qPCR_norm_grouped/"

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

# ---------------------
# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'Station', suffixes = "")

# Preprocess qpcr data
qpcr <- data.frame(Station = pre.qpcr$Station,
                   cast = pre.qpcr$CTD,
                   Depth_m = pre.qpcr$Depth,
                   dms = pre.qpcr$dms,
                   dmspt = pre.qpcr$dmspt,
                   dddP = pre.qpcr$`Average of dddP copies/mL`,
                   dmdA = pre.qpcr$`Average of dmdA copies /mL`,
                   mxaF = pre.qpcr$`Average of mxaF copies/mL`,
                   xoxF = pre.qpcr$`Average of xoxF copies/mL`,
                   Methylobacterium = pre.qpcr$Methylobacterium*100,
                   Rhodobacteraceae = pre.qpcr$Rhodobacteraceae*100,
                   SAR11  = pre.qpcr$SAR11*100,
                   SAR116 = pre.qpcr$SAR116*100,
                   Polaribacter  = pre.qpcr$Polaribacter*100,
                   Thiotrichales  = pre.qpcr$Thiotrichales*100,
                   Methylophylales = pre.qpcr$Methylophylales*100,
                   Oceanospirillales = pre.qpcr$Oceanospirillales*100,
                   Pseudoalteromonas = pre.qpcr$Pseudoalteromonas*100,
                   BA = pre.qpcr$`BA cells/mL`,
                   BP = pre.qpcr$`Bacterial production µgC/L/d`
)
qpcr$Station <- as.numeric(gsub(pattern = "G", replacement = "", qpcr$Station))

# Merge with qpcr
pplot <- merge(x = qpcr, y = pplot, all.x = T, all.y = F, by = c('Station','Depth_m'), suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Remove duplicated rows 
dd <- duplicated(pplot[,c("cast","Depth_m")])
pplot <- pplot[!dd,]

# # Hide data from leg 1b
# pplot[pplot$Station>=400,] <- NA

# Add ratios
pplot$dddP2dmdA <- pplot$dddP/pplot$dmdA                           # dddP/dmdA ratio
pplot$mxaF2dddP <- pplot$mxaF/pplot$dddP                           # mxaF/dddP ratio (methylotrophy)

# ---------------------
# Bin profiles by station categories
df2bin <- pplot
z_class <- cut(df2bin$Depth_m, breaks = c(0,9,21,41,81), labels = c(0,1,2,3))
st_class <- list(owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW")))

# View(pplot[,c("Station","Depth_m","OWD","ANP","dms","dmspt","dddP","dmdA","BA","BP")])

# ---------------------
# Loop on different station classifications

for (sc in "owd_class") { #names(st_class)
  
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
  
  # --------------------------------------------------
  # View some tables, compute some summary stats
  # View(pplot.bin$count[,c("Station","ANP","dms","dmspt","dddP","dmdA","BA","BP")])
}

# --------------------------------------------------
# # Filter out statistics with less than 2 meas.
# pplot.bin$mean[pplot.bin$count<2] <- NA
# pplot.bin$median[pplot.bin$count<2] <- NA

# --------------------------------------------------
# FINAL PLOTTING FOR PAPER

# ---------------------
# Figure with concentrations at all stations

xvarS <- list(ANP = "ANP (-)",
              BA = expression('Bact. abund. (cells mL'^-1*')'),
              BP = expression('Bact. prod. (µg C L'^-1*' d'^-1*')'),
              dddP = expression('dddP (copies mL'^-1*')'),
              dmdA = expression('dmdA (copies mL'^-1*')'),
              dddP2dmdA = "dddP:dmdA",
              mxaF = expression('mxaF (copies mL'^-1*')'),
              xoxF = expression('xoxF (copies mL'^-1*')')
)
yvar <- "Depth_m"
lett <- plet
names(lett) <- names(xvarS)

for (sc in "owd_class") {
  
  if (exportimg) {
    
    png(filename = paste0(opath,"FigS5_qpcr_meth_",sc,".png"), width = 15, height = 9, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')
    
    # Multipanel setup
    m0 <- matrix(data = 0, nrow = 4, ncol = 4)
    mr1 <-  cbind(m0+1,m0+2,m0+3,m0+4)
    m <- rbind(mr1,mr1+4)
    layout(m)
    par(oma = c(1,1,0.5,0.5))
    
    for (xvar in names(xvarS)) {
      
      print(xvar)
      par(mar = c(4,3,2,0.5))
      
      xl <- c(min(c(0,1.1*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))),
              1.1*max(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))
      if (xvar  %in% c("sal","sigt")) {xl[1] <- 0.9*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T)}
      if (xvar  == "ANP") {xl <- rev(xl)}
      print(xl)
      
      plot(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
           y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
           ylim = c(70,0), xlim = xl, pch = 19, col = col[1], cex = 1.9, axes = F, xlab = "", ylab = "")
      box()
      axis(side = 1, cex.axis = 1.1)
      mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
      if (xvar %in% c("ANP","dddP")) {
        axis(side = 2, cex.axis = 1.1)
        mtext(side = 2, "Depth (m)", cex = 0.9, line = 2.5)
      } else {
        axis(side = 2, cex.axis = 1.1, at = seq(0,70,10), labels = rep("",8))
      }
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
      text(x = xl[1]+0.95*(xl[2]-xl[1]), y = 6, labels = lett[xvar], cex = 1.3)
      if (xvar  == "dddP") {
        legend(x = xl[1]+0.4*(xl[2]-xl[1]),
               y = 25,
               pch = 19,
               legend = c("ICE","MIZ","OW"),
               cex = 1.3,
               col = col,
               bg= "white", box.lwd = 0)
        legend(x = xl[1]+0.4*(xl[2]-xl[1]),
               y = 50,
               lty = c(1,3),
               legend = c("medians","means"),
               cex = 1.3,
               col = "black",
               bg= "white", box.lwd = 0)
      }
      
    } # end loop on variables
    dev.off()
  }
}


# ---------------------
# Figure with T5 profiles

xvarS <- list(ANP = "ANP (-)",
              BA = expression('Bact. abund. (cells mL'^-1*')'),
              BP = expression('Bact. prod. (µg C L'^-1*' d'^-1*')'),
              dmspt = "DMSPt (nM)",
              dms = "DMS (nM)",
              dddP = expression('dddP (copies mL'^-1*')'),
              dmdA = expression('dmdA (copies mL'^-1*')'),
              dddP2dmdA = "dddP:dmdA",
              mxaF = expression('mxaF (copies mL'^-1*')'),
              xoxF = expression('xoxF (copies mL'^-1*')')
)
lett <- plet
names(lett) <- names(xvarS)

if (exportfig) {
  
  png(filename = paste0(opath,"FigS6_qpcr_meth_stn_507_519.png"), width = 17, height = 9, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')
  
  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 4)
  mr1 <-  cbind(m0+1,m0+2,m0+3,m0+4,m0+5)
  m <- rbind(mr1,mr1+5)
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  for (xvar in names(xvarS)) {
    
    i507 <- which(pplot$Station==507 & !is.na(pplot[,xvar]))
    i519 <- which(pplot$Station==519 & !is.na(pplot[,xvar]))
    xl <- range(pplot[c(i507,i519),xvar], na.rm = T)
    xl <- c(xl[1]*0.9, xl[2]*1.1)
    if (xvar == "ANP") {xl <- rev(xl)}
    
    print(xvar)
    par(mar = c(4,3,2,0.5))
    plot(pplot[i507,xvar], pplot[i507,"Depth_m"],
         ylim = c(50,0), xlim = xl,
         cex = 0.6, col = col[3],
         type = "l", lwd = 1,
         axes = F, xlab = "", ylab = "")
    lines(pplot[i519,xvar], pplot[i519,"Depth_m"],
          cex = 0.6, col = col[2],
          type = "l", lwd = 1)
    points(pplot[i507,xvar], pplot[i507,"Depth_m"],
           cex = 1, col = col[3], pch = 19)
    points(pplot[i519,xvar], pplot[i519,"Depth_m"],
           cex = 1, col = col[2], pch = 15)
    box()
    axis(side = 1, cex.axis = 1.1)
    mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
    text(x = xl[1]+0.95*(xl[2]-xl[1]), y = 2, labels = lett[xvar], cex = 1.3)
    if (xvar %in% c("ANP","BP")) {
      axis(side = 2, cex.axis = 1.1)
      mtext(side = 2, "Depth (m)", cex = 0.9, line = 2.5)
    } else {
      axis(side = 2, cex.axis = 1.1, at = seq(0,50,10), labels = rep("",6))
    }
    if (xvar == "dmspt") {
      legend(x = xl[1]+0.4*(xl[2]-xl[1]),
             y = 25,
             pch = c(19,15),
             legend = c("st. 507","st. 519"),
             cex = 1.1,
             col = col[c(3,2)],
             bg= "gray95", box.lwd = 0)
    }
  }
  dev.off()
}


# ---------------------------------------
# Correlation analysis for dddP and dmdA
# ra <- cor.test(pplot$dmdA, pplot$dmspt, use = "pairwise", method = "sp")
# print(ra)
# rb <- cor.test(pplot$dddP, pplot$dmspt, use = "pairwise", method = "sp")
# print(rb)
# rc <- cor.test(pplot$dmdA, pplot$dms, use = "pairwise", method = "sp")
# print(rc)
# rd <- cor.test(pplot$dddP, pplot$dms, use = "pairwise", method = "sp")
# print(rd)

# xvars1 <- c("ANP","BA","BP")
# yvarS <- list(dddP = "dddP copies/mL",
#               dmdA = "dmdA copies/mL")
# Xc1 <-
# Xc2 <-
# Xc <- cbind(Xc1,Xc2)

# ---------------------------------------
# Correlation analysis for mxaF, xoxF and DMS
re <- cor.test(pplot$mxaF, pplot$dms, use = "pairwise", method = "sp")
print(re)
rf <- cor.test(pplot$xoxF, pplot$dms, use = "pairwise", method = "sp")
print(rf)
rg <- cor.test(pplot$mxaF+pplot$xoxF, pplot$dms, use = "pairwise", method = "sp")
print(rg)

rp1 <- cor.test(pplot$mxaF, pplot$dddP, use = "pairwise", method = "sp")
print(rp1)
rp2 <- cor.test(pplot$mxaF+pplot$xoxF, pplot$dddP, use = "pairwise", method = "sp")
print(rp2)
rd <- cor.test(pplot$mxaF+pplot$xoxF, pplot$dmdA, use = "pairwise", method = "sp")
print(rd)

