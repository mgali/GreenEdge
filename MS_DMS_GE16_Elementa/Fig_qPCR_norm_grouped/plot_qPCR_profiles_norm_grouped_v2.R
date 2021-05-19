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
                   BP = pre.qpcr$`Bacterial production µgC/L/d`,
                   BP_per_cell = pre.qpcr$`BP per cells`
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
pplot$dms2dmspt <- pplot$dms/pplot$dmspt                           # DMS/DMSPt ratio
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

xvarS <- list(
  ANP = "ANP (-)",
  BA = expression('Bact. abund. (cells mL'^-1*')'),
  BP = expression('Bact. prod. (µg C L'^-1*' d'^-1*')'),
  # BP_per_cell = expression('BP/cell (fg C cell'^-1*' d'^-1*')'),
  # dddP = expression('dddP (copies mL'^-1*')'),
  dddP = expression(paste(italic("dddP")," (copies ",mL^-1*")")),
  dmdA = expression(paste(italic("dmdA")," (copies ",mL^-1*")")),
  dddP2dmdA = expression(paste(italic("dddP"),":",italic(dmdA)))
  # mxaF = expression('mxaF (copies mL'^-1*')'),
  # mxaF2dddP = "mxaF:dddP"
)
yvar <- "Depth_m"
lett <- plet
names(lett) <- names(xvarS)

for (sc in "owd_class") {
  
  if (exportimg) {
    
    png(filename = paste0(opath,"FigS5_qpcr_v2_",sc,".png"), width = 12, height = 9, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'quartz')
    
    # Multipanel setup
    m0 <- matrix(data = 0, nrow = 4, ncol = 4)
    mr1 <-  cbind(m0+1,m0+2,m0+3)
    m <- rbind(mr1,mr1+3)
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
      text(x = xl[1]+0.95*(xl[2]-xl[1]), y = 6, labels = lett[xvar], cex = 1.3)
      if (xvar  == "dddP") {
        legend(x = xl[1]+0.6*(xl[2]-xl[1]),
               y = 40,
               pch = 19,
               legend = c("ICE","MIZ","OW"),
               cex = 1.3,
               col = col,
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
              dmspt = expression("DMSP"[t]*" (nM)"),
              dms = "DMS (nM)",
              BP = expression('Bact. prod. (µg C L'^-1*' d'^-1*')'),
              dddP = expression(paste(italic("dddP")," (copies ",mL^-1*")")),
              dmdA = expression(paste(italic("dmdA")," (copies ",mL^-1*")")),
              dddP2dmdA = expression(paste(italic("dddP"),":",italic(dmdA)))
              # mxaF = expression('mxaF (copies mL'^-1*')'),
              # mxaF2dddP = "mxaF:dddP"
)
lett <- plet
names(lett) <- names(xvarS)

if (exportfig) {
  
  png(filename = paste0(opath,"FigS6_qpcr_v2_stn_507_519.png"), width = 14, height = 9, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'quartz')
  
  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 4)
  mr1 <-  cbind(m0+1,m0+2,m0+3,m0+4)
  m <- rbind(mr1,mr1+4)
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
         cex = 0.6, col = "black",
         type = "l", lwd = 1,
         axes = F, xlab = "", ylab = "")
    lines(pplot[i519,xvar], pplot[i519,"Depth_m"],
          cex = 0.6, col = "black",
          type = "l", lwd = 1)
    points(pplot[i507,xvar], pplot[i507,"Depth_m"],
           cex = 1.5, col = "black", pch = 16)
    points(pplot[i519,xvar], pplot[i519,"Depth_m"],
           cex = 1.5, col = "black", bg = "white", pch = 22)
    box()
    axis(side = 1, cex.axis = 1.1)
    mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
    text(x = xl[1]+0.05*(xl[2]-xl[1]), y = 2, labels = lett[xvar], cex = 1.3)
    if (xvar %in% c("ANP","BP")) {
      axis(side = 2, cex.axis = 1.1)
      mtext(side = 2, "Depth (m)", cex = 0.9, line = 2.5)
    } else {
      axis(side = 2, cex.axis = 1.1, at = seq(0,50,10), labels = rep("",6))
    }
    if (xvar == "BP") {
      legend("bottomright",
             pch = c(16,22),
             legend = c("st. 507","st. 519"),
             cex = 1.5,
             col = "black",
             bg= "gray95", box.lwd = 0)
    }
  }
  dev.off()
}



# ---------------------
# ---------------------
# Reduced figure with T5 profiles potentially for main text

xvarS <- list(dmspt = expression("DMSP"[t]*" (nM)"),
              dms = "DMS (nM)",
              # dms2dmspt = expression('DMS:DMSPt (mol mol'-1*')'),
              dddP = expression(paste(italic("dddP")," (copies ",mL^-1*")")),
              dddP2dmdA = expression(paste(italic("dddP"),":",italic(dmdA)))
)
lett <- plet
names(lett) <- names(xvarS)

if (exportfig) {
  
  png(filename = paste0(opath,"Fig9_qpcr_v2red_stn_507_519.png"), width = 8, height = 9, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'quartz')
  
  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 4)
  mr1 <-  cbind(m0+1,m0+2)
  m <- rbind(mr1,mr1+2)
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  for (xvar in names(xvarS)) {
    
    i507 <- which(pplot$Station==507 & !is.na(pplot[,xvar]))
    i519 <- which(pplot$Station==519 & !is.na(pplot[,xvar]))
    xl <- range(pplot[c(i507,i519),xvar], na.rm = T)
    xl <- c(xl[1]*0.9, xl[2]*1.1)
    if (xvar == "dms2dmspt") {xl <- c(0,0.5)}
    
    print(xvar)
    if (xvar %in% c("dmspt","dms")) {
      par(mar = c(1,2.7,4,0.1))
    } else {
      par(mar = c(4,2.7,1,0.1)) 
    }
    plot(pplot[i507,xvar], pplot[i507,"Depth_m"],
         ylim = c(50,0), xlim = xl,
         cex = 0.6, col = "black",
         type = "l", lwd = 1,
         axes = F, xlab = "", ylab = "")
    lines(pplot[i519,xvar], pplot[i519,"Depth_m"],
          cex = 0.6, col = "black",
          type = "l", lwd = 1)
    points(pplot[i507,xvar], pplot[i507,"Depth_m"],
           cex = 1.5, col = "black", pch = 16)
    points(pplot[i519,xvar], pplot[i519,"Depth_m"],
           cex = 1.5, col = "black", bg = "white", pch = 22)
    box()
    text(x = xl[1]+0.05*(xl[2]-xl[1]), y = 2, labels = lett[xvar], cex = 1.3)
    if (xvar %in% c("dmspt","dddP")) {
      axis(side = 2, cex.axis = 1.2)
      mtext(side = 2, "Depth (m)", cex = 0.9, line = 2.5)
    } else {
      axis(side = 4, cex.axis = 1.2, at = seq(0,50,10), labels = rep("",6))
    }
    if (xvar %in% c("dmspt","dms")) {
      axis(side = 3, cex.axis = 1.2)
      mtext(side = 3, xvarS[[xvar]], cex = 0.9, line = 3)
    } else {
      axis(side = 1, cex.axis = 1.2)
      mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
    }
    if (xvar == "dddP2dmdA") {abline(v = 1, lty = 2)}
    if (xvar == "dddP") {
      legend(x = xl[1]+0.2*(xl[2]-xl[1]),
             y = 35,
             pch = c(16,22),
             legend = c("st. 507 (+3 OWD)","st. 519 (-3 OWD)"),
             cex = 1.3,
             col = "black",
             bg= "white", box.lwd = 0)
    }
  }
  dev.off()
}




# ---------------------------------------
# Quick correlation analysis for dddP and dmdA

# xvars1 <- c("ANP","BA","BP")
# yvarS <- list(dddP = "dddP copies/mL",
#               dmdA = "dmdA copies/mL")
# ra <- cor.test(pplot$dmdA, pplot$dmspt, use = "pairwise", method = "sp")
# print(ra)
# rb <- cor.test(pplot$dddP, pplot$dmspt, use = "pairwise", method = "sp")
# print(rb)
# rc <- cor.test(pplot$dmdA, pplot$dms, use = "pairwise", method = "sp")
# print(rc)
# rd <- cor.test(pplot$dddP, pplot$dms, use = "pairwise", method = "sp")
# print(rd)

# # ---------------------------------------
# # Correlation analysis for mxaF and DMS
# rr <- cor.test(pplot$mxaF, pplot$dms, use = "pairwise", method = "sp")
# print(rr)
