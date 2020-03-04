# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
# prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- F
doexploreplot <- F
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_profiles_norm_grouped/"

showtable <- "owd_class" # either empty, owd_class or sic_class

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

# ---------------------
# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[(!is.na(prof.all$dms) | !is.na(prof.all$dmspt)) & !is.na(prof.all$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(3.93,79.9)
prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)

# Add MIZ classification
icecrit1 <- 0.15
icecrit2 <- 0.70
ICE <- surf.all[,c("SICm2d","SICm1d","SICday")]
icemin <- apply(ICE, 1, min, na.rm = T) # Min-max SIC criterion
icemax <- apply(ICE, 1, max, na.rm = T) # Min-max SIC criterion
icemean <- apply(ICE, 1, min, na.rm = T) # Mean concentration criterion
surf.all$sic_class <- NA
surf.all$sic_class[icemax<icecrit1] <- "OW"
surf.all$sic_class[icemin>icecrit2] <- "ICE"
surf.all$sic_class[icemax>=icecrit1 & icemin<=icecrit2] <- "MIZ"

# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Remove duplicated rows 
# dd <- (duplicated(pplot[,c("dmspt","cast","depth")]) | duplicated(pplot[,c("dms","cast","depth")])) & (!is.na(pplot$tchla) | !is.na(pplot$cpsmooth1)) # Does not work well, too many repeated DMSPt get in
# dd <- (duplicated(pplot[,c("dmspt","cast","depth")]) | duplicated(pplot[,c("dms","cast","depth")])) # Does not work well, too many repeated DMSPt get in
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) | is.na(pplot$dmspt)
pplot <- pplot[!dd,]
View(pplot[,c("stn","cast","depth","dms","dmspt","cpsmooth1","tchla")])

# f_mydiff <- function(x) {y <- as.logical(c(1,diff(x))); y[is.na(y)] <- T; return(y)}
# ddms <- f_mydiff(pplot$dms)
# ddepth <- f_mydiff(pplot$depth)
# dcast <- f_mydiff(pplot$cast)
# pplot <- pplot[ddms&!is.na(pplot$dmspt)&!is.na(pplot$tchla),]

# Hide data from stn 400 (either entire or just surface)
# pplot[pplot$stn<=400,] <- NA
pplot[pplot$stn<=400 & pplot$depth < 5,] <- NA

# Change units of N2 from s-2 to h-1
pplot$N2 <- sqrt(pplot$N2) * 3600

# Calculate photosynthetic and photoprotective carotenoids (Bricaud 2004)
pplot$psc <- rowSums(pplot[,c("fuco","peri","but19_like","hex","hex19_likeSUM")], na.rm = T)
pplot$ppc <- rowSums(pplot[,c("zea","anthera","viola","diadino","diato","allo","tcar")], na.rm = T)
pplot$dd <- rowSums(pplot[,c("diadino","diato")], na.rm = T)
pplot$vaz <- rowSums(pplot[,c("zea","anthera","viola")], na.rm = T)
pplot$tpig <- rowSums(pplot[,seq(59,88,1)], na.rm = T)

# Add ratios
pplot$cp2tchla <- pplot$cpsmooth1/pplot$tchla                           # Cp/tchla ratio
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
pplot$but19like_2_tchla <- pplot$but19_like/pplot$tchla                 # chlc3 to tchla (Phaeocystis proxy?)
pplot$chlc2_2_tchla <- pplot$chlc2group/pplot$tchla                     # chlc2 to tchla (chlc2 widespread pigment)
pplot$but_2_tchla <- pplot$but/pplot$tchla                              # 19-but to tchla (proxy of what? Phaeocystis? other hapto?)
pplot$fuco_2_tchla <- pplot$fuco/pplot$tchla                            # Fucoxanthin to tchla (fuco widespread pigment)
pplot$peri_2_tchla <- pplot$peri/pplot$tchla                            # Peridinin to tchla (peri in dinos)
pplot$chlc3_2_psc <- pplot$chlc3/pplot$psc                              # chlc3 to PSC (Phaeocystis proxy?)
pplot$phaeo2chl <- pplot$phaeo_Tu_ugL/pplot$chla_Tu_ugL                 # Phaeopigments to Chl (Turner)
# pplot$phdaSUM2tchla <- pplot$phbda/pplot$tchla                            # Phaeophorbide a to TChl (HPLC)
pplot$phdaSUM2tchla <- pplot$phdaSUM/pplot$tchla                        # Phaeophorbide a to TChl (HPLC)

# Remove phaeopigments outlier
pplot[pplot$phaeo2chl > 3 & !is.na(pplot$phaeo2chl),c("phaeo_Tu_ugL","chla_Tu_ugL","phaeo2chl")] <- NA

# Fill NA with zeros in pigments frequently below LOD and only if CHl is not NA
pplot[!is.na(pplot$tchla),c("chlc3","but19_like","peri")][is.na(pplot[!is.na(pplot$tchla),c("chlc3","but19_like","peri")])] <- 0

# # Testing some relationships
# plot(pplot$chlc3_2_tchla, pplot$dms2dmspt)
# plot(pplot$chlc3_2_tchla, pplot$dmspt2cp)
# plot(pplot$chlc2_2_tchla, pplot$dmspt2cp)

# ---------------------
# Bin profiles by station categories
df2bin <- pplot
z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3))
st_class <- list(sic_class = pplot$sic_class,
                 # owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,4.5,35), labels = c("ICE","MIZ","OW")))
                 # owd_class = cut(df2bin$OWD, breaks = c(-35,-2.5,5.5,35), labels = c("ICE","MIZ","OW")))
                 # owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW"))) # STANDARD?
                 owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW")))

# ---------------------
# Plot settings. EXPLORATORY
xvarS <- list(dms = "DMS (nM)",
              dmspt = "DMSPt (nM)",
              tchla = "TChla (µg/L)",
              chlc3 = "Chlc3 (µg/L)",
              but19_like = "19-But-like (µg/L)",
              but = "19-But (µg/L)",
              fuco = "Fucoxanthin (µg/L)",
              peri = "Peridinin (µg/L)",
              phbdaSUM = "Phaeophorb_a (µg/L)",
              cpsmooth1 = "Cp (1/m)",
              cp2tchla = "Cp/TChla (m2/mg)",
              dms2dmspt = "DMS/DMSPt",
              dmspt2tchla = "DMSPt/TChla",
              dmspt2cp = "DMSPt/Cp",
              temp = "Temperature (C)",
              sal = "Salinity",
              sigt = "sigma-t (kg/m3)",
              anp = "ANP",
              N2 = "Brunt-Väisälä freq. (1/h)",
              zN2max03 = "N2_max depth (below MLd0.03) (m)",
              zN2max125 = "N2_max depth (below MLd0.125) (m)",
              par_d_p24h_ein_m_2_day_1 = "PAR (µE/m2/d)",
              ppc2psc = "PPC/PSC",
              npp = "PPC/TPig",
              ppc2tchla = "PPC/TChla",
              psc2tchla = "PSC/TChla",
              chlc3_2_tchla = "Chlc3/TChla",
              chlc2_2_tchla = "Chlc2/TChla",
              but19like_2_tchla = "19-But-like/TChla",
              chlc3_2_psc = "Chlc3/PSC",
              but_2_tchla = "19-But/TChla",
              fuco_2_tchla = "Fucoxanthin/TChla",
              peri_2_tchla = "Peridinin/TChla",
              phaeo2chl = "Phaeopigments/Chla (Turner)",
              phdaSUM2tchla = "Phaeophorb_a/TChla (HPLC)",
              dd = "(Dd+Dt)/TChla",
              vaz = "(Vi+Anth+Zea)/TChla")
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
yvar <- "depth"

# ---------------------
# Loop on different station classifications. EXPLORATORY PLOTS

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
  
  for (xvar in names(xvarS)) {
    
    if (doexploreplot) {
      if (exportimg) {png(filename = paste0(opath,paste(sc,xvar,sep = "_"),".png"), width = 6, height = 6, units = 'cm', pointsize = 6, bg = 'white', res = plotres, type = 'cairo')}
      # if (exportimg) {png(filename = paste0(opath,paste(sc,xvar,"m2to3owd",sep = "_"),".png"), width = 6, height = 6, units = 'cm', pointsize = 6, bg = 'white', res = plotres, type = 'cairo')}
      
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
  
  if (showtable == sc) {
    # --------------------------------------------------
    # View some tables, compute some summary stats
    
    # View(pplot[,c("stn","OWD","dms","dmspt")])
    View(pplot.bin$count[,c("stn","OWD","dms","dmspt","cpsmooth1","tchla","N2","anp")])
    # View(pplot.bin$mean[,grep("SIC",names(pplot.bin$mean))]) # equivalent to: View(pplot.bin$mean[,c("SIC_CLASS","SICm2d","SICm1d","SICday")])
    
    # a <- as.matrix(pplot.bin$mean[,c("SICm2d","SICm1d","SICday")])
    # print(mean(a[c(1,2),1]))
    # print(mean(a[c(5,6),1]))
    # b <- as.matrix(pplot.bin$mean[,"OWD"])
    # print(mean(b[c(1,2),1]))
    # print(mean(b[c(5,6),1]))
    
    dmean <- pplot.bin$mean[pplot.bin$mean$Z_CLASS==1,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm","anp")]
    dmin <- pplot.bin$min[pplot.bin$min$Z_CLASS==1,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm","anp")]
    dmax <- pplot.bin$max[pplot.bin$max$Z_CLASS==1,c("SIC_CLASS","mld03","hBD_m","isolume_m_at_0415","Nitracline_m","dbm","anp")]
    
    # dmean <- pplot.bin$mean[pplot.bin$mean$Z_CLASS==0,c("SIC_CLASS","idms_z60","idmspt_z60","icp_z60","iTchla_z60")]
    # dmin <- pplot.bin$min[pplot.bin$min$Z_CLASS==0,c("SIC_CLASS","idms_z60","idmspt_z60","icp_z60","iTchla_z60")]
    # dmax <- pplot.bin$max[pplot.bin$max$Z_CLASS==0,c("SIC_CLASS","idms_z60","idmspt_z60","icp_z60","iTchla_z60")]
    
    View(cbind(dmean,dmin,dmax))
  }
}

# --------------------------------------------------
# Compare classifications based on SIC and OWD

sel <- c(diff(pplot$station),1)
sel <- !is.na(sel) & sel!=0
cl_summary <- cbind(pplot[sel,c("stn","SICm2d","SICm1d","SICday","sic_class","OWD")],st_class$owd_class[sel])

# View(cl_summary)
# print(table(cl_summary$sic_class))
# print(table(cl_summary$`st_class$owd_class[sel]`))

# --------------------------------------------------
# Filter out statistics with less than 3 meas.
pplot.bin$mean[pplot.bin$count<2] <- NA
pplot.bin$median[pplot.bin$count<2] <- NA

# --------------------------------------------------
# FINAL PLOTTING FOR PAPER

# ---------------------
# Figure with concentrations

xvarS <- list(tchla = expression('TChla (µg L'^-1*')'),
              cpsmooth1 = expression('Cp (m'^-1*')'),
              dmspt = "DMSPt (nM)",
              dms = "DMS (nM)",
              chlc3 = expression('Chlc3 (µg L'^-1*')'),
              but19_like = expression('19-But-like (µg L'^-1*')'),
              peri = expression('Peridinin (µg L'^-1*')'),
              # psc = expression('PSC (µg L'^-1*')'), # Choose either photosynthetic carotenoids or phaeophorbide a (below)
              phdaSUM = expression('Phaeophorbide a (µg L'^-1*')'),
              temp = expression('Temperature ('*degree*'C)'),
              N2 = expression('Brunt-Väisälä freq. (h'^-1*')'),
              anp = "ANP (-)",
              par_d_p24h_ein_m_2_day_1 = expression('PAR (µE m'^-2*' d'^-1*')')
)
yvar <- "depth"
lett <- plet
names(lett) <- names(xvarS)

for (sc in "owd_class") {
  
  if (exportimg) {png(filename = paste0(opath,"Fig3_",sc,".png"), width = 17, height = 14, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
  
  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 4)
  mr1 <-  cbind(m0+1,m0+2,m0+3,m0+4)
  m <- rbind(mr1,mr1+4,mr1+8)
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  for (xvar in names(xvarS)) {
    
    print(xvar)
    par(mar = c(4,3,2,0.5))
    
    xl <- c(min(c(0,1.1*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))),
            1.1*max(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))
    if (xvar  %in% c("sal","sigt")) {xl[1] <- 0.9*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T)}
    if (xvar  == "anp") {xl <- rev(xl)}
    print(xl)
    
    plot(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
         y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
         ylim = c(70,0), xlim = xl, pch = 19, col = col[1], cex = 1.9, axes = F, xlab = "", ylab = "")
    box()
    axis(side = 1, cex.axis = 1.1)
    mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
    if (xvar %in% c("tchla","chlc3","temp")) {
      axis(side = 2, cex.axis = 1.1)
      mtext(side = 2, "Depth", cex = 0.9, line = 2.5)
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
    if (xvar  == "temp") {
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
  
  if (exportimg) {dev.off()}
}



# ---------------------
# Figure with ratios

xvarS <- list(dmspt2cp = expression('DMSPt/Cp (µmol m'^-2*')'),
              dmspt2tchla = expression('DMSPt/TChla (nmol µg'^-1*')'),
              cp2tchla = expression('Cp/TChla (m'^2*' mg'^-1*')'),
              dms2dmspt = expression('DMS/DMSPt (mol mol'^-1*')'),
              chlc3_2_tchla = expression('Chlc3/TChla (g g'^-1*')'),
              but19like_2_tchla = expression('19-But-like/TChla (g g'^-1*')'),
              peri_2_tchla = expression('Peridinin/TChla (g g'^-1*')'),
              phdaSUM2tchla = expression('Phaeophorb_a/TChla (g g'^-1*')'),
              psc2tchla = expression('Photosynthetic car./TChla (g g'^-1*')'),
              ppc2tchla = expression('Photoprotective car./TChla (g g'^-1*')'),
              dd = expression('(Dd+Dt)/TChla (g g'^-1*')')
)
yvar <- "depth"
lett <- plet
names(lett) <- names(xvarS)

for (sc in "owd_class") {
  
  if (exportimg) {png(filename = paste0(opath,"Fig5_",sc,".png"), width = 17, height = 14, units = 'cm', pointsize = 8, bg = 'white', res = plotres, type = 'cairo')}
  
  # Multipanel setup
  m0 <- matrix(data = 0, nrow = 4, ncol = 4)
  mr1 <-  cbind(m0+1,m0+2,m0+3,m0+4)
  m <- rbind(mr1,mr1+4,mr1+8)
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  for (xvar in names(xvarS)) {
    
    print(xvar)
    par(mar = c(4,3,2,0.5))
    
    xl <- c(min(c(0,1.1*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))),
            1.1*max(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T))
    if (xvar  %in% c("sal","sigt")) {xl[1] <- 0.9*min(cbind(pplot.bin$mean[,xvar],pplot.bin$median[,xvar]), na.rm = T)}
    if (xvar  == "anp") {xl <- rev(xl)}
    print(xl)
    
    plot(x = pplot.bin$median[pplot.bin$median$SIC_CLASS=="ICE",xvar],
         y = pplot.bin$mean[pplot.bin$median$SIC_CLASS=="ICE",yvar],
         ylim = c(70,0), xlim = xl, pch = 19, col = col[1], cex = 1.9, axes = F, xlab = "", ylab = "")
    box()
    axis(side = 1, cex.axis = 1.1)
    mtext(side = 1, xvarS[[xvar]], cex = 0.9, line = 3)
    if (xvar %in% c("dmspt2cp","chlc3_2_tchla","psc2tchla")) {
      axis(side = 2, cex.axis = 1.1)
      mtext(side = 2, "Depth", cex = 0.9, line = 2.5)
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
    if (xvar  == "psc2tchla") {
      legend(x = xl[1]+0.05*(xl[2]-xl[1]),
             y = 25,
             pch = 19,
             legend = c("ICE","MIZ","OW"),
             cex = 1.3,
             col = col,
             bg= "white", box.lwd = 0)
      legend(x = xl[1]+0.05*(xl[2]-xl[1]),
             y = 50,
             lty = c(1,3),
             legend = c("medians","means"),
             cex = 1.3,
             col = "black",
             bg= "white", box.lwd = 0)
    }
    
  }
  
  if (exportimg) {dev.off()}
}


# ------------------------------------------
# # Checks
# View(pplot.bin$count[,c("tchla","chlc3","but19_like","peri","phdaSUM","chlc3_2_tchla","but19like_2_tchla","peri_2_tchla","phdaSUM2tchla")])
