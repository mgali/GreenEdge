# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(RColorBrewer)
library(classInt) # for function classIntervals
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- F
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"

# Phyto group plotted
pg <- "Phaeo" # Phaeo, diat_cen, diat_pen
xl <- list(diat_pen = "Pennate diatoms (cells/L)",
           diat_cen = "Centric diatoms (cells/L)",
           Phaeo = 'Phaeocystis (cells/L)')

# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Add clustering coefficient as per station by merging with profiles
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove DMS data points where DMSPt or Phaeocystis are missing, and stations where microscopy counts not done
pplot <- pplot[!is.na(pplot$dmspt) & !is.na(pplot$dms),]
pplot <- pplot[pplot$stn >= 418,]

# # If phyto counts duplicated, remove duplicates
# # NOTE: meging DMS(P) and taxonomy by station led to duplicate Phaeo values in stations with DMS in >1 cast
# pplot <- pplot[!duplicated(pplot[[pg]]),]

# Remove duplicated rows
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) & (!is.na(pplot$tchla) & !is.na(pplot$cpsmooth1))
pplot <- pplot[!dd,]

# Add surface and SCM categories
pplot$scm <- 'SCM'
pplot$scm[pplot$depth < 10] <- 'surface'

# ---------------------
# Correlations between DMS(P) and phyto groups
Xall <- pplot[,c("diat_cen","diat_pen","dino_athec","dino_thec","chlor","chrys","dictyo","crypt","eugl","pras","prym","Phaeo","flag","raph","cyan","hetero","choano","cilli")]
Xall$prym_noPhaeo <- Xall$prym - Xall$Phaeo
Xall[Xall == 0] <- NA
Xsel <- Xall[,c("diat_cen","diat_pen","dino_athec","dino_thec","chrys","crypt","Phaeo","flag")]
# boxplot(Xsel, log = "y")
# !!!!!! dictyochales = silicoflagellata, raphidophyceae: filum Heterocontophyta, mostly Heterosigma akashiwo, Euglenophyta, Choanoflagellates !!!!!!

R <- list()
R.r1 <- list()
R.r2 <- list()
R.r3 <- list()
selStats <- c("p.value","estimate")
for (mm in c("spearman","pearson")) {
  for (cc in c("dmspt","dms")) {
    for (zz in c("SCM","surface")) {
      iz <- pplot$scm == zz
      tmp <- lapply(Xsel[iz,], function(x) cor.test(x, pplot[[cc]][iz], use = "pairwise", method = mm))
      tmp <- lapply(tmp, function(x) x <- x[selStats]) # Select desired stats, otherwise different stats in spearman and pearson (latter has confidence interval with length 2 which complicates formatting)
      TMP <- data.table::rbindlist(tmp, use.names = T, fill = F, idcol = "Taxo")
      # TMP <- do.call("rbind", tmp)
      R[[mm]][[cc]][[zz]] <- TMP
    }
    R.r1[[mm]][[cc]] <- data.table::rbindlist(R[[mm]][[cc]], use.names = T, fill = F, idcol = "scm")
  }
  R.r2[[mm]] <- data.table::rbindlist(R.r1[[mm]], use.names = T, fill = F, idcol = "compound")
}
R.r3 <- data.table::rbindlist(R.r2, use.names = T, fill = F, idcol = "method")
# write.csv(x = R.r3, file = paste0(opath,"CorrMat_DMS_P_PhytoCounts.csv"))

# ------------------------------------------
# Plot

# Function to plot significance levels as symbols
f_asterisks <- function(xpos, y.pval){
  aheight <- 0.9
  points(xpos[y.pval < 0.01], rep(aheight,length(xpos))[y.pval < 0.01], pch = 8, col = "darkorange")
  points(xpos[y.pval >= 0.01 & y.pval < 0.05], rep(aheight,length(xpos))[y.pval >= 0.01 & y.pval < 0.05], pch = 4, col = "darkorange")
  points(xpos[y.pval >= 0.05 & y.pval < 0.10], rep(aheight,length(xpos))[y.pval >= 0.05 & y.pval < 0.10], pch = 3, col = "darkorange")
}


# Irradiance color scale
ncol.par <- 9
parpal <- rev(brewer.pal(ncol.par, 'YlGnBu'))
zcol <- pplot$par_d_p24h_ein_m_2_day_1 # Variable used for coloring
zcol[zcol>45] <- 45
zbreaks <- seq(0,45,5)
class <- classIntervals(zcol, ncol.par, style="fixed", fixedBreaks = zbreaks) # fixed breaks
# class <- classIntervals(zcol, ncol.par, style="equal") # fixed breaks
parcol <- findColours(class, parpal)

xbars <- t(R$spearman$dmspt$SCM$Taxo)
xlbars <- list("diat_cen"="Diat_C",
               "diat_pen"="Diat_P",
               "dino_athec"="Dino_A",
               "dino_thec"="Dino_T",
               "chrys"="Chryso",
               "crypt"="Crypto",
               "Phaeo"="Phaeocystis",
               "flag"="Flag_other")
xpos <- seq(1, length(xbars))*3-1.5

for (mm in c("spearman","pearson")) {
  if (exportimg) {png(filename = paste0(opath,"Fig6_phytoCounts_dms_dmspt_corr_",mm,".png"), width = 16, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}
  
  # ---------------------
  # Multipanel setup
  m1_2 <- rbind(matrix(data = 1, nrow = 5, ncol = 9), matrix(data = 2, nrow = 6, ncol = 9))
  m3_4 <- rbind(matrix(data = 3, nrow = 5, ncol = 5), matrix(data = 4, nrow = 6, ncol = 5))
  m <- cbind(m1_2, m3_4, matrix(data = 5, nrow = 11, ncol = 2))
  layout(m)
  par(oma = c(1,1,0.5,0.5))
  
  # ---------------------
  # a) DMSPt vs. Phyto counts correlations
  par(mar = c(3,4,1,1))
  Ybars <- t(cbind(R[[mm]][["dmspt"]][["surface"]][["estimate"]], R[[mm]][["dmspt"]][["SCM"]][["estimate"]]))
  barplot(height = Ybars,
          names.arg = rep("", length(xbars)),
          beside = T,
          legend.text = c("surface","SCM"),
          args.legend = list(x = "topleft", bty = "n"),
          col = parpal[c(7,1)],
          border = F,
          ylab = "Spearman's correlation with DMSPt",
          ylim = c(min(0,min(Ybars, na.mm = T)-0.1),1),
          cex.axis = 0.9,
          cex.names = 0.95,
          las = 1)
  f_asterisks(xpos, t(R[[mm]][["dmspt"]][["surface"]][["p.value"]]))
  f_asterisks(xpos+1, t(R[[mm]][["dmspt"]][["SCM"]][["p.value"]]))
  abline(h = seq(0,1,0.2), col = "gray70", lwd = 0.5, lty = 3)
  mtext(text = "a)", line = 0, adj = 0.05, cex = 0.9)
  
  # ---------------------
  # c) DMS vs. Phyto counts correlations
  par(mar = c(5,4,0,1))
  Ybars <- t(cbind(R[[mm]][["dms"]][["surface"]][["estimate"]], R[[mm]][["dms"]][["SCM"]][["estimate"]]))
  barplot(height = Ybars,
          names.arg = xlbars[xbars],
          beside = T,
          col = parpal[c(7,1)],
          border = F,
          ylab = "Spearman's correlation with DMS",
          ylim = c(min(0,min(Ybars, na.rm = T)-0.1),1),
          cex.axis = 0.9,
          cex.names = 0.95,
          las = 2)
  f_asterisks(xpos, t(R[[mm]][["dms"]][["surface"]][["p.value"]]))
  f_asterisks(xpos+1, t(R[[mm]][["dms"]][["SCM"]][["p.value"]]))
  abline(h = seq(-0.5,1,0.5), col = "gray70", lwd = 0.5, lty = 3)
  mtext(text = "c)", line = 0, adj = 0.05, cex = 0.9)
  
  # ---------------------
  # b) DMSPt vs. Phaeocystis
  par(mar = c(2,5,1,1))
  plot(x = pplot[[pg]], y = pplot$dmspt, log = "xy",
       type = "p",
       col = parcol,
       pch = 19,
       cex = 1.1,
       cex.lab = 0.9,
       cex.axis = 0.9,
       xlab = "",
       ylab = "DMSPt (nM)",
       xlim = c(1e3,2e7),
       ylim = c(10,550),
       axes = F
  )
  points(x = pplot[[pg]][pplot$scm=='surface'], y = pplot$dmspt[pplot$scm=='surface'],
         pch = 1,
         lwd = 1,
         col = "gray",
         cex = 1.6)
  points(x = pplot[[pg]][pplot$scm=='SCM'], y = pplot$dmspt[pplot$scm=='SCM'],
         pch = 1,
         lwd = 1,
         col = "black",
         cex = 1.6)
  axis(1, labels = F, cex.axis = 0.9)
  axis(2, labels = T, cex.axis = 0.9)
  # box()
  mtext(text = "b)", line = 0, adj = 0.05, cex = 0.9)
  
  # ---------------------
  # d) DMS vs. Phaeocystis
  par(mar = c(5,5,1,1))
  plot(x = pplot[[pg]], y = pplot$dms, log = "xy",
       type = "p",
       col = parcol,
       pch = 19,
       cex = 1.1,
       cex.lab = 0.9,
       cex.axis = 0.9,
       xaxp = c(10^4, 10^7, 1), # Uncomment only if (pg=="Phaeo") 
       # yaxp = c(1, 100, 2),
       xlab = xl[[pg]],
       ylab = "DMS (nM)",
       xlim = c(1e3,2e7),
       ylim = c(0.5,100),
       axes = F
  )
  points(x = pplot[[pg]][pplot$scm=='surface'], y = pplot$dms[pplot$scm=='surface'],
         pch = 1,
         lwd = 1,
         col = "gray",
         cex = 1.6)
  points(x = pplot[[pg]][pplot$scm=='SCM'], y = pplot$dms[pplot$scm=='SCM'],
         pch = 1,
         lwd = 1,
         col = "black",
         cex = 1.6)
  axis(1, labels = T, cex.axis = 0.9)
  axis(2, labels = T, cex.axis = 0.9)
  # box()
  mtext(text = "d)", line = 0, adj = 0.05, cex = 0.9)
  
  # ---------------------
  # e) Common color bar showing PAR
  par(mar = c(9,2,7,6)) # May cause margins error depending on figure window size
  # par(mar = c(4,4,0,0))
  cbar <- zbreaks[1:length(zbreaks)-1]
  image(t(cbar), cbar,
        col = parpal,
        cex = 1,
        cex.main = 1,
        cex.axis = 0.9,
        cex.lab = 0.9,
        xaxt = "n", yaxt = "n")
  axis(4, cex.axis=1, mgp = c(0, .5, 0), at = zbreaks-2.5, labels = zbreaks)
  mtext(text = expression(paste('daily PAR (mol photons ',m^-2,' ',d^-1,')')), side = 4, cex = 0.7, line = 2.5)
  
  if (exportimg) {dev.off()}
}

