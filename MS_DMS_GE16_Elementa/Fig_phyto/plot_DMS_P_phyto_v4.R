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
exportimg <- T
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
f_asterisks <- function(xpos, ypos, y.pval, tcol){
  points(xpos[y.pval < 0.01], rep(ypos,length(xpos))[y.pval < 0.01], pch = 8, col = tcol, cex = 0.8, lwd = 0.5)
  points(xpos[y.pval >= 0.01 & y.pval < 0.05], rep(ypos,length(xpos))[y.pval >= 0.01 & y.pval < 0.05], pch = 4, col = tcol, cex = 0.8, lwd = 0.5)
  points(xpos[y.pval >= 0.05 & y.pval < 0.10], rep(ypos,length(xpos))[y.pval >= 0.05 & y.pval < 0.10], pch = 3, col = tcol, cex = 0.8, lwd = 0.5)
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
xpos <- (seq(1, length(xbars)))*6.1 - 1.7

mm <- "spearman"
mb <- "pearson"
if (exportimg) {png(filename = paste0(opath,"Fig6_phytoCounts_dms_dmspt_corr_v4.png"), width = 16, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

# ---------------------
# Multipanel setup
m1_2 <- rbind(matrix(data = 1, nrow = 4, ncol = 9), matrix(data = 2, nrow = 4, ncol = 9))
m3_4 <- rbind(matrix(data = 3, nrow = 4, ncol = 5), matrix(data = 4, nrow = 4, ncol = 5))
m <- cbind(m1_2, m3_4, matrix(data = 5, nrow = 8, ncol = 2))
layout(m)
par(oma = c(1,1,0.5,0.5))

# ---------------------
# a) DMSPt vs. Phyto counts correlations
par(mar = c(4,4,1,1))
Ybars <- t(cbind(R[[mm]][["dmspt"]][["surface"]][["estimate"]], R[[mm]][["dmspt"]][["SCM"]][["estimate"]]))
barplot(height = Ybars,
        names.arg = xlbars[xbars],
        beside = T,
        space = c(0.2,4.58),
        width = 0.9,
        legend.text = c("surface","SCM"),
        args.legend = list(x = "topleft", bty = "n"),
        col = parpal[c(6,1)],
        border = parpal[c(6,1)],
        ylab = "Correlation coefficient, DMSPt",
        ylim = c(min(0,min(Ybars, na.mm = T)-0.1),1.13),
        cex.axis = 0.9,
        cex.names = 0.8,
        las = 1)
# abline(h = seq(0,1,.2), lwd = 0.5, lty = 3, col = "gray")
# barplot(height = Ybars,
#         names.arg = xlbars[xbars],
#         beside = T,
#         add = T,
#         space = c(0.2,4.58),
#         width = 0.9,
#         legend.text = c("surface","SCM"),
#         args.legend = list(x = "topleft", bty = "n"),
#         col = parpal[c(6,1)],
#         border = parpal[c(6,1)],
#         ylab = "Correlation coefficient, DMSPt",
#         ylim = c(min(0,min(Ybars, na.mm = T)-0.1),1.13),
#         cex.axis = 0.9,
#         cex.names = 0.8,
#         las = 1)
f_asterisks(xpos, 1.05, t(R[[mm]][["dmspt"]][["surface"]][["p.value"]]), "black")
f_asterisks(xpos+1, 1.05, t(R[[mm]][["dmspt"]][["SCM"]][["p.value"]]), "black")
# Add pearson corr
Ybars <- t(cbind(R[[mb]][["dmspt"]][["surface"]][["estimate"]], R[[mb]][["dmspt"]][["SCM"]][["estimate"]]))
barplot(height = Ybars,
        names.arg = NULL,
        beside = T,
        add = T,
        space = c(1.1,9.1),
        width = 0.5,
        col = "white",
        border = parpal[c(6,1)],
        ylim = c(min(0,min(Ybars, na.mm = T)-0.1),1.13),
        cex.axis = 0.9,
        cex.names = 0.8,
        las = 1)
f_asterisks(xpos+0.5, 1.1, t(R[[mb]][["dmspt"]][["surface"]][["p.value"]]), "darkgray")
f_asterisks(xpos+1.5, 1.1, t(R[[mb]][["dmspt"]][["SCM"]][["p.value"]]), "darkgray")
abline(h = 0, lwd = 1, lty = 1, col = "black")

# ---------------------
# b) DMS vs. Phyto counts correlations
par(mar = c(4,4,1,1))
Ybars <- t(cbind(R[[mm]][["dms"]][["surface"]][["estimate"]], R[[mm]][["dms"]][["SCM"]][["estimate"]]))
barplot(height = Ybars,
        names.arg = xlbars[xbars],
        beside = T,
        space = c(0.2,4.58),
        width = 0.9,
        col = parpal[c(6,1)],
        border = parpal[c(6,1)],
        ylab = "Correlation coefficient, DMS",
        ylim = c(min(0,min(Ybars, na.rm = T)),1.13),
        cex.axis = 0.9,
        cex.names = 0.8,
        las = 1)
# abline(h = seq(-0.5,1,.25), lwd = 0.5, lty = 3, col = "gray")
# barplot(height = Ybars,
#         names.arg = xlbars[xbars],
#         beside = T,
#         add = T,
#         space = c(0.2,4.58),
#         width = 0.9,
#         col = parpal[c(6,1)],
#         border = parpal[c(6,1)],
#         ylab = "Correlation coefficient, DMS",
#         ylim = c(min(0,min(Ybars, na.rm = T)),1.13),
#         cex.axis = 0.9,
#         cex.names = 0.8,
#         las = 1)
f_asterisks(xpos, 1.05, t(R[[mm]][["dms"]][["surface"]][["p.value"]]), "black")
f_asterisks(xpos+1, 1.05, t(R[[mm]][["dms"]][["SCM"]][["p.value"]]), "black")
# Add pearson corr
Ybars <- t(cbind(R[[mb]][["dms"]][["surface"]][["estimate"]], R[[mb]][["dms"]][["SCM"]][["estimate"]]))
barplot(height = Ybars,
        names.arg = NULL,
        beside = T,
        add = T,
        space = c(1.1,9.1),
        width = 0.5,
        col = "white",
        border = parpal[c(6,1)],
        ylim = c(min(0,min(Ybars, na.mm = T)-0.2),1.13),
        cex.axis = 0.9,
        cex.names = 0.8,
        las = 1)
f_asterisks(xpos+0.5, 1.1, t(R[[mb]][["dms"]][["surface"]][["p.value"]]), "darkgray")
f_asterisks(xpos+1.5, 1.1, t(R[[mb]][["dms"]][["SCM"]][["p.value"]]), "darkgray")
abline(h = 0, lwd = 1, lty = 1, col = "black")

# ---------------------
# c) DMSPt vs. Phaeocystis
par(mar = c(4,5,1,1))
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
axis(2, labels = T, cex.axis = 0.9, las = 1)
# box()

# ---------------------
# d) DMS vs. Phaeocystis
par(mar = c(5,5,0,1))
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
axis(2, labels = T, cex.axis = 0.9, las = 1)
# box()

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
axis(4, cex.axis=1, mgp = c(0, .8, 0), at = zbreaks-2.5, labels = zbreaks, las = 1)
mtext(text = expression(paste('daily PAR (mol photons ',m^-2,' ',d^-1,')')), side = 4, cex = 0.7, line = 2.5)

if (exportimg) {dev.off()}



# ----
# Testing
# Ybars <- t(cbind(R[[mm]][["dmspt"]][["surface"]][["estimate"]], R[[mb]][["dmspt"]][["surface"]][["estimate"]],
#                  R[[mm]][["dmspt"]][["SCM"]][["estimate"]], R[[mb]][["dmspt"]][["SCM"]][["estimate"]]))
# barplot(height = Ybars,
#         names.arg = xlbars[xbars],
#         beside = T,
#         space = c(0.1,3),
#         legend.text = c("surface","rp","SCM","rp"),
#         args.legend = list(x = "topleft", bty = "n"),
#         col = parpal[c(7,NA,1,NA)],
#         border = parpal[c(7,7,1,1)],
#         ylab = "Correlation coefficient, DMSPt",
#         ylim = c(min(0,min(Ybars, na.mm = T)-0.1),1.13),
#         cex.axis = 0.9,
#         cex.names = 0.8,
#         las = 1)