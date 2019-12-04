# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(classInt)
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"

# Phyto group plotted
pg <- "Phaeo" # Phaeo, diat_cen, diat_pen
xl <- list(diat_pen = "Pennate diatoms (cells/mL)",
           diat_cen = "Centric diatoms (cells/mL)",
           Phaeo = 'Phaeocystis (cells/mL)')


# Add clustering coefficient as per station by merging with profiles
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove DMS data points where DMSPt or Phaeocystis are missing
pplot <- pplot[!is.na(pplot$dmspt) & !is.na(pplot$dms_consens_cf68),]

# If phyto counts duplicated, remove duplicates
# NOTE: meging DMS(P) and taxonomy by station led to duplicate Phaeo values in stations with DMS in >1 cast
pplot <- pplot[!duplicated(pplot[[pg]]),]

# Add surface and SCM categories
pplot$scm <- 'SCM'
pplot$scm[pplot$depth <= 10] <- 'surface'

# Irradiance color scale
# col.par <- colorRamp(brewer.pal(11, 'Spectral'))

ncol.par <- 9
parpal <- rev(brewer.pal(ncol.par, 'YlGnBu'))
zcol <- pplot$par_d_p24h_ein_m_2_day_1 # Variable used for coloring
zcol[zcol>45] <- 45
zbreaks <- seq(0,45,5)
class <- classIntervals(zcol, ncol.par, style="fixed", fixedBreaks = zbreaks) # fixed breaks
# class <- classIntervals(zcol, ncol.par, style="equal") # fixed breaks
parcol <- findColours(class, parpal)


# pal.par <- colorRamp(rev(brewer.pal(9, 'YlGnBu')))
# col.par <- pal.par(zcol / max(zcol, na.rm = T))

# ORGANIZE PLOTTING IN LOOP IF SEVERAL PLOTS ARE GOING TO BE MADE. For (j in c(dms,dmspt)), for (k in c(phyt1,2,...)), rows are j, cols are k

# Converting phyto counts to biomass would make a much stronger case for DMSP contribution to S and C cycling
# Lines of C_DMSP:C_tot could be added to plot panels

# ---------------------
# Correlations

rs.sur <- cor(pplot[[pg]][pplot$scm=='surface'],pplot$dmspt[pplot$scm=='surface'], use = "pairwise", method = "spearman")
rs.scm <- cor(pplot[[pg]][pplot$scm!='surface'],pplot$dmspt[pplot$scm!='surface'], use = "pairwise", method = "spearman")
ps.sur <- cor.test(pplot[[pg]][pplot$scm=='surface'],pplot$dmspt[pplot$scm=='surface'], use = "pairwise", method = "spearman")
ps.scm <- cor.test(pplot[[pg]][pplot$scm!='surface'],pplot$dmspt[pplot$scm!='surface'], use = "pairwise", method = "spearman")


# ---------------------
# Plot

if (exportimg) {png(filename = paste0(opath,"dms_dmspt_",pg,"_counts.png"), width = 6, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

m1_2 <- rbind(matrix(data = 1, nrow = 4, ncol = 5), matrix(data = 2, nrow = 4, ncol = 5))
m <- cbind(m1_2, matrix(data = 3, nrow = 8, ncol = 2))
layout(m)
par(oma = c(1,1,0.5,0.5))

# ---------------------
# DMSPt vs. Phaeocystis
par(mar = c(3,4,1,0))
plot(x = pplot[[pg]], y = pplot$dmspt, log = "xy",
     type = "p",
     col = parcol,
     pch = 19,
     cex = 1.2,
     cex.axis = 0.9,
     xaxt = "n",
     xlab = "",
     ylab = "DMSPt (nM)",
     ylim = c(10,550)
     )
if (pg=="Phaeo") {axis(side = 1, at = 10^(4:7), cex.axis = 1, labels = NULL)}
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

# ---------------------
# DMS vs. Phaeocystis
par(mar = c(4,4,0,0))
plot(x = pplot[[pg]], y = pplot$dms_consens_cf68, log = "xy",
     type = "p",
     col = parcol,
     pch = 19,
     cex = 1.2,
     cex.axis = 0.9,
     xaxp = c(10^4, 10^7, 1), # Uncomment only if (pg=="Phaeo") 
     # yaxp = c(1, 100, 2),
     xlab = xl[[pg]],
     ylab = "DMS (nM)",
     ylim = c(0.5,100)
     )
points(x = pplot[[pg]][pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
       pch = 1,
       lwd = 1,
       col = "gray",
       cex = 1.6)
points(x = pplot[[pg]][pplot$scm=='SCM'], y = pplot$dms_consens_cf68[pplot$scm=='SCM'],
       pch = 1,
       lwd = 1,
       col = "black",
       cex = 1.6)

# ---------------------
# Common color bar showing PAR
par(mar = c(9,2,7,4))
cbar <- zbreaks[1:length(zbreaks)-1]
image(t(cbar), cbar,
      col = parpal,
      cex = 1,
      cex.main = 1,
      cex.axis = 0.9,
      xaxt = "n", yaxt = "n")
axis(4, cex.axis=1, mgp = c(0, .5, 0), at = zbreaks-2.5, labels = zbreaks)
mtext(text = expression(paste('daily PAR (mol photons ',m^-2,' ',d^-1,')')), side = 4, cex = 0.7, line = 2.5)

if (exportimg) {dev.off()}


# # Correlations
# 
# print(cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "spearman"))
# print(cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "pearson"))
# print((cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "pearson"))^2)

       