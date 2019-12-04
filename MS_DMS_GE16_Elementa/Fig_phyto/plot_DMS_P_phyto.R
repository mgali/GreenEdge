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

# Add clustering coefficient as per station by merging with profiles
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove DMS data points where DMSPt or Phaeocystis are missing
pplot <- pplot[!is.na(pplot$dmspt) & !is.na(pplot$dms_consens_cf68) &!is.na(pplot$Phaeo),]

# If phyto counts duplicated, remove duplicates
# NOTE: meging DMS(P) and taxonomy by station led to duplicate Phaeo values in stations with DMS in >1 cast
pplot <- pplot[!duplicated(pplot$Phaeo),]

# Add surface and SCM categories
pplot$scm <- 'SCM'
pplot$scm[pplot$depth <= 10] <- 'surface'

# Irradiance color scale
# col.par <- colorRamp(brewer.pal(11, 'Spectral'))

ncol.par <- 9
parpal <- rev(brewer.pal(ncol.par, 'YlGnBu'))
zcol <- pplot$par_d_p24h_ein_m_2_day_1 # Variable used for coloring
zcol[zcol>45] <- 45
class <- classIntervals(zcol, ncol.par, style="fixed", fixedBreaks = seq(0,45,5)) # fixed breaks
# class <- classIntervals(zcol, ncol.par, style="equal") # fixed breaks
parcol <- findColours(class, parpal)


# pal.par <- colorRamp(rev(brewer.pal(9, 'YlGnBu')))
# col.par <- pal.par(zcol / max(zcol, na.rm = T))

# ORGANIZE PLOTTING IN LOOP IF SEVERAL PLOTS ARE GOING TO BE MADE. For (j in c(dms,dmspt)), for (k in c(phyt1,2,...)), rows are j, cols are k

# Converting phyto counts to biomass would make a much stronger case for DMSP contribution to S and C cycling
# Lines of C_DMSP:C_tot could be added to plot panels

# ---------------------
# DMSPt vs. Phaeocystis

rs.sur <- cor(pplot$Phaeo[pplot$scm=='surface'],pplot$dmspt[pplot$scm=='surface'], use = "pairwise", method = "spearman")
rs.scm <- cor(pplot$Phaeo[pplot$scm!='surface'],pplot$dmspt[pplot$scm!='surface'], use = "pairwise", method = "spearman")
ps.sur <- cor.test(pplot$Phaeo[pplot$scm=='surface'],pplot$dmspt[pplot$scm=='surface'], use = "pairwise", method = "spearman")
ps.scm <- cor.test(pplot$Phaeo[pplot$scm!='surface'],pplot$dmspt[pplot$scm!='surface'], use = "pairwise", method = "spearman")

# Put corr stats in structure and construct string for plot title, highlighting significance with asterisk

if (exportimg) {png(filename = paste0(opath,"dmsp_phaeo_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$Phaeo, y = pplot$dmspt, log = "xy",
     type = "p", pch = 21, bg = "gray", cex = 1.5,
     xlab = "Phaeocystis (cells/mL)",
     ylab = "DMSPt (nM)",
     # xlim = c(5e3,2e7),
     ylim = c(10,550))
points(x = pplot$Phaeo[pplot$scm=='surface'], y = pplot$dmspt[pplot$scm=='surface'],
     pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMSPt vs. unidentified flagellates

if (exportimg) {png(filename = paste0(opath,"dmsp_flag_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

rs.sur <- cor(pplot$flag[pplot$scm=='surface'],pplot$dmspt[pplot$scm=='surface'], use = "pairwise", method = "spearman")
rs.scm <- cor(pplot$flag[pplot$scm!='surface'],pplot$dmspt[pplot$scm!='surface'], use = "pairwise", method = "spearman")

plot(x = pplot$flag, y = pplot$dmspt, log = "xy",
     type = "p", pch = 21, bg = "gray", cex = 1.5,
     xlab = "Unidentified flagellates (cells/mL)",
     ylab = "DMSPt (nM)",
     # xlim = c(5e3,2e7),
     ylim = c(10,550),
     main = "corr")
points(x = pplot$flag[pplot$scm=='surface'], y = pplot$dmspt[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMSPt vs. centric diatoms

if (exportimg) {png(filename = paste0(opath,"dmsp_diatcen_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$diat_cen, y = pplot$dmspt, log = "xy",
     type = "p", pch = 21, bg = "gray", cex = 1.5,
     xlab = "Centric diatoms (cells/mL)",
     ylab = "DMSPt (nM)",
     ylim = c(10,550))
points(x = pplot$diat_cen[pplot$scm=='surface'], y = pplot$dmspt[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMSPt vs. pennate diatoms

if (exportimg) {png(filename = paste0(opath,"dmsp_diatpen_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$diat_pen, y = pplot$dmspt, log = "xy",
     type = "p", pch = 21, bg = "gray", cex = 1.5,
     xlab = "Pennate diatoms (cells/mL)",
     ylab = "DMSPt (nM)",
     ylim = c(10,550))
points(x = pplot$diat_pen[pplot$scm=='surface'], y = pplot$dmspt[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMS vs. Phaeocystis

if (exportimg) {png(filename = paste0(opath,"dms_phaeo_counts.png"), width = 6, height = 5, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

m <- cbind(matrix(data = 1, nrow = 5, ncol = 5), matrix(data = 2, nrow = 5, ncol = 1))
layout(m)
par(oma = c(1,1,0.5,0.5))

par(mar = c(2,4,2,0))
plot(x = pplot$Phaeo, y = pplot$dms_consens_cf68, log = "xy",
     type = "p",
     col = parcol,
     pch = 19,
     cex = 1.5,
     # bg = "gray",
     xlab = "Phaeocystis (cells/mL)",
     ylab = "DMS (nM)",
     ylim = c(0.5,100)
     )
# points(x = pplot$Phaeo[pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
#        pch = 21, bg = "white", cex = 1.5)
# Adding z variable in color scale. Could show PAR for DMS
# points(x = pplot$Phaeo[pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
#        pch = 19,
#        # col = parcol,
#        cex = 1.5) #col.par(1)

par(mar = c(2,1,2,2))
cbar <- seq(0, max(zcol), length.out = ncol.par)
image(t(cbar), cbar, col = plotc, cex.main=1, xaxt = "n", yaxt = "n", main = expression(paste(E_d)))

if (exportimg) {dev.off()}

# ---------------------
# DMS vs. unidentified flagellates

if (exportimg) {png(filename = paste0(opath,"dms_flag_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$flag, y = pplot$dms_consens_cf68, log = "xy",
     type = "p", pch = 21, cex = 1.5,
     bg = "gray",
     xlab = "Unidentified flagellates (cells/mL)",
     ylab = "DMS (nM)",
     ylim = c(0.5,100)
)
points(x = pplot$flag[pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMS vs. centric diatoms

if (exportimg) {png(filename = paste0(opath,"dms_diatcen_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$diat_cen, y = pplot$dms_consens_cf68, log = "xy",
     type = "p", pch = 21, cex = 1.5,
     bg = "gray",
     xlab = "Centric diatoms (cells/mL)",
     ylab = "DMS (nM)",
     ylim = c(0.5,100)
)
points(x = pplot$diat_cen[pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}

# ---------------------
# DMS vs. centric diatoms

if (exportimg) {png(filename = paste0(opath,"dms_diatpen_counts.png"), width = 5, height = 5, units = 'cm', pointsize = 5, bg = 'white', res = 600, type = 'cairo')}

plot(x = pplot$diat_pen, y = pplot$dms_consens_cf68, log = "xy",
     type = "p", pch = 21, cex = 1.5,
     bg = "gray",
     xlab = "Pennate diatoms (cells/mL)",
     ylab = "DMS (nM)",
     ylim = c(0.5,100)
)
points(x = pplot$diat_pen[pplot$scm=='surface'], y = pplot$dms_consens_cf68[pplot$scm=='surface'],
       pch = 21, bg = "white", cex = 1.5)

if (exportimg) {dev.off()}



# # Correlations
# 
# print(cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "spearman"))
# print(cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "pearson"))
# print((cor(pplot$xvar,pplot$yvar, use = "pairwise", method = "pearson"))^2)

       