# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(tidyverse)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)
surf.all <- surf.all[,c('stn','AW_ArW_clustering_coefficient')]
colnames(surf.all) <- c('stn','CC')

# Restrict depth range for plotting
prof.plots <- prof.all[prof.all$depth < 60 & !is.na(prof.all$depth),]

# Create factor variables from continuous ones for plotting
prof.plots$Depth <- factor(sapply(floor(prof.plots$depth/10)*10, function(x) paste0(x,'_',x+10)))
prof.plots$ANP <- prof.plots$anp

# Add clustering coefficient as per station
prof.plots <- merge(x = prof.plots, y = surf.all, all.x = T, all.y = F, by = 'stn')

# Create labels dictionary? Not possible if using xvar and yvar as colnames


# Rename variables and create labels
# prof.plots$xvar <- prof.plots[,'tchla']; xl <- 'TChla (µg/L)'; xn <- 'tchla'
# prof.plots$xvar <- prof.plots[,'cpsmooth1']; xl <- 'Cp (1/m)'; xn <- 'Cp'
prof.plots$xvar <- prof.plots[,'dmspt']; xl <- 'DMSPt (nM)'; xn <- 'dmspt'
prof.plots$yvar <- prof.plots[,'dms_consens_cf68']; yl <- 'DMS (nM)'; yn <- 'dms'
# prof.plots$yvar <- prof.plots[,'dmspt']; yl <- 'DMSPt (nM)'; yn <- 'dmspt'
bb <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500)

# ----------------
# Plot xvar (chl, cp, DMSPt) vs. yvar (DMS, DMSPt) by depth and ANP
p <- ggplot(prof.plots, aes(x = xvar, y = yvar, shape = Depth, colour = ANP)) + geom_point(size = 3)
p + scale_color_gradient(low="red", high="blue") +
  scale_x_continuous(name = xl, breaks = bb, trans = 'log10') + # labels, limits
  scale_y_continuous(name = yl, breaks = bb, trans = 'log10') # labels, limits

# Superimpose reference lines for y/x ratios
# Consider adding light (including ice effect) in some way

# Save
fname <- paste('fig',xn,yn,'z','ANP',sep = '_')
ggsave(paste0(genpath,pdir,fname,'.png'), width = 8, height = 6, units = 'cm', dpi = 300, scale = 2)

# ----------------
# Plot xvar (chl, cp, DMSPt) vs. yvar (DMS, DMSPt) by depth and CC
p <- ggplot(prof.plots, aes(x = xvar, y = yvar, shape = Depth, colour = CC)) + geom_point(size = 3)
p + scale_color_gradient(low="blue", high="red") +
  scale_x_continuous(name = xl, breaks = bb, trans = 'log10') + # labels, limits
  scale_y_continuous(name = yl, breaks = bb, trans = 'log10') # labels, limits

# Superimpose reference lines for y/x ratios
# Consider adding light (including ice effect) in some way

# Save
fname <- paste('fig',xn,yn,'z','CC',sep = '_')
ggsave(paste0(genpath,pdir,fname,'.png'), width = 8, height = 6, units = 'cm', dpi = 300, scale = 2)


# Correlations

print(cor(prof.plots$xvar,prof.plots$yvar, use = "pairwise", method = "spearman"))
print(cor(prof.plots$xvar,prof.plots$yvar, use = "pairwise", method = "pearson"))
print((cor(prof.plots$xvar,prof.plots$yvar, use = "pairwise", method = "pearson"))^2)

       