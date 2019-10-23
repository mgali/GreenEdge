# MAKE PLOTS OF DMS AND DMSP VS. OTHER RELEVANT X VARIABLES

library(tidyverse)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
View(prof.all)

# Restrict depth range for plotting
prof.plots <- prof.all[prof.all$depth < 60 & !is.na(prof.all$depth),]

# Create factor variables from continuous ones for plotting
prof.plots$Depth <- factor(sapply(floor(prof.plots$depth/10)*10, function(x) paste0(x,'_',x+10)))
prof.plots$ANP <- prof.plots$anp

# Create labels dictionary? Not possible if using xvar and yvar as colnames


# Rename variables and create labels
# prof.plots$xvar <- prof.plots[,'chla']; xl <- 'Chla (µg/L)'
# prof.plots$xvar <- prof.plots[,'cpsmooth1']; xl <- 'Cp (1/m)'
prof.plots$xvar <- prof.plots[,'dmspt']; xl <- 'DMSPt (nM)'
prof.plots$yvar <- prof.plots[,'dms_consens_cf68']; yl <- 'DMS (nM)'
# prof.plots$yvar <- prof.plots[,'dmspt']; yl <- 'DMSPt (nM)'
bb <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500)

# Plot xvar (chl, cp, DMSPt) vs. yvar (DMS, DMSPt) by depth and ANP
p <- ggplot(prof.plots, aes(x = xvar, y = yvar, shape = Depth, colour = ANP)) + geom_point(size = 3)
p + scale_color_gradient(low="blue", high="red") +
  # xlab('Chla (µg/L)') +
  # ylab('DMS (nM)') +
  # coord_trans(x = 'log10', y = 'log10') +
  scale_x_continuous(name = xl, breaks = bb, trans = 'log10') + # labels, limits
  scale_y_continuous(name = yl, breaks = bb, trans = 'log10') # labels, limits

# Superimpose reference lines for y/x ratios

# Consider adding light (including ice effect) in some way

# Save
fname <- paste('fig','DMSPt','DMS','z','ANP',sep = '_')
ggsave(paste0(genpath,pdir,fname,'.png'))
       