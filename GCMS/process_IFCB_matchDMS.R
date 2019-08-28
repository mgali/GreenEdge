# ANALYZE PHYTOPLANKTON COMMUNITY FROM IMAGING FLOW CYTOBOT (IFCB) DATA
# ONLY MEASUREMENTS MATCHING DMS(P) DURING GREEN EDGE CRUISE LEG 1B

library(dplyr)

# Set general data path
genpath <- '~/Desktop/GreenEdge/GCMS/'
ifcbpath <- '~/Desktop/GreenEdge/IFCB/'

# ------------------------------------------------------------------------
# Load DMS-CTD, pigments, taxo. Merge (rbind, no shared coluns). Sort by cast/depth
prof.ctd.dms <- read.csv(paste0(genpath,'GE2016.profiles.DATA_CTD.DMS.csv'))
prof.taxo <- read.csv(paste0(genpath,'GE2016.profiles.TAXO.csv'))
prof.pigm <- read.csv(paste0(genpath,'GE2016.profiles.PIGM.csv'))
prof.all <- cbind(prof.ctd.dms, prof.taxo, prof.pigm)
prof.all <- prof.all[order(prof.all$cast, prof.all$depth),]

# ------------------------------------------------------------------------
# Load IFCB logsheet. Rename cast and depth columns for matching. Sort by cast/depth
ifcb.log <- read.csv(paste0(ifcbpath,'IFCB_logsheet_amundsen2016.csv'))
nn <- colnames(ifcb.log)
nn[nn == 'Cast'] <- 'cast'
nn[nn == 'Depth_.m.'] <- 'depth'
colnames(ifcb.log) <- nn
ifcb.log <- ifcb.log[order(ifcb.log$cast, ifcb.log$depth),]

# ------------------------------------------------------------------------
# Loop through casts and depths and match to IFCB
cast <- unique(prof.all$cast)
for (jc in cast) {
  
  # Subset prof.all data and find IFCB file match (merge)
  tmp.cast <- subset.data.frame(prof.all[prof.all$cast == jc,], select = c(cast,depth))
  tmp.cast$depth[tmp.cast$depth == 0.7] <- 0  # Set 0.7 to 0 m depth for exact matching of IFCB
  tmp.cast <- merge(x = tmp.cast, y = ifcb.log, by = c('cast','depth'), all.x = T, all.y = F)
  # tmp.cast <- dplyr::left_join(x = tmp.cast, ifcb.log) # equivalent
  
  # Loop through ifcb filenames in cast (match unique names, then expand for repeated depths)
  cfile <- unique(tmp.cast$File_Name)
  for (jf in cfile) {
    
    if (!is.na(jf)) {
      fpath <- paste0(ifcbpath,'post_processing/export_144_20190111_1439/processed_counts/',
                      jf,'.csv')
      ifcb.data <- read.csv(fpath, header = T)
      print(ifcb.data)
    }
  }
  
  # print(tmp.cast)
  
}