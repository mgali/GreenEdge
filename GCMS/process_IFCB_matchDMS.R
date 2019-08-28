# ANALYZE PHYTOPLANKTON COMMUNITY FROM IMAGING FLOW CYTOBOT (IFCB) DATA
# ONLY MEASUREMENTS MATCHING DMS(P) DURING GREEN EDGE CRUISE LEG 1B
#
# Note: instead of preallocating, fill empty list in loop (all items are equal length vectors)
# and then reshape to df using do.call("rbind, mylist)
# With this procedure, col names of matched ifcb variables directly set in f_get_ifcb4DMS1.R
#
# Martí Galí August 2019

library(dplyr)
source('~/Desktop/GreenEdge/GCMS/f_get_ifcb4dms1.R')

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
prof.all[is.na(prof.all)] <- NA

# ------------------------------------------------------------------------
# Load IFCB logsheet. Rename cast and depth columns for matching. Sort by cast/depth
ifcb.log <- read.csv(paste0(ifcbpath,'IFCB_logsheet_amundsen2016.csv'))
nn <- colnames(ifcb.log)
nn[nn == 'Cast'] <- 'cast'
nn[nn == 'Depth_.m.'] <- 'depth'
colnames(ifcb.log) <- nn
ifcb.log <- ifcb.log[order(ifcb.log$cast, ifcb.log$depth),]

# ------------------------------------------------------------------------
# Create empty list that will hold all cast-depth ifcb data
tmp.cd <- list()

# ------------------------------------------------------------------------
# Loop through casts and depths and match to IFCB
cast <- unique(prof.all$cast)
for (jc in cast) {
  
  # Subset prof.all data [and find IFCB file match (merge): avoided because of variable n of files per depth]
  tmp.cast <- subset.data.frame(prof.all[prof.all$cast == jc,], select = c(cast,depth))
  tmp.cast$depth[tmp.cast$depth == 0.7] <- 0  # Set 0.7 to 0 m depth for exact matching of IFCB
  # tmp.cast <- merge(x = tmp.cast, y = ifcb.log, by = c('cast','depth'), all.x = T, all.y = F)
  # tmp.cast <- dplyr::left_join(x = tmp.cast, ifcb.log) # equivalent to merge
  
  # Loop through depths in cast (match unique names, then expand for repeated depths)
  cdepth <- unique(tmp.cast$depth)
  for (jd in cdepth) {
    
    # Create empty list that will hold all samples for a given depth
    tmp.cdfiles <- list()
    
    # Search all filenames
    cdfiles <- unique(ifcb.log$File_Name[ifcb.log$depth == jd & ifcb.log$cast == jc])
    cdfiles <- cdfiles[!is.na(cdfiles)]
    if (!is_empty(cdfiles)) {
      
      for (jf in cdfiles){
        fpath <- paste0(ifcbpath,'post_processing/export_144_20190111_1439/processed_counts/',
                        jf,'.csv')
        if (file.exists(fpath)) {
          ifcb.data <- read.csv(fpath, header = T)
          jf_data <- f_get_ifcb4dms1(ifcb.data)
          jf_data$cast <- jc
          jf_data$depth <- jd
          tmp.cdfiles[[ as.character(jf) ]] <- jf_data
        } # test file exists
      } # loop on file names within depth
      
      # Convert to df and average by columns
      tmp.cdfiles <- do.call("rbind",tmp.cdfiles) # elegant. Line below does the same without getting the col names, and less elegantly
      # tmp.cdfiles <- data.frame(matrix(unlist(tmp.cdfiles), nrow=length(tmp.cdfiles), byrow=T))
      tmp.cd[[ paste(jc,jd,sep = "_") ]] <- colMeans(tmp.cdfiles, na.rm = T)
      
    } # test non-empty list of files
  } # loop on depths within cast
} # loop on casts

match_ifcb <- do.call("rbind",tmp.cd)

# ------------------------------------------------------------------------
# Match prof.all and ifcb
prof.all <- merge(x = prof.all, y = match_ifcb, by = c('cast','depth'), all.x = T, all.y = F)
View(prof.all)

