# ANALYZE PHYTOPLANKTON COMMUNITY FROM IMAGING FLOW CYTOBOT (IFCB) DATA
# ONLY MEASUREMENTS MATCHING DMS(P) DURING GREEN EDGE CRUISE LEG 1B

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

# ------------------------------------------------------------------------
# Load IFCB logsheet. Rename cast and depth columns for matching. Sort by cast/depth
ifcb.log <- read.csv(paste0(ifcbpath,'IFCB_logsheet_amundsen2016.csv'))
nn <- colnames(ifcb.log)
nn[nn == 'Cast'] <- 'cast'
nn[nn == 'Depth_.m.'] <- 'depth'
colnames(ifcb.log) <- nn
ifcb.log <- ifcb.log[order(ifcb.log$cast, ifcb.log$depth),]

# ------------------------------------------------------------------------
# Define names of ifcb variables exported in next loop
# IMPORTANT: update this vector to match categories in f_get_ifcb4dms1.R, if changed
names_ifcb <- c("detritus_mg_L","diat_pelagic_mg_L","melosira_mg_L","prym_mg_L","phaeocystis_mg_L",
                "prym_clumped_mg_L","dino_mg_L","badfocus_mg_L","phaeocystis_n_mL")

# Preallocate ifcb data matching prof.all
match_ifcb <- data.frame(NA, dim = c(dim(prof.all[1]),length(names_ifcb)))
colnames(match_ifcb) <- names_ifcb

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
    
    # Preallocate here dataframe of dim [length(cdepth) * length(jf_data)]
    
    # Search all filenames
    cdfiles <- unique(ifcb.log$File_Name[ifcb.log$depth == jd])
    cdfiles <- cdfiles[!is.na(cdfiles)]
    if (!is_empty(cdfiles)) {
      
      # Preallocate here dataframe of dim [length(cdfiles) * length(jf_data)]
      tmp.cdfiles <- list()
      
      for (jf in cdfiles){
        fpath <- paste0(ifcbpath,'post_processing/export_144_20190111_1439/processed_counts/',
                        jf,'.csv')
        if (file.exists(fpath)) {
          ifcb.data <- read.csv(fpath, header = T)
          jf_data <- f_get_ifcb4dms1(ifcb.data)
          tmp.cdfiles[[ as.character(jf) ]] <- jf_data
        } # test file exists
      } # loop on file names within depth
      
      # Convert to df and sum columns
      tmp.cdfiles <- do.call("rbind",tmp.cdfiles) # elegant. Line below does the same without getting the col names, and less elegantly
      # tmp.cdfiles <- data.frame(matrix(unlist(tmp.cdfiles), nrow=length(tmp.cdfiles), byrow=T))
      sum.cdfiles <- colSums(tmp.cdfiles)
      
    } # test non-empty list of files
  } # loop on depths within cast
  
  # print(tmp.cast)
  
} # loop on casts