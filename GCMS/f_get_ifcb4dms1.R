# FUNCTION TO EXTRACT DESIRED INFORMATION FROM IFCB FILES FROM GREENEDGES
# Input: ifcb sample file as formatted by Philippe-Israel Morin
#
# Output: data frame with as many columns as taxonomic levels targetted for data aggregation,
# plus Phaeocystis counts. Currently length(out) = 10.
#
# Martí Galí August 2019

f_get_ifcb4dms1 <- function(y) {
  
  # Define taxonomic groups to be extracted
  g <- list()
  g$detritus_mg_L <- 'detritus (not-living)'
  g$diat_pelagic_mg_L <- c('centric (Bacillariophyta)',
                      'Thalassiosira (Mediophyceae)',
                      'Chaetoceros (Mediophyceae)',
                      'pennate (Bacillariophyta)',
                      'Thalassiosira nordenskioeldii (Thalassiosira)',
                      'Porosira (Mediophyceae)',
                      'Eucampia (Mediophyceae)',
                      'Bacillariophyta (Ochrophyta)',
                      'Attheya (Mediophyceae)',
                      'Cylindrotheca closterium (Cylindrotheca)',
                      'Pseudo-nitzschia (Bacillariophyceae)',
                      'Nitzschia frigida (Nitzschia)',
                      'chain (centric)',
                      'chain (pennate)',
                      'Entomoneis (Bacillariophyceae)',
                      'Polarella glacialis (Polarella)')
  g$melo_mg_L <- 'Melosira (Melosirids)'
  g$prym_mg_L <- 'Prymnesiophyceae (Haptophyta)'
  g$phaeo_mg_L <- 'Phaeocystis (Phaeocystaceae)'
  g$prym_clumped_mg_L <- 'clumped Prymnesiophyceae'
  g$dino_mg_L <- 'Dinophyceae (Holodinophyta)'
  g$badfocus_mg_L <- 'badfocus (artefact)'
  
  out <- lapply(g, function(x) sum(y$mgC_per_L[y$object_annotation_category %in% x], na.rm = T))
  out <- as.data.frame(t(unlist(out)))
  out$tot_mg_L <- sum(y$mgC_per_L, na.rm = T)
  out$phaeo_n_mL <- sum(y$validated_counts_per_ml[y$object_annotation_category %in% 'Phaeocystis (Phaeocystaceae)'], na.rm = T)
  return(out)
  
}